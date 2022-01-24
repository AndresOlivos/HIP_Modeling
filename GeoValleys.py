# ---------------------------------------------------------------------------------------------------------------------#
# Tool: GeoValleys       # Toolbox: HIP Toolbox
# Description:
# Inputs: gdb workspace, drainage lines, slope raster (in %), channel depth field
# Parameters: cost multiplier for cost distance analysis, resolution of valley model (pixel size in meters), PPF
# Outputs: valley polygons
# ---------------------------------------------------------------------------------------------------------------------#
# BY J. Andres Olivos (Oregon State University)
# email: jaolivosh@gmail.com                                       # alternative contact: ivan.arismendi@oregonstate.edu
# UPDATED: July 2021
# Citation: Olivos, J.A., Arismendi, I., Flitcroft, R., Penaluna, B., Firman, J. (In Prep.). Modeling the intrinsic
# potential of riverscapes to sustain animal invasions. Methods in Ecology and Evolution (?).
# -------------------------------------------------------------------------------------------------------------------- #

from arcpy import *
from arcpy.da import *
from arcpy.management import *
from arcpy.sa import *

########################################################################################################################
################################################### INPUT PARAMETERS ###################################################
########################################################################################################################

env.workspace = GetParameterAsText(0) # Input workspace geodatabase
lines = GetParameterAsText(1) # Input drainage lines
slope = GetParameterAsText(2) # Input slope surface(in %)
capacity = GetParameterAsText(3) # Input field with line's capacity (default is Depth)
cost_multiplier = GetParameterAsText(4) # Input cost multiplier (default is 0.0025)
valleys = GetParameterAsText(5) # Output valleys
valley_res = GetParameterAsText(6) # Pixel size for valley delineation (default is slope raster)
env.parallelProcessingFactor = GetParameterAsText(7) # Input parallel processing factor

########################################################################################################################
################################################### INPUT VALIDATION ###################################################
########################################################################################################################

# Check version of ArcGIS Pro (2.8 or above)
CheckOutExtension("Spatial") # Check out ArcGIS Spatial Analyst extension license

########################################################################################################################
############################################### SECONDARY PARAMETERS ###################################################
########################################################################################################################

# SET WORKSPACE AND ENVIRONMENT
env.overwriteOutput = True
env.extent = slope # All spatial processes are limited to the extent of the slope raster
env.snapRaster = slope # Set Snap Raster environment -- all rasters saved will be snapped to slope raster
env.cellSize = valley_res

## VALLEY CONFINEMENT
# Correct slope values
out_raster = Raster(slope) + 0.0001
out_raster.save("slope_corr")

## Delineate valleys
out_allocation_raster = PathAllocation(lines, "slope_corr", None, None, "BINARY 1 45", None, "BINARY 1 -30 30", None, None, "HydroID", None, None, cost_multiplier, None, None, capacity, '');
out_allocation_raster.save("valleys_ras")

## Vectorize valleys
conversion.RasterToPolygon("valleys_ras", valleys, "NO_SIMPLIFY", "Value", "SINGLE_OUTER_PART", None)
AlterField(valleys, "gridcode", "HydroID", "HydroID", "LONG", 4, "NULLABLE", "DO_NOT_CLEAR")

# Join valley area to lines
DeleteField(lines, "Shape_Area") # First make sure Shape_Area is not a field yet
JoinField(lines, "HydroID", valleys, "HydroID", "Shape_Area")

## Calculate valley width and VWI
AddFields(lines, "VallWidth FLOAT # # # #;VWI FLOAT # # # #")
fields = ['Shape_Area', 'Shape_Length', 'ChanWidth', 'VallWidth', 'VWI']
with UpdateCursor(lines, fields) as cursor:
    for row in cursor:
        if row[0] is not None: # If valley area is not null:
            row[3] = row[0]/row[1] # valley width = area/length

        if row[2] is not None and row[3] is not None: # If chann and valley widths are not null
            row[4] = row[3]/row[2] # VWI = vall/chann

        elif row[2] is None and row[3] is not None: # If width is null, but valley is true
            row[2] = row[3] # Then channel and valley are equal
            row[4] = 1 # and VWI = 1

        if row[3] is not None and row[2] > row[3]: # If valley is smaller than channel
            row[3] = row[2] # Then valley is equal to channel

        if row[4] is None or row[4] < 1: # And if VWI is less than 1
            row[4] = 1 # Then, VWI is 1

        cursor.updateRow(row)


## Add relevant layers to current map


SetProgressorPosition()
AddMessage("Processing completed (GeoValleys).")

########################################################################################################################
################################################ CLEAN WORKSPACE #######################################################
########################################################################################################################

## Delete secondary data
# Files
secondary_files = ["slope_corr", "valleys_ras"]
secondary_files = ";".join(x for x in secondary_files)
Delete(secondary_files)

# Fields
DeleteField(lines, "Shape_Area")
DeleteField(valleys, "Id")

# Delete all user-defined Python objects
for element in dir():
    if element[0:2] != "__":
        del globals()[element]
del element

