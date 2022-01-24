
# ---------------------------------------------------------------------------------------------------------------------#
# Tool: GeoBasins       # Toolbox: HIP Toolbox
# Description:
# Inputs:
# Outputs
# BY J. ANDRES OLIVOS (Oregon State University)
# email: jaolivosh@gmail.com
# UPDATED: July 2021
# Citation: Olivos, J.A., Arismendi, I., Flitcroft, R., Penaluna, B., Firman, J. (In Prep.). Modeling the physical
# potential of riverscapes to bear animal invasions. Methods in Ecology and Evolution (?).
#
# -------------------------------------------------------------------------------------------------------------------- #

import arcpy, numpy, os
from arcpy import *
from arcpy.da import *
from arcpy.management import *
from arcpy.sa import *

########################################################################################################################
################################################### INPUT PARAMETERS ###################################################
########################################################################################################################

env.workspace = GetParameterAsText(0) # Input workspace geodatabase
fdir = GetParameterAsText(1) # Input flow direction
facc = GetParameterAsText(2) # Input flow accumulation
order_ras = GetParameterAsText(3) # Input stream order raster
min_basin_area = GetParameterAsText(4) # Input minimum basin area
segmentation_order = GetParameterAsText(5) # Input order of segmentation (by Shreve)
basins = GetParameterAsText(6) # Input basins
sub_basins = GetParameterAsText(7) # Output subbasins

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
env.extent = fdir # Spatial processing is limited to the extent of the DEM
env.snapRaster = fdir  # Set Snap Raster environment -- all rasters saved will be snapped to original DEM
env.cellSize = "MINOF"

## Define names of secondary data
basin_statement = "Shape_Area < {0}".format(str(float(min_basin_area) * 1000000))
subbasin_statement = "Value > {0}".format(segmentation_order)

########################################################################################################################
################################################ BASIN IDs         #####################################################
########################################################################################################################

# Delete basins below the minimum area threshold
tiny_basins = SelectLayerByAttribute(basins, "NEW_SELECTION", basin_statement, None)
DeleteRows(tiny_basins)

# Assign Basin_ID (four digit code)
AddField(basins, "Basin_ID", "LONG")
currentID = 999 # starts Basin_IDs at 1000

with UpdateCursor(basins, ['Basin_ID']) as cursor:
    for row in cursor:
        currentID = currentID + 1
        row[0] = currentID
        cursor.updateRow(row)

########################################################################################################################
############################################# SUB-BASIN DELINEATION ####################################################
########################################################################################################################

order_lines = Con(order_ras, 1, None, subbasin_statement)
lines_wOrder = StreamToFeature(order_lines, fdir, "lines_wOrder", "NO_SIMPLIFY")
pour_points = FeatureVerticesToPoints(lines_wOrder, "pour_points", "BOTH_ENDS")
snapped_pour_points = SnapPourPoint(pour_points, facc, 25, "OBJECTID")

## Delineate subbasins
subbasins_ras = Watershed(fdir, snapped_pour_points, "Value")


# Vectorize subbasins
sub_basins_temp = conversion.RasterToPolygon(subbasins_ras, "sub_basins_temp", "NO_SIMPLIFY", "Value", "SINGLE_OUTER_PART", None)

# Unite temporary subbasins (only at big basins) with basins to obtain sub_basins
subbasins_to_unite = "{0};{1}".format(str(basins),str(sub_basins_temp))
analysis.Union(subbasins_to_unite, sub_basins) # Try Identity instead, to remove inconsistencies at ridgelines

# Assign parent Basin_ID to sub_basins
analysis.SpatialJoin(sub_basins, basins, sub_basins_temp, "#", "#", '#', "WITHIN")
CopyFeatures(sub_basins_temp, sub_basins) # Goes around the creation of a new layer by Spatial Join tool

# Assign SubBasin_ID (six digit code)
AddField(sub_basins, "SubBasin_ID", "LONG")
currentID = 9 # Starts SubBasin_ID at 100010

with UpdateCursor(sub_basins, ['Basin_ID', 'SubBasin_ID']) as cursor:
    for row in cursor:
        currentID = currentID + 1
        row[1] = str(row[0]) + str(currentID)
        cursor.updateRow(row)


## Delete intermediate data
# unnecessary fields
DeleteField(basins, "Id; gridcode")
DeleteField(sub_basins, "Id; gridcode")

# unnecessary files
Delete(r"lines_wOrder;pour_points;sub_basins_temp")

## Add relevant layers to current map
#aprx = mp.ArcGISProject("CURRENT")
#m = aprx.listMaps("Map")[0]
#m.addDataFromPath(sub_basins)
#m.addDataFromPath(basins)

SetProgressorPosition()
AddMessage("Processing completed (GeoBasins).")

# Delete all user-defined Python objects
for element in dir():
    if element[0:2] != "__":
        del globals()[element]
del element

