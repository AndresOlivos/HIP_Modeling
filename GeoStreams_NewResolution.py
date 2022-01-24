
# GeoStreams - New scale - Sub tool for modeling additional resolutions
# GEOHYDROLOGICAL MODELING PROTOCOL FOR ARCGIS PRO 2.8.x
# BY J. ANDRES OLIVOS (Oregon State University)
# email: jaolivosh@gmail.com
# UPDATED: July 2021
# Citation: Olivos, J.A., Arismendi, I., Flitcroft, R., Penaluna, B., Firman, J. (In Prep.). Modeling the physical
# potential of riverscapes to bear animal invasions. Methods in Ecology and Evolution (?).

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
fill = GetParameterAsText(1) # Input filled DEM
fdir = GetParameterAsText(2) # Input flow direction
lines = GetParameterAsText(3) # Input lines
hydro_regions = GetParameterAsText(4) # Input polygon with hydro-regions and regression coefficients
reach_resolution = GetParameterAsText(5) # Input preferred reach resolution (in meters)
new_lines = GetParameterAsText(6) # Output lines with new resolution
new_points = GetParameterAsText(7) # Output new reading points associated to new lines

########################################################################################################################
################################################### INPUT VALIDATION ###################################################
########################################################################################################################

# Check version of ArcGIS Pro (2.8 or above)
CheckOutExtension("Spatial") # Check out ArcGIS Spatial Analyst extension license

# Check spatial references (projected and equal area)
# Check linear and vertical units (in meters)
# Snapped rasters and resampled rasters
# Check hydro_regions has necessary fields:
#('RegionID', 'First_B', 'DrainArea_B', 'CatchPrecip_B', 'CatchSlope_B')

########################################################################################################################
############################################### SECONDARY PARAMETERS ###################################################
########################################################################################################################

# SET WORKSPACE AND ENVIRONMENT
env.overwriteOutput = True

## DEFINE OBJECTS
cell_size_result = GetRasterProperties(fill, "CELLSIZEX")# raster resolution in meters
cell_size = float(cell_size_result.getOutput(0))
cell_area = cell_size**2 # area (square meters) per pixel
pix_sqk = 1000000 / cell_area # pixels per square kilometer
dist_statement = str(reach_resolution) + " Meters"

## Define names of secondary data
lines_ras = "lines_ras"
facc_lines = "facc_lines"
facc_precip_lines = "facc_precip_lines"
facc_slope_lines = "facc_slope_lines"
order_ras = "order_lines"
split_points = "split_points"
temp_points = "read_points_temp"
lines_dissolved = "lines_dissolved"
lines_or = "lines_or"

########################################################################################################################
################################################# SPLIT NEW LINES # ####################################################
########################################################################################################################

# Update progressor
SetProgressorPosition()
SetProgressorLabel("Generating new drainage network...")

## Create and segment new network
# Rasterize input lines
AddField(lines, "Single_ID", "SHORT")
with UpdateCursor(lines, ["Single_ID"]) as cursor:
    for row in cursor:
        row[0] = 1
        cursor.updateRow(row)
conversion.PolylineToRaster(lines, "Single_ID", lines_ras, "MAXIMUM_LENGTH", "NONE", fdir, "DO_NOT_BUILD")

# Vectorize
StreamToFeature(lines_ras, fdir, lines_or, "NO_SIMPLIFY")
GeneratePointsAlongLines(lines_or, split_points, "DISTANCE", dist_statement, None, "END_POINTS")
new_lines = SplitLineAtPoint(lines_or, split_points, new_lines, "1 Meters")

# Assign unique IDs (HydroIDs):
AddField(new_lines, "HydroID", "LONG")
with UpdateCursor(new_lines, ["HydroID","OBJECTID"]) as cursor:
    for row in cursor:
        row[0] = row[1]
        cursor.updateRow(row)

# Create reading points
FeatureToPoint(new_lines, new_points, "INSIDE")

########################################################################################################################
################################################# DISCHARGE MODELING ###################################################
########################################################################################################################

# Update progressor
SetProgressorPosition()
SetProgressorLabel("Re-calculating network discharge...")

# READ RAW FLOW ACC
AddSurfaceInformation(new_points, facc_lines, "Z")
AlterField(new_points, "Z", 'DrainArea', 'DrainArea', "DOUBLE", 8, "NULLABLE", "DO_NOT_CLEAR")

# READ FLOW ACC WITH PRECIPITATION
AddSurfaceInformation(new_points, facc_precip_lines, "Z")
AlterField(new_points, "Z", 'CatchPrecip', 'CatchPrecip', "DOUBLE", 8, "NULLABLE", "DO_NOT_CLEAR")

# READ FLOW ACC WITH SLOPE
AddSurfaceInformation(new_points, facc_slope_lines, "Z")
AlterField(new_points, "Z", 'CatchSlope', 'CatchSlope', "DOUBLE", 8, "NULLABLE", "DO_NOT_CLEAR")

# READ HYDRO REGION DATA
analysis.SpatialJoin(new_points, hydro_regions, temp_points)
CopyFeatures(temp_points, new_points, '', None, None, None)

## CALCULATE: 'DrainArea','CatchPrecip','CatchSlope', and 'MeanFlow'
AddField(new_points, "MeanFlow", "FLOAT")
fields = ['DrainArea','CatchPrecip','CatchSlope','MeanFlow','RegionID', 'a', 'b', 'c', 'd']

with UpdateCursor(new_points, fields) as cursor:
    for row in cursor:
        # Test if there are null values in the table
        input_values = [row[0],row[1],row[2],row[4],row[5],row[6],row[7],row[8]]
        if None not in input_values:
            # Calculate catchment area
            row[0] = row[0]/pix_sqk

            # Calculate average catchment precipitation
            row[1] = (row[1]/1000)/(row[0]*pix_sqk)

            # Calculate average catchment slope
            row[2] = row[2]/(row[0]*pix_sqk)

            # Calculate mean annual flow
            row[3] = row[5] * row[0] ** row[6] * row[1] ** row[7] * row[2] ** row[8]

        cursor.updateRow(row)

# ASSIGN STREAM ORDER (STRAHLER)
AddSurfaceInformation(new_points, order_ras, "Z")
AlterField(new_points, "Z", 'StreamOrder', 'StreamOrder', "DOUBLE", 8, "NULLABLE", "DO_NOT_CLEAR")

## JOIN DATA BACK TO STREAMS
JoinField(new_lines, "HydroID", new_points, "HydroID", "DrainArea;CatchPrecip;CatchSlope;MeanFlow")

########################################################################################################################
############################################### GEOMETRICAL MODELING ###################################################
########################################################################################################################

# Update progressor
SetProgressorPosition()
SetProgressorLabel("Modeling in-stream geophysical conditions...")

## GRADIENT
AddSurfaceInformation(new_lines, fill, "AVG_SLOPE", "BILINEAR", None, 1, 0, '')

## CHANNEL WIDTH AND DEPTH
# Add necessary fields
AddFields(new_lines, "Depth FLOAT # # # #;ChanWidth FLOAT # # # #;Velocity FLOAT # # # #")

## Calculate bankfull depth and width based on regional regressions (Andreadis et al.)
fields = ['DrainArea','MeanFlow','Depth','ChanWidth','Avg_Slope','Velocity']

with UpdateCursor(new_lines, fields) as cursor:
    for row in cursor:
        ## Calculate width and depth
        if row[1] is not None:
            row[2] = 0.27 * row[1] ** 0.3 # bankfull depth by Andreadis et al. 2013
            row[3] = 7.2 * row[1] ** 0.5 # channel width by Andreadis et al. 2013
        else:
            row[2] = 0
            row[3] = 0
            row[5] = 0

        ## Calculate velocity
        if row[2] > 0 < row[3]:
            # Get roughness index (r)
            if row[2] is not None and row[3] is not None:
                r = (row[2] * row[3]) / (row[2] * 2 + row[3])

            # Get n coefficient
            if row[4] <= 8:
                n = 0.05
            elif row[4] > 8 and row[3] < 30:
                n = 0.03
            elif row[4] > 8 and row[3] > 30:
                n = 0.025

            row[5] = (r ** 0.66 * (row[4] / 100) ** 0.5) / n # Velocity according to Manning's equation (in meters per second)
        else:
            row[5] = 0
        cursor.updateRow(row)

########################################################################################################################
################################################ CLEAN WORKSPACE #######################################################
########################################################################################################################

# Update progressor
SetProgressorPosition()
SetProgressorLabel("Deleting secondary files and fields...")

## Delete unneccesary fields
DeleteField(lines, "Single_ID")
DeleteField(new_lines, "arcid;from_node;grid_code;to_node;Single_ID")
DeleteField(new_points, "Join_Count;TARGET_FID;arcid;from_node;grid_code;to_node;ORID_FID;a;b;c;d")

## Delete secondary data
secondary_files = [temp_points, lines_dissolved, lines_or]
secondary_files = ";".join(x for x in secondary_files)
Delete(secondary_files)

## Add relevant layers to current map
#aprx = mp.ArcGISProject("CURRENT")
#m = aprx.listMaps("Map")[0]
#m.addDataFromPath(new_lines)

SetProgressorPosition()
AddMessage("Processing completed (GeoStreams - New resolution).")

# Delete all user-defined Python objects
for element in dir():
    if element[0:2] != "__":
        del globals()[element]
del element