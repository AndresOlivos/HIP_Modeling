
# GEOSTREAMS 1.0
# DRAINAGE NETWORK MODELING FOR ARCGIS PRO 2.8.x
# BY J. ANDRES OLIVOS (Oregon State University)
# email: jaolivosh@gmail.com
# UPDATED: July 2021
# Citation: Olivos, J.A., Arismendi, I., Flitcroft, R., Penaluna, B., Firman, J. (In Prep.). Modeling the physical
# potential of riverscapes to bear animal invasions. Methods in Ecology and Evolution (?).

# -------------------------------------------------------------------------------------------------------------------- #

import arcpy, numpy, os, datetime
from arcpy import *
from arcpy.da import *
from arcpy.management import *
from arcpy.sa import *

########################################################################################################################
################################################### INPUT PARAMETERS ###################################################
########################################################################################################################

env.workspace = GetParameterAsText(0) # Input workspace geodatabase
dem = GetParameterAsText(1) # Input DEM
precip = GetParameterAsText(2) # Input precipitation raster (in cc)
hydro_regions = GetParameterAsText(3) # Input polygon with hydro-regions and regression coefficients
drain_thres = GetParameterAsText(4) # Input drainage area threshold for initial delineation of drainage lines
flow_thres = GetParameterAsText(5) # Input minimum flow of interest
reach_resolution = GetParameterAsText(6) # Input preferred reach resolution (in meters)
fill = GetParameterAsText(7) # Output fill
fdir = GetParameterAsText(8) # Input deterministic flow direction for watershed analysis (when?)
facc = GetParameterAsText(9) # Output flow accumulation raster (raw)
slope = GetParameterAsText(10) # Output slope (in %)
lines = GetParameterAsText(11) # Output vector lines
read_points = GetParameterAsText(12) # Output reading points
facc_lines = GetParameterAsText(13) # Output facc lines
facc_precip_lines = GetParameterAsText(14) # Output precip lines
facc_slope_lines = GetParameterAsText(15) # Output slope lines
order_lines = GetParameterAsText(16) # Output stream order raster
env.parallelProcessingFactor = GetParameterAsText(17) # Number of parallel processes available (int - eg 6 or percent - eg 50%)

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
env.extent = dem # All spatial processes are limited to the extent of the DEM
env.snapRaster = dem  # Set Snap Raster environment -- all rasters saved will be snapped to original DEM
env.cellSize = "MINOF"

## DEFINE OBJECTS
cell_size_result = GetRasterProperties(dem, "CELLSIZEX")# raster resolution in meters
cell_size = float(cell_size_result.getOutput(0))
cell_area = cell_size**2 # area (square meters) per pixel
pix_sqk = 1000000 / cell_area # pixels per square kilometer
pixel_thres = float(drain_thres) * pix_sqk # drain area threshold to number of pixels
dist_statement = str(reach_resolution) + " Meters"

## Define names of secondary data
facc_precip = "facc_precip"
facc_slope = "facc_slope"
lines_ras = "lines_ras"
lines_or = "lines_or"
split_points = "split_points"
temp_points = "read_points_temp"

########################################################################################################################
################################################ DRAINAGE MODELING #####################################################
########################################################################################################################
# Start time
start = datetime.datetime.now()

# Set progressor
processing_steps = 17
SetProgressor("step", "Identifying and filling sinks in the DEM...", 0, processing_steps, 1)

## Identify and fill sinks in DEM
out_surface_raster = Fill(dem, None)
out_surface_raster.save(fill)

# Update progressor
SetProgressorPosition()
SetProgressorLabel("Modeling flow accumulation (unweighted)...")

## Run raw flow accumulations
out_accumulation_raster = FlowAccumulation(fdir, None, "FLOAT", "D8")
out_accumulation_raster.save(facc)

# Extract lines
out_raster = ExtractByAttributes(facc, "VALUE > {}".format(pixel_thres))
out_raster.save(facc_lines)

# Update progressor
SetProgressorPosition()
SetProgressorLabel("Modeling flow accumulation (precipitation weighted)...")

## Run flow accumulation weighted with mean annual precipitation
# First extract and snap rasters
out_raster = ExtractByMask(precip, fill)
out_accumulation_raster = FlowAccumulation(fdir, precip, "FLOAT", "D8")
out_accumulation_raster.save(facc_precip)

# Extract lines
# Extract by mask using facc_lines to obtain facc_precip_lines
out_raster = ExtractByMask(facc_precip, facc_lines)
out_raster.save(facc_precip_lines)

# Update progressor
SetProgressorPosition()
SetProgressorLabel("Generating slope raster...")


## Run flow accumulation weighted with slope
# Generate slope surface (%)
out_raster = Slope(fill, "PERCENT_RISE", 1, "PLANAR", "METER")
out_raster.save(slope)

# Update progressor
SetProgressorPosition()
SetProgressorLabel("Modeling flow accumulation (slope weighted)...")

# Run flow acc with slope
out_accumulation_raster = FlowAccumulation(fdir, slope, "FLOAT", "D8")
out_accumulation_raster.save(facc_slope)

# Extract lines
# Extract by mask using facc_lines to obtain facc_slope_lines
out_raster = ExtractByMask(facc_slope, facc_lines)
out_raster.save(facc_slope_lines)

# End time
end = datetime.datetime.now()
duration = end - start

AddMessage("Done with flow simulations (Duration:   {} minutes)".format(str(round(duration.total_seconds()/60,1))))

########################################################################################################################
################################################# VECTORIZE NETWORK ####################################################
########################################################################################################################


# Update progressor
SetProgressorPosition()
SetProgressorLabel("Generating drainage network...")

# Con lines
out_raster = Con(facc_lines, 1, None, "VALUE > 0")
out_raster.save(lines_ras)

# Vectorize
StreamToFeature(lines_ras, fdir, lines_or, "NO_SIMPLIFY")

# Split lines based user-specified distance
edit.FlipLine(lines_or)

# Generate points at specified distance
GeneratePointsAlongLines(lines_or, split_points, "DISTANCE", dist_statement, None, "END_POINTS")

# Split lines based on those points
SplitLineAtPoint(lines_or, split_points, lines, "1 Meters")

# Assign unique IDs (HydroIDs):
AddField(lines, "HydroID", "LONG")
with UpdateCursor(lines, ["HydroID","OBJECTID"]) as cursor:
    for row in cursor:
        row[0] = row[1]
        cursor.updateRow(row)

# Create reading points
FeatureToPoint(lines, read_points, "INSIDE")

# End time
end = datetime.datetime.now()
duration = end - start

AddMessage("Done with network vectorization (Duration {} minutes)".format(str(round(duration.total_seconds()/60,1))))

########################################################################################################################
################################################# DISCHARGE MODELING ###################################################
########################################################################################################################
# Start time
start = datetime.datetime.now()

# Update progressor
SetProgressorPosition()
SetProgressorLabel("Assigning stream order...")

# ASSIGN STREAM ORDER (SHREVE)
out_ras = StreamOrder(facc_lines, fdir, "SHREVE")
out_ras.save(order_lines)

# Update progressor
SetProgressorPosition()
SetProgressorLabel("Reading flow accumulation outputs...")

## Add flow data to points
# Accumulated pixels
AddSurfaceInformation(read_points, facc_lines, "Z")
AlterField(read_points, "Z", 'DrainArea', 'DrainArea', "DOUBLE", 8, "NULLABLE", "DO_NOT_CLEAR")
SetProgressorPosition()

# Accumulated pixels weighted with precipitation
AddSurfaceInformation(read_points, facc_precip_lines, "Z")
AlterField(read_points, "Z", 'CatchPrecip', 'CatchPrecip', "DOUBLE", 8, "NULLABLE", "DO_NOT_CLEAR")
SetProgressorPosition()

# Accumulated pixels weighted with slope
AddSurfaceInformation(read_points, facc_slope_lines, "Z")
AlterField(read_points, "Z", 'CatchSlope', 'CatchSlope', "DOUBLE", 8, "NULLABLE", "DO_NOT_CLEAR")

# Update progressor
SetProgressorPosition()
SetProgressorLabel("Joining hydrological regions' data...")

# READ HYDRO REGION DATA
analysis.SpatialJoin(read_points, hydro_regions, temp_points)
CopyFeatures(temp_points, read_points, '', None, None, None)

# Update progressor
SetProgressorPosition()
SetProgressorLabel("Fitting discharge regression model...")

## CALCULATE: 'DrainArea','CatchPrecip','CatchSlope', and 'MeanFlow'
AddField(read_points, "MeanFlow", "FLOAT")
fields = ['DrainArea','CatchPrecip','CatchSlope','MeanFlow','RegionID', 'a', 'b', 'c', 'd']

with UpdateCursor(read_points, fields) as cursor:
    for row in cursor:
        # Calculate catchment area
        if row[0] is not None:
            row[0] = row[0]/pix_sqk
        else:
            row[0] = 0

        # Calculate average catchment precipitation
        if row[1] is not None:
            row[1] = (row[1]/1000)/(row[0]*pix_sqk)
        else:
            row[1] = 0

        # Calculate average catchment slope
        if row[2] is not None:
            row[2] = row[2]/(row[0]*pix_sqk)
        else:
            row[2] = 0

        # Calculate mean annual flow
        if row[4] is not None: # Only if we got regional coefficients
            row[3] = row[5] * row[0] ** row[6] * row[1] ** row[7] * row[2] ** row[8]
        cursor.updateRow(row)

# Update progressor
SetProgressorPosition()
SetProgressorLabel("Transferring data from reading points to lines...")

## JOIN DATA BACK TO STREAMS
JoinField(lines, "HydroID", read_points, "HydroID", "DrainArea;CatchPrecip;CatchSlope;MeanFlow;StreamOrder")

# Update progressor
SetProgressorPosition()
SetProgressorLabel("Removing lines below discharge thresholds...")

## DELETE ROWS BELOW FLOW THRESHOLD (from lines and points)
subflow_lines = SelectLayerByAttribute(lines, "NEW_SELECTION", "MeanFlow < {}".format(flow_thres), None)
DeleteRows(subflow_lines)
subflow_points = SelectLayerByAttribute(read_points, "NEW_SELECTION", "MeanFlow < {}".format(flow_thres), None)
DeleteRows(subflow_points)
SetProgressorPosition()

# End time
end = datetime.datetime.now()
duration = end - start

AddMessage("Done with discharge model (Duration:   {} minutes)".format(str(round(duration.total_seconds()/60,1))))

########################################################################################################################
############################################### GEOMETRICAL MODELING ###################################################
########################################################################################################################
# Start time
start = datetime.datetime.now()

# Update progressor
SetProgressorPosition()
SetProgressorLabel("Calculating line gradients...")

## GRADIENT
AddSurfaceInformation(lines, fill, "AVG_SLOPE", "BILINEAR", None, 1, 0, '')

# Update progressor
SetProgressorPosition()
SetProgressorLabel("Modeling channel width, depth, and water velocity...")

## CHANNEL WIDTH AND DEPTH
# Add necessary fields
AddFields(lines, "Depth FLOAT # # # #;ChanWidth FLOAT # # # #;Velocity FLOAT # # # #")

## Calculate bankfull depth and width based on regional regressions (Andreadis et al.)
fields = ['DrainArea','MeanFlow','Depth','ChanWidth','Avg_Slope','Velocity']

with UpdateCursor(lines, fields) as cursor:
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

SetProgressorPosition()

# End time
end = datetime.datetime.now()
duration = end - start

AddMessage("Done with geometrical model (Duration:   {} minutes)".format(str(round(duration.total_seconds()/60,1))))

########################################################################################################################
################################################ CLEAN WORKSPACE #######################################################
########################################################################################################################

# Update progressor
SetProgressorLabel("Deleting secondary files and fields...")


## Add relevant layers to current map
# Get session info
#aprx = mp.ArcGISProject("CURRENT")
#m = aprx.listMaps("Map")[0]

# List relevant layers
#layers = [lines, read_points,facc,fdir,order_lines,slope,facc_lines,facc_precip_lines,facc_slope_lines]

# Add each of them to the session
#m.addDataFromPath(lines)
#m.addDataFromPath(read_points)

## Delete secondary data
# Delete unnecessary fields
DeleteField(lines, "arcid;from_node;grid_code;to_node")
DeleteField(read_points, "Join_Count;TARGET_FID;arcid;from_node;grid_code;to_node;ORID_FID;a;b;c;d")

# Delete unnecessary files
secondary_files = [facc_precip,facc_slope, lines_ras, lines_or, temp_points, split_points]
secondary_files = ";".join(x for x in secondary_files)
Delete(secondary_files)

AddMessage("Process completed (GeoStreams).")

# Delete all user-defined Python objects
for element in dir():
    if element[0:2] != "__":
        del globals()[element]
del element