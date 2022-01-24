# ---------------------------------------------------------------------------------------------------------------------#
# Tool: GeoAnadromy      # Toolbox: HIP Toolbox
# Description: calculates the maximum gradient occurring between every reach and the ocean
# Inputs: drainage lines, reading points, documented barriers, flow direction, lines segmented at different scale (opt.)
# Parameters: Maximum flow included in model (excludes mainstem artifacts), range of gradients to model, PPF
# Outputs: updated drainage lines with user-defined field name
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
################################################## INPUT PARAMETERS ####################################################
########################################################################################################################

lines = GetParameterAsText(0) # Input lines dataset
read_points = GetParameterAsText(1) # Input reading_points
barriers = GetParameterAsText(2) # Input documented barriers
fdir = GetParameterAsText(3) # Input flow direction raster (D8)
max_flow = GetParameterAsText(4) # Max. flow to model for exclusion of mainstem artifacts (Default is 500 cms)
slope_low = GetParameterAsText(5) # Lower bound of interest (X where Y starts decreasing from 1 to 0)
slope_high = GetParameterAsText(6) # Upper bound of interest (X where Y is 0)
field_name = GetParameterAsText(7) # Output field with maximum gradient downstream
old_lines = GetParameterAsText(8)  # Optional input for joining connectivity to original line resolution
env.parallelProcessingFactor = GetParameterAsText(9) # Input parallel processing factor

########################################################################################################################
############################################### SECONDARY PARAMETERS ###################################################
########################################################################################################################
# SET WORKSPACE AND ENVIRONMENT
env.overwriteOutput = True
env.extent = fdir # All spatial processes are limited to the extent of the DEM
env.snapRaster = fdir  # Set Snap Raster environment -- all rasters saved will be snapped to original DEM
env.cellSize = "MINOF"

breaks = [x for x in range(int(slope_low),(int(slope_high)+1),1)] # Range of gradients of interest for anadromous migrations
mgd_fields = [(field_name + '_' + str(x)) for x in breaks]
x_fc = 'x_RP' # Temporary name for reaches above breaks
fields = ["Barrier"] + mgd_fields + [field_name]
break_n = len(breaks)
mgd_rows = [i+1 for i in range(break_n)]

# Clean old secondary fields if present
old_fields = ";".join(x for x in mgd_fields)
DeleteField(read_points, old_fields)

########################################################################################################################
############################################### CONNECTIVITY MODELING ##################################################
########################################################################################################################

# Start by modeling connectivity based on documented barriers
out_raster = Watershed(fdir, barriers)
AddSurfaceInformation(read_points, out_raster, "Z", "BILINEAR", None, 1, 0, '')
DeleteField(read_points, "Barrier")
AlterField(read_points, "Z", "Barrier", "Barrier", "DOUBLE", 8, "NULLABLE", "DO_NOT_CLEAR")

# Set progressor
processing_steps = len(breaks) + 1
SetProgressor("step", "Calculating maximum gradients downstream...", 0, processing_steps, 1)

# Join slope data from lines to read points
DeleteField(read_points, "Avg_Slope")
JoinField(read_points, "HydroID", lines, "HydroID", "Avg_Slope")

# Select all reaches above each slope break, run the Watershed tool and read the results
for x in breaks:
    statement = 'MeanFlow < {0} And Avg_Slope > {1}'.format(max_flow, str(x))
    analysis.Select(read_points, x_fc, statement)
    output_name = 'Area_MGD_' + str(x)
    CalculateField(x_fc, "grad_code", x, "PYTHON3", '', "SHORT", "NO_ENFORCE_DOMAINS")
    out_raster = Watershed(fdir, x_fc, "grad_code")
    AddSurfaceInformation(read_points, out_raster, "Z", "BILINEAR", None, 1, 0, '')
    current_field = field_name + "_" + str(x)
    AlterField(read_points, "Z", current_field, current_field, "DOUBLE", 8, "NULLABLE", "DO_NOT_CLEAR")

    # Update progressor
    SetProgressorPosition()
    AddMessage("Done with knickpoints above {0}%...".format(str(x)))

## Delete previous existing field containung Maximum Gradients Downstream at the same scale and recreate it
DeleteField(read_points, field_name)
AddField(read_points, field_name, "SHORT", None, None, None, '', "NULLABLE", "NON_REQUIRED", '')

SetProgressorPosition()
AddMessage("Updating field {}".format(field_name))

# Update cursor with a field collecting barrier and MGD data (eg. MGD_1K):
with UpdateCursor(read_points, fields) as cursor:
    for row in cursor:
        # if barrier field is not null, convert to 1 (Barrier present downstream)
        if row[0] is not None:
            row[0] = 1
        else: # Barrier absent downstream
            row[0] = 0

        # Collect gradient data
        GDs = [0]
        for i in mgd_rows:
            GDs = GDs + [row[i]]
        GDs = [0 if x is None else x for x in GDs]  ## Change nulls for zeros
        row[-1] = max(GDs)
        cursor.updateRow(row)

## Join MGD data back to respective lines
DeleteField(lines, field_name)
DeleteField(lines, "Barrier")
JoinField(lines, "HydroID", read_points, "HydroID", "Barrier;{}".format(field_name))

## Optional, share MGD fields between two scales:
if len(old_lines) > 0:
    # From the previous scale (old_lines) to the current (lines)
    analysis.SpatialJoin(lines, old_lines, "temp_new_lines", "#", "#", "#", "SHARE_A_LINE_SEGMENT_WITH")
    new_lines_fields = ListFields("temp_new_lines")
    for field in new_lines_fields:
        if field.name[-2:] == "_1":
            DeleteField("temp_new_lines", field.name)


    # From the current scale (lines) to the previous scale (current_lines)
    analysis.SpatialJoin(old_lines, lines, "temp_old_lines", "#", "#", "#", "SHARE_A_LINE_SEGMENT_WITH")
    old_lines_fields = ListFields("temp_old_lines")
    for field in old_lines_fields:
        if field.name[-2:] == "_1":
            DeleteField("temp_old_lines", field.name)

    CopyFeatures("temp_new_lines", lines, '', None, None, None)
    CopyFeatures("temp_old_lines", old_lines, '', None, None, None)

SetProgressorPosition()
AddMessage("Processing completed (GeoAnadromy).")

########################################################################################################################
################################################ CLEAN WORKSPACE #######################################################
########################################################################################################################

## Delete secondary data
# Files
secondary_files = [x_fc, "temp_new_lines", "temp_old_lines"]
secondary_files = ";".join(x for x in secondary_files)
Delete(secondary_files)

# Fields
secondary_fields = mgd_fields + ["Join_Count", "TARGET_FID"]
secondary_fields = ";".join(x for x in secondary_fields)
DeleteField(lines, secondary_fields)
DeleteField(read_points, secondary_fields)

# Delete all user-defined Python objects
for element in dir():
    if element[0:2] != "__":
        del globals()[element]
del element

