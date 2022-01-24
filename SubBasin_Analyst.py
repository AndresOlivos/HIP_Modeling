# ---------------------------------------------------------------------------------------------------------------------#
# Tool: SubBasin Analyst     # Toolbox: HIP Toolbox
# Description:
# Inputs:
# Parameters:
# Outputs:
# ---------------------------------------------------------------------------------------------------------------------#
# BY J. Andres Olivos (Oregon State University)
# email: jaolivosh@gmail.com                                       # alternative contact: ivan.arismendi@oregonstate.edu
# UPDATED: July 2021
# Citation: Olivos, J.A., Arismendi, I., Flitcroft, R., Penaluna, B., Firman, J. (In Prep.). Modeling the intrinsic
# potential of riverscapes to sustain animal invasions. Methods in Ecology and Evolution (?).
# -------------------------------------------------------------------------------------------------------------------- #

import arcpy, numpy, os, statistics
from arcpy import *
from arcpy.da import *
from arcpy.management import *
from arcpy.sa import *

########################################################################################################################
################################################### INPUT PARAMETERS ###################################################
########################################################################################################################

env.workspace = GetParameterAsText(0) # Input workspace geodatabase
sub_basins = GetParameterAsText(1) # Input sub-basins
predefined_variables = GetParameterAsText(2) # Input list with predefined variables of interest (e.g., hab_complement) and output field names
fc_table = GetParameterAsText(3) # Input table with feature classes, input fields, statistics of interest, and output field names
raster_table = GetParameterAsText(4) # Input table of rasters, statistics of interest, and output field names

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
env.extent = sub_basins # Spatial processing is limited to the extent of the input sub_basins
env.cellSize = "MINOF"

## Define secondary data
## DEFINE OBJECTS
subbasin_points = "subbasin_points"

## Convert input tables to readable lists
if len(raster_table) > 0:
    raster_variables = raster_table.split(";")
    raster_variables = (x.split(" ") for x in raster_variables)
    raster_variables = [(x[0],x[1],x[2]) for x in raster_variables]

fc_variables = fc_table.split(";")
fc_variables = (x.split(" ") for x in fc_variables)

## Add all required fields
# First test that they don't exist in table and delete them if they are
existing_fields = ListFields(sub_basins)
all_fields = []

if len(predefined_variables) > 0 == len(raster_table) and len(fc_table) == 0:
    pred_fields = [i.split(" ",1)[1] for i in predefined_variables.split(";")]
    all_fields = pred_fields

if len(raster_table) > 0 and len(fc_table) == 0:
    raster_fields = [x[2] for x in raster_variables]
    all_fields = raster_fields

if len(fc_table) > 0 and len(raster_table) == 0 == len(predefined_variables):
    fc_fields = [x[1] for x in fc_variables]
    all_fields = fc_fields

if len(fc_table) > 0 < len(raster_table) and len(predefined_variables) > 0:
    all_fields = fc_fields + raster_fields + pred_fields

if len(all_fields) != 0:
    for field in all_fields:
        if field in existing_fields:
            DeleteField(sub_basins, field)



# Always add fields: (...)

# Go through input tables and identify variables of interest


########################################################################################################################
############################################## PREDEFINED VARIABLES ####################################################
########################################################################################################################


# Sub-basin (area in sqk):
CalculateGeometryAttributes(sub_basins, "Area_Sqk AREA", '', "SQUARE_KILOMETERS", None, "SAME_AS_INPUT")

# Sub-basin perimeter

#
for i in predefined_variables.split(";"):
    lines = i.split(" ",1)[0]
    variable = i.split(" ",1)[1]

    ## COHO SALMON HABITAT COMPLEMENTARITY
    if variable == "'Hab. complementarity - Coho salmon'":
        # Spawning lines
        lines_coho_s = analysis.Select(lines, "lines_coho_s", "HIP_Coho_S >= .75")
        AddFields(lines_coho_s, "Coho_S_Km FLOAT # # # #")
        fields = ["Shape_Length", "Coho_S_Km"]
        with UpdateCursor(lines_coho_s, fields) as cursor:
            for row in cursor:
                row[1] = row[0]/1000
                cursor.updateRow(row)

        # Rearing lines
        lines_coho_r = analysis.Select(lines, "lines_coho_r", "HIP_Coho_R >= .75")
        AddFields(lines_coho_r, "Coho_R_Km FLOAT # # # #")
        fields = ["Shape_Length", "Coho_R_Km"]
        with UpdateCursor(lines_coho_r, fields) as cursor:
            for row in cursor:
                row[1] = row[0]/1000
                cursor.updateRow(row)

        # Merge spawning and rearing lines
        coho_lines = Merge("{0};{1}".format(lines_coho_s,lines_coho_r),"coho_lines")

        # Join to sub-basins
        analysis.SpatialJoin(sub_basins, coho_lines, "sub_basins_temp", "#", "#", '#', "CONTAINS")
        CopyFeatures("sub_basins_temp", sub_basins, '', None, None, None)

        ## Calculate Habitat Complementarity Indexes (HCI):
        # Add fields for indexes (HCI_Coho_S, HCI_Coho_R, HCI_Coho):
        AddFields(sub_basins, "HCI_Coho_S SHORT # # # #; HCI_Coho_R SHORT ####;HCI_Coho SHORT ####")
        fields = ['Coho_S_Km','Coho_R_Km','HCI_Coho_S','HCI_Coho_R','HCI_Coho']

        ## Calculate the mean and standard deviation of the abundance of rearing habitats:
        # Collect all values in dataset
        with SearchCursor(sub_basins,fields) as cursor:
            rear_kms = []
            for row in cursor:
                if row[1] is not None:
                    rear_kms = rear_kms + [row[1]]

        # If we have enough subbasins get the mean and SD:
        if len(rear_kms) > 1:
            mean_rear_kms = statistics.mean(rear_kms) # Mean
            sd_rear_kms = statistics.stdev(rear_kms) # SD

        # If not, raise an error message:
        else:
            raise ValueError("Not enough sub-basins for the calculation of habitat complementarity scores")


        with UpdateCursor(sub_basins,fields) as cursor:
            for row in cursor:
                # Spawning index
                if row[0] is None or row[0] < 0.1:
                    row[2] = 1
                elif row[0] >= 0.1:
                    row[2] = 2

                # Rearing index
                if row[1] is None or row[1] < 0.1:
                    row[3] = 0
                elif row[1] < (mean_rear_kms - sd_rear_kms):
                    row[3] = 1
                elif (mean_rear_kms - sd_rear_kms) < row[1] < (mean_rear_kms + sd_rear_kms):
                    row[3] = 2
                elif row[1] > (mean_rear_kms + sd_rear_kms):
                    row[3] = 3

                # Final HCI Index
                row[4] = row[2] * row[3]
                cursor.updateRow(row)

        # Delete junk fields (REPLACE THIS LINE FOR THE USE OF FIELD MAPPINGS WHEN DOING SJs)
        DeleteField(sub_basins,"Avg_Slope;Barrier;Basin_ID;Basin_ID_1;CatchPrecip;CatchSlope;ChanWidth;Conn_Chinook;Conn_Coho;CW_Beaver;Depth;DrainArea;FID_basins;FID_sub_basins_temp;Flow_Chinook_R;Flow_Coho_R;Flow_Coho_S;Grad_Beaver;Grad_Chinook_R;Grad_Chinook_S;Grad_Coho_R;Grad_Coho_S;gridcode_1;gridcode_12;HIP_Beaver;HIP_Chinook_R;HIP_Chinook_S;HIP_Coho_R;HIP_Coho_S;HydroID;Id_1;Id_12;Join_Count;Join_Count_1;MeanFlow;MGD_100m;MGD_1K;Shape_Area_1;Shape_Length_1;Shape_Length_12;Single_ID;StreamOrder;TARGET_FID;TARGET_FID_1;VallWidth;Velocity;VW_Beaver;VWI;VWI_Chinook_R;VWI_Chinook_S;VWI_Coho_R;VWI_Coho_S;Width_Chinook_S")

    ## CHINOOK SALMON HABITAT COMPLEMENTARITY
    if variable == "'Hab. complementarity - Chinook salmon'":
        # Spawning lines
        lines_chinook_s = analysis.Select(lines, "lines_chinook_s", "HIP_Chinook_S >= .75")
        AddFields(lines_chinook_s, "Chinook_S_Km FLOAT # # # #")
        fields = ["Shape_Length", "Chinook_S_Km"]
        with UpdateCursor(lines_chinook_s, fields) as cursor:
            for row in cursor:
                row[1] = row[0]/1000
                cursor.updateRow(row)

        # Rearing lines
        lines_chinook_r = analysis.Select(lines, "lines_chinook_r", "HIP_Chinook_R >= .75")
        AddFields(lines_chinook_r, "Chinook_R_Km FLOAT # # # #")
        fields = ["Shape_Length", "Chinook_R_Km"]
        with UpdateCursor(lines_chinook_r, fields) as cursor:
            for row in cursor:
                row[1] = row[0]/1000
                cursor.updateRow(row)

        # Merge spawning and rearing lines
        chinook_lines = Merge("{0};{1}".format(lines_chinook_s,lines_chinook_r),"chinook_lines")

        # Join to sub-basins
        analysis.SpatialJoin(sub_basins, chinook_lines, "sub_basins_temp", "#", "#", '#', "CONTAINS")
        CopyFeatures("sub_basins_temp", sub_basins, '', None, None, None)
        
        ## Calculate Habitat Complementarity Indexes (HCI):
        # Add fields for indexes (HCI_Chinook_S, HCI_Chinook_R, HCI_Chinook):
        AddFields(sub_basins, "HCI_Chinook_S SHORT # # # #; HCI_Chinook_R SHORT ####;HCI_Chinook SHORT ####")
        fields = ['Chinook_S_Km','Chinook_R_Km','HCI_Chinook_S','HCI_Chinook_R','HCI_Chinook']

        ## Calculate the mean and standard deviation of the abundance of rearing habitats:
        # Collect all values in dataset
        with SearchCursor(sub_basins,fields) as cursor:
            rear_kms = []
            for row in cursor:
                if row[1] is not None:
                    rear_kms = rear_kms + [row[1]]

        # If we have enough subbasins get the mean and SD:
        if len(rear_kms) > 1:
            mean_rear_kms = statistics.mean(rear_kms)  # Mean
            sd_rear_kms = statistics.stdev(rear_kms)  # SD

        # If not, raise an error message:
        else:
            raise ValueError("Not enough sub-basins for the calculation of habitat complementarity scores")

        with UpdateCursor(sub_basins,fields) as cursor:
            for row in cursor:
                # Spawning index
                if row[0] is None or row[0] < 0.1:
                    row[2] = 1
                elif row[0] >= 0.1:
                    row[2] = 2

                # Rearing index
                if row[1] is None or row[1] < 0.1:
                    row[3] = 0
                elif row[1] < (mean_rear_kms - sd_rear_kms):
                    row[3] = 1
                elif (mean_rear_kms - sd_rear_kms) < row[1] < (mean_rear_kms + sd_rear_kms):
                    row[3] = 2
                elif row[1] > (mean_rear_kms + sd_rear_kms):
                    row[3] = 3

                # Final HCI Index
                row[4] = row[2] * row[3]
                cursor.updateRow(row)

        # Delete junk fields (REPLACE THIS LINE FOR THE USE OF FIELD MAPPINGS WHEN DOING SJs)
        DeleteField(sub_basins,"Avg_Slope;Barrier;CatchPrecip;CatchSlope;ChanWidth;Conn_Chinook;Conn_Coho;CW_Beaver;Depth;DrainArea;FID_basins;FID_sub_basins_temp;Flow_Chinook_R;Flow_Coho_R;Flow_Coho_S;Grad_Beaver;Grad_Chinook_R;Grad_Chinook_S;Grad_Coho_R;Grad_Coho_S;gridcode_1;gridcode_12;HIP_Beaver;HIP_Chinook_R;HIP_Chinook_S;HIP_Coho_R;HIP_Coho_S;HydroID;Id_1;Id_12;Join_Count;Join_Count_1;MeanFlow;MGD_100m;MGD_1K;Shape_Area_1;Shape_Length_1;Shape_Length_12;Single_ID;StreamOrder;TARGET_FID;TARGET_FID_1;VallWidth;Velocity;VW_Beaver;VWI;VWI_Chinook_R;VWI_Chinook_S;VWI_Coho_R;VWI_Coho_S;Width_Chinook_S")

########################################################################################################################
########################################## SUB-BASIN CHARACTERIZATION ##################################################
################################################ FEATURE CLASSES #######################################################
########################################################################################################################

## Total habitats (connected and suitable stream km.)

## Habitat complementation indexes

## Stream density (line km/subbasin area)

## Road density (road km/subbasin area)

##


########################################################################################################################
########################################## SUB-BASIN CHARACTERIZATION ##################################################
################################################ RASTER CLASSES #######################################################
########################################################################################################################

if len(raster_table) != 0:                                      ## Run this only if there are inputs in the raster table
    FeatureToPoint(sub_basins, subbasin_points, "INSIDE")       # Create subbasin reading point
    raster_fields = ";".join(x for x in raster_fields)
    for ras,stat,field in raster_variables:                 # Iterate through rasters and statistics of interest
        if stat == "PercentageArea":                            # If user wants percentage area:
            AddMessage("checkpoint")
            cell_area = float(GetRasterProperties(ras, "CELLSIZEX").getOutput(0)) ** 2 / 1000000           # pixel area in square kilometers
            current_zs = ZonalStatistics(sub_basins, "SubBasin_ID", ras, "SUM")           # Zonal statistics summing all pixels
            ExtractMultiValuesToPoints(subbasin_points, "{0} {1}".format(current_zs,field))   # Read zonal statistics raster
            with UpdateCursor(subbasin_points,[field,'Area_Sqk']) as cursor:
                for row in cursor:
                    if row[0] is not None:
                        row[0] = row[0] * cell_area / row[1] * 100 # area percentage
                    else:
                        row[0] = 0
                    cursor.updateRow(row)

        else:
            AddMessage("No statistics provided")

    # Join fields back to sub-basins
    JoinField(sub_basins, "SubBasin_ID", subbasin_points, "SubBasin_ID", raster_fields)

########################################################################################################################
################################################### CLEAN WORKSPACE ####################################################
########################################################################################################################

## Delete intermediate data
# unnecessary fields
#DeleteField(basins, "Id; gridcode")
#DeleteField(sub_basins, "Id; gridcode")

# unnecessary files
Delete(r"sub_basins_temp")

## Add relevant layers to current map

# Delete all user-defined Python objects
for element in dir():
    if element[0:2] != "__":
        del globals()[element]
del element

