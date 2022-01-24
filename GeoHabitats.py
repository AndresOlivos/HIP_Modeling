# ---------------------------------------------------------------------------------------------------------------------#
# Tool: GeoHabitats       # Toolbox: HIP Toolbox
# Description: classifies drainage lines based on geophysical variables and user-defined species/life-stage of interest
# Inputs: drainage lines (reaches) with required attributes, and species/life-stage(s) of interest
# Parameters: N/A
# Outputs: updated drainage lines
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

lines = GetParameterAsText(0) # Streams generated with GeoStreams
sps_list = GetParameterAsText(1) # Species and life-stages of interest
sps = sps_list.split(";")

# Base fields
fields = ['Avg_Slope','MeanFlow','Velocity', 'Depth', 'ChanWidth', 'VallWidth', 'VWI', 'MGD_100m', 'MGD_1K', 'Barrier']
    #ROWS:      0          1         2         3            4            5        6        7          8          9

if "'Atlantic salmon - Spawning'" in sps:
    # Add fields for salar spawning
    AddFields(lines,"Conn_Salar FLOAT # # # #;Depth_Salar_S FLOAT # # # #;Vel_Salar_S FLOAT # # # #;VWI_Salar_S FLOAT # # # #;HIP_Salar_S FLOAT # # # #")
    fields_salar_s = fields + ['Conn_Salar','Depth_Salar_S', 'Vel_Salar_S', 'VWI_Salar_S', 'HIP_Salar_S']
                             #      10          11              12              13              14
    with UpdateCursor(lines, fields_salar_s) as cursor:
        for row in cursor:
            ## CONNECTIVITY SALAR
            # Barriers
            if row[9] is None or row[9] == 0:
                barrier_score = 1
            else:
                barrier_score = 0

            # MGD_1K
            if row[8] is None or row[8] <= 10:
                mgd1K_score = 1
            elif 10 < row[8] <= 14:
                mgd1K_score = 3.5 - 0.25 * row[8]
            elif row[8] > 14:
                mgd1K_score = 0

            # MGD_100m
            if row[7] is None or row[7] <= 16:
                mgd100m_score = 1
            elif 16 < row[7] <= 20:
                mgd100m_score = 5 - 0.25 * row[7]
            elif row[7] > 20:
                mgd100m_score = 0

            # Geometric mean of the three connectivity variables
            row[10] = (mgd1K_score * mgd100m_score * barrier_score) ** (1 / 3)

            ## DEPTH SALAR SPAWNING
            if row[3] is None or row[3] <= 0.15:
                row[11] = 0
            elif 0.15 <= row[3] < 0.3:
                row[11] = -1 + 6.67 * row[3]
            elif 0.3 <= row[3] < 0.6:
                row[11] = 1
            elif 0.6 <= row[3] < 0.9:
                row[11] = 2 - 1.67 * row[3]
            elif row[3] >= 0.9:
                row[11] = 0

            ## VELOCITY SALAR SPAWNING
            if row[2] is None or row[2] <= 0.15:
                row[12] = 0.1
            elif 0.15 <= row[2] < 0.4:
                row[12] = -0.44 + 3.6 * row[2]
            elif 0.4 <= row[2] < 1.05:
                row[12] = 1
            elif 1.05 <= row[2] < 1.35:
                row[12] =  3.8 - 2.67 * row[2]
            elif row[2] >= 1.35:
                row[12] = 0.2

            ## VWI SALAR SPAWNING
            if row[6] is None or row[6] < 4:
                row[13] = 0.1
            elif 4 <= row[6] < 20:
                row[13] = -0.125 + 0.06 * row[6]
            elif 20 <= row[6] < 40:
                row[13] = 1
            elif 40 <= row[6] < 60:
                row[13] = 2 - 0.025 * row[6]
            elif row[6] >= 60:
                row[13] = 0.5

            ## HIP SALAR SPAWNING
            row[14] = min([row[10],abs(row[11] * row[12] * row[13]) ** (1 / 3)]) #abs() fixes inprecise coeff. (negative scores tending to 0)

            if row[14] > 1: # this fixes values above 1
                row[14] = 1
            cursor.updateRow(row)

if "'Atlantic salmon - Rearing'" in sps:
    # Add fields Depth_Salar_R, Vel_Salar_R, VWI_Salar_R, HIP_Salar_R, Conn_Salar?
    fields_salar_r = fields + ['Depth_Salar_R', 'Vel_Salar_R', 'VWI_Salar_R', 'HIP_Salar_R']

if "'North American beaver'" in sps:
    # Add fields 'Grad_Beaver', 'CW_Beaver', 'VW_Beaver', 'HIP_Beaver'
    AddFields(lines,
              "Grad_Beaver FLOAT # # # #;CW_Beaver FLOAT # # # #;VW_Beaver FLOAT # # # #;HIP_Beaver FLOAT # # # #")
    fields_beaver = fields + ['Grad_Beaver', 'CW_Beaver', 'VW_Beaver', 'HIP_Beaver']
    #                               10              11          12          13
    with UpdateCursor(lines, fields_beaver) as cursor:
        for row in cursor:
            ## GRADIENT BEAVER
            if row[0] is None or row[0] < 3:
                row[10] = 1
            elif 3 <= row[0] < 10:
                row[10] = 1.428571429 + -0.142857143 * row[0]
            elif row[0] >= 10:
                row[10] = 0

            ## CHANNEL WIDTH BEAVER
            if row[4] is None:
                row[11] = 0
            elif row[4] <= 7:
                row[11] = 1
            elif 7 < row[4] <= 24:
                row[11] = 1.411764706 + -0.058823529 * row[4]
            elif row[4] > 24:
                row[11] = 0

            ## VALLEY WIDTH BEAVER
            if row[5] is None or row[5] <= 10:
                row[12] = 0
            elif 10 < row[5] < 25:
                row[12] = -0.666666667 + 0.066666667 * row[5]
            elif row[5] >= 25:
                row[12] = 1

            ## HIP BEAVER
            row[13] = abs(row[10] * row[11] * row[12]) ** (1 / 3)

            if row[13] > 1:  # this fixes values above 1
                row[13] = 1

            cursor.updateRow(row)

if "'Coho salmon - Spawning'" in sps:
    # Add fields Conn_Coho, Flow_Coho_S, Grad_Coho_S, VWI_Coho_S, Spawn_Coho_S
    AddFields(lines,
                               "Conn_Coho FLOAT # # # #;Flow_Coho_S FLOAT # # # #;Grad_Coho_S FLOAT # # # #;VWI_Coho_S FLOAT # # # #;HIP_Coho_S FLOAT # # # #")

    fields_coho_s = fields + ['Conn_Coho','Flow_Coho_S', 'Grad_Coho_S', 'VWI_Coho_S', 'HIP_Coho_S']
                             #      10          11              12              13          14             15

    with UpdateCursor(lines, fields_coho_s) as cursor:
        for row in cursor:
            ## CONNECTIVITY COHO
            # Barriers
            if row[9] is None or row[9] == 0:
                barrier_score = 1
            else:
                barrier_score = 0

            # MGD_1K
            if row[8] is None or row[8] <= 8:
                mgd1K_score = 1
            elif 8 < row[8] <= 12:
                mgd1K_score = 3 - 0.25 * row[8]
            elif row[8] > 12:
                mgd1K_score = 0

            # MGD_100m
            if row[7] is None or row[7] <= 12:
                mgd100m_score = 1
            elif 12 < row[7] <= 16:
                mgd100m_score = 4 - 0.25 * row[7]
            elif row[7] > 16:
                mgd100m_score = 0

            # Geometric mean of the three connectivity variables
            row[10] = (mgd1K_score * mgd100m_score * barrier_score) ** (1 / 3)

            ## FLOW COHO SPAWNING
            if row[1] is None or row[1] < 0.01:
                row[11] = 0
            elif 0.01 <= row[1] < 0.05:
                row[11] = -0.25 + 25 * row[1]
            elif 0.05 <= row[1] < 4.5:
                row[11] = 1
            elif 4.5 <= row[1] < 9:
                row[11] = 1.4 - 0.09 * row[1]
            elif row[1] >= 9:
                row[11] = 0.6

            ## GRADIENT COHO SPAWNING
            if row[0] is None:
                row[12] = 0
            elif row[0] <= 6.5:
                row[12] = 1 - 0.15 * row[0]
            elif row[0] > 6.5:
                row[12] = 0

            ## VWI COHO SPAWNING
            if row[6] is None or row[6] <= 2:
                row[13] = 0
            elif 2 < row[6] < 5:
                row[13] = -0.67 + 0.33 * row[6]
            elif 5 <= row[6] < 21:
                row[13] = 1
            elif 21 <= row[6] < 40:
                row[13] = 1.55 - 0.03 * row[6]
            elif row[6] >= 40:
                row[13] = 0.5

            ## HIP COHO SPAWNING
            row[14] = min([row[10],abs(row[11] * row[12] * row[13]) ** (1 / 3)])
            if row[14] > 1: # this fixes values above 1
                row[14] = 1

            cursor.updateRow(row)

if "'Coho salmon - Rearing'" in sps:
    AddFields(lines,
              "Conn_Coho FLOAT # # # #;Flow_Coho_R FLOAT # # # #;Grad_Coho_R FLOAT # # # #;VWI_Coho_R FLOAT # # # #;HIP_Coho_R FLOAT # # # #")

    fields_coho_s = fields + ['Conn_Coho', 'Flow_Coho_R', 'Grad_Coho_R', 'VWI_Coho_R', 'HIP_Coho_R']
    #                             10          11              12              13          14

    with UpdateCursor(lines, fields_coho_s) as cursor:
        for row in cursor:
            ## CONNECTIVITY COHO
            # Barriers
            if row[9] is None or row[9] == 0:
                barrier_score = 1
            else:
                barrier_score = 0

            # MGD_1K
            if row[8] is None or row[8] <= 8:
                mgd1K_score = 1
            elif 8 < row[8] <= 12:
                mgd1K_score = 3 - 0.25 * row[8]
            elif row[8] > 12:
                mgd1K_score = 0

            # MGD_100m
            if row[7] is None or row[7] <= 12:
                mgd100m_score = 1
            elif 12 < row[7] <= 16:
                mgd100m_score = 4 - 0.25 * row[7]
            elif row[7] > 16:
                mgd100m_score = 0

            # Geometric mean of the three connectivity variables
            row[10] = (mgd1K_score * mgd100m_score * barrier_score) ** (1 / 3)

            ## FLOW COHO REARING
            if row[1] is None or row[1] < 0.01:
                row[11] = 0
            elif 0.01 <= row[1] < 0.06:
                row[11] = -0.2 + 20 * row[1]
            elif 0.06 <= row[1] < 21.24:
                row[11] = 1
            elif 21.24 <= row[1] < 76.45:
                row[11] = 1.192356457 - 0.00905633 * row[1]
            elif row[1] >= 76.45:
                row[11] = 0.5

            ## GRADIENT COHO REARING
            if row[0] is None:
                row[12] = 0
            elif row[0] <= 5:
                row[12] = 1 - 0.2 * row[0]
            elif row[0] > 5:                # Override of envelope aggregation operator
                row[12] = 0

            ## VWI COHO REARING
            if row[6] is None or row[6] <= 5.06:
                row[13] = 0.25
            elif 5.06 < row[6] < 8.86:
                row[13] = -0.748684211 + 0.197368421 * row[6]
            elif row[6] >= 8.86:
                row[13] = 1

            ## HIP COHO REARING
            row[14] = min([row[10], abs(row[11] * row[12] * row[13]) ** (1 / 3)])
            if row[14] > 1:  # this fixes values above 1
                row[14] = 1
            cursor.updateRow(row)

if "'Chinook salmon - Spawning'" in sps:
    # Add fields Flow_Chinook_R, Grad_Chinook_R, VWI_Chinook_R, HIP_Chinook_R
    AddFields(lines,
              "Conn_Chinook FLOAT # # # #;Width_Chinook_S FLOAT # # # #;Grad_Chinook_S FLOAT # # # #;VWI_Chinook_S FLOAT # # # #;HIP_Chinook_S FLOAT # # # #")

    fields_chinook_s = fields + ['Conn_Chinook', 'Width_Chinook_S', 'Grad_Chinook_S', 'VWI_Chinook_S', 'HIP_Chinook_S']
    #                                10                 11              12              13                  14

    with UpdateCursor(lines, fields_chinook_s) as cursor:
        for row in cursor:
            ## CONNECTIVITY CHINOOK
            # Barriers
            if row[9] is None or row[9] == 0:
                barrier_score = 1
            else:
                barrier_score = 0

            # MGD_100m
            if row[7] is None or row[7] <= 14:
                mgd100m_score = 1
            elif 14 < row[7] <= 18:
                mgd100m_score = 4.5 - 0.25 * row[7]
            elif row[7] > 18:
                mgd100m_score = 0

            # MGD_1K
            if row[8] is None or row[8] <= 8:
                mgd1K_score = 1
            elif 8 < row[8] <= 12:
                mgd1K_score = 3 - 0.25 * row[8]
            elif row[8] > 12:
                mgd1K_score = 0

            # Geometric mean of the three connectivity variables
            row[10] = (mgd1K_score * mgd100m_score * barrier_score) ** (1 / 3)

            ## CHANNEL WIDTH CHINOOK SPAWNING
            if row[4] is None or row[4] < 3.7:
                row[11] = 0
            elif 3.7 <= row[4] < 5.7:
                row[11] = -1.85 + 0.5 * row[1]
            elif row[4] >= 5.7:
                row[11] = 1

            ## GRADIENT CHINOOK SPAWNING
            if row[0] is None:
                row[12] = 0
            elif row[0] < 2:
                row[12] = 1
            elif 2 < row[0] <= 4:
                row[12] = 2 - 0.5 * row[0]
            elif row[0] > 4:
                row[12] = 0

            ## VWI CHINOOK SPAWNING
            if row[6] is None or row[6] <= 1:
                row[13] = 0.75
            elif 1 < row[6] < 8.87:
                row[13] = 0.718233799 + 0.031766201 * row[6]
            elif row[6] >= 8.87:
                row[13] = 1

            ## HIP CHINOOK SPAWNING
            row[14] = min([row[10], abs(row[11] * row[12] * row[13]) ** (1 / 3)])
            if row[14] > 1:  # this fixes values above 1
                row[14] = 1
            cursor.updateRow(row)

if "'Chinook salmon - Rearing'" in sps:
    # Add fields Flow_Chinook_R, Grad_Chinook_R, VWI_Chinook_R, HIP_Chinook_R
    AddFields(lines,
              "Conn_Chinook FLOAT # # # #;Flow_Chinook_R FLOAT # # # #;Grad_Chinook_R FLOAT # # # #;VWI_Chinook_R FLOAT # # # #;HIP_Chinook_R FLOAT # # # #")

    fields_chinook_r = fields + ['Conn_Chinook', 'Flow_Chinook_R', 'Grad_Chinook_R', 'VWI_Chinook_R', 'HIP_Chinook_R']
    #                                10                 11              12              13                  14

    with UpdateCursor(lines, fields_chinook_r) as cursor:
        for row in cursor:
            ## CONNECTIVITY CHINOOK
            # Barriers
            if row[9] is None or row[9] == 0:
                barrier_score = 1
            else:
                barrier_score = 0

            # MGD_100m
            if row[7] is None or row[7] <= 14:
                mgd100m_score = 1
            elif 14 < row[7] <= 18:
                mgd100m_score = 4.5 - 0.25 * row[7]
            elif row[7] > 18:
                mgd100m_score = 0

            # MGD_1K
            if row[8] is None or row[8] <= 8:
                mgd1K_score = 1
            elif 8 < row[8] <= 12:
                mgd1K_score = 3 - 0.25 * row[8]
            elif row[8] > 12:
                mgd1K_score = 0

            # Geometric mean of the three connectivity variables
            row[10] = (mgd1K_score * mgd100m_score * barrier_score) ** (1 / 3)

            ## FLOW CHINOOK REARING
            if row[1] is None or row[1] < 0.1:
                row[11] = 0
            elif 0.1 <= row[1] < 1.5:
                row[11] = -0.071428571 + 0.714285714 * row[1]
            elif row[1] >= 1.5:
                row[11] = 1

            ## GRADIENT CHINOOK REARING
            if row[0] is None:
                row[12] = 0
            elif row[0] < 1.75:
                row[12] = 1
            elif 1.75 < row[0] <= 4:                            # Override of envelope aggregator operator
                row[12] = 1.777777778-0.444444444 * row[0]
            elif row[0] > 4:
                row[12] = 0

            ## VWI CHINOOK REARING
            if row[6] is None or row[6] <= 1:
                row[13] = 0.8
            elif 1 < row[6] < 4:
                row[13] = 0.733333333 + 0.066666667 * row[6]
            elif row[6] >= 4:
                row[13] = 1

            ## HIP CHINOOK REARING
            row[14] = min([row[10], abs(row[11] * row[12] * row[13]) ** (1 / 3)])
            if row[14] > 1:  # this fixes values above 1
                row[14] = 1
            cursor.updateRow(row)

if "'Steelhead/Rainbow trout - Spawning'" in sps:
    # Add fields Conn_Steelhead, CWidth_Steelhead_S, Grad_Steelhead_S, VWI_Steelhead_S, Spawn_Steelhead_S
    AddFields(lines,
                               "Conn_Steelhead FLOAT # # # #;CWidth_Steelhead_S FLOAT # # # #;Grad_Steelhead_S FLOAT # # # #;VWI_Steelhead_S FLOAT # # # #;HIP_Steelhead_S FLOAT # # # #")

    fields_coho_s = fields + ['Conn_Steelhead','CWidth_Steelhead_S', 'Grad_Steelhead_S', 'VWI_Steelhead_S', 'HIP_Steelhead_S']
                             #      10                   11              12                      13              14

    with UpdateCursor(lines, fields_coho_s) as cursor:
        for row in cursor:
            ## CONNECTIVITY STEELHEAD
            # Barriers
            if row[9] is None or row[9] == 0:
                barrier_score = 1
            else:
                barrier_score = 0

            # MGD_1K
            if row[8] is None or row[8] <= 10:
                mgd1K_score = 1
            elif 10 < row[8] <= 14:
                mgd1K_score = 3.5 - 0.25 * row[8]
            elif row[8] > 14:
                mgd1K_score = 0

            # MGD_100m
            if row[7] is None or row[7] <= 16:
                mgd100m_score = 1
            elif 16 < row[7] <= 20:
                mgd100m_score = 5 - 0.25 * row[7]
            elif row[7] > 20:
                mgd100m_score = 0

            # Geometric mean of the three connectivity variables
            row[10] = (mgd1K_score * mgd100m_score * barrier_score) ** (1 / 3)

            ## BANKFULL WIDTH STEELHEAD
            if row[4] is None or row[4] < 2.2:
                row[11] = 0
            elif 2.2 <= row[4] < 3.8:
                row[11] = -1.375 + 0.625 * row[4]
            elif 3.8 <= row[4] < 25:
                row[11] = 1
            elif 25 <= row[4] < 50:
                row[11] = 1.9 - 0.036 * row[4]
            elif row[4] >= 50:
                row[11] = 0.1

            ## GRADIENT STEELHEAD
            if row[0] is None or row[0] == 0:
                row[12] = 0
            elif 0 < row[0] <= 0.5:
                row[12] = 2 * row[0]
            elif 0.5 < row[0] <= 4:
                row[12] = 1
            elif 4 < row[0] <= 8:
                row[12] = 2 - 0.25 * row[0]
            elif row[0] > 8:
                row[12] = 0

            ## VWI STEELHEAD
            if row[6] is None or row[6] <= 2:
                row[13] = 0
            elif 2 < row[6] < 5:
                row[13] = -0.67 + 0.33 * row[6]
            elif 5 <= row[6] < 21:
                row[13] = 1
            elif 21 <= row[6] < 40:
                row[13] = 1.55 - 0.03 * row[6]
            elif row[6] >= 40:
                row[13] = 0.5

            ## HIP STEELHEAD
            row[14] = min([row[10],abs(row[11] * row[12] * row[13]) ** (1 / 3)])

            if row[14] > 1: # this fixes values above 1 (due to approximations in classifications)
                row[14] = 1
            cursor.updateRow(row)

if "'Steelhead/Rainbow trout - Rearing'" in sps:
    # Add fields Flow_Steelhead_R, Grad_Steelhead_R, VWI_Steelhead_R, HIP_Steelhead_R
    fields = fields + ['Flow_Steelhead_R', 'Grad_Steelhead_R', 'VWI_Steelhead_R', 'HIP_Steelhead_R']

SetProgressorPosition()
AddMessage("Processing completed (GeoHabitats).")

########################################################################################################################
################################################ CLEAN WORKSPACE #######################################################
########################################################################################################################

# Delete all user-defined Python objects
for element in dir():
    if element[0:2] != "__":
        del globals()[element]
del element


