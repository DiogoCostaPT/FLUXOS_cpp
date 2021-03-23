# FLUXOS-OVERLAND (supporting scripts to examine the model outputs)
# Copyright (C) 2019-2021 Diogo Costa
#
# This file is part of FLUXOS-OVERLAND
#
# For more information see: http://www.ral.ucar.edu/projects/summa
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


from datetime import datetime
import sys
import os
import numpy as np
from joblib import Parallel, delayed
import multiprocessing
from tqdm import tqdm
# my libraries
import googleearth_klmgen as geklm
import data_management as dm
import cross_section_extract as cse
import vtk_generator as vtkgen


'''
#########################################################################################
#########################################################################################
### General Model Settings (user to provide input)
#########################################################################################
#########################################################################################
'''

# Provide directory with all the simulations to examine (or leave sim_bactch blank and define resultdir_list below
run_batch_flag = True  # TRUE: reads simulations from sim_batch_dir folder
                       # FALSE: only reads simulations listed in resultdir_list_use
'''
sim_batch_dir = '/media/dcosta/data/megasync/my_server/fluxos/TESTS_VARIA/2'
resultdir_list_select = [ # list here the simulations to examine (if you want to use this, leave
        '/media/dcosta/data/megasync/my_server/fluxos/batch_1_selected_paper_CRHM/t_36_paper_crhm/Results/',
        '/media/dcosta/data/megasync/my_server/fluxos/batch_1_selected_paper_CRHM/t_49_paper_crhm/Results/',
        '/media/dcosta/data/megasync/my_server/fluxos/batch_1_selected_paper_CRHM/t_65_paper_crhm/Results/',
    ]
'''

sim_batch_dir = '/media/dcosta/data/megasync/my_server/fluxos/Simulations/SD_Kevin_2021/4_force_pond_1_new_with_IC'
resultdir_list_select = [ # list here the simulations to examine (if you want to use this, leave
        '/media/dcosta/data/megasync/my_server/fluxos/Simulations/SD_Kevin_2021/4_force_pond_1_new_with_IC/Results' \
        ##'/media/dcosta/data/megasync/my_server/fluxos/batch_1_selected_paper_CRHM/t_49_paper_crhm/Results/', \
        ##'/media/dcosta/data/megasync/my_server/fluxos/batch_1_selected_paper_CRHM/t_65_paper_crhm/Results/', \
    ]

# DEM file
#dempath = '/media/dcosta/data/megasync/my_server/fluxos/TESTS_VARIA/2/dem_clip_SRTM_resample_500m.asc'
dempath = '/media/dcosta/data/megasync/my_server/fluxos/Simulations/SD_Kevin_2021/4_force_pond_1_new_with_IC/Rosa_2m.asc'

# Cross section shapefile (for cross-section data extraction)
#CrossSecLine_path_shapefile = \
#            '/media/dcosta/data/megasync/ec_main/models/fluxos/support/STC_data_pre-processing/00_Cross_Sections/MS9.shp'

CrossSecLine_path_shapefile = '/media/dcosta/data/megasync/ec_main/models/fluxos/support/Kevin_data_stDenis_SmithCreek/StDenis/cross_section_exit.shp'

# Define model configuration
TimeStrgStart = datetime(2011, 3, 31, 0, 0, 0) # start date (yyyy,mm,dd,HH,MM,SS)
Tinitial = 0 # start file name (numeric because it corresponds to time in seconds)
dxy = 3 # regular grid size (m)
runlag = 32400 # lag factor to account for the fact that we the model is being forced using streamflow time series

# If GOOGLE EARTH is selected, add the information below
# Coordinates for St Mary's (Janina)
coords = [43.341240, # N
          43.333936, # S
          -81.134184, # E
          -81.142515, # W
          ]
rotation = 0  # in degrees
trimrow = np.linspace(0,90,num=91).astype(int)
trimcol = np.linspace(0,0,num=0).astype(int)
resolImage = 200  # resolution of the images in Google Earth (in dpi)
var_1_graphymax = 0.12
mapoverlay_opaqueness = 50

'''
#########################################################################################
#########################################################################################
### CALCULATIONS (DON'T CANGE)
#########################################################################################
#########################################################################################
'''

# Here the user will be asked to select an option
simType = input("Options:\n"
                "# Analyse a cross-section (type 'cs')\n"
                "# GoogleEarth results overlay ('gm')\n"
                "# VTK generator ('vtk)\n Answer: ")

# Read geo to get nx and ny
geomatrix = np.loadtxt(dempath, dtype='f', delimiter=' ', skiprows=6)
ny = geomatrix.shape[0]
nx = geomatrix.shape[1]

# Read coordinates (x,y) of corner
XLL_YLL_CORNER_AND_CELL_SIZE_DEM = np.loadtxt(dempath, skiprows=2, max_rows=3, usecols=1, dtype=float)

# Identify the simulation folder to analyse
if run_batch_flag == True:
    resultdir_list_raw = [x[0] for x in os.walk(sim_batch_dir)]
    resultdir_list = []
    for i in range(0,len(resultdir_list_raw)):
        resultdir_list_candidate = resultdir_list_raw[i]
        if (resultdir_list_candidate[len(resultdir_list_candidate)-7:len(resultdir_list_candidate)]=='Results'):
            resultdir_list.append(resultdir_list_candidate + '/')
else:
    resultdir_list = resultdir_list_select

########################################################
# SECTION 1: if simType == 'cs' -> cross section analysis
#########################################################
if (simType == 'cs'):

    simType = input("Options:\n# Examine Flow (f)\n# Water Quality (wq)\n# Soil Quality (sq)\n Answer: ")

    if (simType == 'f'):
        var_col_1 = 6  # 3-h, 6-qx, 9-C, 10 - soil mass, 11 - fn_1, 12 - fe_1
        var_col_2 = 7  # 3-h, 6-qx, 9-C, 10 - soil mass, 11 - fn_1, 12 - fe_1
    elif (simType == 'wq'):
        var_col_1 = 9
        var_col_2 = 3
    elif (simType == 'sq'):
        var_col_1 = 10
        var_col_2 = 0
    else:
        sys.exit("Error: Not a valid entry (only accepts 'f', 'wq' or 'sq'")

    Overwrite_csfile = input(
        "Overwrite *.out files if existent?\n# Yes (y)\n# No (n), append new time steps to existing *.out file \n Answer: ")

    # Extract Cross-Section points from the shapefile
    shp_nodes_CS = cse.lineCSshapefile(
        CrossSecLine_path_shapefile)  # get nodes of the shapefile with the desired cross-section
    xy_CS = cse.InterpCSpoints_from_y_CS(
        shp_nodes_CS)  # interpolate a continuous line between the end nodes of "shp_nodes_CS"

    for sim in range(0, len(resultdir_list)):
        # Extract simulation name
        resultdir = resultdir_list[sim]
        try:
            resfiles_list = [x for x in os.listdir(resultdir) if x not in 'vtk' if
                             x not in 'cs']  # listing the files in folder
            print("    " + resultdir + "(found)")
        except:
            print("    " + resultdir + " (NOT found; simulation skipped)")
            continue
        try:
            os.mkdir(resultdir + "cs/")
            print("     cs directory: CREATED")
        except:
            print("     cs directory: ALREADY EXISTS")

        # read cs file if it exists
        simname = resultdir + "cs/"
        filenam = simType + ".out"
        output_dir = simname + filenam

        if (Overwrite_csfile != 'y'):
            try:
                output_dir_file = [x for x in os.listdir(simname) if x in filenam]
                csfile = open(output_dir)
                crosecval_list_exist = csfile.readlines()
                crosecval_exist = np.array([float(x) for x in crosecval_list_exist[0].split(",")])
                for i in range(1,len(crosecval_list_exist)):
                    filcontent_cs_i = np.array([float(x) for x in crosecval_list_exist[i].split(",")])
                    crosecval_exist = np.vstack((crosecval_exist, filcontent_cs_i))
                crosecval_exist_time = np.transpose(crosecval_exist)[0][:]
                csfile.close()
                # remove existing time steps in cs file from the list of output files to load
                for x in range(0, len(crosecval_exist_time)):
                    resfiles_list.remove(str(int(crosecval_exist_time[x])) + ".txt")
                file_cs_exists = True
            except:
                file_cs_exists = False
        else:
            file_cs_exists = False

        if (len(resfiles_list) == 0):
            continue

        # Extract the results over the cross-section
        crosecval = cse.csextract(simType,resultdir,resfiles_list, xy_CS, XLL_YLL_CORNER_AND_CELL_SIZE_DEM,sim, var_col_1, var_col_2, nx, ny, dxy)

        # stack existig and new crosecval
        if file_cs_exists:
            crosecval = np.vstack((crosecval_exist, crosecval))

        crosecval = crosecval[np.argsort(crosecval[:, 0])][:] # sort results (parallel returns it unsorted)

        dm.savereslt(output_dir,crosecval)

########################################################
# SECTION 2: if simType == 'im' -> animated inundation map in Google Earth
#########################################################
elif (simType == 'gm'):

    for sim in range(0, len(resultdir_list)):

        simname = dm.getsimname(resultdir_list[sim], '', simType)

        # Google Maps
        var_col_1 = 3
        geklm.google_eart_animation(resultdir_list[sim], simname, var_col_1, TimeStrgStart, Tinitial, nx, ny,
                                    dxy, coords, rotation, resolImage, var_1_graphymax, mapoverlay_opaqueness,
                                    trimrow,trimcol)


########################################################
# SECTION 2: if simType == 'vtk' -> save results in VTK format
#########################################################
elif (simType == 'vtk'):
    print ("Running VTK generator...")
    for sim in range(0, len(resultdir_list)):
        # Extract simulation name
        resultdir = resultdir_list[sim]

        try:
            resfiles_list = [x for x in os.listdir(resultdir) if x not in 'vtk' if x not in 'cs']
            print("    " + resultdir + "(found)")
        except:
            print("    " + resultdir + " (NOT found; simulation skipped)")
            continue

        try:
            os.mkdir(resultdir + "vtk/")
            print("     vtk directory: CREATE")
        except:
            print("     vtk directory: ALREADY EXISTS")
        num_cores = multiprocessing.cpu_count()
        Parallel(n_jobs=num_cores)(delayed(vtkgen.vtk_generator)(simType, resultdir, resfiles_list[res_i], dempath, nx, ny, dxy) for res_i in tqdm(range(0,len(resfiles_list))))

else:
    sys.exit("Error: Not a valid entry (only accepts 'cs', 'im' or 'vtk'")