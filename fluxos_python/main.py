# -------------------------------------------
# -------------------------------------------
#  MAIN FUNCTION
# -------------------------------------------
# -------------------------------------------

from datetime import datetime
# my libraries
import googleearth_klmgen as geklm
import data_management as dm
import cross_section_extract as cse
import graphing_functions as grph
import sys
import vtk_generator as vtkgen
import os
import numpy as np
from joblib import Parallel, delayed
import multiprocessing
from tqdm import tqdm

'''
General Model Settings
'''

#sim_batch_dir = '/media/dcosta/data/megasync/my_server/fluxos/Janina_batch_1/Essex_1/'
#dempath = '/media/dcosta/data/megasync/my_server/fluxos/Janina_batch_1/Essex_1/Essex_DEM_ascii'

#sim_batch_dir = '/media/dcosta/data/megasync/my_server/fluxos/Janina_batch_1/St_Marys_2/'
#dempath = '/media/dcosta/data/megasync/my_server/fluxos/Janina_batch_1/St_Marys_2/St_Marys_DEM_corrected.asc'

#sim_batch_dir = '/media/dcosta/data/megasync/my_server/fluxos/batch_1_select_paper/'
#sim_batch_dir = '/media/dcosta/data/megasync/my_server/fluxos/batch_1_selected_paper_additional/'
#sim_batch_dir = '/media/dcosta/data/megasync/my_server/fluxos/batch_1_selected_paper_additional_noIC/'
#sim_batch_dir = '/media/dcosta/data/megasync/my_server/fluxos/batch_1_selected_paper_additional_previousFLUXOSversion/'
#sim_batch_dir = '/media/dcosta/data/megasync/my_server/fluxos/batch_1_selected_paper_additional_previousFLUXOSversion/'
#sim_batch_dir = '/media/dcosta/data/megasync/my_server/fluxos/batch_1_selected_paper_additional_previousFLUXOSversion_graham/'
#sim_batch_dir = '/media/dcosta/data/megasync/my_server/fluxos/batch_1_selected_paper_additional_noIC_graham/'

#dempath = '/media/dcosta/data/megasync/my_server/fluxos/batch_1_select_paper/t_36_paper/ersi_grid_dem_alldomain_basinwalls'

sim_batch_dir = '/media/dcosta/data/megasync/my_server/fluxos/TESTS_VARIA/2'
dempath = '/media/dcosta/data/megasync/my_server/fluxos/TESTS_VARIA/2/dem_clip_SRTM_resample_500m.asc'


try:
    resultdir_list_raw = [x[0] for x in os.walk(sim_batch_dir)]
    resultdir_list = []
    for i in range(0,len(resultdir_list_raw)):
        resultdir_list_candidate = resultdir_list_raw[i]
        if (resultdir_list_candidate[len(resultdir_list_candidate)-7:len(resultdir_list_candidate)]=='Results'):
            resultdir_list.append(resultdir_list_candidate + '/')


except:
    # resultdir_list = [
    #                '/media/dcosta/DATADRIVE1/fluxos_tests/SIMULATIONS_sync/batch_1/t_36/Results/',
    #                '/media/dcosta/DATADRIVE1/fluxos_tests/SIMULATIONS_sync/batch_1/t_49/Results/',
    #                '/media/dcosta/DATADRIVE1/fluxos_tests/SIMULATIONS_sync/batch_1/t_65/Results/',
    #                ]
    resultdir_list = [
        '/media/dcosta/data/megasync/my_server/fluxos/batch_1_selected_paper_CRHM/t_36_paper_crhm/Results/',
        '/media/dcosta/data/megasync/my_server/fluxos/batch_1_selected_paper_CRHM/t_49_paper_crhm/Results/',
        '/media/dcosta/data/megasync/my_server/fluxos/batch_1_selected_paper_CRHM/t_65_paper_crhm/Results/',
    ]

TimeStrgStart = datetime(2011, 3, 31, 0, 0, 0)
Tinitial = 0
Timee = 18000 #1468800
t_step_read = 3600

dxy = 3
runlag = 32400 # lag factor to account for the fact that we the model is being forced using streamflow time series

# Coordinates for STC (3m resolution) - GOOGLE EARTH
#coords = [49.339205,  # N
#          49.317788,  # S
#          -98.342849,  # E
#          -98.402657]  # W
#rotation = 0  # in degrees
#nanrow = []
#nancol = []

# Coordinates for Essex (Janina)
#coords = [42.135017,  # N
#          42.128802,  # S
#          -82.790974,  # E
#          -82.795857]  # W
#rotation = 0  # in degrees
#nanrow = []
#nancol = []

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
Choose the functions to run
'''
# Sim results analysis
Map3d_Dem_and_sim_plotly = False # plot 3d maps of dem and sim results in plotly (interactive)
GoogleEarthSim = True # export kml to google earth
Quiver_flowpaths = False

simType = input("Options:\n# Analyse a cross-section (type 'cs')\n# Inundation map ('im')\n# VTK generator ('vtk)\n Answer: ")

# Read geo to get nx and ny


geomatrix = np.loadtxt(dempath, dtype='i', delimiter=' ', skiprows=6)
ny = geomatrix.shape[0]
nx = geomatrix.shape[1]


if (simType == 'cs'):

    simType = input("Options:\n# Examine Flow (f)\n# Water Quality (wq)\n# Soil Quality (sq)\n Answer: ")

    #yearselect =  int(input("Simulation year (STC): "))

    if (simType == 'f'):
        CrossSecLine_path_shapefile = '/media/dcosta/data/megasync/ec_main/models/fluxos/support/STC_data_pre-processing/00_Cross_Sections/MS9.shp'
        #if yearselect==2009:
        #    obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/1_Compiled_for_FLUXOS_validation/1_2009_Compiled/Streamflow_MS9C_2009_trimmed_for_simulation.csv'
        #elif yearselect==2010:
        #    obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/1_Compiled_for_FLUXOS_validation/2_2010_Compiled/Streamflow_MS9C_2010_trimmed_to_simulation.csv'
        #elif yearselect==2011:
        #    obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/0_Obs_used_in_first_2011_tests/Snowmelt_Runoff_MS9_2011_justflow.csv'  # this is the one
        var_col_1 = 6  # 3-h, 6-qx, 10-C, 11 - soil mass, 12 - fn_1, 13 - fe_1
        var_col_2 = 7  # 3-h, 6-qx, 10-C, 11 - soil mass, 12 - fn_1, 13 - fe_1
    elif (simType == 'wq'):
        CrossSecLine_path_shapefile = '/media/dcosta/data/megasync/ec_main/models/fluxos/support/STC_data_pre-processing/00_Cross_Sections/MS_lake.shp'
        #if yearselect == 2009:
        #   obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/1_Compiled_for_FLUXOS_validation/2_2009_Compiled/Streamflow_WQ_MS12_2009.csv'
        #elif yearselect == 2010:
        #    obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/1_Compiled_for_FLUXOS_validation/2_2010_Compiled/Streamflow_WQ_MS12_2010.csv'
        #elif yearselect == 2011:
        #    obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/1_Compiled_for_FLUXOS_validation/3_2011_Compiled/Streamflow_WQ_MS12_2011_NO3event.csv'
        #    # obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/1_Compiled_for_FLUXOS_validation/3_2011_Compiled/Streamflow_WQ_MS9C_2011.csv'
        var_col_1 = 10  # 3-h, 6-qx, 10-C, 11-soil
        var_col_2 = 3  # 7-
    elif (simType == 'sq'):
        CrossSecLine_path_shapefile = '/media/dcosta/data/megasync/ec_main/models/fluxos/support/STC_data_pre-processing/00_Cross_Sections/MS_lake.shp'
        #obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/1_Compiled_for_FLUXOS_validation/3_2011_Compiled/Streamflow_WQ_MS12_2011_NO3event.csv'
        ## obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/1_Compiled_for_FLUXOS_validation/3_2011_Compiled/Streamflow_WQ_MS9C_2011.csv'
        var_col_1 = 11  # 3-h, 6-qx, 10-C, 11-soil
        var_col_2 = 0  # 7-
    else:
        sys.exit("Error: Not a valid entry (only accepts 'f', 'wq' or 'sq'")


    '''
    Call of functions (do not change)
    '''

    # START....

    Overwrite_csfile = input("Overwrite *.out files if existent?\n# Yes (y)\n# No (n), append new time steps to existing *.out file \n Answer: ")

    # Extract Cross-Section points from the shapefile
    shp_nodes_CS = cse.lineCSshapefile(CrossSecLine_path_shapefile)  # get nodes of the shapefile with the desired cross-section
    xy_CS = cse.InterpCSpoints_from_y_CS(shp_nodes_CS)  # interpolate a continuous line between the end nodes of "shp_nodes_CS"

    for sim in range(0, len(resultdir_list)):
        # Extract simulation name
        resultdir = resultdir_list[sim]
        try:
            resfiles_list = [x for x in os.listdir(resultdir) if x not in 'vtk' if x not in 'cs'] # listing the files in folder
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
            #simname = dm.getsimname(resultdir,obsPath,simType)
            #simname = dm.getsimname_2(resultdir,simType)

        if (len(resfiles_list) == 0):
            continue

        # Extract the results over the cross-section
        crosecval = cse.csextract(simType,resultdir,resfiles_list, xy_CS, sim, var_col_1, var_col_2, nx, ny, dxy)
        #crosecval_time = resfiles_list
        #crosecval_time = [s.strip('.txt') for s in crosecval_time]
        #crosecval_time = np.array(list(map(float, crosecval_time)))
        #crosecval = np.transpose(np.vstack((crosecval_time, np.transpose(crosecval))))

        # stack existig and new crosecval
        if file_cs_exists:
            crosecval = np.vstack((crosecval_exist, crosecval))

        crosecval = crosecval[np.argsort(crosecval[:, 0])][:] # sort results (parallel returns it unsorted)

        dm.savereslt(output_dir,crosecval)

    # Extract Observations
    #time_col = 0  # time column
    #val_col = 1  # obs column
    #obsval = dm.obsextract(obsPath, time_col, val_col)

    # Plot cross-section data
    #grph.plotCSvals(crosecval, obsval, simname,resultdir_list,simType,runlag)


elif (simType == 'im'):

    for sim in range(0, len(resultdir_list)):

        simname = dm.getsimname(resultdir_list[sim], '', simType)

        simType = input("Options:\n# kml for Google maps ('gm')\n# Scatter Plot 3D ('sp')\n# Quiver flowpaths ('qf')")

        # Google Maps
        if (simType == 'gm'):
            var_col_1 = 3
            geklm.google_eart_animation(resultdir_list[sim], simname, var_col_1, TimeStrgStart, Tinitial, nx, ny,
                                        dxy, coords, rotation, resolImage, var_1_graphymax, mapoverlay_opaqueness,
                                        trimrow,trimcol)

        # plot 3D dem and sim
        elif (simType == 'sp'):
            grph.scatter3d_pltly(simname,resultdir_list[sim],simnum,dempath,nx,ny) #  Good
            #scatter3d_matplotlib() # probablyt not working anymore

        # Quiver
        elif (simType == 'qf'):
            grph.quivergen(simname,resultdir_list[sim],simnum,dempath,nx,ny,dxy)

        else:
            sys.exit("Error: Not a valid entry")

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