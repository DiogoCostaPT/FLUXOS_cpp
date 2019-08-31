
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
from joblib import Parallel, delayed
import multiprocessing
from tqdm import tqdm

'''
General Model Settings
'''

dempath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/DEM_ASCII/model_geo.csv'

resultdir_list = [
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_36/Results/',
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_40/Results/',
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_41/Results/',
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_43/Results/',
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_44/Results/',
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_45/Results/',
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_46/Results/',
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_47/Results/',
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_48/Results/',
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_49/Results/',
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_50/Results/',
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_51/Results/',
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_52/Results/',
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_53/Results/',
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_54/Results/',
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_55/Results/',
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_56/Results/',
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_57/Results/',
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_58/Results/',
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_59/Results/',
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/t_60/Results/',
                '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/l_1/Results/',
]

TimeStrgStart = datetime(2011, 3, 31, 0, 0, 0)
Tinitial = 0
Timee = 18000 #1468800
t_step_read = 3600
nx = 722
ny = 1034
dxy = 3
runlag = 32400 # lag factor to account for the fact that we the model is being forced using streamflow time series

# Coordinates for STC (3m resolution) - GOOGLE EARTH
coords = [49.339205,  # N
          49.317788,  # S
          -98.342849,  # E
          -98.402657]  # W

resolImage = 200  # resolution of the images in Google Earth (in dpi)
var_1_graphymax = 5000
mapoverlay_opaqueness = 100

'''
Choose the functions to run
'''
# Sim results analysis
Map3d_Dem_and_sim_plotly = False # plot 3d maps of dem and sim results in plotly (interactive)
GoogleEarthSim = False # export kml to google earth
Quiver_flowpaths = True

simType = input("Options:\n# Analyse a cross-section (type 'cs')\n# Inundation map ('im')\n# VTK generator ('vtk)\n")

if (simType == 'cs'):

    simType = input("Options:\n# Examine Flow (f)\n# Water Quality (wq)\n# Soil Quality (sq)\n")

    #yearselect =  int(input("Simulation year (STC): "))

    if (simType == 'f'):
        CrossSecLine_path_shapefile = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/00_Cross_Sections/MS9.shp'
        #if yearselect==2009:
        #    obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/1_Compiled_for_FLUXOS_validation/1_2009_Compiled/Streamflow_MS9C_2009_trimmed_for_simulation.csv'
        #elif yearselect==2010:
        #    obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/1_Compiled_for_FLUXOS_validation/2_2010_Compiled/Streamflow_MS9C_2010_trimmed_to_simulation.csv'
        #elif yearselect==2011:
        #    obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/0_Obs_used_in_first_2011_tests/Snowmelt_Runoff_MS9_2011_justflow.csv'  # this is the one
        var_col_1 = 6  # 3-h, 6-qx, 10-C, 11 - soil mass, 12 - fn_1, 13 - fe_1
        var_col_2 = 7  # 3-h, 6-qx, 10-C, 11 - soil mass, 12 - fn_1, 13 - fe_1
    elif (simType == 'wq'):
        #CrossSecLine_path_shapefile = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/00_Cross_Sections/MS_lake.shp'
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
        #CrossSecLine_path_shapefile = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/00_Cross_Sections/MS_lake.shp'
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

    # Extract Cross-Section points from the shapefile
    shp_nodes_CS = cse.lineCSshapefile(CrossSecLine_path_shapefile)  # get nodes of the shapefile with the desired cross-section
    xy_CS = cse.InterpCSpoints_from_y_CS(shp_nodes_CS)  # interpolate a continuous line between the end nodes of "shp_nodes_CS"

    for sim in range(0, len(resultdir_list)):
        # Extract simulation name
        resultdir = resultdir_list[sim]
        resfiles_list = os.listdir(resultdir)
        #simname = dm.getsimname(resultdir,obsPath,simType)
        simname = dm.getsimname_2(resultdir,simType)

        # Extract the results over the cross-section
        crosecval = cse.csextract(simType,resultdir,resfiles_list, xy_CS, sim, var_col_1, var_col_2, nx, ny, dxy)

        dm.savereslt(simname,crosecval)

    # Extract Observations
    #time_col = 0  # time column
    #val_col = 1  # obs column
    #obsval = dm.obsextract(obsPath, time_col, val_col)

    # Plot cross-section data
    #grph.plotCSvals(crosecval, obsval, simname,resultdir_list,simType,runlag)


elif (simType == 'im'):

    simname = dm.getsimname(resultdir_list, '', simType)

    simType = input("Options:\n# kml for Google maps ('gm')\n# Scatter Plot 3D ('sp')\n# Quiver flowpaths ('qf')")

    # Google Maps
    if (simType == 'gm'):
        var_col_1 = 3
        geklm.google_eart_animation(resultdir_list, simname, var_col_1, TimeStrgStart, Tinitial, Timee, t_step_read, nx, ny,
                                    dxy, coords, resolImage, var_1_graphymax, mapoverlay_opaqueness)

    # plot 3D dem and sim
    elif (simType == 'sp'):
        grph.scatter3d_pltly(simname,resultdir_list,simnum,dempath,nx,ny) #  Good
        #scatter3d_matplotlib() # probablyt not working anymore

    # Quiver
    elif (simType == 'qf'):
        grph.quivergen(simname,resultdir_list,simnum,dempath,nx,ny,dxy)

    else:
        sys.exit("Error: Not a valid entry")

elif (simType == 'vtk'):
    print ("Running VTK generator...")
    for sim in range(0, len(resultdir_list)):
        # Extract simulation name
        resultdir = resultdir_list[sim]
        resfiles_list = os.listdir(resultdir)
        print("    " + resultdir)
        try:
            os.mkdir(resultdir + "vtk/")
            print("     vtk directory: CREATE")
        except:
            print("     vtk directory: ALREADY EXISTS")
        num_cores = multiprocessing.cpu_count()
        Parallel(n_jobs=num_cores)(delayed(vtkgen.vtk_generator)(simType, resultdir, resfiles_list[res_i], dempath, nx, ny, dxy) for res_i in tqdm(range(0,len(resfiles_list))))

else:
    sys.exit("Error: Not a valid entry (only accepts 'cs', 'im' or 'vtk'")