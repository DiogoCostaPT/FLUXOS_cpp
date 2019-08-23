
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

'''
General Model Settings
'''

dempath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/DEM_ASCII/model_geo.csv'

#resultdir = '/media/dcosta/DATADRIVE1/fluxos_tests/local/test_wintra_11/Results/'
#resultdir = '/media/dcosta/DATADRIVE1/fluxos_tests/local/test_wintra_14/Results/'
resultdir = '/media/dcosta/DATADRIVE1/fluxos_tests/local/test_wintra_15_re-test_of_test_5_with_a_commit_earlier_-Seg_fault_with_crhono../Results/'
#resultdir = '/media/dcosta/DATADRIVE1/fluxos_tests/local/test_wintra_16/Results/'
resultdir = '/media/dcosta/DATADRIVE1/fluxos_tests/local/test_wintra_14/Results/'
resultdir = '/media/dcosta/DATADRIVE1/fluxos_tests/local/test_wintra_18/Results/'
resultdir = '/media/dcosta/DATADRIVE1/fluxos_tests/local/test_wintra_19/'
#resultdir = '/media/dcosta/DATADRIVE1/fluxos_tests/local/test_wintra_20/Results/'
resultdir = '/media/dcosta/DATADRIVE1/fluxos_tests/local/test_wintra_23_before_swapping_commit_2bf78954f3db3b07a6f3480fb515c987ee546b04/'
#resultdir = '/media/dcosta/DATADRIVE1/fluxos_tests/local/test_wintra_22/Results/'
resultdir = '/media/dcosta/DATADRIVE1/fluxos_tests/local/test_wintra_11/Results/'
resultdir = '/media/dcosta/DATADRIVE1/fluxos_tests/local/test_wintra_28/Results/'
resultdir = '/media/dcosta/DATADRIVE1/fluxos_tests/local/test_wintra_31/Results/'
resultdir = '/media/dcosta/DATADRIVE1/fluxos_tests/local/test_wintra_32/Results/'
#resultdir = '/media/dcosta/DATADRIVE1/fluxos_tests/local/test_wintra_33/Results/'
resultdir = '/media/dcosta/DATADRIVE1/fluxos_tests/local/test_wintra_34/Results/'

resultdir = '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/test_wintra_35/Results/'  # Ks = 0.01
resultdir = '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/test_wintra_36/Results/'  # Ks = 0.005
#resultdir = '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/test_wintra_37/Results/'  # Ks = 0.015
#resultdir = '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/test_wintra_38/Results/'  # Ks = 0.005 - wintra par changed
#resultdir = '/media/dcosta/DATADRIVE1/fluxos_tests/local/STC/test_wintra_39/Results/'  # Ks = 0.002

TimeStrgStart = datetime(2011, 3, 31, 0, 0, 0)
Tinitial = 0
Timee = 1468800
t_step_read = 3600
simnum = 972000  # simulation to plot map
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

    if (simType == 'f'):
        CrossSecLine_path_shapefile = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/00_Cross_Sections/MS9.shp'
        obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/0_Obs_used_in_first_2011_tests/Snowmelt_Runoff_MS9_2011_justflow.csv'  # this is the one
        var_col_1 = 6  # 3-h, 6-qx, 10-C, 11 - soil mass, 12 - fn_1, 13 - fe_1
        var_col_2 = 7  # 3-h, 6-qx, 10-C, 11 - soil mass, 12 - fn_1, 13 - fe_1
    elif (simType == 'wq'):
        CrossSecLine_path_shapefile = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/00_Cross_Sections/MS_lake.shp'
        obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/1_Compiled_for_FLUXOS_validation/3_2011_Compiled/Streamflow_WQ_MS12_2011_NO3event.csv'
        # obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/1_Compiled_for_FLUXOS_validation/3_2011_Compiled/Streamflow_WQ_MS9C_2011.csv'
        var_col_1 = 10  # 3-h, 6-qx, 10-C, 11-soil
        var_col_2 = 3  # 7-
    elif (simType == 'sq'):
        CrossSecLine_path_shapefile = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/00_Cross_Sections/MS_lake.shp'
        obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/1_Compiled_for_FLUXOS_validation/3_2011_Compiled/Streamflow_WQ_MS12_2011_NO3event.csv'
        # obsPath = '/media/dcosta/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/0_Obs/1_Compiled_for_FLUXOS_validation/3_2011_Compiled/Streamflow_WQ_MS9C_2011.csv'
        var_col_1 = 11  # 3-h, 6-qx, 10-C, 11-soil
        var_col_2 = 0  # 7-
    else:
        sys.exit("Error: Not a valid entry (only accepts 'f', 'wq' or 'sq'")


    '''
    Call of functions (do not change)
    '''

    # Extract simulation name
    simname = dm.getsimname(resultdir,obsPath,simType)

    # analyse the cross-section results
    # Extract Cross-Section points from the shapefile
    shp_nodes_CS = cse.lineCSshapefile(CrossSecLine_path_shapefile)  # get nodes of the shapefile with the desired cross-section
    xy_CS = cse.InterpCSpoints_from_y_CS(shp_nodes_CS)  # interpolate a continuous line between the end nodes of "shp_nodes_CS"

    # Extract the results over the cross-section
    crosecval = cse.csextract(simType,resultdir, xy_CS, Tinitial, Timee, t_step_read, var_col_1, var_col_2, nx, ny, dxy)

    # Extract Observations
    time_col = 0  # time column
    val_col = 1  # obs column
    obsval = dm.obsextract(obsPath, time_col, val_col)

    # Plot cross-section data
    grph.plotCSvals(crosecval, obsval, simname,resultdir,simType,runlag)


elif (simType == 'im'):

    simname = dm.getsimname(resultdir, '', simType)

    simType = input("Options:\n# kml for Google maps ('gm')\n# Scatter Plot 3D ('sp')\n# Quiver flowpaths ('qf')")

    # Google Maps
    if (simType == 'gm'):
        var_col_1 = 3
        geklm.google_eart_animation(resultdir, simname, var_col_1, TimeStrgStart, Tinitial, Timee, t_step_read, nx, ny,
                                    dxy, coords, resolImage, var_1_graphymax, mapoverlay_opaqueness)

    # plot 3D dem and sim
    elif (simType == 'sp'):
        grph.scatter3d_pltly(simname,resultdir,simnum,dempath,nx,ny) #  Good
        #scatter3d_matplotlib() # probablyt not working anymore

    # Quiver
    elif (simType == 'qf'):
        grph.quivergen(simname,resultdir,simnum,dempath,nx,ny,dxy)

    else:
        sys.exit("Error: Not a valid entry")

elif (simType == 'vtk'):
    simname = dm.getsimname(resultdir, '', simType)
    vtkgen.vtk_generator(simname, resultdir, simnum, dempath, nx, ny, dxy)

else:
    sys.exit("Error: Not a valid entry (only accepts 'cs', 'im' or 'vtk'")