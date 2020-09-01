
# import libraries
import geopandas as gpd
import os
from joblib import Parallel, delayed
import multiprocessing
from tqdm import tqdm
from scipy import interpolate
import numpy as np
# my libraries
import data_management as dm


# Extract the values for the cross-section (for all time steps) - parallel
def csextract(simType,resultdir,resfiles_list, xy_CS,  sim, var_col_1, var_col_2, nx, ny, dxy):

    # extracting the cross-section values
    #ntimstp = round((Timee - Tinitial) / t_step_read)+1  # number of time steps
    ntimstp = len(resfiles_list) # number of outputfiles
    #timevec = np.linspace(Tinitial, Timee, num=ntimstp)  # time of simulation
    #timevec = timevec.astype(int)

    # Get angle of the CS with the vertical and horizontal planes - change direction to guarantee always b1 < b2 (xy_CS_cor)
    geom_CS, xy_CS_cor = getanglesCS(xy_CS)

    num_cores = round(multiprocessing.cpu_count()*2/3)
    crosecvals = np.vstack(Parallel(n_jobs=num_cores)(delayed(Extract_File_Res)(simType,resultdir,resfiles_list,t_int, xy_CS_cor, geom_CS, nx, ny, sim, var_col_1, var_col_2) for t_int in tqdm(range(0, ntimstp))))

    return crosecvals

# Extract the values for the cross-section (each time step)
def Extract_File_Res(simType,resultdir,resfilepath_all, t_int,xy_CS_cor, geom_CS, nx, ny, sim, var_col_1, var_col_2):  # Loop over the result files

    # generate the file name string
    crosecvals_t = np.zeros(((len(xy_CS_cor)) + 1))  # + 1 because the first column is time
    #timei = timevec[t]
    refile = resfilepath_all[t_int-1]
    resfilepath = resultdir + refile

    refile_time = refile[0:len(refile)-4]
    try:
        crosecval_time = int(refile_time)

        # open result file
        if os.path.exists(resfilepath):
            with open(resfilepath, 'r') as fid:  # open the result file x
                # read the result file x
                try:
                    dataraw = np.genfromtxt(resfilepath, delimiter=',')

                    if (simType == 'sq'): # SQ case (soil concentrations)
                        # Var 1 only
                        xyz_columndata = dm.xyz_extract_z_column(dataraw, 0, 1, var_col_1,0)  # extract relevant column
                        xyz_matrix_var_1 = dm.xyz_to_matrix(xyz_columndata, nx, ny)  # convert into matrix
                    else:
                        xyz_columndata = dm.xyz_extract_z_column(dataraw, 0, 1, var_col_1, var_col_2)  # extract relevant column
                        xyz_matrix_var_1 = dm.xyz_to_matrix(xyz_columndata[:, [0, 1, 2]], nx, ny)  # convert into matrix (var 1)
                        xyz_matrix_var_2 = dm.xyz_to_matrix(xyz_columndata[:, [0, 1, 3]], nx, ny)  # convert into matrix (var 2)

                    fid.close()

                    # Get geometry of the cross-section to calculate parallel and ortogonal flow (with respect to the cross-section)
                    a = geom_CS[0]
                    a = a.astype(int)
                    b = geom_CS[1]
                    b = b.astype(int)
                    alpha_h = geom_CS[2]

                    # extract the values of the cross section
                    #crosecvals_t[0] = timei

                    for segi in range(0, len(xy_CS_cor)):
                        xi = xy_CS_cor[segi, 0].astype(int)
                        yi = xy_CS_cor[segi, 1].astype(int)
                        if simType == 'sq':  # Water levels or concentrations
                            crosecvals_t[segi+1] = xyz_matrix_var_1[yi, xi]  # cross-section values
                        elif simType == 'wq':
                            crosecvals_t[segi + 1] = xyz_matrix_var_1[yi, xi]#/xyz_matrix_var_2[yi, xi]  # cross-section values
                        else:  # calculation of the ortogonal flow to the cross section
                            p_flow = xyz_matrix_var_1[yi, xi]
                            q_flow = xyz_matrix_var_2[yi, xi]
                            # different conditions (see notes)
                            if a == 0:  # perfectly horizontal cross-section
                                crosecvals_t[segi+1] = q_flow
                            elif b == 0:  # perfectly vertical cross-section
                                crosecvals_t[segi+1] = 0
                            else:  # general case
                             # crosecvals_t[segi + 1] = q_flow * np.cos(alpha_h) + p_flow * np.sin(alpha_h)
                             #crosecvals_t[segi + 1] = max(q_flow * np.cos(alpha_h) + p_flow * np.sin(alpha_h),0)
                                crosecvals_t[segi + 1] = max(q_flow * np.cos(alpha_h), 0) + max(p_flow * np.sin(alpha_h), 0)
                                #crosecvals_t[segi+1] = q_flow * np.cos(alpha_h) #+ abs(p_flow * np.sin(alpha_h))
                               # crosecvals_t[segi + 1] = abs(q_flow * np.cos(alpha_h))
                                #crosecvals_t[segi + 1] = abs(p_flow * np.sin(alpha_h))
                                #crosecvals_t[segi + 1] = q_flow + p_flow
                except:
                    print("Cannot open file:" + resfilepath + "or exception in function <Extract_File_Res>")

        else:
            print("File does not exist: " + resfilepath)
        crosecvals_t = np.hstack((crosecval_time, crosecvals_t))

    except:
        print("Invalid result file found: " + refile_time)
        crosecvals_t = np.hstack((0,crosecvals_t))

    return crosecvals_t

#  Get geometry of cross-section (for computation of perpendicular flow)
def getanglesCS(xy_CS):

    # Change direction of cross-section points to make sure that b1 < b2
    if xy_CS[0, 1] > xy_CS[len(xy_CS)-1, 1]:
        xy_CS_cor = np.flipud(xy_CS)
    else:
        xy_CS_cor = xy_CS

    a1 = xy_CS_cor[0,1]
    a2 = xy_CS_cor[len(xy_CS_cor)-1, 1]
    b1 = xy_CS_cor[0,0]
    b2 = xy_CS_cor[len(xy_CS_cor)-1, 0]

    a = abs(a1 - a2)
    b = abs(b1 - b2)
    # ratio = a / b  #  not needed

    alpha_h = np.arctan(a/b)

    geom_CS = np.vstack((a, b, alpha_h.astype(float)))

    return geom_CS, xy_CS_cor

#  Extract nodes from shapefile
def lineCSshapefile(polygons_shp_file):

    lineCS = gpd.read_file(polygons_shp_file)  # read shapefile containing the cross-section
    points_bounds = lineCS.unary_union.boundary.minimum_rotated_rectangle.xy  # extract minimum rectangle
    shp_nodes_CS = np.vstack((points_bounds[0], points_bounds[1]))

    return shp_nodes_CS

#  Get the xy of all points in the cross-section from the cross-section shapefile nodes (interpolation)
def InterpCSpoints_from_y_CS(shp_nodes_CS):



    x_start = int(round(shp_nodes_CS[0, 0]))
    x_end = int(round(shp_nodes_CS[0, 1]))
    y_start = int(round(shp_nodes_CS[1, 0]))
    y_end = int(round(shp_nodes_CS[1, 1]))

    # Fix the xy_CS_x (all values between the CS end nodes) - step is 1 (all grid cells)
    numcells = abs(x_start - x_end) + 1
    xy_CS_x = np.linspace(x_start, x_end, num = numcells)
    xy_CS_x = xy_CS_x.astype(int)

    # Interpolate the corresponding y values in the grid
    f = interpolate.interp1d([x_start, x_end], [y_start, y_end])
    xy_CS_y = f(xy_CS_x)
    xy_CS_y = xy_CS_y.astype(int)

    # The final CS array
    xy_CS = np.vstack((xy_CS_x,xy_CS_y)).T

    return xy_CS



