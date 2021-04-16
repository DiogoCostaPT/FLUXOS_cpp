# FLUXOS-OVERLAND (supporting scripts to examine the model outputs)
# Copyright (C) 2019-2021 ECCC
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

# import libraries
import geopandas as gpd
import pandas as pd
import os
from joblib import Parallel, delayed
import multiprocessing
from tqdm import tqdm
from scipy import interpolate
import numpy as np
# my libraries
import data_management as dm


# Extract the values for the cross-section (for all time steps) - parallel
def csextract(
        simType,
        resultdir,
        resfiles_list,
        xy_CS,
        XLL_YLL_CORNER_AND_CELL_SIZE_DEM,
        sim,
        nx,
        ny,
        dxy):

    # extracting the cross-section values
    #ntimstp = round((Timee - Tinitial) / t_step_read)+1  # number of time steps
    ntimstp = len(resfiles_list) # number of outputfiles
    #timevec = np.linspace(Tinitial, Timee, num=ntimstp)  # time of simulation
    #timevec = timevec.astype(int)

    # Get angle of the CS with the vertical and horizontal planes - change direction to guarantee always b1 < b2 (xy_CS_cor)
    geom_CS, xy_CS_cor = getanglesCS(xy_CS)

    num_cores = round(multiprocessing.cpu_count()*2/3)
    crosecvals = np.vstack(Parallel(n_jobs=num_cores)(
        delayed(Extract_File_Res)(
            simType,
            resultdir,
            resfiles_list,
            t_int,
            xy_CS_cor,
            XLL_YLL_CORNER_AND_CELL_SIZE_DEM,
            geom_CS,
            nx,
            ny,
            sim) for t_int in tqdm(range(0, ntimstp))))

    return crosecvals

# Extract the values for the cross-section (each time step)
def Extract_File_Res(
        simType,
        resultdir,
        resfilepath_all,
        t_int,xy_CS_cor,
        XLL_YLL_CORNER_AND_CELL_SIZE_DEM,
        geom_CS,
        nx,
        ny,
        sim):  # Loop over the result files

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

                try:
                    # read the result file x
                    dataraw = np.genfromtxt(resfilepath, delimiter=',', skip_header=1)
                    header = pd.read_csv(resfilepath, header=0, nrows=1).columns.tolist()

                    # Coordinates (from Results file)
                    # x-dir
                    coord_x_index = header.index('Xcoord')
                    coord_x = dm.xyz_extract_z_column(dataraw, 0, 1, coord_x_index, 0)  # extract relevant column
                    coord_x = np.round(coord_x)
                    #coord_x_matrix = dm.xyz_to_matrix(coord_x, nx, ny)  # convert into matrix
                    # y-dir
                    coord_y_index = header.index('Ycoord')
                    coord_y = dm.xyz_extract_z_column(dataraw, 0, 1, coord_y_index, 0)  # extract relevant column
                    coord_y = np.round(coord_y)
                    #coord_y_matrix = dm.xyz_to_matrix(coord_y, nx, ny)  # convert into matrix

                    # Vars
                    if simType == 'sq': # SQ case (soil concentrations)
                        # Var 1 only
                        var_col_1 = header.index('soil_mass [g]')
                        xyz_columndata = dm.xyz_extract_z_column(dataraw, 0, 1, var_col_1,0)  # extract relevant column
                        xyz_matrix_var_1 = dm.xyz_to_matrix(xyz_columndata, nx, ny)  # convert into matrix
                    else:
                        if simType == 'f':
                            var_col_1 = header.index('qx * dxy [m3/sec]')
                            var_col_2 = header.index('qy * dxy [m3/sec]')
                        elif (simType == 'wq'):
                            var_col_2 = header.index('conc_SW [mg/l]')
                            var_col_1 = header.index('h [m]')

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

                        # Get shapefile X and Y coordinates
                        xi_CCshapefile = xy_CS_cor[segi, 0].astype(float)
                        yi_CCshapefile = xy_CS_cor[segi, 1].astype(float)

                        xi_loc = np.where(coord_x[:,2] == xi_CCshapefile)
                        yi_loc = np.where(coord_y[:,2] == yi_CCshapefile)

                        xi_yi_row_loc = np.intersect1d(xi_loc,yi_loc)

                        # check if empty (if empty skip)
                        if xi_yi_row_loc:

                            xi = coord_y[xi_yi_row_loc,0]
                            yi = coord_y[xi_yi_row_loc,1]

                            # Reset to origin (0,0) to identify location in FLUXOS output files
                            #xi = round((xi - XLL_YLL_CORNER_AND_CELL_SIZE_DEM[0])/XLL_YLL_CORNER_AND_CELL_SIZE_DEM[2])
                            #yi = round((yi - XLL_YLL_CORNER_AND_CELL_SIZE_DEM[1])/XLL_YLL_CORNER_AND_CELL_SIZE_DEM[2])

                            # convert from float into int
                            xi = xi.astype(int)
                            yi = ny - yi.astype(int)

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
                                    crosecvals_t[segi + 1] = max(q_flow * np.cos(alpha_h), 0) + max(p_flow * np.sin(alpha_h), 0)
                        else:
                            # Element in cross-section is dry (could not find the element in Results files)
                            #print("Cross-Section element is dry (" + str(int(xi_CCshapefile)) + "," + str(
                            #    int(yi_CCshapefile)) + ") - SKIPPED")
                            crosecvals_t[segi + 1] = 0
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



