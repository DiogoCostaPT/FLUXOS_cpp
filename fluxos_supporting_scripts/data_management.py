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
import numpy as np
import re

# Extract the simulation name
def getsimname(resultdir,obsPath,simType):

    # namstart = resultdir.find("/")
    nameend = resultdir.find('Results')

    # Csectiom_nstart = obsPath.find('MS')
    # Csectiom_nend = obsPath.find('20')

    resultdir_i = resultdir[0:nameend-1]

    res = [i.start() for i in re.finditer("/", resultdir_i)]
    namstart = res[len(res) - 1]

    simname = resultdir_i[namstart+1:len(resultdir_i)]

    # locbrack = simname.find('/')
    # simname = simname[0:locbrack] + '_' + simname[locbrack+1:len(simname)] + '_' + simType

    return simname

def getsimname_2(resultdir,simType):

    namstart = resultdir.find('/t_')
    nameend = resultdir.find('Results/')

    simname = [resultdir[namstart+1:nameend - 1] + '_' + simType]

    return simname

#  generated xyz_data by extracting  x, y and z values from the selected columns
def xyz_extract_z_column(table_data, x_loc, y_loc, var1_loc, var2_loc):
    # var_col - collum that we want to extract (it assummes that x and y are in columns 1 and 2, respectively
    if var2_loc == 0:
        xyz_columndata_list = [table_data[:, x_loc], table_data[:, y_loc], table_data[:, var1_loc]]
    else:
        xyz_columndata_list = [table_data[:, x_loc], table_data[:, y_loc], table_data[:, var1_loc], table_data[:, var2_loc]]

    xyz_columndata = np.array(xyz_columndata_list).T
    return xyz_columndata


#  Convert xyz into matrix
def xyz_to_matrix(xyz_columndata, nx, ny):
    # xyzdata = xyz data
    # matrix dimensions: nx * ny
    # assumes that following order: x, y and z in columns 1, 2 and 3, respectively

    var = np.zeros((ny,nx))
    # loop to extract each point of the cross-section from the result file and put in the varall_1 matrix
    for row in range(1, len(xyz_columndata)):

        xi = xyz_columndata[row, 1] - 1
        yi = xyz_columndata[row, 0] - 1
        var_i = xyz_columndata[row, 2]
        var[int(yi), int(xi)] = var_i

    #var[var == 0] = 'nan'
    return var

# Extract relevant observation data
def obsextract(obsPath,time_col,val_col):

    with open(obsPath, 'r') as fid:  # open the result file x
        dataraw = np.genfromtxt(obsPath, delimiter=',')
        obsval = np.vstack((dataraw[:, time_col], dataraw[:, val_col]))

        obsval_noNAN = obsval[:, ~np.all(np.isnan(obsval), axis=0)]
        obsval_noNAN = obsval_noNAN.T # put in columns

    return obsval_noNAN

def savereslt(output_dir,crosecval):
    np.savetxt(output_dir, crosecval, delimiter=',')

