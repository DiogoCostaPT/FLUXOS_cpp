
# import libraries
import numpy as np

# Extract the simulation name
def getsimname(resultdir,obsPath,simType):

    namstart = resultdir.find('local/')
    nameend = resultdir.find('Results/')

    Csectiom_nstart = obsPath.find('MS')
    Csectiom_nend = obsPath.find('20')

    simname = resultdir[namstart:nameend - 1] + '_' + obsPath[Csectiom_nstart:Csectiom_nstart+4]

    locbrack = simname.find('/')
    simname = simname[0:locbrack] + '_' + simname[locbrack+1:len(simname)] + '_' + simType

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

    var = np.zeros((nx, ny))
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