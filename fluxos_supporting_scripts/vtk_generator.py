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

import vtktools
import pandas as pd
import data_management as dm

def vtk_generator(simType,resultdir,resfile_i,dempath,nx,ny,dxy):

    #simname = dm.getsimname_2(resultdir, simType)

    vtk_writer = vtktools.VTK_XML_Serial_Unstructured()

    simpath = resultdir + resfile_i
    #simpath = 'model_geo.txt'
    xyz_columndata_all = pd.read_csv(simpath,error_bad_lines=False)
    xyz_columndata_all = xyz_columndata_all.values
    xyz_columndata = dm.xyz_extract_z_column(xyz_columndata_all, 0, 1, 4, 5)  # extract relevant column (x,y,z,h)

    x = xyz_columndata[:,0]
    y = xyz_columndata[:,1]
    z = xyz_columndata[:,2]
    h = xyz_columndata[:,3]

    xyz_columndata_2 = dm.xyz_extract_z_column(xyz_columndata_all, 8, 9,11,12)  # extract relevant column (ux,uy,qx,qy)
    ux = xyz_columndata_2[:, 0] * 1000 / (dxy*dxy)   # m3/s -> l/s
    uy = xyz_columndata_2[:, 1] * 1000 / (dxy*dxy)  # m3/s -> l/s
    qx = xyz_columndata_2[:, 0]  # m3/s
    qy = xyz_columndata_2[:, 1]  # m3/s
    conc_sw = xyz_columndata_2[:, 2]
    conc_soil = xyz_columndata_2[:, 3]

    xyz_columndata_4 = dm.xyz_extract_z_column(xyz_columndata_all, 15, 0, 0, 0)  # extract relevant column (ux,uy,qx,qy)
    timeconnect_h = xyz_columndata_4[:, 0]

    outputnam = resultdir + "vtk/" + resfile_i[0:len(resfile_i)-4]

    #vtk_writer.snapshot(outputnam + ".vtu", x, y, z, h, ux, uy, qx, qy, conc_sw, conc_soil)
    vtk_writer.snapshot(outputnam + ".vtu", x, y, z, h, ux, uy, qy, qx, conc_sw, conc_soil,timeconnect_h)
    #vtk_writer.snapshot(outputnam + ".vtu", x, y, h, ux, uy, qy, qx, conc_soil)
    vtk_writer.writePVD(outputnam + ".pvd")