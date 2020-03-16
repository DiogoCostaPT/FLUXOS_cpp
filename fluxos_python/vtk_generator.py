
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
    xyz_columndata = dm.xyz_extract_z_column(xyz_columndata_all, 0, 1, 2, 3)  # extract relevant column (x,y,z,h)

    x = xyz_columndata[:,0]
    y = xyz_columndata[:,1]
    z = xyz_columndata[:,2]
    h = xyz_columndata[:,3]

    xyz_columndata_2 = dm.xyz_extract_z_column(xyz_columndata_all, 4, 5, 12, 13)  # extract relevant column (ux,uy,qx,qy)
    ux = xyz_columndata_2[:, 2] * 1000 / (dxy*dxy)   # m3/s -> l/s
    uy = xyz_columndata_2[:, 3] * 1000 / (dxy*dxy)  # m3/s -> l/s
    qx = xyz_columndata_2[:, 2]  # m3/s
    qy = xyz_columndata_2[:, 3]  # m3/s

    xyz_columndata_3 = dm.xyz_extract_z_column(xyz_columndata_all, 8, 9, 10, 11)  # extract relevant column (ux,uy,qx,qy)
    umag = xyz_columndata_3[:, 0]
    us = xyz_columndata_3[:, 1]
    conc_sw = xyz_columndata_3[:, 2]
    conc_soil = xyz_columndata_3[:, 3]

    xyz_columndata_4 = dm.xyz_extract_z_column(xyz_columndata_all, 12, 13, 14, 14)  # extract relevant column (ux,uy,qx,qy)
    timeconnect_h = xyz_columndata_4[:, 3]

    outputnam = resultdir + "vtk/" + resfile_i[0:len(resfile_i)-4]

    #vtk_writer.snapshot(outputnam + ".vtu", x, y, z, h, ux, uy, qx, qy, conc_sw, conc_soil)
    vtk_writer.snapshot(outputnam + ".vtu", x, y, z, h, ux, uy, qy, qx, conc_sw, conc_soil,timeconnect_h)
    vtk_writer.writePVD(outputnam + ".pvd")