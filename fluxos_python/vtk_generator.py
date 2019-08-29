
import vtktools
import pandas as pd
import data_management as dm

def vtk_generator(simname,resultdir,simnum,dempath,nx,ny,dxy):

    vtk_writer = vtktools.VTK_XML_Serial_Unstructured()

    simpath = resultdir + str(simnum) + '.txt'
    #simpath = 'model_geo.txt'
    xyz_columndata_all = pd.read_csv(simpath)
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

    outputnam = "vtk/" + simname + "_" + str(simnum)

    #vtk_writer.snapshot(outputnam + ".vtu", x, y, z, h, ux, uy, qx, qy, conc_sw, conc_soil)
    vtk_writer.snapshot(outputnam + ".vtu", x, y, z, h, ux, uy, qy, qx, conc_sw, conc_soil)
    vtk_writer.writePVD(outputnam + ".pvd")