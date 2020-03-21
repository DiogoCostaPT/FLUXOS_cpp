
# create png image
def plotPNGplotly(googlefolder,simname,xyz_matrix_var,nx,ny,dxy,timei,resolImage,var_1_graphymax):

    # import plotly
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    import numpy as np
    from matplotlib.colors import BoundaryNorm
    from matplotlib.ticker import MaxNLocator

    # Make data.
    x = np.linspace(0, nx*dxy, nx)
    y = np.linspace(0, ny*dxy, ny)
    x, y = np.meshgrid(x, y)
    z = xyz_matrix_var

    z[z <= 0] = np.nan
    #z[z >= 6] = np.nan

    # x and y are bounds, so z should be the value *inside* those bounds.
    # Therefore, remove the last value from the z array.
    # z = z[:-1, :-1]
    levels = MaxNLocator(nbins=50).tick_values(0, var_1_graphymax)

    # pick the desired colormap, sensible levels, and define a normalization
    # instance which takes data values and translates those into levels.
    cmap = plt.get_cmap('jet')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    fig, (ax0) = plt.subplots(nrows=1)

    im = ax0.pcolormesh(x, y, z, cmap=cmap, norm=norm)
    clb = fig.colorbar(im, ax=ax0)
    clb.set_label('Water level [m]')

    plt.axis('off')

    fig.tight_layout()

    simpngname = googlefolder + '/' + simname + str(timei) + '.png'
    plt.savefig(simpngname, transparent=True, dpi=resolImage, bbox_inches='tight') #

    return simpngname


# PNG image creator for klm (Google Earth)
def pgncreator(resultdir, googlefolder,simname, timevec, t, var_col,nx,ny,dxy,resolImage,var_1_graphymax):
    import os
    import numpy as np
    import data_management as dm

    # generate the file name string
    timei = timevec[t]
    resfilepath = resultdir + str(timei) + ".txt"

    # open result file
    if os.path.exists(resfilepath):
        with open(resfilepath, 'r') as fid:  # open the result file x
        # read the result file x
          try:
            dataraw = np.genfromtxt(resfilepath, delimiter=',')

            xyz_columndata = dm.xyz_extract_z_column(dataraw, 0, 1, var_col,0)  # extract relevant column
            xyz_matrix_var = dm.xyz_to_matrix(xyz_columndata[:, [0, 1, 2]], ny,nx)  # convert into matrix (var 1)

            fid.close()

            simpngname = plotPNGplotly(googlefolder,simname,xyz_matrix_var, nx, ny, dxy, timei,resolImage,var_1_graphymax)

          except:
            print("Cannot open file or problem parsing it:" + resfilepath)
            simpngname = ''

    else:
        print("File does not exist: " + resfilepath)
        simpngname = ''

    return simpngname

#% Google Earth KLM generator
def google_eart_animation(resultdir,simname,var_col,TimeStrgStart,Tinitial,nx,ny,dxy,coords,resolImage,var_1_graphymax,mapoverlay_opaqueness):

    from joblib import Parallel, delayed
    import multiprocessing
    from tqdm import tqdm
    import os
    from datetime import timedelta, datetime
    import numpy as np
    from os import listdir
    from os.path import isfile, join

    googlefolder = resultdir + 'GOOGLE_EARTH'
    if not os.path.exists(googlefolder):
        os.makedirs(googlefolder)

    resultfiles = [f for f in listdir(resultdir) if isfile(join(resultdir, f))]
    resultfiles_notxt = list(map(lambda x: x.replace('.txt',''),resultfiles))
    for i in range(0, len(resultfiles_notxt)):
        resultfiles_notxt[i] = int(resultfiles_notxt[i])

    resultfiles_notxt.sort()

    num_cores = multiprocessing.cpu_count()

    pgn_filepaths = np.vstack(Parallel(n_jobs=num_cores)(
        delayed(pgncreator)(resultdir, googlefolder, simname, resultfiles_notxt, t, var_col,nx,ny,dxy,resolImage,var_1_graphymax) for t in
        tqdm(range(0, len(resultfiles_notxt)))))

    cwd = os.getcwd()

    # write klm
    kmlfullname = googlefolder + '/' + simname + '.kml'
    with open(kmlfullname, 'w') as fid:
        fid.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fid.write('<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">\n')
        fid.write('<name>' + kmlfullname + '</name>\n')
        fid.write('<Document id="1">\n')
        fid.write('<open>1</open>\n')
        fid.write('<Style>\n')
        fid.write('<ListStyle >\n')
        fid.write('<listItemType>check</listItemType>\n')
        fid.write('<bgColor>00ffffff</bgColor >\n')
        fid.write('<maxSnippetLines>2</maxSnippetLines>\n')
        fid.write('</ListStyle>\n')
        fid.write('</Style>\n')

        for row in range(1, len(pgn_filepaths)):
            gn_filepath_i = pgn_filepaths[row]
            gn_filepath_i = gn_filepath_i[0:1][0]

            # Extract name simulations
            str1 = gn_filepath_i[::-1]
            slashfind = str1.find("/")
            simname = str1[0:slashfind]
            simname = simname[::-1]

            fid.write('<GroundOverlay >\n')
            fid.write('<name>'+ simname +'</name>\n')
            fid.write('<visibility>0</visibility>\n')
            fid.write('<TimeSpan >\n')

            # Get time as string
            t_step_read = resultfiles_notxt[1] - resultfiles_notxt[0]
            timediff = row * t_step_read
            timestr_start = TimeStrgStart + timedelta(hours=timediff)
            timestr_end = timestr_start + timedelta(hours=1)

            timestr_start_GEform = datetime.strftime(timestr_start,'%Y-%m-%dT%XZ') # Google maps format: YYYY-MM-DDThh:mm:ssZ
            timestr_end_GEform = datetime.strftime(timestr_end, '%Y-%m-%dT%XZ')

            fid.write('<begin>' + timestr_start_GEform + '</begin>\n')
            fid.write('<end>' + timestr_end_GEform + '</end>\n')
            fid.write('</TimeSpan>\n')
            fid.write('<color>' + str(mapoverlay_opaqueness) + 'ffffff</color>\n')
            fid.write('<Icon>\n')
            fid.write('<href>' + simname + '</href>\n')
            fid.write('<viewBoundScale>0.75</viewBoundScale>\n')
            fid.write('</Icon>\n')
            fid.write('<LatLonBox>\n')
            fid.write('<north>' + str(coords[0]) + '</north>\n')
            fid.write('<south>' + str(coords[1]) + '</south>\n')
            fid.write('<east>' + str(coords[2]) + '</east>\n')
            fid.write('<west>' + str(coords[3]) + '</west>\n')
            fid.write('</LatLonBox>\n')
            fid.write('</GroundOverlay>\n')

        fid.write('</Document>\n')
        fid.write('</kml>\n')

    fid.close()