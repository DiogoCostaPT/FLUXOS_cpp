

import plotly as py
import plotly.graph_objs as go
#plotly.tools.set_credentials_file(username='dcosta', api_key='m1o7WiVD2ylaQTwQfnFS') # if I want to plot online
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# my libraries
import data_management as dm

#  Plot scatter
def scatter3d_matplotlib():
    '''
    ==============
    3D scatterplot
    ==============

    Demonstration of a basic scatterplot in 3D.
    '''

    # Fixing random state for reproducibility
    np.random.seed(19680801)

    # path and extract the results
    resfilepath = "/media/DATADRIVE1/fluxos_tests/plato/stc_54_pl/Results/SW/200_SW.csv"
    reader = np.genfromtxt(resfilepath, delimiter=',')

    # For each set of style and range settings, plot n random points in the box
    # defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
    xs = reader[:, 0]
    ys = reader[:, 1]
    zs = reader[:, 3]

    colors = ["#049FBB", "#5DF10A", "#BF5B3F", "#A95058", "#9C6BE1", "#6CC43C", "#7B3001", "#B3D407",
              "#97CE66"]  # hex colours

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xs, ys, c=[colors[i] for i in zs])

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()


def scatter3d_pltly(simname,resultdir,simnum,dempath,nx,ny):


    # DEM: Read data from a csv
    xyz_columndata = pd.read_csv(dempath)
    xyz_columndata = xyz_columndata.values
    xyz_matrix_dem = dm.xyz_to_matrix(xyz_columndata, nx, ny)  # convert into matrix
    xyz_matrix_dem = pd.DataFrame.from_records(xyz_matrix_dem)

    simpath = resultdir + str(simnum) + '_SW.csv'
    xyz_columndata_all = pd.read_csv(simpath)
    xyz_columndata_all = xyz_columndata_all.values
    xyz_columndata = dm.xyz_extract_z_column(xyz_columndata_all, 0, 1, 2, 0)  # extract relevant column
    xyz_columndata_colour = dm.xyz_extract_z_column(xyz_columndata_all, 0, 1, 3, 0)

    xyz_matrix_var1 = dm.xyz_to_matrix(xyz_columndata, nx, ny)  # convert into matrix
    xyz_matrix_var1[xyz_matrix_var1 < 10] = 0
    xyz_matrix_var1 = pd.DataFrame.from_records(xyz_matrix_var1)

    data = [
        go.Surface(z=xyz_matrix_dem.as_matrix(), showscale=False, colorscale='Greys'),
        go.Surface(z=xyz_matrix_var1.as_matrix(), showscale=True, surfacecolor=xyz_columndata_colour, opacity=0.65, cmin=0, cmax=1.5),
    ]

    layout1 = go.Layout(
        title='Simulation Results (South Tabacco Creek)',
        autosize=False,
        scene=dict(
            xaxis=dict(
                nticks=4, range=[0, nx], ),
            yaxis=dict(
                nticks=4, range=[0, ny], ),
            zaxis=dict(
                nticks=4, range=[0, max(xyz_matrix_dem)], ), ),
        width=1500,
        height=1000,
        margin=dict(
            l=150,
            r=150,
            b=65,
            t=90
        )

    )

    figname = 'DEM_sim_3D_Results/' + simname + '_CS_Map.html'
    #data = [trace1, trace2]
    layout = layout1
    fig = go.Figure(data=data, layout=layout)
    py.offline.plot(fig, filename=figname)
