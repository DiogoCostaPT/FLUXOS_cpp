

import plotly as py
import plotly.graph_objs as go
#plotly.tools.set_credentials_file(username='dcosta', api_key='m1o7WiVD2ylaQTwQfnFS') # if I want to plot online
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# my libraries
import data_management as dm


# Plot CS extracted values
def plotCSvals(crosecval,obsval,simname,resultdir,simType,runlag):

    import plotly as py
    from plotly import tools
    import plotly.graph_objs as go
    import pandas as pd

    # Simulations
    time_CS = (crosecval[:, 0]  - runlag) / 3600 # sec to hour (for lag time between forcing and response)
    values_CS = crosecval[:, 1:len(crosecval.T)-1]
    values_CS_orig = values_CS

    if(simType=='wq'): # water levels and concentrations (average value) - remove zeros
        values_CS[values_CS==0] = np.nan
        ymodel = np.nanmean(values_CS, axis=1)
    elif(simType=='sq'):
        values_CS[values_CS == 0] = np.nan
        ymodel = np.nanmean(values_CS, axis=1)
    elif(simType=='f'): # flow
        ymodel = np.sum(values_CS, axis=1)


    # Observations
    obsval_adj = obsval
    obsval_adj_val = obsval_adj[:,1]
    #obsval_adj_val = obsval_adj[:,1] / 1000 / 3600 /24 * 210290 * 9
    time_temp = ((obsval_adj[:,0] - obsval_adj[0,0])) * 24
    #time_temp = (obsval_adj[:, 0] - obsval_adj[1, 0])


    # Create a trace
    trace_simfluxos = go.Scatter(
        x = time_CS,
        y = ymodel,
        name = 'Sim',
        line = dict(
                color='rgb(183, 32, 32)'
        )
    )

    trace_obs = go.Scatter(
        x = time_temp,
        y = obsval_adj_val,
        name = 'Obs',
        mode = 'markers',
        marker = dict(
            size=9,
            color='rgba(0,0,0,0)',
            line=dict(
                width=2,
                color='rgb(150, 150, 150)'
            )
        )
    )

    trace_surf = go.Heatmap(
        z=values_CS_orig.T,
        x=time_CS,
        y=np.linspace(1, len(values_CS_orig[1,:]), num = len(values_CS_orig[1,:])),
        opacity=0.5,
        colorbar=dict(
            title=simType,
            titleside='top',
            tickmode='array',
            ticktext=['Hot', 'Mild', 'Cool'],
            ticks='outside',
            lenmode = 'fraction',len=0.75
        )
    )

    fig = tools.make_subplots(rows=2, cols=1,subplot_titles=('Flow in each cell (m3/s)', 'Total flow (m3/s)'))
    fig.append_trace(trace_surf, 1, 1)
    fig.append_trace(trace_obs, 2, 1)
    fig.append_trace(trace_simfluxos, 2, 1)

    fig['layout']['xaxis1'].update(title='Time [hours]', showgrid=True)
    fig['layout']['xaxis2'].update(title='Time [hours]', showgrid=True)
    fig['layout']['yaxis1'].update(title='Cross-Section', showgrid=True)
    fig['layout']['yaxis2'].update(title='Flow [m3/s]', showgrid=True)

    titleplot = 'Simulation Results (Cross-Section): ' + resultdir
    figname = 'CS_Results/' + simname + '.html'
    fig['layout'].update(height=1000, width=1500, title=titleplot)
    py.offline.plot(fig, filename=figname)


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


# Scatter plot 3D
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


# Quiver for flowpaths
def quivergen(simname,resultdir,simnum,dempath,nx,ny,dxy):

    import plotly.figure_factory as ff

    X1, Y1 = np.meshgrid(np.arange(0, 2 * np.pi, .2), np.arange(0, 2 * np.pi, .2))
    U1 = np.cos(X1)
    V1 = np.sin(Y1)
    M1 = np.hypot(U1, V1)

    X, Y = np.meshgrid(np.arange(dxy/2, ny*dxy+dxy/2, dxy), np.arange(dxy/2, nx*dxy+dxy/2, dxy))
    simpath = resultdir + str(simnum) + '.txt'
    xyz_columndata_all = pd.read_csv(simpath)
    xyz_columndata_all = xyz_columndata_all.values
    xyz_columndata = dm.xyz_extract_z_column(xyz_columndata_all, 0, 1, 6, 7)  # extract relevant column
    #xyz_columndata_colour = dm.xyz_extract_z_column(xyz_columndata_all, 0, 1, 3, 0)
    U = dm.xyz_to_matrix(xyz_columndata[:,[0,1,2]], nx, ny)  # convert into matrix
    V = dm.xyz_to_matrix(xyz_columndata[:,[0,1,3]], nx, ny)  # convert into matrix
    #M = (U**2 + V**2)**(1/2)
    M = np.hypot(U, V)

    fig3, ax3 = plt.subplots()
    ax3.set_title("pivot='tip'; scales with x view")

    Q = ax3.quiver(X, Y, U, V, M,scale=1 / 0.1)
    qk = ax3.quiverkey(Q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E',
                       coordinates='figure')
    #ax3.scatter(X, Y, color='0.5', s=1)

    plt.show()
    plt.savefig(['Quiver_plots/' + simname + '.png'])


    #fig = ff.create_quiver(X, Y, U, V)
    #py.offline.plot(fig, filename='Quiver Plot Example')