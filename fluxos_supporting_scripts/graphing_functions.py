

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
       # values_CS[values_CS < 0] = 0
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

