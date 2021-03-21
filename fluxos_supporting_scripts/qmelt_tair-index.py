
'''
Calculation of the Tair index for Qmelt Calculations
'''

import numpy as np
import scipy
from scipy import interpolate
from sklearn import datasets, linear_model
from plotly import tools
import plotly as py
import plotly.graph_objs as go
# my libs
import data_management as dm

AreaCatch = 2050000

Qmelt_file = '/media/DATADRIVE1/fluxos_tests/0_Obs/0_Obs_used_in_first_2011_tests/Snowmelt_Runoff_MS9_2011_justflow.csv' #  from observations of runoff
Tair_file = '/media/DATADRIVE1/MegaSync/FLUXOS/STC_data_pre-processing/Qmelt_estimation/Qmelt_Tair_relationship/Tair_STC_2011.csv'

# Read both files
time_col = 0
val_col = 1
Qmelt_obs = dm.obsextract(Qmelt_file, time_col, val_col)
Tair_obs = dm.obsextract(Tair_file, time_col, val_col)

Qmelt_obs[:,1] = Qmelt_obs[:,1]/ AreaCatch * 1000 * 3600 * 24  # conversion for FLUXOS: from m3/s -> mm/day

# interpolate Tair into the time step of Qmelt
f = interpolate.interp1d(Tair_obs[:, 0], Tair_obs[:, 1])

Tair_obs_interp = f(Qmelt_obs[:, 0])
Tair_obs_new = np.vstack((Qmelt_obs[:, 0],Tair_obs_interp)).T


'''Fits a linear fit of the form mx+b to the data'''
x_raw = Tair_obs_new[:, 1]
y_raw = Qmelt_obs[:, 1]
x = x_raw[x_raw>=0]
y = y_raw[x_raw>=0]

fitfunc = lambda params, x: params[0] * x    #create fitting function of form mx (goes through origin)
errfunc = lambda p, x, y: fitfunc(p, x) - y              #create error function for least squares fit

init_a = 0.5                            #find initial value for a (gradient)
init_p = np.array((init_a))  #bundle initial values in initial parameters

#calculate best fitting parameters (i.e. m and b) using the error function
p1, success = scipy.optimize.leastsq(errfunc, init_p.copy(), args = (x, y))
f = fitfunc(p1, x)          #create a fit with those parameters

Qmelt_Tair = x
Qmelt_predcurve = f
# Plotting

trace0 = go.Scatter( # Tair (original)
    x = Tair_obs[:,0],
    y = Tair_obs[:,1],
    mode = 'markers',
    name = 'Tair (original)'
)

trace1 = go.Scatter( # Qmelt (original)
    x = Qmelt_obs[:,0],
    y = Qmelt_obs[:,1],
    mode = 'markers',
    name = 'Qmelt (original)'
)


trace2 = go.Scatter( # Tair -> interpolated to Qmelt time step
    x = Tair_obs_new[:,0],
    y = Tair_obs_new[:,1],
    mode = 'lines+markers',
    name = 'Tair (interp to Qmelt timestep)'
)

trace3 = go.Scatter( # Tair Vs Qmelt
    x = Tair_obs_new[:, 1],
    y = Qmelt_obs[:, 1],
    mode = 'markers',
    name = 'Tair Vs Qmelt (obs)'
)

trace4 = go.Scatter(
    x = Qmelt_Tair,
    y = Qmelt_predcurve,
    mode = 'lines',
    name = 'Tair->Qmelt (model)'
)


fig = tools.make_subplots(rows=2, cols=1)

fig.append_trace(trace0, 1, 1)
fig.append_trace(trace1, 1, 1)
fig.append_trace(trace2, 1, 1)
fig.append_trace(trace3, 2, 1)
fig.append_trace(trace4, 2, 1)

fig['layout']['xaxis1'].update(title='Time')
fig['layout']['xaxis2'].update(title='Tair (degree celsius)')
fig['layout']['yaxis1'].update(title='Tair (degree celsius)')
fig['layout']['yaxis2'].update(title='Qmelt (mm/day)')

fig['layout'].update(height=1000, width=1500, title='i <3 annotations and subplots')

py.offline.plot(fig, filename='Qmelt-Tair_model.html')

