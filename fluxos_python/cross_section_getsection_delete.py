
# Extract values from Cross-Section

import numpy as np
import tkinter as tk
import re
import os
import csv
from tkinter import filedialog
root = tk.Tk()
root.withdraw()
root.update()
from io import StringIO
import plotly.plotly as py
import plotly.graph_objs as go

address_DEM = "E:/fluxos_tests/plato/stc_54_pl/data_input_STC/2dmb/model_geo.txt"
address_dain = "E:/fluxos_tests/plato/stc_54_pl/dain.txt"
filepathhResults = "E:/fluxos_tests/plato/stc_54_pl/Results/SW/200_SW.csv" # as reference to choose cross-section

# Get the print step and hdry
#address_dain = filedialog.askopenfilename(title = "Locate the 2dmb dain.txt") \ to ask the user for the file
with open(address_dain,'r') as fid:
    stopLine = 0
    while stopLine == 1:
        line = fid.readline().rstrip()
        print_step_flag = line.find('print_step')
        if print_step_flag != -1:
            print_step_cell = re.findall('\d.\d+', line)
            hdry_dain = float(print_step_cell[0])
            stopLine = stopLine + 1
fid.close()

# Get basin info: dx, dy, nx, ny
#address_DEM = filedialog.askopenfilename(title = "Locate 2dmb model_geo.txt") \ to ask the user for the file
with open(address_DEM,'r') as fid:
    line = fid.readline().rstrip()
    line = fid.readline().rstrip()
    linenum = re.findall('\d*\.?\d+', line)
    dxySW = linenum[0]
    nx = linenum[1]
    line = fid.readline().rstrip()
    linenum = re.findall('\d*\.?\d+', line)
    ny = linenum[1]
    line = fid.readline().rstrip()
    line = fid.readline().rstrip()
    reader = csv.reader(fid,delimiter=',') # much slower - but needed because there are field that aren't numeric
    DEMdata_raw = []
    for row in reader:
        try:
            dem_i = [float(x) for x in row[0:3]]
            DEMdata_raw.append(dem_i)
        except:
            pass # the last row is \\\\ so needs to be jumped
fid.close()

# Read Result to support cross-section delineation
with open(filepathhResults,'r') as fid:
    line = fid.readline().rstrip()
    reader = csv.reader(fid, quoting=csv.QUOTE_NONNUMERIC) # can only have numeric values - but it's much much faster
    simSW = []
    for row in reader:
        simSW.append(row[0:3])
    fid.close()


# Plotting
x = [row[0] for row in simSW]
y = [row[1] for row in simSW]
z = [row[2] for row in simSW]

x, y, z = np.random.multivariate_normal(np.array([0,0,0]), np.eye(3), 400).transpose()

trace1 = go.Scatter3d(
    x=x,
    y=y,
    z=z,
    mode='markers',
    marker=dict(
        size=12,
        color=z,                # set color to an array/list of desired values
        colorscale='Viridis',   # choose a colorscale
        opacity=0.8
    )
)

data = [trace1]
layout = go.Layout(
    margin=dict(
        l=0,
        r=0,
        b=0,
        t=0
    )
)
fig = go.Figure(data=data, layout=layout)
py.iplot(fig, filename='3d-scatter-colorscale')


