import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# read data from file
data = np.loadtxt('C:\\tmp\\FVMc\\results_graph_ex3.csv', delimiter=',')
# Order the data following the first column without destroy the vector so that it can be plotted as a heatmap 
data = data[data[:,2].argsort()] # First sort doesn't need to be stable.
data= data[data[:,1].argsort(kind='mergesort')]
data= data[data[:,0].argsort(kind='mergesort')]

# extract x, y, and temperature values
x = data[:, 0]
y = data[:, 1]
temp = data[:, 2]

# define grid
xi = np.linspace(min(x), max(x), 30)
yi = np.linspace(min(y), max(y), 30)
xi, yi = np.meshgrid(xi, yi)

# interpolate temperature values onto grid
zi = griddata((x, y), temp, (xi, yi), method='linear') # method='

# plot the heatmap
fig, ax = plt.subplots()
im = ax.imshow(zi, cmap='hot', 
               interpolation='none',
               extent =[x.min(), x.max(), y.min(), y.max()])
ax.invert_yaxis()  # invert y-axis
fig.colorbar(im)

# set x and y axis units to meters
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')

# rescale x and y axis limits to physical values in meters
x_ticks = np.linspace(min(x), max(x), 5)
y_ticks = np.linspace(min(y), max(y), 5)
ax.set_xticks(x_ticks)
ax.set_yticks(y_ticks)
ax.set_xticklabels(['{:.2f}'.format(x_tick) for x_tick in x_ticks], rotation=45)
ax.set_yticklabels(['{:.2f}'.format(y_tick) for y_tick in y_ticks], rotation=45)

# adjust tick spacing
fig.subplots_adjust(bottom=0.15)

plt.show()