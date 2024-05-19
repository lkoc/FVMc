# Create a python program that opens a csv file (x,y,T) and reads the data and plot it as a heat map
# hte name of the file is given as an argument in the command line
# the program should be able to read the data and plot it as a heat map
# the program should be able to save the plot as a png file
# the program should be able to print the plot on the screen

# import libraries
import numpy as np
import matplotlib.pyplot as plt
import sys

# read the file name from the command line
filename = "C:\\tmp\FVMc\\results_graph_ex3.csv" # sys.argv[1]

# read the data from the file
data = np.loadtxt(filename,delimiter=',')
x = data[:,0]
y = data[:,1]
T = data[:,2]

# plot the data as a heat map for x,y as coordinates and T as the color
# no scatter, heatap style
plt.figure()
plt.scatter(x,y,c=T,marker='s')
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Heat Map')
plt.savefig('C:\\tmp\FVMc\\results_graph_ex3.png')
plt.show()

# end of program


