import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# read data from file
data = np.loadtxt('results_graph_ex3.csv', delimiter=',')
# Order the data following the first column without destroy the vector so that it can be plotted as a heatmap 
#data = data[data[:,2].argsort()] # First sort doesn't need to be stable.
#data= data[data[:,1].argsort(kind='mergesort')]
#data= data[data[:,0].argsort(kind='mergesort')]

# extract x, y, and temperature values
x = data[:, 0]
y = data[:, 1]
temp = data[:, 2]

# define grid
xi = np.linspace(min(x), max(x), 100)
yi = np.linspace(min(y), max(y), 100)
xi, yi = np.meshgrid(xi, yi)

fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(15, 5))
ax1.set_title('Seaborn kdeplot')
sns.kdeplot(x=x, y=y, cmap='Blues', fill=True, thresh=0.02, ax=ax1)

ax2.set_title('Scatter plot with high transparency')
ax2.scatter(x, y, color='blue', alpha=0.02)

plt.tight_layout()
plt.show()