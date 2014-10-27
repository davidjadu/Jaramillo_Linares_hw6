from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import sys

filename = sys.argv[1]

data = np.loadtxt(filename);

fig = plt.figure()
ax=plt.axes()
ax.plot(data[:,1],data[:,2])
ax.set_xlabel('$x[m]$')
ax.set_ylabel('$y[m]$')
plt.savefig(filename[:-4] +'.pdf', format = 'pdf', transparent = True)
plt.close


fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
ax.plot(data[:,1], data[:,2], data[:,3])
ax.set_xlabel('$x[m]$')
ax.set_ylabel('$y[m]$')
ax.set_zlabel('$z[m]$')
plt.savefig(filename[:-4] +'3d.pdf', format = 'pdf', transparent = True)
plt.show()
plt.close

