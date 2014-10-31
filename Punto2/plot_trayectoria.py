from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import sys

filename = sys.argv[1]

data = np.loadtxt(filename);

fig = plt.figure()
ax=plt.axes()
ax.plot(data[:,1],data[:,2], color ='r')
ax.set_xlabel('$x[R_t]$')
ax.set_xlim(-2.2,2.2)
ax.set_ylim(-2.2,2.2)
ax.set_ylabel('$y[R_t]$')
plt.savefig(filename[:-4] +'.pdf', format = 'pdf', transparent = True)
plt.close

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
theta = np.linspace(0,2*np.pi,100);
phi = np.linspace(0,np.pi, 100);
x = np.outer(np.cos(theta), np.sin(phi))
y = np.outer(np.sin(theta), np.sin(phi))
z = np.outer(np.ones(100), np.cos(phi))
ax.plot_surface(x,y,z, rstride=3,cstride=3, color ='b')
ax.plot(data[:,1], data[:,2], data[:,3], color = 'r')
ax.set_xlabel('$x[R_t]$')
ax.set_ylabel('$y[R_t]$')
ax.set_zlabel('$z[R_t]$')
ax.set_xlim(-2.2,2.2)
ax.set_ylim(-2.2,2.2)
ax.set_zlim(-2.2,2.2)
plt.savefig(filename[:-4] +'3d.pdf', format = 'pdf', transparent = True)
plt.close

