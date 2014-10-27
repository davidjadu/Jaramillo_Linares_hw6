from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import sys

filename = sys.argv[1]

data = np.loadtxt(filename);

fig = plt.figure()
ax=plt.axes()
ax.plot(data[:,1],data[:,2])
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.savefig(filename[:-4] +'.pdf', format = 'pdf', transparent = True)
plt.close


fig = plt.figure()
bx = fig.add_subplot(111, projection = '3d')
bx.plot_wireframe(data[:,1],data[:,2], data[:,3])
bx.set_xlabel('x')
bx.set_ylabel('y')
bx.set_zlabel('z')
plt.savefig(filename[:-4] +'3d.pdf', format = 'pdf', transparent = True)
plt.close

