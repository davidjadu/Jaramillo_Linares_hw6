import matplotlib.pyplot as plt
import numpy as np
import sys

filename = sys.argv[1]

data = np.loadtxt(filename);

fig = plt.figure()
ax=plt.axes()

ax.plot(data[:,1],data[:,2])
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
plt.savefig(filename[:-4] +'.pdf', format = 'pdf', transparent = True)
plt.close
