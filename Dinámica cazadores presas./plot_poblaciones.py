import matplotlib.pyplot as plt
import numpy as np
import sys

filename = sys.argv[1]

data = np.loadtxt(filename);

fig = plt.figure()
ax=plt.axes()

ax.plot(data[:,0],data[:,1], label="$x$")
ax.plot(data[:,0],data[:,2], label="$y$")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels)
ax.set_xlabel('$t$')
ax.set_ylabel('Population')
plt.savefig(filename[:-4] +'.pdf', format = 'pdf', transparent = True)
plt.close
