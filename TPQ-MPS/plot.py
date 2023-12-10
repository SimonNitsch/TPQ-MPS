import sys
import os 
import numpy
import matplotlib.pyplot as plt


arguments = sys.argv
path = arguments[1]
x_max = arguments[2]

vector = np.fromfile(path, dtype=np.float64)
length = vector.shape[0] / 2
vector = vector.reshape((2,length))

x = np.linspace(0,x_max,length)


plt.figure()
plt.plot(x,vector[:,0],"r")
plt.plot(x,vector[:,0]+vector[:,1],"r-.")
plt.plot(x,vector[:,0]-vector[:,1],"r-.")
plt.show()



