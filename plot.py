import numpy as np
import matplotlib.pyplot as plt

file = "KHeE.txt"
gsefile = "KHeG.txt"
N = 48

vector = np.loadtxt(file,delimiter=",") / N
gse = np.loadtxt(gsefile,delimiter=",") / N

ons = 100

x1 = np.linspace(0,5,100)
x2 = np.linspace(5,15,100)
x3 = np.linspace(15,35,100)
y1 = vector[1:ons+1,:]
y2 = vector[ons+2:2*ons+2,:]
y3 = vector[2*ons+3:3*ons+3,:]
plt.figure()
plt.subplot(1,3,1)

plt.plot(x1,y1[:,0],"y")
plt.plot([0,5],[gse,gse],"r:")
plt.plot(x1,y1[:,0]+y1[:,1],"y-.")
plt.plot(x1,y1[:,0]-y1[:,1],"y-.")

plt.ylabel("Energy")
plt.xlabel("Inverse Temperature")
plt.legend(["Calculated Energy","DMRG Ground State Energy"])


plt.subplot(1,3,2)
x1 = np.hstack((x1,x2))
y1 = np.vstack((y1,y2))

plt.plot(x1,y1[:,0],"y")
plt.plot([0,15],[gse,gse],"r:")
plt.plot(x1,y1[:,0]+y1[:,1],"y-.")
plt.plot(x1,y1[:,0]-y1[:,1],"y-.")
plt.plot([5,5],[gse,0],":",color="orange")

plt.ylabel("Energy")
plt.xlabel("Inverse Temperature")
plt.legend(["Calculated Energy","DMRG Ground State Energy"])


plt.subplot(1,3,3)
x1 = np.hstack((x1,x3))
y1 = np.vstack((y1,y3))

plt.plot(x1,y1[:,0],"y")
plt.plot([0,35],[gse,gse],"r:")
plt.plot(x1,y1[:,0]+y1[:,1],"y-.")
plt.plot(x1,y1[:,0]-y1[:,1],"y-.")
plt.plot([5,5],[gse,0],":",color="orange")
plt.plot([15,15],[gse,0],":",color="orange")

plt.ylabel("Energy")
plt.xlabel("Inverse Temperature")
plt.legend(["Calculated Energy","DMRG Ground State Energy"])
plt.show()
