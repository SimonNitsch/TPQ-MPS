import numpy as np
import matplotlib.pyplot as plt

txtname = "Case4874.txt"
foldername = "/Users/simon/Documents/TPQ-MPS/SomeResults/ExactEnergies/"
particles = 12


filename = foldername + txtname
a = np.loadtxt(filename,usecols=3)
x2 = np.logspace(-1.5,3,200)
y2 = np.zeros(x2.size)

for i in range(x2.size):
    aexp = np.exp(-x2[i]*a)
    aexp /= np.sum(aexp)
    y2[i] = x2[i]**2 * (np.sum(a**2*aexp) - np.sum(a*aexp)**2)
    #y2[i] = np.sum(a*aexp)

x2 = 1 / x2
y2 /= particles

plt.figure(figsize=[12,6])
plt.plot(x2,y2)
plt.xlabel("T")
plt.ylabel("C")
plt.xscale("log")
plt.legend(["Exact Diagonalization"])
plt.show()



'''

plt.figure(figsize=[12,6])
plt.plot(x,yc[:,0],"r")
plt.plot(x2,y2)
plt.legend(["TDVP TPQ-MPS","Exact Diagonalisation"])
plt.plot(x,yc[:,0]+yc[:,1],"r-.")
plt.plot(x,yc[:,0]-yc[:,1],"r-.")
plt.xlabel("T")
plt.ylabel("C")
plt.xscale("log")



plt.figure(figsize=[12,6])
plt.plot(x,yc[:,0],"r")
plt.plot(x2,y2)
plt.legend(["tanTRG","Exact Diagonalisation"])
plt.plot(x,yc[:,0]+yc[:,1],"r-.")
plt.plot(x,yc[:,0]-yc[:,1],"r-.")
plt.xlabel("T")
plt.ylabel("C")
plt.xscale("log")

'''

