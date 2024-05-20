import numpy as np
import matplotlib.pyplot as plt

txtname = "Case5075.txt"
foldername = "/Users/simon/Documents/TPQ-MPS/SomeResults/ExactEnergies/"
particles = 12


filename = foldername + txtname
a = np.loadtxt(filename,usecols=3)
x2 = np.logspace(-2,2,200)
y2 = np.zeros(x2.size)

for i in range(x2.size):
    aexp = np.exp(-x2[i]*a)
    aexp /= np.sum(aexp)
    y2[i] = x2[i]**2 * (np.sum(a**2*aexp) - np.sum(a*aexp)**2)
    #y2[i] = np.sum(a*aexp)

x2 = np.flip(x2)
y2 /= particles

plt.figure()
plt.plot(x2,y2)
plt.xscale("log")
plt.show()



