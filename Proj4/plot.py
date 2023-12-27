import numpy as np
import matplotlib.pyplot as plt

file = "Heisenberg_1D"

vector = np.fromfile(file, dtype=np.float64)
length = int(np.round(vector.shape[0] / 2))
vector = vector.reshape((length,2)) / (20)

x = np.linspace(0,10,100)
plt.figure()

plt.plot(x,vector[:,0],"y")
plt.plot([0,10],[-0.4438,-0.4438],"r:")
plt.plot(x,vector[:,0]+vector[:,1],"y-.")
plt.plot(x,vector[:,0]-vector[:,1],"y-.")

plt.ylabel("Energy")
plt.xlabel("Inverse Temperature")
plt.legend(["Calculated Energy","Ground State Energy"])
plt.show()


