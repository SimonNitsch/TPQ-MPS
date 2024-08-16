import numpy as np
import matplotlib.pyplot as plt
import os


particles = 16
filename = "Issp/4x4S1-01m/"
filename2 = "Classic"
K = 1.
    
os.chdir("../SomeResults")   





x = np.loadtxt(filename+"xdata.txt")

YE = np.loadtxt(filename+"E.txt",delimiter=",")
G = np.loadtxt(filename2+"_GSE.txt") / particles
print(G)

Curie = lambda T, J : 1 / (4*T + J)
Dimerz = lambda T, J : 1/2/T * np.exp(J/4/T) / (np.exp(J/4/T) - np.exp(-J/4/T))
Dimerx = lambda T, J : 1/J * np.tanh(J/4/T)
    
    
CS = True
W = True
M = True
Chi = True
try:
    YC = np.loadtxt(filename+"C.txt",delimiter=",")
    YS = np.loadtxt(filename+"S.txt",delimiter=",")
except:
    CS = False
try:
    YW = np.loadtxt(filename+"W.txt",delimiter=",")
    yw = []
except:
    W = False
try: 
    YMx = np.loadtxt(filename+"Mx.txt",delimiter=",")
    YMy = np.loadtxt(filename+"My.txt",delimiter=",")
    YMz = np.loadtxt(filename+"Mz.txt",delimiter=",")
    YMx2 = np.loadtxt(filename+"Mx2.txt",delimiter=",")
    YMy2 = np.loadtxt(filename+"My2.txt",delimiter=",")
    YMz2 = np.loadtxt(filename+"Mz2.txt",delimiter=",")
except:
    M = False
try:
    YCHIx = np.loadtxt(filename+"Chix.txt",delimiter=",")
    YCHIy = np.loadtxt(filename+"Chiy.txt",delimiter=",")
    YCHIz = np.loadtxt(filename+"Chiz.txt",delimiter=",")
except:
    Chi = False

        
        
        
xb = 1/x
ye = YE / particles
plt.figure(figsize=[12,6])
plt.plot(x,ye[:,0],"b")
plt.plot([x[-1],x[0]],[G,G],"r--")
plt.legend(["Energy","Ground State Energy"])
plt.plot(x,ye[:,0]+ye[:,1],"b-.")
plt.plot(x,ye[:,0]-ye[:,1],"b-.")
plt.xscale("log")
plt.xlabel("T")
plt.show()
    
Emin = np.min(ye[:,0])
sum = 0
for i in range(xb.size-1):
    sum += (xb[i+1]-xb[i])/2 *(ye[i,0]+ye[i+1,0]-2*Emin)
print(sum)
    
    
    
if CS:
    yc = YC / particles
    ys = YS / particles
    plt.figure(figsize=[12,6])
    plt.plot(x,yc[:,0],"r")
    plt.plot(x,ys[:,0],"g")
    plt.legend(["Heat Capacity", "Entropy Integral"])
    plt.plot(x,yc[:,0]+yc[:,1],"r-.")
    plt.plot(x,yc[:,0]-yc[:,1],"r-.")
    plt.plot(x,ys[:,0]+ys[:,1],"g-.")
    plt.plot(x,ys[:,0]-ys[:,1],"g-.")
    plt.xscale("log")
    plt.xlabel("T")
    plt.show()
        
if W:
    yw = YW
    plt.figure(figsize=[12,6])
    plt.plot(x,yw[:,0],"y")
    plt.plot(x,yw[:,0]+yw[:,1],"y-.")
    plt.plot(x,yw[:,0]-yw[:,1],"y-.")
    plt.xscale("log")
    plt.legend(["Plaquette Flux"])
    plt.xlabel("T")
    plt.show()

if M:
    ymx = YMx / particles
    ymy = YMy / particles
    ymz = YMz / particles
    ymx2 = YMx2 / particles
    ymy2 = YMy2 / particles
    ymz2 = YMz2 / particles

    plt.figure(figsize=[12,6])
    plt.plot(x,ymx[:,0],"c")
    plt.plot(x,ymy[:,0],"y")
    plt.plot(x,ymz[:,0],"m")
    plt.legend([r"$M_x$", r"$M_y$", r"$M_z$"])
    plt.plot(x,ymx[:,0]+ymx[:,1],"c-.")
    plt.plot(x,ymx[:,0]-ymx[:,1],"c-.")
    plt.plot(x,ymy[:,0]+ymy[:,1],"y-.")
    plt.plot(x,ymy[:,0]-ymy[:,1],"y-.")
    plt.plot(x,ymz[:,0]+ymz[:,1],"m-.")
    plt.plot(x,ymz[:,0]-ymz[:,1],"m-.")
    plt.xscale("log")
    plt.xlabel("T")
    plt.show()

    plt.figure(figsize=[12,6])
    plt.plot(x,ymx2[:,0],"c")
    plt.plot(x,ymy2[:,0],"y")
    plt.plot(x,ymz2[:,0],"m")
    plt.legend([r"$M_x^2$", r"$M_y^2$", r"$M_z^2$"])
    plt.plot(x,ymx2[:,0]+ymx2[:,1],"c-.")
    plt.plot(x,ymx2[:,0]-ymx2[:,1],"c-.")
    plt.plot(x,ymy2[:,0]+ymy2[:,1],"y-.")
    plt.plot(x,ymy2[:,0]-ymy2[:,1],"y-.")
    plt.plot(x,ymz2[:,0]+ymz2[:,1],"m-.")
    plt.plot(x,ymz2[:,0]-ymz2[:,1],"m-.")
    plt.xscale("log")
    plt.xlabel("T")
    plt.show()

        
if Chi:
    ycx = YCHIx / particles
    ycy = YCHIy / particles
    ycz = YCHIz / particles
    plt.figure(figsize=[12,6])
    plt.plot(x,ycx[:,0],"c")
    plt.plot(x,ycy[:,0],"y")
    plt.plot(x,ycz[:,0],"m")
    plt.legend([r"$\chi^x$", r"$\chi^y$", r"$\chi^z$"])
    plt.plot(x,ycx[:,0]+ycx[:,1],"c-.")
    plt.plot(x,ycx[:,0]-ycx[:,1],"c-.")
    plt.plot(x,ycy[:,0]+ycy[:,1],"y-.")
    plt.plot(x,ycy[:,0]-ycy[:,1],"y-.")
    plt.plot(x,ycz[:,0]+ycz[:,1],"m-.")
    plt.plot(x,ycz[:,0]-ycz[:,1],"m-.")
    plt.plot(x,Curie(x,K),"k:")
    plt.xlabel("T")
    #plt.plot(x,Dimerx(x,K),"c:")
    #plt.plot(x,Dimerz(x,K),"m:")
    plt.xscale("log")
    plt.ylim([0,np.max(ycz[:,0])*1.5])
    plt.show()
        





