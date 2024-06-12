import numpy as np
import matplotlib.pyplot as plt
import os


    
intervals = np.array([1,9,25,25,25,25,40,50,50,50,50,50,100])
particles = 24
filename = "MagnetTest/Dataproto/"
filename2 = "Classic"
K = 1.33
    
os.chdir("../SomeResults")   






x = []
sum_i = 0
YE = np.loadtxt(filename+"E.txt",delimiter=",")
G = np.loadtxt(filename2+"_GSE.txt") / particles
print(G)

Curie = lambda T, J : 1 / (4*T + J)
Dimerz = lambda T, J : 1/2/T * np.exp(J/4/T) / (np.exp(J/4/T) - np.exp(-J/4/T))
Dimerx = lambda T, J : 1/J * np.tanh(J/4/T)
    
interval_steps = int((YE.shape[0]-1) / intervals.shape[0])
print(interval_steps)
    
CS = True
W = True
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
    YCHIx = np.loadtxt(filename+"Chix.txt",delimiter=",")
    YCHIy = np.loadtxt(filename+"Chiy.txt",delimiter=",")
    YCHIz = np.loadtxt(filename+"Chiz.txt",delimiter=",")
except:
    Chi = False
        
        
              
        
for i in range(intervals.size):
    x.append(np.linspace(sum_i,sum_i+intervals[i],interval_steps,endpoint=False)**(-1))
    sum_i += intervals[i]
        
x.append(np.array([1/sum_i]))
        

x = np.hstack(x)[1:]
xb = 1/x
ye = YE[1:,:] / particles
plt.figure()
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
    yc = YC[1:,:] / particles
    ys = YS[1:,:] / particles
    plt.figure()
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
    yw = YW[1:,:]
    plt.figure()
    plt.plot(x,yw[:,0],"y")
    plt.plot(x,yw[:,0]+yw[:,1],"y-.")
    plt.plot(x,yw[:,0]-yw[:,1],"y-.")
    plt.xscale("log")
    plt.legend(["Plaquette Flux"])
    plt.show()
        
if Chi:
    ycx = YCHIx[1:,:] / particles
    ycy = YCHIy[1:,:] / particles
    ycz = YCHIz[1:,:] / particles
    plt.figure()
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
    plt.plot(x,Dimerx(x,K),"c:")
    plt.plot(x,Dimerz(x,K),"m:")
    plt.xscale("log")
    plt.ylim([0,20])
    plt.show()
        




        
        
 


plt.figure(figsize=[10,10])

ymax = 4
xmax = 2

for i in range(ymax):
    for j in range(xmax):
        xh = j * 2 * np.cos(np.pi/6) + np.ceil(i/2) * np.cos(np.pi/6)
        yh = np.ceil(i/2) * np.sin(np.pi/6) + np.floor(i/2)

        if i%2==1:
            plt.plot([xh,xh-np.cos(np.pi/6)],[yh,yh-np.sin(np.pi/6)],"r")
            if j!=(xmax-1):
                plt.plot([xh,xh+np.cos(np.pi/6)],[yh,yh-np.sin(np.pi/6)],"b")
            if i!=(ymax-1):
                plt.plot([xh,xh],[yh,yh+1],"g")


        if i==0:
            plt.plot([xh,xh],[yh,yh-0.5],"g:")
        if i==ymax-1:
            plt.plot([xh,xh],[yh,yh+0.5],"g:")
        if i==0 and j ==0:
            plt.plot([xh,xh-0.5],[yh,yh],"k:")
        if i==ymax-1 and j==xmax-1:
            plt.plot([xh,xh+0.5],[yh,yh],"k:")
        if j==0 and i%2==0:
            plt.plot([xh,xh-0.5*np.cos(np.pi/6)],[yh,yh+0.5*np.sin(np.pi/6)],"b:")
        if j==xmax-1 and i%2==1:
            plt.plot([xh,xh+0.5*np.cos(np.pi/6)],[yh,yh-0.5*np.sin(np.pi/6)],"b:")

for i in range(ymax):
    for j in range(xmax):
        xh = j * 2 * np.cos(np.pi/6) + np.ceil(i/2) * np.cos(np.pi/6)
        yh = np.ceil(i/2) * np.sin(np.pi/6) + np.floor(i/2)
        plt.plot(xh,yh,"k.")
                
            
plt.axis("equal")
plt.show()



plt.figure(figsize=[10,10])

ymax = 6
xmax = 4



'''
os.chdir("../SomeResults")

file = "1_less"



x1 = np.linspace(0,5,100)
x2 = np.linspace(5,15,100)
x3 = np.linspace(15,35,100)
x4 = np.linspace(35,55,100)
x5 = np.linspace(55,75,100)
x6 = np.linspace(75,100,100)

x = np.hstack((x1,x2,x3,x4,x5,x6)) **-1

YE = np.loadtxt(file+"_E.txt",delimiter=",")
YC = np.loadtxt(file+"_C.txt",delimiter=",")
YS = np.loadtxt(file+"_S.txt",delimiter=",")


YE = np.vstack((YE[1:101,:],YE[102:202,:],YE[203:303,:],YE[304:404,:],YE[405:505,:],YE[506:606,:])) / 24
YC = np.vstack((YC[1:101],YC[102:202],YC[203:303],YC[304:404],YC[405:505],YC[506:606])) / 24
YS = np.vstack((YS[0:100],YS[101:201],YS[202:302],YS[303:403],YS[404:504],YS[505:605])) / 24


plt.figure()
plt.plot(x,YE[:,0],"b")
plt.plot(x,YE[:,0]+YE[:,1],"b-.")
plt.plot(x,YE[:,0]-YE[:,1],"b-.")
plt.xscale("log")
plt.show()
plt.plot()
plt.plot(x,YC[:,0],"r")
plt.plot(x,YS[:,0],"g")
plt.legend(["Heat Capacity","Entropy Integral"])
plt.plot(x,YC[:,0]+YC[:,1],"r-.")
plt.plot(x,YC[:,0]-YC[:,1],"r-.")
plt.plot(x,YS[:,0]+YS[:,1],"g-.")
plt.plot(x,YS[:,0]-YS[:,1],"g-.")
plt.xscale("log")
plt.show()



'''