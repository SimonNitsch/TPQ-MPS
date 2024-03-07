import numpy as np
import matplotlib.pyplot as plt
import os



def Plot(filename, filename2, intervals, particles):
    x = []
    sum_i = 0
    YE = np.loadtxt(filename+"_E.txt",delimiter=",")
    G = np.loadtxt(filename2+"_GSE.txt") / particles
    print(G)
    
    interval_steps = int((YE.shape[0]-1) / intervals.shape[0])
    print(interval_steps)
    
    CS = True
    W = True
    try:
        YC = np.loadtxt(filename+"_C.txt",delimiter=",")
        YS = np.loadtxt(filename+"_S.txt",delimiter=",")
        YS = np.insert(YS,0,np.array([0,0]),0)
        yc = []
        ys = []
    except:
        CS = False
    try:
        YW = np.loadtxt(filename+"_W.txt",delimiter=",")
        yw = []
    except:
        W = False
        
    ye = []
        
              
        
    for i in range(intervals.size):
        x.append(np.linspace(sum_i,sum_i+intervals[i],interval_steps)**(-1))
        sum_i += intervals[i]
        
        ye.append((YE[i*interval_steps:(i+1)*interval_steps,:])/particles)
        if CS:
            yc.append((YC[i*interval_steps:(i+1)*interval_steps,:])/particles)
            ys.append((YS[i*interval_steps:(i+1)*interval_steps,:])/particles)
        if W:
            yw.append(YW[i*interval_steps:(i+1)*interval_steps,:])
        

    x = np.hstack(x)
    xb = 1/x
    ye = np.vstack(ye)
    plt.figure()
    plt.plot(x,ye[:,0],"b")
    plt.plot([x[-1],x[1]],[G,G],"r--")
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
        yc = np.vstack(yc)
        ys = np.vstack(ys)
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
        yw = np.vstack(yw)
        plt.figure()
        plt.plot(x,yw[:,0],"y")
        plt.plot(x,yw[:,0]+yw[:,1],"y-.")
        plt.plot(x,yw[:,0]-yw[:,1],"y-.")
        plt.xscale("log")
        plt.legend(["Plaquette Flux"])
        plt.show()
        
        


if __name__=="__main__":
    
    intervals = np.array([100])
    particles = 48
    name = "Check5TwoSiteBigTest"
    name2 = "ResutFolder2/43GSE"
    
    
    
    os.chdir("../SomeResults")    
    Plot(name,name2,intervals,particles)


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