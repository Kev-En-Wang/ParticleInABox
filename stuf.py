#Made by Zheng Yi Wang for CHEM 350, February 12th 2019
#Defining some constants and importing some libraries
#we assume that the potential of the infinite square well is 0
import numpy as np
import math
import matplotlib.pyplot as plt
pi=3.1415
hbar=1/pi
m=1
a=10#length of the total well
af=1#length of the finite well
ai=a-af#length of the infinite well
Ei=64/(8*ai) #energy of the infinite sqaure well

def Energyfinite(Vo,af): #function to find energy of the finite well
    ######### arrays for the intercepts
    symint=[]
    asymint=[]
    symbound=[]
    asymbound=[]
    #lists for graphs
    sym=[] #for even parity
    asym=[] #for odd parity
    blist=[]
    zlist=[]
    #########
    qo=1.0 #this is used to find the intercepts
    q=1.0
    p=1.0
    po=1.0
    E=-Vo+0.00001 #used to scan over an energy
    z=(af/hbar)*(2*m*(Vo+E))**0.5 #refer to discussion
    zo=(af/hbar)*(2*m*Vo)**0.5 #taken from Griffiths
    ######
    ######This finds all the bound states within the finite well
    while zo>=z: #while the particle is bound
        if q==-qo: #this is also part of finding when the b value intersects
            symint.append(z)
            qo=q
        elif p==-po:
            asymint.append(z)
            po=p
    ################################
        z=(af/hbar)*(2*m*(Vo+E))**0.5#calculates the current z value
        zo=(af/hbar)*(2*m*Vo)**0.5
        b=((zo**2-z**2)**0.5) #calculates the current b value
        sym.append(z*(math.tan(z)))
        asym.append(-z/(math.tan(z)))
        blist.append(b)
        zlist.append(z)
        p=(np.sign(b+z/(math.tan(z))))#sees if it intersected or not
        q=(np.sign(b-z*(math.tan(z))))
        E=E+0.0001 #this is to make sure we have a nice smooth curve
    ###################
    plt.plot(zlist,sym)
    plt.plot(zlist,asym)
    plt.plot(zlist,blist)
    plt.ylabel("function value")
    plt.xlabel("z value")
    plt.ylim(0,zo)
    plt.xlim(zlist[0],z)
    plt.show()
    #this is for showing the plot of the finite state
    #####################
    #this filters all the intercepts so the first one is counted only
    for i in range(0,len(symint)):
        if i%2==0:
            symbound.append(symint[i])
        else:
            continue
    for o in range(0,len(asymint)):
        if o%2==0:
            asymbound.append(asymint[o])
        else:
            continue
    i=o=0 #resets the variable
        ##################################
    Ef=0 #energy of the finite well
    for i in range(0,len(symbound)-1):
        Ef=Ef+(2*hbar*(symbound[i])**2)/(m*af**2)
    if len(symbound)==1:
        Ef=Ef+(2*hbar*(symbound[0])**2)/(m*af**2)
    for o in range(0,len(asymbound)-1):
        Ef=Ef+(2*hbar*(asymbound[o])**2)/(m*af**2)
    if len(asymbound)==1:
        Ef=Ef+(2*hbar*(asymbound[0])**2)/(m*af**2)
        #sums up all the energy from the states in the finite well
    return (Ef)            
##########################
for x in range (1,99):   ##When I change the potential to match
    af=1                ##
    ai=a-af             ##
    Vo=0.01*x*Ei        ##if the whole system was an infinite square well
    Ef=Energyfinite(Vo,af)#####################
    print("This is the finite energy:",Ef)   ##
    print("This is the infinite energy:",Ei) ##
    print("This is the ratio:",Ef/Ei)        ##
###############################################
for y in range(1,30):   ##When I make the well as wide as the whole 
    Vo=0.1*Ei           ##well
    af=y*0.01           ##
    ai=a-af             #length of the infinite well
    Ei=64/(8*ai)        #energy of the infinite sqaure well
    Ef=Energyfinite(Vo,af)#####################
    print("This is the finite energy:",Ef)   ##
    print("This is the infinite energy:",Ei) ##
    print("This is the ratio:",Ef/Ei)        ##
###############################################