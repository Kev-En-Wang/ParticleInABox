# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 22:31:36 2019

@author: Kevin
"""
import math
import matplotlib.pyplot as plt
h=6.626*10**-34
T1=300
T2=1000
T3=50
c=3*10**10
k=1.38*10**-23
B=0.037
sigma=2
Plist1=[]
Jlist1=[]
Plist2=[]
Jlist2=[]
Plist3=[]
Jlist3=[]
for J in range (400):
    Jlist1.append(J)
    P=math.exp(-h*B*J*c*(J+1)/(k*T1))*(J*2+1)/((k*T1)/(sigma*h*c*B))
    Plist1.append(P)

for J in range (400):
    Jlist2.append(J)
    P=math.exp(-h*B*J*c*(J+1)/(k*T2))*(J*2+1)/((k*T2)/(sigma*h*c*B))
    Plist2.append(P)

for J in range (400):
    Jlist3.append(J)
    P=math.exp(-h*B*J*c*(J+1)/(k*T3))*(J*2+1)/((k*T3)/(sigma*h*c*B))
    Plist3.append(P)
    
plt.plot (Jlist1,Plist1)
plt.plot (Jlist2,Plist2)
plt.plot (Jlist3,Plist3)
plt.xlabel("J value")
plt.ylabel("probability distribution")

plt.show()
