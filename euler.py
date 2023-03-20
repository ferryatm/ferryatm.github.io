import matplotlib.pyplot as plt
import numpy as np
import math

#Define parameters
h1=9*1e-2 #Step size
x1=np.arange(0,10,h1) #Numerical grid
y1=np.zeros(len(x1))
z1=np.zeros(len(x1))
y1[0]=90 #Initial condition
z1[0]=0 #Initial condition

#Explicit Euler Method
def fa(x):return math.exp(-0.16*x)*(5.04429*
                                        math.sin(2.85472*x)+90*math.cos(2.85472*x))
def dydx(x,y,z):return z
def ddyddx(x,y,z):return (-0.32*z)-(8.175*y)
def mp(dydx,ddyddx,x,y,z,h):
    k1=h*dydx(x,y,z)
    l1=h*ddyddx(x,y,z)
    yi=y+k1
    zi=z+l1
    return yi,zi

X1=[0]
Y1MP=[y1[0]]
YA=[fa(0)]

for i in range(0,4):
    yi,zi=rk4(dydx,ddyddx,x1[i],y1[i],z1[i],h1)
    ya=fa(x1[i])
    y1[i+1]=yi
    z1[i+1]=zi

for i in range(3,len(x1)-1):
    z1[i+1]=z1[i]+((h1/24)*(55*ddyddx(x1[i],y1[i],z1[i])-59*
                     ddyddx(x1[i-1],y1[i-1],z1[i-1])+37*ddyddx(x1[i-2],y1[i-2],z1[i-2])-9*ddyddx(x1[i-3],y1[i-3],z1[i-3])))
    y1[i+1]=y1[i]+((h1/24)*(55*z1[i]-59*
                     z1[i-1]+37*z1[i-2]-9*z1[i-3]))
    ya=fa(x1[i+1])
    X1.append(x1[i+1])
    YA.append(ya)
    YAB4.append(y1[i+1])

#Plotting Graphic
plt.figure(figsize=(12,8),dpi=75)
plt.plot(Y1,YA,'b',label='Analytical Method')
plt.plot(X1,Y1MP,'r',label='E for h=0.09')
plt.xlabel('t (s)')
plt.ylabel('theta (θ)')
plt.legend()
plt.show()
