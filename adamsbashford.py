import matplotlib.pyplot as plt
import numpy as np
import math

#Define parameters
h1=45*1e-3 #Step size
x1=np.arange(0,10,h1) #Numerical grid
y1=np.zeros(len(x1))
z1=np.zeros(len(x1))
y1[0]=90 #Initial condition
z1[0]=0 #Initial condition

#Explicit Adams Bashford Orde 4 Method
def fa(x):return math.exp(-0.16*x)*(5.04429*
                                        math.sin(2.85472*x)+90*math.cos(2.85472*x))
def dydx(x,y,z):return z
def ddyddx(x,y,z):return (-0.32*z)-(8.175*y)
def rk4(dydx,ddyddx,x,y,z,h):
    k1=h*dydx(x,y,z)
    l1=h*ddyddx(x,y,z)
    k2=h*dydx(x+h/2.0,y+k1/2.0,z+l1/2.0)
    l2=h*ddyddx(x+h/2.0,y+k1/2.0,z+l1/2.0)
    k3=h*dydx(x+h/2.0,y+k2/2.0,z+l2/2.0)
    l3=h*ddyddx(x+h/2.0,y+k2/2.0,z+l2/2.0)
    k4=h*dydx(x+h,y+k3,z+l3)
    l4=h*ddyddx(x+h,y+k3,z+l3)
    yi=y+(k1+2.0*k2+2.0*k3+k4)/6.0
    zi=z+(l1+2.0*l2+2.0*l3+l4)/6.0
    return yi,zi

X1=[0]
YAB4=[y1[0]]
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
plt.plot(X1,YA,'bo',label='Analytical Method')
plt.plot(X1,YAB4,'ro',label='AB4 for h=0.09')
plt.xlabel('t (s)')
plt.ylabel('theta (θ)')
plt.legend()
plt.show()
