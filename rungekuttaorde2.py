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

h2=45*1e-3 #Step size
x2=np.arange(0,10,h2) #Numerical grid
y2=np.zeros(len(x2))
z2=np.zeros(len(x2))
y2[0]=90
z2[0]=0

#Explicit Runge-Kutta Orde 2 Method
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
Y1RK4=[y1[0]]
X2=[0]
Y2RK4=[y2[0]]
YA=[fa(0)]

for i in range(len(x1)-1):
    X1.append(x1[i])
    yi,zi=rk4(dydx,ddyddx,x1[i],y1[i],z1[i],h1)
    ya=fa(x1[i])
    y1[i+1]=yi
    z1[i+1]=zi
    Y1RK4.append(yi)

for i in range(len(x2)-1):
    X2.append(x2[i])
    yi,zi=rk4(dydx,ddyddx,x2[i],y2[i],z2[i],h2)
    ya=fa(x2[i])
    y2[i+1]=yi
    z2[i+1]=zi
    Y2RK4.append(yi)
    YA.append(ya)

#Plotting Graphic
plt.figure(figsize=(12,8),dpi=75)
plt.plot(X2,YA,'b',label='Analytical Method')
plt.plot(X1,Y1RK4,'r',label='RK for h=0.09')
plt.plot(X2,Y2RK4,'k',label='RK for h=0.045')
plt.xlabel('t (s)')
plt.ylabel('theta (θ)')
plt.legend()
plt.show()