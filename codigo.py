import math
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np
from sympy import sin,cos

#codigo branch 2 ...""
# VARIABLES
alfa = 1
a = 1
beta = 1
c1 = 1 #debe ser mayor a 0
vr = 1 #solo se puede dar velocidad lineal y angular
wr = math.radians(50)
gamma = 1
n = 1
c2 = 1
a = 2

#POSICIONES
#xe = np.array([5,7,9,10,7])
#ye = np.array([1,2,3,4])
num_datos = 50
t = np.linspace(0,2*np.pi,num_datos) #por defecto son 50 valores 
r = 10
xe = np.zeros(num_datos)
ye = np.zeros(num_datos)
fi = np.zeros(num_datos)
xe = r*np.cos(t)
ye = r*np.sin(t)
v  = np.zeros(xe.shape)
w  = np.zeros(xe.shape)
#print(xe.shape)
#plt.plot(xe,ye)
#plt.show()
dxdt = np.zeros(num_datos)
print(xe[1])
print(ye[1])

def dynamicsStateSpace(x,t):
    fi = (math.pi/2)*((x[1]*vr)/(1+(x[1]*vr)**2))
    v = vr*cos(x[2])+c1*x[0]
    w = (1/(1+(fi*x[0]*vr)))*((gamma*n*x[1]*vr)+(wr)+(fi*(vr**2*sin(x[2])+x[1]*a))+(c2*gamma*(x[2]+fi)))
    dxdt = [w*x[1]-v+vr*cos(x[2]), -w*x[0]+vr*sin(x[2]),wr-w]
    return dxdt

initialState = np.array([-1,-1,-1])
simulationTime = np.linspace(0, 50, 200)
# generate the state-space trajectory
solutionState = odeint(dynamicsStateSpace, initialState, simulationTime)

plt.figure()
fig,axes = plt.subplots()
plt.title('Señales ', fontsize=14)
plt.xlabel('Tiempo')
plt.ylabel('velocidad')
axes.plot(simulationTime,solutionState[:, 0],'--',color ='r')
axes.plot(simulationTime,solutionState[:, 1],'--')
#axes.plot(simulationTime,solutionState[:, 2],'--')
#axes.plot(solutionState[:, 0],solutionState[:, 1],'--')
plt.savefig('señales.png', dpi=600)

plt.figure()
fig,axes = plt.subplots()
plt.title('Plano de fase ', fontsize=14)
plt.xlabel('x[0]')
plt.ylabel('x[1]')
axes.plot(solutionState[:, 0],solutionState[:, 1],'--')
plt.savefig('plano de fase.png', dpi=600)
