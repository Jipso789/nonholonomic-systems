import math
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np
from sympy import sin,cos


num_datos = 50
t = np.linspace(0,2*np.pi,num_datos) #por defecto son 50 valores 
xd = 1
yd = 1
tetha_d = math.radians(45)
fi_d = math.radians(45)
k1 = 1
k2 = 1
vd = 0.1
wd = vd/math.sqrt(xd**2+yd**2)
d = 0.3 #DISTANCIA ENTRE LAS RUEDAS
#print(xe.shape)
#plt.plot(xe,ye)
#plt.show()
dxdt = np.zeros(num_datos)
print(wd)
#print(ye[1])

def sec(x): #ESTA FUNCION CALCULA LA SECANTE
    return 1/math.cos(x)


def dynamicsStateSpace(x,t):
    E1 = x[1]-yd-math.tan(x[2]+x[3])*(x[0]-xd)
    E2 = math.tan(tetha_d+fi_d)-math.tan(x[2]+x[3])-(math.sin(x[3])+(sec(x[2]+x[3])**3)*(x[0]-xd))/(d)
    E3 = (math.sin(fi_d)*((sec(tetha_d+fi_d))**3)-math.sin(x[3]*sec(x[2]+x[3])**3))/(d)
    E4 = x[1]-xd
    zi2= (math.sin(x[3])*(sec(x[2]+fi_d))**3)
    ui1d = vd*math.cos(tetha_d+fi_d)
    ui2d = (((3*math.sin(tetha_d)*math.sin(tetha_d+fi_d)*(sec(tetha_d+fi_d))**4+math.cos(fi_d)*sec(tetha_d+fi_d)**3)*wd)/(d))+((3*math.sin(tetha_d+fi_d)*(math.sin(fi_d)**2)*(sec(tetha_d+fi_d)**4)*vd)/(d**2))
    ui2 = ui2d + 2*ui1d*E2-k1*(E1+E3)
    ui1 = ui1d + (ui2*E2-k2*E4)/(1-2*E1*zi2-zi2*E3)
    v = (ui1)/(math.cos(x[2]+x[3]))
    w = (ui2*(d**2)-3*math.sin(x[2]+x[3])*math.sin(x[3])**2*(sec(x[2]+x[3])**4)*v)/(d**2*(3*math.sin(x[3]))*math.sin(x[2]+x[3])*((sec(x[2]+x[3]))**4)+math.cos(x[3])*(sec(x[2]+x[3]))**3)
    #X[0] ES POSICION X, X[1] POSICION EN Y, X[2] TETHA, X[3] FI
    dxdt = [math.cos(x[2])*v, math.sin(x[2])*v,(math.tan(x[3])/d)*v,w]
    return dxdt

initialState = np.array([0,0,0,0])
simulationTime = np.linspace(0,10,2000)
# generate the state-space trajectory
solutionState = odeint(dynamicsStateSpace, initialState, simulationTime)

plt.figure()
fig,axes = plt.subplots()
plt.title('Señales ', fontsize=14)
plt.xlabel('Posicion x')
plt.ylabel('Posicion y')
axes.plot(simulationTime,solutionState[:, 0],'--',color ='r')
axes.plot(simulationTime,solutionState[:, 1],'--')
#axes.plot(simulationTime,solutionState[:, 2],'--')
#axes.plot(solutionState[:, 0],solutionState[:, 1],'--')
#axes.plot(solutionState[:, 0],solutionState[:, 1],'--')
plt.savefig('señales_control.png', dpi=600)

#plt.figure()
#fig,axes = plt.subplots()
#plt.title('Plano de fase ', fontsize=14)
#plt.xlabel('x[0]')
#plt.ylabel('x[1]')
#axes.plot(solutionState[:, 0],solutionState[:, 1],'--')
#plt.savefig('plano de fase.png', dpi=600)
