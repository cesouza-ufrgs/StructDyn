#!/usr/bin/env python3
'''
Este e um script para simular um sistema de 2 GDL, sujeito a movimento de base.

'''


import matplotlib.pyplot as plt

import numpy as np
from numpy.linalg import inv, eig
from numpy import pi, sin, cos

from scipy.integrate import odeint

from numpy.fft import fftn


def basemove(t,Y,v,L,d12):


    ti=1
    tf=1.1


    y1 = 0
    vy1= 0
    y2 = 0
    vy2= 0

    omega = 2.*pi*v/L

    # diferenca de tempo para o inicio da
    dt12  = d12/v

    if t>ti and  t <= tf:

        y1 = Y*sin(omega*(t-ti))

        # derivadas

        vy1 = omega*Y*cos(omega*(t-ti))

        print(ti,tf,t,y1,vy1)


    if t>ti+dt12 and  t <= tf+dt12:

        # o y2 precisa descontar a chegada no obstaculo
        y2 = Y*sin(omega*(t-ti-dt12))


        vy2 = omega*Y*cos(omega*(t-ti-dt12))


    print(ti,tf,t,y1,y2,vy1,vy2)
    return y1,y2,vy1,vy2



# function that returns dy/dt
def compute_dqdt(q,t, A,B,u,k1,k2,c1,c2,Y,v,L,d12):

    dqdt = np.empty((4,))


    y1,y2,vy1,vy2= basemove(t,Y,v,L,d12)


    u[0] =  k1*y1   +k2*y2   +c1*vy1   +c2*vy2
    u[1] = -k1*y1*d1+k2*y2*d2-c1*vy1*d1+c2*vy2*d2


    dqdt = np.dot(A , q)+np.dot(B , u)

    return dqdt


# compute fft returns half signal
def ffthalfamp(t, y, T, N):

    yf   = fftn(y)
    #yf1   = np.fftn(q[:,1])

    amp = np.abs(yf) # get amplitude spectrum

    amp = (2/amp.size)*amp[0:amp.size//2]

    freq = np.linspace(0.0, 1.0 /(2.0*T) , N//2) # get freq axis

    return freq, amp


#-----------------------------------


#Constantes do problema
M = 1.5
r = .125
J= M*r*r

k1=100
k2=1000

c1 = 10
c2 = 10


d1=.12
d2=.13

d12 = d1+d2

# entrada do sistema - movimento de base
Y = 0.05
vel = 5
L = 1



#Condições iniciais
x10=0.0
x20=0.0
v10=0.0
v20=0.0

#matrizes do sistema
ngl = 2
M=np.array([[M, 0],
            [0, J]])

K=np.array([[k1+k2,      k2*d2-k1*d1],
            [k2*d2-k1*d1, k1*d1**2+k2*d2*2]])


C=np.array([[c1+c2, c2*d2-c1*d1    ],
            [c2*d2-c1*d1, c1*d1**2+c2*d2**2]])

print('matrix M\n',M)
print('matrix K\n',K)

q0 = np.empty((4,))
#q0 = np.zeros((2*ngl))

q0[0] = x10;
q0[1] = x20;
q0[2] = v10;
q0[3] = v20;


#Resolvendo o problema de autovalores e autovetores

D, V=eig(np.dot(inv(M),K))
#Obtendo as frequências
omega1=np.sqrt(D[1])
omega2=np.sqrt(D[0])


print( 'frequencies')
print( '1', omega1,' f1:',omega1/2/pi,' Hz')
print( '2', omega2,' f2:',omega2/2/pi,' Hz')


print(V)
#Ajustando os autovetores

v1 = V[:,0]
v2 = V[:,1]

v1 = v1/v1[0]
v2 = v2/v2[0]


print('autovetores normalizados\n')
print(v1)
print(v2)





#----------------------------------------
# simulando no espaço de estados

# cria o modelo no espaco de estados

A=np.zeros((2*ngl,2*ngl))

A[0:ngl, ngl:2*ngl]   = np.identity(ngl)
A[ngl:2*ngl,0:ngl]     = np.dot(-inv(M),K)
A[ngl:2*ngl,ngl:2*ngl] = np.dot(-inv(M),C)

print('matrix A\n',A)

B=np.zeros((2*ngl,ngl))

B[ngl:2*ngl,0:ngl] =  inv(M)
print('matrix B\n',B)


u=np.empty((ngl,))
u[0] = 0.
u[1] = 0.


print('q0',q0)
print('u',u)

#----------------------------------------
# time integration with runge kutta

# time points
ttotal = 5
T = 0.001
N = int(np.floor(ttotal/T))
t = np.linspace(0.0, N*T, N)

q = odeint(compute_dqdt,q0,t,rtol=1e-9, args=(A,B,u,k1,k2,c1,c2,Y,vel,L,d12 ))

xcg = q[:,0]
theta = q[:,1]
vcg = q[:,2]
vtheta = q[:,3]

x1 = xcg-d1*theta
x2 = xcg+d2*theta
v1 = vcg-d1*vtheta
v2 = vcg+d2*vtheta

#----------------------------------------
#----------------------------------------
# plotando tudo
fig, axs = plt.subplots(3)

axs[0].title.set_text('time signal')
axs[0].plot(t,x1,label='x1')
axs[0].plot(t,x2,label='x2')
axs[0].set(ylabel='displ $[m]')
axs[0].set_xlim([0,4])
axs[0].legend()

axs[1].plot(t,v1,label='v1')
axs[1].plot(t,v2,label='v2')
axs[1].set(ylabel='vel. [m/s]')
axs[1].set_xlim([0,4])
axs[1].legend()


#----------------------------------------
#-- plot the y function
#----------------------------------------
funy = np.vectorize(basemove,otypes=[float,float,float,float])

yt=funy(t,Y,vel,L,d12)


axs[2].plot(t,yt[0] , label='y1')
axs[2].plot(t,yt[1] , label='y2')
axs[2].set(ylabel='y(t) [m]')
axs[2].set(xlabel='time [s]')
axs[2].set_xlim([0,4])
axs[2].legend()


plt.tight_layout()


#----- FFT
freq, yf  = ffthalfamp(t, yt[0], T, N)
freq, yf0 = ffthalfamp(t, q[:,0], T, N)
freq, yf1 = ffthalfamp(t, q[:,1], T, N)

fig2, axs2 = plt.subplots(3)
axs2[0].title.set_text('frequency domain')
axs2[0].plot(freq, yf1)
axs2[0].set(ylabel='$A_2(f)$[m/s]')
axs2[1].plot(freq, yf0)
axs2[1].set(ylabel='$A_1(f)$[m/s]')
axs2[2].plot(freq, yf)
axs2[2].set(ylabel='$A_f(f)$[m/s]')
axs2[2].set(xlabel='f [Hz]')

plt.tight_layout()
plt.show()


#plt.plot(T, x1h(T)+x1p(T), label='xt(t)')
#plt.plot(T, x2h(T)+x2p(T), label='θt(t)')
#plt.legend()
#plt.xlabel('Tempo [s]')
#plt.ylabel('Posição [m/rad]')
#plt.grid()
#plt.show()
