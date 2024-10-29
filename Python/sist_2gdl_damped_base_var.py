#!/usr/bin/env python3


import matplotlib.pyplot as plt

import numpy as np
from numpy.linalg import inv, eig
from numpy import pi, sin, cos

from scipy.integrate import odeint

from numpy.fft import fftn


def basemove(t,Y,omega):
   
    y=Y*np.sin(omega*t)
    
    return y

# defines the base movement    
def compQ(t, k1, c1, Y, omega):
        
    Q = k1*Y*np.sin(omega*t) + c1*omega*Y*np.cos(omega*t)
    
    return Q
 
    
# function that returns dy/dt
def compute_dqdt(q,t, A,B,u,k1,c1,Y,omega):
   
    dqdt = np.empty((4,))
        
    u[0] = compQ(t, k1, c1, Y, omega)  
        
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
m =5
m1=m
m2=10*m
k1=4000
k2=8000

# entrada do sistema - movimento de base
Y = 0.05
f = 15
omega = 2*pi*f

c1 =  0
c2 = 100

#Condições iniciais
x10=0.0
x20=0.0
v10=0.0
v20=0.0

#matrizes do sistema
ngl = 2
M=np.array([[m1, 0],
            [0, m2]])

K=np.array([[k1+k2,   -k2],
            [-k2, k2]])


C=np.array([[c1, 0    ],
            [0, c2]])

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

'''   
simulando no espaço de estados
'''
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
N = 2048
T = 0.005
t = np.linspace(0.0, N*T, N)

q = odeint(compute_dqdt,q0,t, args=(A,B,u,k1,c1,Y,omega, ))


#----------------------------------------
#----------------------------------------
# plot 
fig, axs = plt.subplots(3)

axs[0].title.set_text('time signal')
axs[0].plot(t,q[:,1])
#axs[0].plot(t,q[:,1])
axs[0].set(ylabel='$x_2$[m/s]')

axs[1].plot(t,q[:,0])
#axs[1].plot(t,q[:,3])
axs[1].set(ylabel='$x_1$[m/s]')


#----------------------------------------
#-- plot the y function
#----------------------------------------
yt = np.vectorize(basemove)
axs[2].plot(t, yt(t,Y,omega))
axs[2].set(ylabel='y(t) [m/s]')
axs[2].set(xlabel='time')


plt.tight_layout()


#----- FFT
freq, yf  = ffthalfamp(t, yt(t,Y,omega), T, N)
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
