import numpy as np
import matplotlib.pyplot as plt
import os
from numpy import cos,arccos,sin,arcsin,tan,arctan,exp,pi

#MATH OR NUMPY test speed.
#On peut tester performance si T est un dictionnaire T = {"ss" : ..., "sp" : ...,"ps" : ...,"pp" : ...} où ... = complex
def a(n_i,n_t,th_i):
    a_p=(n_t*cos(th_i)-n_i/n_t*(np.sqrt(n_t**2-n_i**2*sin(th_i)**2)))/\
    (n_i/n_t*np.sqrt(n_t**2-n_i**2*sin(th_i)**2) + n_t*cos(th_i))
    a_s=(n_i*cos(th_i)-np.sqrt(n_t**2-n_i**2*sin(th_i)**2))/\
    (n_i*cos(th_i) + np.sqrt(n_t**2-n_i**2*sin(th_i)**2))
    return a_p, a_s

def Jonex_Matrix(th_i, th_r, phi):
    T = np.zeros(4,dtype=complex) # Tss Tsp Tps Tpp  
    T[0]=a(n_i,n_t,th_i)[1]cos(et_i)*cos(et_r)+a(n_i,n_t,th_i)[0]sin(et_i)*sin(et_r)
    T[1]=-a(n_i,n_t,th_i)[1]sin(et_i)*cos(et_r)+a(n_i,n_t,th_i)[0]cos(et_i)*sin(et_r)
    T[2]=-a(n_i,n_t,th_i)[1]cos(et_i)*sin(et_r)+a(n_i,n_t,th_i)[0]sin(et_i)*cos(et_r)
    T[3]=a(n_i,n_t,th_i)[1]sin(et_i)*sin(et_r)+a(n_i,n_t,th_i)[0]cos(et_i)*cos(et_r)
    return T

def Mueller(T):
    M = np.ones((4,4),dtype=complex)
    M[0,0] = np.absolute(T[0])**2 + np.absolute(T[1])**2 + np.absolute(T[2])**2 + np.absolute(T[3])**2
    M[0,1] = np.absolute(T[0])**2 + np.absolute(T[1])**2 - np.absolute(T[2])**2 - np.absolute(T[3])**2
    M[0,2] = T[0]*np.conj(T[2]) + np.conj(T[0]*np.conj(T[2])) + T[1]*np.conj(T[3]) + np.conj(T[1]*np.conj(T[3]))
    M[0,3] = 1j*(T[2]*np.conj(T[0]) - np.conj(T[2]*np.conj(T[0]))) + 1j*(T[3]*np.conj(T[1]) + np.conj(T[3]*np.conj(T[1])))

    M[1,0] = np.absolute(T[0])**2 - np.absolute(T[1])**2 + np.absolute(T[2])**2 - np.absolute(T[3])**2
    M[1,1] = np.absolute(T[0])**2 + np.absolute(T[1])**2 - np.absolute(T[2])**2 + np.absolute(T[3])**2
    M[1,2] = T[0]*np.conj(T[2]) + np.conj(T[0]*np.conj(T[2])) - T[1]*np.conj(T[3]) - np.conj(T[1]*np.conj(T[3]))
    M[1,3] = 1j*(T[2]*np.conj(T[0]) - np.conj(T[2]*np.conj(T[0]))) - 1j*(T[3]*np.conj(T[1]) + np.conj(T[3]*np.conj(T[1])))

    M[2,0] = T[0]*np.conj(T[1]) + np.conj(T[0]*np.conj(T[1])) + T[2]*np.conj(T[3]) + np.conj(T[2]*np.conj(T[3]))
    M[2,1] = T[0]*np.conj(T[1]) + np.conj(T[0]*np.conj(T[1])) - T[2]*np.conj(T[3]) - np.conj(T[2]*np.conj(T[3]))
    M[2,2] = T[0]*np.conj(T[3]) + np.conj(T[0]*np.conj(T[3])) + T[2]*np.conj(T[1]) + np.conj(T[2]*np.conj(T[1]))
    M[2,3] = 1j*(T[2]*np.conj(T[1]) - np.conj(T[2]*np.conj(T[1]))) - 1j*(T[0]*np.conj(T[3]) - np.conj(T[0]*np.conj(T[3])))

    M[3,0] = 1j*(T[0]*np.conj(T[1]) - np.conj(T[0]*np.conj(T[1]))) + 1j*(T[2]*np.conj(T[3]) - np.conj(T[2]*np.conj(T[3])))
    M[3,1] = 1j*(T[0]*np.conj(T[1]) - np.conj(T[0]*np.conj(T[1]))) - 1j*(T[2]*np.conj(T[3]) - np.conj(T[2]*np.conj(T[3])))
    M[3,2] = 1j*(T[0]*np.conj(T[3]) - np.conj(T[2]*np.conj(T[1]))) + 1j*(T[0]*np.conj(T[3]) - np.conj(T[0]*np.conj(T[3])))
    M[3,3] = 1j*(T[0]*np.conj(T[1]) - np.conj(T[0]*np.conj(T[1]))) - 1j*(T[2]*np.conj(T[3]) + np.conj(T[2]*np.conj(T[3])))
    return M/2
 
 def beta(th_i, th_r, phi):
    co_b=cos(th_i)*cos(th_r)+sin(th_i)*sin(th_r)*cos(phi)
    return arccos(co_b)/2

def eta(th_i,th_r,b):
    co_et_i = ((cos(th_i)+cos(th_r))/(2*cos(b))-cos(th_i)*cos(b))/(sin(th_i)*sin(b))
    co_et_r = ((cos(th_i)+cos(th_r))/(2*cos(b))-cos(th_r)*cos(b))/(sin(th_r)*sin(b))
    return arccos(co_et_i),arccos(co_et_r)

def theta(th_i,th_r,b):
    co = (cos(th_i)+cos(th_r))/(2*cos(b))
    return arccos(co)
    
def f(th_i,th_r,phi):
    #f = np.ones((4,4),dtype=complex)
    b = beta(th_i, th_r, phi)
    th = theta(th_i,th_r,b)
    f = exp(-(tan(th)**2)/(2*sigma**2))\
    /(2*pi*4*(sigma**2)*(cos(th)**4)*cos(th_r)*cos(th_i))\
    *Mueller(Jonex_Matrix(th_i, th_r, phi))
    return f
N = 90
th_r = np.linspace(0,np.pi/2,N) # 0 - 90 degrés
phi = np.linspace(0,np.pi/2,N) # 0 - 90 degrés
sigma = 0.15
n_t = 1.57
n_i = 1
th_i = np.pi/3 # 60 degrés
Z = f(th_i,th_r,phi)*cos(th_r)
