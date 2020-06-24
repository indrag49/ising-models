import random
import numpy as np
import matplotlib.pyplot as plt
import math
import random

L=2**4
N=L**2
iterations=1000
J=1
H=0
Kb=1
Temp=[]
t=0.5
while t<=5:
    Temp+=[t, ]
    t+=0.1

Energy=[]
Magnetization=[]
Susceptibility=[]
Specific_heat=[]
Tm=[]

phi=np.zeros((L, L), dtype=int)
for i in range(L):
    for j in range(L):
        r=random.uniform(0, 1)
        if r<=0.5:
            phi[i, j]=1
        else:
            phi[i, j]=-1

phi1=phi
# implementing the boundary conditions
p=[0]*L
m=[0]*L
for i in range(L):
    p[i]=i+1
    m[i]=i-1
p[-1]=0
m[0]=L-1
    
for T in Temp:
    Beta=1/(Kb*T)
    phi=phi1

    # metropolis algorithm
    E=[]
    M=[]
    chi=[]
    Cv=[]
    time=[]

    for k in range(iterations):
        e=0
        e_sq=0
        mag=0
        mag_sq=0
        for i in range(L):
            for j in range(L):
                en=(J/2)*phi[i, j]*(phi[p[i],j]+phi[m[i],j]+phi[i,p[j]]+phi[i,m[j]])+H*phi[i,j]
                e+=en
                e_sq+=en**2
                mag+=phi[i,j]
                mag_sq+=phi[i,j]**2

        E+=[-e/N, ]
        M+=[mag/N, ]
        chi+=[(mag_sq/N-(mag/N)**2)*Beta, ]
        Cv+=[(e_sq/N-(e/N)**2)*Beta/T, ]
        time+=[k, ]
        for i in range(L):
            for j in range(L):
                deltaE=2*J*phi[i,j]*(phi[p[i],j]+phi[m[i],j]+phi[i,p[j]]+phi[i,m[j]])+2*H*phi[i,j]
                if deltaE>0:
                    omega=math.exp(-Beta*deltaE)
                    r=random.uniform(0, 1)
                    if r<=omega:
                        phi[i,j]=-phi[i,j]
                else:
                    phi[i,j]=-phi[i,j]


    Tm+=[T*Kb, ]
    Energy+=[sum(E)/iterations, ]
    Magnetization+=[abs(sum(M))/iterations, ]
    Susceptibility+=[sum(chi)/iterations, ]
    Specific_heat+=[sum(Cv)/iterations, ]

plt.plot(Tm, Energy, 'k*', ms=4)
plt.ylabel('<E>')
plt.xlabel('Kb*T')
plt.show()

plt.plot(Tm, Magnetization, 'k*', ms=4)
plt.ylabel('<|M|>')
plt.xlabel('Kb*T')
plt.show()

plt.plot(Tm, Susceptibility, 'k*-', ms=4)
plt.ylabel('$<\chi>$')
plt.xlabel('Kb*T')
plt.show()

plt.plot(Tm, Specific_heat, 'k*-', ms=4)
plt.ylabel('$<C_v>$')
plt.xlabel('Kb*T')
plt.show()
    
