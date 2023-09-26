#!/usr/bin/env python
# coding: utf-8

# In[11]:


# Implemented reverse Euler technique, requiring linear algebra to solve a set of linear equations on every timestep
get_ipython().run_line_magic('matplotlib', 'notebook')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import datetime
from scipy import sparse
from scipy.sparse.linalg import spsolve

starttime = datetime.datetime.now()

def animate(i):
    MetalRod.set_ydata(Temperature[i,:])
    return MetalRod

#Implement PDE solver
def BackwardEuler(u_prev,A):
    #Using linear algebra to solve simultaneous equations to find next temperature values(Central difference scheme)
    u = spsolve(A,u_prev)
    
    #Boundary conditions (Heat bath by default/Closed system below)
    """
    u[0] = u[1]
    u[-1]=u[-2]
    """
    return u

#Define space and time
TotalTime = 0.5
Steps = 1000
TimeStep = TotalTime/Steps
Length = 1
Granularity = 500
LengthStep = Length/Granularity
Position = np.arange(0,Length,LengthStep)
Temperature = np.zeros((Steps,Granularity))
delta = TimeStep/(LengthStep**2) #Mesh Fourier number

#Inital temperature conditions (Heat bath/Closed system/Temperature discontinuity)#
"""
Temperature[:,0] = 10
Temperature[:,-1] = 10
"""
stdev = 0.1
mean = Position[int(np.floor(Granularity/2))]
Temperature[0,:] = (10/stdev**0.5)*np.exp(-(Position-mean)**2/(2*stdev**2))
"""
Temperature[0,int(np.floor(Granularity/3)):int(np.floor(2*Granularity/3))] = 10
"""


#Precalculate sparse matrix for solving linear equations 
#1 = lower, 2=middle, 3=upper
diag1 = np.zeros(Granularity-1)
diag2 = np.zeros(Granularity)
diag3 = np.zeros(Granularity-1)

diag1[:-1] = -delta
diag2[0] = diag2[-1] = 1
diag2[1:-1] = 1+2*delta
diag3[1:] = -delta

A = sparse.diags(diagonals = [diag1,diag2,diag3], offsets = [-1,0,1],format="csr")

np.savetxt(str(0)+"OpenGaussianData.txt",Temperature[0,:])
#Time step and use solver on each step (Backward euler)
for i in range(1,Steps):
    Temperature[i,:] = BackwardEuler(Temperature[i-1,:],A)
    np.savetxt(str(i)+"OpenGaussianData.txt",Temperature[i,:])
    
"""
#Plot animated graph
f0,ax = plt.subplots(figsize = (10,6))
MetalRod = plt.plot(Position,Temperature[0,:],color='red')[0]
ani = animation.FuncAnimation(fig=f0, func=animate, frames=int(Steps), interval=2)
ax.set(xlim=[0, Length], ylim=[0, 10], xlabel='Position', ylabel='Temperature')
plt.show()

#Save as gif

savefile = r"GaussianTemperature.gif"
pillowwriter = animation.PillowWriter(fps=20)
ani.save(savefile, writer=pillowwriter)
"""




#Approximates error based on energy loss of system (closed system boundary conditions only)
print(np.sum(Temperature[0,:])/Granularity-np.sum(Temperature[-1,:])/Granularity)
#Outputs time taken to solve and plot
print(datetime.datetime.now()-starttime)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




