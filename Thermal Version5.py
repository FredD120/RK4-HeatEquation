#!/usr/bin/env python
# coding: utf-8

# In[5]:


# Visualising the data in a nicer way using heatmaps, and setup a way to save animations as gifs
get_ipython().run_line_magic('matplotlib', 'notebook')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import datetime
from scipy import sparse
from scipy.sparse.linalg import spsolve
import seaborn as sns

starttime = datetime.datetime.now()

def animate(i):
    """ 
    data = np.ones(3,len(Temperature[i,:]))
    for j in range(0,10):
        data[j,:] = Temperature[i,:]
    sns.heatmap(data, vmax=.8, square=True, cbar=False)
    """
    sns.heatmap(((Temperature[i*10,:],Temperature[i*10,:])), vmax=.8, square=True, cbar=False)

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
TotalTime = 0.6
Steps = 2000
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

stdev = 0.1
mean = Position[int(np.floor(Granularity/2))]
Temperature[0,:] = (1/stdev**0.5)*np.exp(-(Position-mean)**2/(2*stdev**2))
"""
Temperature[0,int(np.floor(Granularity/3)):int(np.floor(2*Granularity/3))] = 10



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


#Time step and use solver on each step (Backward euler)
for i in range(1,Steps):
    Temperature[i,:] = BackwardEuler(Temperature[i-1,:],A)
  

#Plot heatmap over time
f0 = plt.figure(figsize = (10,6))
sns.heatmap(Temperature, vmax=1, square=False, cbar=False)


#Plot animated graph
f0 = plt.figure(figsize = (10,6))
ani = animation.FuncAnimation(fig=f0, func=animate, frames=int(Steps/10), interval = 2)


#Save animation as a gif
savefile = r"ThermalBathTemperatureHeatmap.gif"
pillowwriter = animation.PillowWriter(fps=50)
ani.save(savefile, writer=pillowwriter)


plt.show()

#Approximates error based on energy loss of system (closed system boundary conditions only)
#print((np.sum(Temperature[0,:])-np.sum(Temperature[-1,:]))/(np.sum(Temperature[0,:])))

#Outputs time taken to solve and plot
print(datetime.datetime.now()-starttime)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




