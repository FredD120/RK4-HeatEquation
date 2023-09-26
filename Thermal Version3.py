#!/usr/bin/env python
# coding: utf-8

# In[2]:


# Added proper descriptors, vectorised code to avoid for loops
get_ipython().run_line_magic('matplotlib', 'notebook')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import datetime

starttime = datetime.datetime.now()

def animate(i):
    MetalRod.set_ydata(Temperature[i,:])
    return MetalRod

#Implement PDE solver
def ForwardEuler(u,d):
    #Heat flux is difference between cell temp and adjacent cell temps (Central difference scheme)
    u[1:-1] = u[1:-1] + (u[2:] -2*u[1:-1] + u[:-2])*d
    #Boundary conditions (Heat bath/Closed system)
    """
    u[0] = 10
    u[-1] = 10
    """
    u[0] = u[1]
    u[-1]=u[-2]
    
    return u

#Define space and time
TotalTime = 1
Steps = 50000
TimeStep = TotalTime/Steps
Length = 1
Granularity = 150
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
Temperature[0,:] = (1/stdev**0.5)*np.exp(-(Position-mean)**2/(2*stdev**2))
"""
Temperature[0,int(np.floor(Granularity/3)):int(np.floor(2*Granularity/3))] = 10
"""

#Step through time (Forward euler)
for i in range(1,Steps):
    Temperature[i,:] = ForwardEuler(Temperature[i-1,:],delta)
  
#Plot animated graph
f0,ax = plt.subplots(figsize = (10,6))
MetalRod = plt.plot(Position,Temperature[0,:],color='red')[0]
ani = animation.FuncAnimation(fig=f0, func=animate, frames=Steps, interval=2)
ax.set(xlim=[0, Length], ylim=[0, 10], xlabel='Position', ylabel='Temperature')
plt.show()

print(np.sum(Temperature[0,:])/Granularity-np.sum(Temperature[-1,:])/Granularity)
print(datetime.datetime.now()-starttime)


# In[8]:


f0,ax = plt.subplots(figsize = (10,6))
plt.plot(Position,np.random.normal(1,0.1,Granularity),color='red')
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:




