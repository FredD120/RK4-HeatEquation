#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('matplotlib', 'notebook')
import numpy as np
from scipy import integrate 
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def animate(i):
    MetalRod.set_ydata(Temperature[i,:])
    return MetalRod

def Euler(f,u,dt,dx):
    return u + f(u)*dt/(dx**2)

def f(u):
    v = np.zeros(len(u))
    for i in range(1,len(u)-1):
        v[i] = u[i+1] -2*u[i] + u[i-1]
    return v
    

TotalTime = 1
Steps = 50000
TimeStep = TotalTime/Steps
Length = 1
Granularity = 100
LengthStep = Length/Granularity
Temperature = np.zeros((Steps,Granularity))
#InitialTemp = np.arange(Length,0,-LengthStep)
Position = np.arange(0,Length,LengthStep)
Temperature[:,0] = 10
Temperature[:,-1] = 10

for i in range(1,Steps):
    Temperature[i,:] = Euler(f,Temperature[i-1,:],TimeStep,LengthStep)
    

f0,ax = plt.subplots(figsize = (10,6))
MetalRod = plt.plot(Position,Temperature[0,:],color='red')[0]
ani = animation.FuncAnimation(fig=f0, func=animate, frames=Steps, interval=2)
ax.set(xlim=[0, Length], ylim=[0, 10], xlabel='Position', ylabel='Temperature')
plt.show()
        


# In[ ]:




