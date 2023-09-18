#!/usr/bin/env python
# coding: utf-8

# In[17]:


get_ipython().run_line_magic('matplotlib', 'notebook')
import numpy as np
from scipy import integrate 
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def RK4(f,t,y,h):
    K1 = f(t,y)
    K2 = f(t+h/2,y+h*K1/2)
    K3 = f(t+h/2,y+h*K2/2)
    K4 = f(t+h,y+h*K3)
    return y+(h/6)*(K1+2*K2+2*K3+K4)

def dx(t,x):
    return 5*t


TotalTime = 5
Steps = 1000
TimeStep = TotalTime/Steps
Position = np.zeros(Steps)
Time = np.arange(0,TotalTime,TimeStep)
xInitial = 1
Position[0] = xInitial


for i in range(1,Steps):
    Position[i] = RK4(dx,Time[i],Position[i-1],TimeStep)
    
f0 = plt.plot(Time,Position)
plt.xlabel('Time')
plt.ylabel('Position')
plt.show()


# In[ ]:




