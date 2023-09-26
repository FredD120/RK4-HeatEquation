#!/usr/bin/env python
# coding: utf-8

# In[6]:


#2D Temperature calculation
get_ipython().run_line_magic('matplotlib', 'notebook')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import datetime
from scipy import sparse
from scipy.sparse.linalg import spsolve
import seaborn as sns

#Measure time taken
starttime = datetime.datetime.now()

#File name to save off data
FileName = "ModelTest2D.txt"

#Setup number of time and space steps
Time = 0.3
Steps = 1001
SideLength = 1
Granularity = 501
TimeStep = Time/Steps
LengthStep = SideLength/Granularity
delta = TimeStep/(LengthStep**2) #Mesh Fourier number
Nx = Ny = Granularity
N = Nx*Ny

#Set up matrix to solve linear alebra
diag1 = np.zeros(N-Nx)
diag2 = np.zeros(N-1)
diag3 = np.ones(N)
diag4 = np.zeros(N-1)
diag5 = np.zeros(N-Nx)

for i in range(Nx+1,N-Nx-1):
    if 0<i%Nx<Nx-1:
        diag3[i] = 1+4*delta
        diag2[i-1] = -delta
        diag4[i] = -delta
        diag1[i-Nx] = -delta
        diag5[i] = -delta

A = sparse.diags(diagonals = [diag1,diag2,diag3,diag4,diag5], offsets = [-Nx,-1,0,1,Nx],format="csr")


#Set up initial temperature state
CurrentTemp = np.zeros((Nx,Ny))

for i in range(Nx//2-50,Nx//2+51):
    CurrentTemp[Nx//2-50:Nx//2+51,i] = 25


#Backward Euler solver, called every time step
def BackwardEuler(u_prev,A):
    u = np.zeros((Nx,Ny))
    u.flat = spsolve(A,u_prev.flat) 
    return u

for i in range(0,Steps):
    #Save every n(=4) frames
    if i%4 == 0:
        np.savetxt(str(i)+FileName,CurrentTemp)
    CurrentTemp = BackwardEuler(CurrentTemp,A)
print(datetime.datetime.now()-starttime)


# In[ ]:




