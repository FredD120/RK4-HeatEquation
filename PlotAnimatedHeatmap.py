#!/usr/bin/env python
# coding: utf-8

# In[26]:


#Make animated heatmap gif from files
get_ipython().run_line_magic('matplotlib', 'notebook')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import datetime
from scipy import sparse
from scipy.sparse.linalg import spsolve
import seaborn as sns
import os

starttime = datetime.datetime.now()

#File name to save off data
FileName = "ModelTest2D.txt"

#Load data on each frame and animate
data = np.loadtxt(str(0)+FileName)
def animate(i): 
    print(i)
    data = np.loadtxt(str(i*4)+FileName)
    hm = sns.heatmap(data, vmax=.8, square=True, cbar=False)
    

#Setup heatmap and animation
fig = plt.figure()
anim = animation.FuncAnimation(fig, animate,  frames =200, interval=1)
hm = sns.heatmap(data, vmax=.8, square=True, cbar=False)

#Save animation as gif
savefile = r"Large2DTemp.gif"
pillowwriter = animation.PillowWriter(fps=40)
anim.save(savefile, writer=pillowwriter)


plt.show()
print(datetime.datetime.now()-starttime)


# In[ ]:


import os
for i in range(0,999):
    if os.path.exists(str(i)+"OpenGaussianData.txt"):
        os.remove(str(i)+"OpenGaussianData.txt")
    else:
        print("file not found")


# In[24]:


import matplotlib.pyplot as plt
import seaborn as sns

FileName = "ModelTest2D.txt"
fig = plt.figure()
data = np.loadtxt(str(160*4)+FileName)
hm = sns.heatmap(data, vmax=.8, square=True, cbar=False)


# In[ ]:




