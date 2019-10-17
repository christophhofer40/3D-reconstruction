
# coding: utf-8

# In[1]:


from maptools import autotune
#import matplotlib.pyplot as plt
#import tifffile
import numpy as np
import sys
from PIL import Image

# In[4]:

args=sys.argv
print(args)
im=autotune.Imaging()
arglist=np.zeros(8) 
arglist[:len(args)]=args

# In[25]:


ar=im.image_grabber(imsize=6,frame_parameters={'pixeltime':-1,'fov':args[0],'size_pixels':(256,256)},
                    aberrations={'EHTFocus': args[1], 'C12_a': args[2], 'C12_b': args[3], 'C21_a': args[4],
                                 'C21_b': args[5], 'C23_a': args[6], 'C23_b': args[7]}
,relative_aberrations=False)

# In[26]:


#get_ipython().run_line_magic('matplotlib', 'auto')
#plt.matshow(ar)
#plt.show()

