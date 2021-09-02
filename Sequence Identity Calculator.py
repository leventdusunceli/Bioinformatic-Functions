#!/usr/bin/env python
# coding: utf-8

# In[4]:


def Percent_Identity (seq1,seq2):
    
    same_aas = 0 
    
    for i in range(len(seq1)): 
        
        if seq1[i] == seq2[i]:
            
            same_aas += 1
            
    return (same_aas*100)/len(seq1) 


# 

# In[3]:





# In[ ]:




