#!/usr/bin/env python
# coding: utf-8

# In[5]:


def motif_finder(string,substring):
    
    print(substring,"motif found between nucleotides:")
    
    index = 0 
    while index < len(string):
        
        index = string.find(substring,index)
        
        if index == -1: 
            break
        
        print  (index + 1 , "and", index + 1 + len(substring))    #(+1) because we want to know the actual start position of the substring, not the index 
        index +=1


# In[ ]:




