#!/usr/bin/env python
# coding: utf-8

# In[2]:


def transcriber(string): 
    
    string = string.replace("T","U")
    
    return string


# In[1]:


def complementer(DNA_strand):
    
    aa_dict = {"G":"C", "C":"G","T":"A","A":"T"}
    complementary_strand = ""
    
    for char in DNA_strand: 
        complementary_strand += aa_dict[char]
        
    return complementary_strand


# In[3]:


def reverse_complementer(DNA):
    
    reversed_strand= DNA[::-1]    #this line reverses the given DNA Sequence
    
    reversed_complement = complementer(reversed_strand)
    
    print(reversed_complement)


# In[ ]:




