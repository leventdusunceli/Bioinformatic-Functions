#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#Counts the number of each nucleotide in a given DNA sequence


def nucleotide_counter(string):
    
    d={}
    
    liste = ["A","C","T","G"]

    
    for i in range(len(liste)):
        d[liste[i]] = string.count(liste[i])
    
    print (d)

