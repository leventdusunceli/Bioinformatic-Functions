#!/usr/bin/env python
# coding: utf-8

# ## Function for a DNA sequence as string

# In[4]:


def GC_Content(st):    
    total_GC = 0 
    
    for char in st:
        
        if char == "G" or char == "C":
            total_GC+=1
    print ("Length of DNA Seq:" , len(st))  
    print ("Total GC Number:",total_GC)
    print("GC Content Percentage:", total_GC*100/len(st))     


# ## If DNA sequences are contained in FASTA file

# In[3]:


from Bio import SeqIO

for seq_record in SeqIO.parse("fastafile.txt", "fasta"):
    a =  seq_record.id +  "-->" + str(GC_Content(seq_record.seq))   #we used GC_content function here
    print( "\n", a)


# In[ ]:




