#!/usr/bin/env python
# coding: utf-8

# ### First we need to define aminoacids corresponding to different codons 

# In[1]:


aa_dict = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
    "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}


# ### Now onto function translating RNA to protein sequence

# In[2]:


def RNA_Translator(RNA_seq):

    #First we need to split the RNA sequnce into codons 
    
    codons = []
    n = 3 

    for i in range(0,len(RNA_seq),n):
        codons.append(RNA_seq[i:i+n])
    
    #This point on the translation process commences
    
    protein_sequence = ""
    
    for i in range(len(codons)):
        
        if aa_dict[codons[i]] == "STOP":
            continue
        
        elif codons[i] in aa_dict:
            protein_sequence +=(aa_dict[codons[i]])
            
    return(protein_sequence)


## Alternatively, you can define 2 functions, one for codon separation and one for translation.  
## This method would be particularly useful if you aim to do other manupulations with just the codons of the sequence. 

# In[3]:


def codon_separator (RNA_seq):
    
    codons = []
    n = 3 

    for i in range(0,len(RNA_seq),n):
        codons.append(RNA_seq[i:i+n])

    return codons 


# In[4]:


def RNA_Translator2 (RNA_seq):
    
    codons2 = codon_separator(RNA_seq)
    
    protein_sequence = ""
    
    for i in range(len(codons2)):
        
        if aa_dict[codons2[i]] == "STOP":
            continue
        
        elif codons2[i] in aa_dict:
            protein_sequence +=(aa_dict[codons2[i]])
            
    return(protein_sequence)
