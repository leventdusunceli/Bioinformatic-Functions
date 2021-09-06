#!/usr/bin/env python
# coding: utf-8

# ## Finding the longest shared sequence  

# 

# In[ ]:


def shared_motif_finder(fasta_file):
    
    #In order to use data in fasta file, we need SeqIO package from Biopython
    
    from Bio import SeqIO
    
    fasta_seqs = list(SeqIO.parse(fasta_file,'fasta'))  #This will insert each entry in the fasta file into a list
    
    #Following lines are for being able to insert sequences into loops 
    sequences = []

    for i in range(len(fasta_seqs)):
        sequences.append(fasta_seqs[i].seq)
    
    
    sorted_sequences = sorted(sequences,key=len)  
    shortest_seq =sorted_sequences[0]  
    other_seqs = sorted_sequences[1:]
    motif= ''
    
    
    for i in range(len(shortest_seq)):        #starting from the very first residue of the shortest sequence
        for j in range(i,len(shortest_seq)):  #check if there's a match
            m = shortest_seq[i:j+1]
            found = False
            for b in other_seqs:
                if m in b:
                    found = True
                else:
                    found = False 
                    break

            if found and len(m)>len(motif):
                motif = m
    
    return(motif)          
                

