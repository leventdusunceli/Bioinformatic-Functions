
#First, let's define 2 functions for transcription and translation 

def transcriber(string): 
    
    string = string.replace("T","U")
    
    return string

def RNA_Translator(RNA_seq):

    #First we need to split the RNA sequnce into codons 
    aa_dict = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L","UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S", "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
    "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W", "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L", "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P","CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R","AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M","ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K","AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R","GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A","GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E", "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
    
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

#Now, let's define the splicer function.  Functions defined above will be used in this function. 

def rna_splicer (filename):
  
  #We'll use BioPython to parse the given fasta file

  !install Bio 
  from Bio import SeqIO
  
  #parsing each fasta formatted sequence with SeqIO and appending just the sequences into a list 

  seq_list =[]  
  for seq_record in SeqIO.parse(filename,"fasta"): 
    seq_list.append(str(seq_record.seq))

  #Full length DNA will be in index 0 of seq_list because:
  #seqIO.parse parses each fasta in the document in given order 
  #since the first fasta in our txt file is full length DNA 
  #It'll be the first element to be appended to seq_list, making it index 0 
  #Sequences in the rest of the indexes are introns 

  #now we'll remove introns from the full length DNA, aka perform splicing 

  for i in range(len(seq_list)):

    if i == 0 : 
      continue 

    else: 
      seq_list[0] = seq_list[0].replace(seq_list[i],"") 
      
  mature_DNA = seq_list[0]  

  #We'll perform transcription and translation with functions defined previously 
  
  rna_string = transcriber(mature_DNA)
  protein_sequence = RNA_Translator(rna_string)

  return protein_sequence
