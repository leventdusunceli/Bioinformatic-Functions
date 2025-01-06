import pandas as pd


class BioinformaticFunctions:

    def return_hello(self):
        return "hello"

    def transcriber(self, string):
        """Changes a DNA string to an RNA strinng
        Args:
            string (str): Input DNA string

        Returns:
            str: RNA string
        """
        transcribed_string = string.replace("T", "U")

        return transcribed_string

    def complementer(self, string):
        """Creates the complementary string of a given DNA string

        Args:
            string (str): Input DNA string

        Returns:
            str: Complementary DNA string
        """
        aa_dict = {"G": "C", "C": "G", "T": "A", "A": "T"}
        complementary_strand = ""

        for char in string:
            complementary_strand += aa_dict[char]

        return complementary_strand

    def reverse_complementer(self, string):
        """Creates a reverse complement of a givenn DNA string

        Args:
            string (str): Input DNA string

        Returns:
            _type_: _description_
        """

        reversed_strand = string[::-1]  # this line reverses the given DNA Sequence

        reversed_complement = self.complementer(reversed_strand)

        return reversed_complement


    def GC_content(self,string):    
        """Function for calculating GC content of a DNA string 

        Args:
            string (str): Input DNA string 

        Returns: 

        """
        total_GC = 0 
        
        for char in string:
            
            if char == "G" or char == "C":
                total_GC+=1
        gc_content = total_GC*100/len(string)

        return f"GC Percentage Content: {gc_content}"
    
    def motif_finder(self,string,substring):
    
    print(substring,"motif found between nucleotides:")
    
    index = 0 
    while index < len(string):
        
        index = string.find(substring,index)
        
        if index == -1: 
            break
        
        print  (index + 1 , "and", index + 1 + len(substring))    #(+1) because we want to know the actual start position of the substring, not the index 
        index +=1