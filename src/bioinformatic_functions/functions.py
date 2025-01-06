import pandas as pd
from Bio import SeqIO
from urllib.request import urlopen
import re


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

    def GC_content(self, string):
        """Function for calculating GC content of a DNA string

        Args:
            string (str): Input DNA string

        Returns:

        """
        total_GC = 0

        for char in string:

            if char == "G" or char == "C":
                total_GC += 1
        gc_content = total_GC * 100 / len(string)

        return f"GC Percentage Content: {gc_content}"

    def motif_finder(self, string, substring):
        """Finding a motif in a longer genomic string

        Args:
            string (str): Query string
            substring (str): Substring being searched

        Returns:
            List: locations where the substring was found
        """

        print(substring, "motif found between nucleotides:")

        substr_locations = []

        index = 0
        while index < len(string):

            index = string.find(substring, index)

            if index == -1:
                break
            start_pos = index + 1
            end_pos = (
                index + 1 + len(substring)
            )  # (+1) because we want to know the actual start position of the substring, not the index
            loc = f"{start_pos}:{end_pos}"
            substr_locations.append(loc)

            index += 1

        return substr_locations

    def shared_motif_finder(self, fasta_file):
        """Function for finding the longest shared motif between strings within a fasta file

        Args:
            fasta_file (fasta): Fasta file containing multiple sequences

        Returns:
            str: the longest shared motif sequence
        """

        # In order to use data in fasta file, we need SeqIO package from Biopython

        fasta_seqs = list(
            SeqIO.parse(fasta_file, "fasta")
        )  # This will insert each entry in the fasta file into a list

        # Following lines are for being able to insert sequences into loops
        sequences = []

        for i in range(len(fasta_seqs)):
            sequences.append(fasta_seqs[i].seq)

        sorted_sequences = sorted(sequences, key=len)
        shortest_seq = sorted_sequences[0]
        other_seqs = sorted_sequences[1:]
        motif = ""

        for i in range(
            len(shortest_seq)
        ):  # starting from the very first residue of the shortest sequence
            for j in range(i, len(shortest_seq)):  # check if there's a match
                m = shortest_seq[i : j + 1]
                found = False
                for b in other_seqs:
                    if m in b:
                        found = True
                    else:
                        found = False
                        break

                if found and len(m) > len(motif):
                    motif = m

        return motif

    def motif_search(self, access_ID):
        """for finding N-glycosylation motif positions in a protein sequence if the UniProt access IDs are given.

        Args:
            access_ID (str): UniProt ID of protein

        Returns:
            list: positions where N-glycosylation motifs can be found
        """

        # pulling fasta file of given access ID
        URL = "http://www.uniprot.org/uniprot/" + access_ID + ".fasta"
        data = urlopen(URL)
        fasta = data.read().decode("utf-8")

        # writing fasta file into a txt file
        file_name = "seq_file1.txt"
        with open("seq_file1.txt", "w") as txt_file:
            txt_file.write(fasta)

        # parsing the protein sequence from fasta format
        for seq_record in SeqIO.parse(file_name, "fasta"):
            protein_seq = str(seq_record.seq)

        # searching for n-glycosylation motif
        # one can use this function for searching other motifs by changing the
        # regular expression in the following lines
        motifs = []

        motifs.append(re.findall("(?=(N[^P][ST][^P]))", protein_seq))

        # we added ?= infront of our regular expression denoting the n-glycosylation
        # motif in order to cover overlapping sequences
        # ?= is a lookahead assertion

        # re.findall returns found motifs as a list of strings inside a list
        # i.e [['NFSD', 'NSSN', 'NWTE', 'NLSK', 'NISA']]
        # following for loop adds each motif as a single element into a list
        # i.e ["NFSD","NSSN","NWTE","NLSK","NISA"]

        motifs2 = []

        for i in range(len(motifs)):
            for j in range(len(motifs[i])):
                motifs2.append(motifs[i][j])

        # following for loops return the positions of motifs inside the protein seq

        indexes = []
        positions = []

        for i in motifs2:
            indexes.append(protein_seq.index(i))

        for i in indexes:
            positions.append(i + 1)

        return access_ID, positions
