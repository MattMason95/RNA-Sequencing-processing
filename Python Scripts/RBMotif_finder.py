import pandas as pd
import numpy as np

class RBMotif_finder:
    '''
        RBmotif_finder is a function for analysing submitted fasta sequences for the presence of RNA-binding factor motif sites.
        
    '''
    def __init__(self, fasta_path, threshold=0.8, measurement='hamming'):
        self.fasta_path = fasta_path
        self.threshold = threshold
        self.measurement = measurement   

        valid_measurements = ['hamming']
        if self.measurement not in valid_measurements:
            raise Exception(f"Invalid parameter. {measurement} was passed to the measurement parameter. Only {valid_measurements} are permitted.")
        
        if isinstance(self.threshold, float):
            pass
        else:
            raise Exception(f"Invalid data type. Only integers can be passed to threshold parameter.")
            
    
    def fasta_parser(self):
        data = pd.read_csv(self.fasta_path,sep='\t',header=None)
        return data
        

    def do_something(self):
        parsed_data = self.fasta_parser()
        print(parsed_data)
    # For each K-mer, scan through the fasta sequence with a read window of length K. How will this be calculated? 
        # e.g. if sequence is 5 and k-mer is 3, there are 3 overlapping read windows; for a sequence length of 6, there are 4. What formula dictates this relationship?
        #   
        


        # Need to compute Hamming distance between strings - probably best to write a separate function for this

