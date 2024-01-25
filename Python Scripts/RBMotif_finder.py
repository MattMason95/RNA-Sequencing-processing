import pandas as pd
import numpy as np

class RBMotif_finder:
    '''
        RBmotif_finder is a function for analysing submitted fasta sequences for the presence of RNA-binding factor motif sites.
        
    '''
    def __init__(self, fasta_path, motif_path, threshold=0.8, measurement='hamming'):
        self.fasta_path = fasta_path
        self.threshold = threshold
        self.motif_path = motif_path
        self.measurement = measurement   

        valid_measurements = ['hamming']
        if self.measurement not in valid_measurements:
            raise Exception(f"Invalid parameter. {measurement} was passed to the measurement parameter. Only {valid_measurements} are permitted.")
        
        if isinstance(self.threshold, float):
            pass
        else:
            raise Exception(f"Invalid data type. Only integers can be passed to threshold parameter.")       
    # ~~~~~~~~~~~~~~~~~~~~~~
    def load_fasta(self):
        fasta = pd.read_csv(self.fasta_path,sep='\t',header=None)
        labels = fasta.iloc[::2]
        sequences = fasta.iloc[1::2]

        if len(labels) != len(sequences):
            raise Exception('Mismatch between number of extracted labels and extracted sequences.')
        
        return labels, fasta
    # ~~~~~~~~~~~~~~~~~~~~~~
    def load_motifs(self):
        motifs = pd.read_csv(self.motif_path,header=None)
        motif_labels = motifs.iloc[:,0]
        motif_seqs = motifs.iloc[:,1]
        
        if len(motif_labels) != len(motif_seqs):
            raise Exception('Mismatch between number of extracted labels and extracted sequences.')
        return motif_label, motif_seq

    def hamming_distance(self):
    
    def compute_distances(self):
        labels,sequences = self.load_fasta()
        motif_labels,motif_seqs = self.load_motifs()

        print(f'Loaded {len(labels)} FASTA sequences.')
        print(f'Scanning with {len(motifs_labels)} RNA binding motifs.')
        print('---------------------------------------------------------------')
        print(f'Calculating distance with {self.measurement}')
        
        if self.measurement == 'hamming':
            for label, sequence in zip(labels,sequences):
                for motif_label, motif_seq in zip(motif_labels, motif_seqs):
                    
                    motif_length = len(motif_seq
            
            
    
    



    # For each K-mer, scan through the fasta sequence with a read window of length K. How will this be calculated? 
        # e.g. if sequence is 5 and k-mer is 3, there are 3 overlapping read windows; for a sequence length of 6, there are 4. What formula dictates this relationship?
        #   
        


        # Need to compute Hamming distance between strings - probably best to write a separate function for this

