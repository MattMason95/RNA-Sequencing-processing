import pandas as pd
import numpy as np
import itertools
from itertools import islice

class RBMotif_finder:
    '''
        RBmotif_finder is a function for analysing submitted fasta sequences for the presence of RNA-binding factor motif sites.
        
    '''
    def __init__(self, fasta_path, motif_path, background_path=None, threshold=0.8, measurement='hamming'):
        self.fasta_path = fasta_path
        self.motif_path = motif_path
        self.background_path = background_path
        self.threshold = threshold
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
        ## Extract FASTA data from the path provided
        fasta = pd.read_csv(self.fasta_path,sep='\t',header=None)
        labels = fasta.iloc[::2][0].to_list()
        sequences = fasta.iloc[1::2][0].to_list()

        if len(labels) != len(sequences):
            raise Exception('Mismatch between number of extracted FASTA labels and sequences.')
        
        return labels, sequences
    # ~~~~~~~~~~~~~~~~~~~~~~
    def load_motifs(self):
        ## Extract motif data from the path provided
        motifs = pd.read_csv(self.motif_path,header=None)
        motif_labels = motifs.iloc[:,0]
        motif_seqs = motifs.iloc[:,1]
        
        if len(motif_labels) != len(motif_seqs):
            raise Exception('Mismatch between number of extracted motif labels and sequences.')
        return motif_labels, motif_seqs

    def hamming_distance(self,chain1,chain2):
        self.chain1 = chain1
        self.chain2 = chain2

        return sum(c1 != c2 for c1, c2 in zip(chain1,chain2))
    
    def window(self,sequence,n):
        '''
        Taken from Itertools documentation: 
        Returns a sliding window (of width n) over data from the iterable.
        s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...
        '''
        it = iter(sequence)
        result = tuple(islice(it, n))
        if len(result) == n:
            yield result    
        for elem in it:
            result = result[1:] + (elem,)
            yield result

    def background_model(self):
        '''
        In order to evaluate the statistical significance of a detected consensus motif, the consensus scores must be compared against a background model.
        The background model represents a selection of randomly chosen regulatory regions.
        Z-scores
        '''


    def compute_distances(self):
        labels,sequences = self.load_fasta()
        motif_labels,motif_seqs = self.load_motifs()

        print(f'Loaded {len(labels)} FASTA sequences.')
        print(f'Scanning with {len(motif_labels)} RNA binding motifs.')
        print('---------------------------------------------------------------')
        print(f'Calculating distance with {self.measurement} metric.')
        print(sequences)

        if self.measurement == 'hamming':
            ## Utilising the Hamming distance metric for comparing FASTA and motif sequences
            ## Prepare output container
            
            output = pd.DataFrame(columns=['seq_label','motif','hamming_distances','hamming_max'])

            for fasta_label, fasta_seq in zip(labels,sequences):
                ## Prepare output containers
                

                for motif_label, motif_seq in zip(motif_labels, motif_seqs):
                    motif_length = len(motif_seq)
                    sequence_slices = [(motif_length-self.hamming_distance(x,motif_seq))/motif_length for x in self.window(fasta_seq, motif_length)]
                    
                    # hamming_dist = self.hamming_distance(sequence_slices,motif_seq)
                    output.loc[f'{fasta_label}_{motif_label}'] = [motif_label,motif_seq,sequence_slices,max(sequence_slices)]
                    print(motif_label,motif_seq,sequence_slices)
        
               
        return output
                    
    ## To do: 
        ## Build up background model for statistical validation
        ## Build up the weighted_rank function 
        ## Try with some real FASTA sequences and real motifs 
        
    ##def weighted_rank(self):
        ''' 
        The weighted rank function accounts for the propensity of suboptimal motifs to cluster around the significant motif.
        
        '''

    
