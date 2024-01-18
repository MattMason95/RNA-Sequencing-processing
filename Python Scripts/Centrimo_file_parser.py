def data_parser(textfile,test_split):
    '''
    The objective of this function is to iterate through a slice of a text file that summarises motif localisations.
    There are three data points that need to be extracted:
    - Genomic coordinate {Chromosome:rangeStart-rangeEnd}
    - P-value {float}
    - Motif information [Repeated for each motif {indexStart_[strand+ID(p_value)_indexEnd]}
    
    To achieve this, we will iterate through each line of the summary file and use Regex patterns to extract the information into arrays for reassembly, summary statistics, and plotting. 
    '''
    ## Library imports
    import regex as re
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Pre-Processing
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ## Find the line corresponding to the Combined Block Diagram summary
    offset = []
    pattern='Combined block diagrams:'

    for idx, line in enumerate(textfile):
        if pattern in line:
            offset = idx
    
    ## Truncate input file to summary - i.e. data after the block diagram summary
    summary = textfile[offset+2:]
    
    ## For testing purposes, test_split can be used to truncate input file size
    if test_split:
        summary = summary[:test_split]
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    ## Instantiate all regex patterns
    coordinate_pattern = re.compile(r"([0-9]+)+:([0-9]+)-([0-9]+)", re.IGNORECASE)
    stat_pattern = re.compile(r"(\d*)_\[(\S\d*)\((\d*\.\d*e-\d*)\)")
    
    ## Prepare coordinates output container
    outputs = []
    
    ## Prepare temporary containers for retaining indices
    temp_idx = None
    temp_chromosome = None
    temp_start = None
    temp_end = None
    
    for idx,line in enumerate(summary):
        ## First check the line for coordinate patterns - these are essentially indices, because of the text wrapping that has occured in this text file
        if coordinate_pattern.match(line):
            ## For each line search for and extract matches to the coordinate pattern
            coord_search = coordinate_pattern.search(line)
            chromosome = coord_search.group(1).strip()
            start = coord_search.group(2).strip()
            end = coord_search.group(3).strip()
            
            ## Store coordinates in both the output container and the temporary containers
            # coordinates.append([chromosome,start,end])
            temp_idx = idx
            temp_chromosome = chromosome
            temp_start = start
            temp_end = end
            
            ## Search the same line for and extract matches to the stats pattern
            ## As some lines will contain information for multiple motifs, using regex findall and iterating through the pattern matches will yield all results
            stat_search = stat_pattern.findall(line)
            
            ## Iterate through the matches from find all and append a new line for each (with the same genomic index)
            for match in stat_search:
                motif_start = match[0].strip()
                motif_id = match[1].strip()
                motif_stats = match[2].strip()
                try:
                    motif_start = int(motif_start)
                    outputs.append([f'seq_{idx}',chromosome,start,end,motif_start,motif_id,motif_stats])
                except:
                    print('Skipping empty array')
                    continue 
              
        ## If the line is missing a coordinate pattern, check whether stats information is present - if so, the information on the line is attributed to the previous coordinate index
        elif stat_pattern.findall(line): #'stat_pattern.match(line):
            ## Check for stats information
            stat_search = stat_pattern.findall(line)
            
            ## Iterate through the matches from find all and append a new line for each (with the same genomic index)
            for match in stat_search:
                motif_start = match[0].strip()
                motif_id = match[1].strip()
                motif_stats = match[2].strip()
                try: 
                    motif_start = int(motif_start)
                    outputs.append([f'seq_{temp_idx}',temp_chromosome,temp_start,temp_end,motif_start,motif_id,motif_stats])
                except:
                    print('Skipping empty array')
                    continue
                    
        else:
            ## In the event that the line continues no coordinate or stat information, skip line
            continue
  
    ## Prepare output dataframe
    df = pd.DataFrame(outputs,columns=['seq_id','chromosome','sequence_start','sequence_end','motif_start','motif_id','motif_pvalue'])
    df = df.astype({'seq_id':'object','chromosome': 'object', 'sequence_start': 'int', 'sequence_end': 'int','motif_start':'int','motif_id':'str','motif_pvalue':'str'})
    
    ## Add additional columns/edit existing
    df['sequence_length'] = df['sequence_end'] - df['sequence_start']
    df['strand'] = np.where(df['motif_id'].str[0] == '+', '+','-')
    df['motif_id'] = df['motif_id'].str[1]
    
    ## Changing order for sake of appearance
    df = df[['seq_id','chromosome','sequence_start','sequence_end','sequence_length','motif_start','motif_id','strand','motif_pvalue']]
    
    ## Produce normalised sequence lengths and positions
    df['normalised_motif_start'] = df['motif_start']/df['sequence_length']
    df['normalised_motif_start'] = np.where(df['strand'] == '+',df['normalised_motif_start'],1-df['normalised_motif_start'])
    ## For negative strand motif information, fraction is subtracted from 1 to give the same fractional distance from 5' orientation
    return df
