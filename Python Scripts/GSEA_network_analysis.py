# Geneset enrichment analysis (GSEA) yields an abubdance of data that is often redundant in nature, due to the overlapping nature of the biological genesets assessed. 
# The aim of this tool is to compress the informational complexity of the data produced by GSEA by building a relational graph with NetworkX.
# In this graph, Nodes represent genesets and Edges represent the Jaccard index of shared information between the two genesets. 
# Subgraphs within the network will be evaluated with basic NLP to summarise the semantic content of each subgraph (i.e. the biological theme of genesets) 
# Author: Matthew Mason
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import libraries 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import pandas as pd 
import numpy as np
import seaborn as sns
import matplotlib as mp
import matplotlib.pyplot as plt 
import palettable as pal 
import tarfile
import re
import os
import itertools
import networkx as nx
import regex as re

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def Jaccard(moduleNames,geneLists):
  '''
  This function will be nested within the fileAccessor function to provide the capability to calculate Jaccard indices for comparing genesets. 
  The basic principle of the Jaccard index is a calculation of the percentage of overlapping elements within two lists.
  The biological assumption is that genesets with overlapping membership will be related to the same biological process, and will thus reduce redundancy/ 
  '''
  ## Build master and output containers
  masterList = list(zip(moduleNames,geneLists))
  
  node1List = []
  node2List = []
  jaccard = []
  

  ## Use itertools.combinations to iteratively compare all possible pairs in the masterlist 
  for x,y in itertools.combinations(masterList, 2):
    ## Fetch names of compared genesets (node1 and node2) and append to output containers  
    node1 = str(x[0])
    node2 = str(y[0])
    node1List.append(node1)
    node2List.append(node2)

    ## Fetch genes associated with compared genesets and assess these for overlap 
    s1 = set(x[1])
    s2 = set(y[1])
    ## set.intersetion() = Overlap; set.union() = Total unique 
    statistic = float(len(s1.intersection(s2)) / len(s1.union(s2)))
    jaccard.append(statistic)

  ## Zip together output containers for export
  jaccardOutput = list(zip(node1List,node2List,jaccard))  
  
  return jaccardOutput

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def edbFileParser(edbData,label,mode='full'):
  '''
  This function will take the full array of EDB data and extract the geneset identifier, normalised enrichment scores, and FDR statistics using Regex pattern recognition
  May build up Full versus Basic modes to modify the extracted specific information. 
  '''
  modes = ['full'] #'basic'
  if mode not in modes:
    raise ValueError("Invalid mode. Expected one of: %s" % modes)

  if mode == 'full':      
    ## Build outputs containers
    genesetList = []
    nesList = []
    fdrList = []
    labels = []

    ## Iterate through the EDB 
    for i in range(len(edbData)):
      for j in range(len(edbData[i])):
        string = edbData[i][j]
        if 'DTG RANKED_LIST' in string:
          genesetParser = re.search(r'(?<=gene_sets.gmt#).*?(?=" ES=")', string).group(0)
          nesParser = float(re.search(r'(?<=NES=").*?(?=" NP=")', string).group(0))
          fdrParser = float(re.search(r'(?<=FDR=").*?(?=" FWER=")', string).group(0))

          ## If the absolute Normalised Enrichment Score is less than 1, skip the geneset
          if abs(nesParser) < 1:
              continue
          ## Otherwise, append all of the extracted information 
          else:
              genesetList.append(genesetParser)
              nesList.append(nesParser)
              fdrList.append(fdrParser)
              labels.append(label)
        else:
          continue
                  
      edbParserOutput = list(zip(genesetList,nesList,fdrList,labels))
      return edbParserOutput
 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
def significantGenesets(zippedList,threshold,unzip=False):
  '''
  This is a simple function to apply a statistical threshold to the input data and, when requested (unzip=True) to parse output the elements of the zipped input list.
  '''
  ## Valid conditions for unzip
  unzips = [True, False]
  if unzip not in unzips:
    raise ValueError('Unzip requires a boolean operator: True or False.')
        
  ## This loop goes through the numerical data for the output Jaccard list and applies a threshold filter - retain only those <= the statistical threshold
  filteredGenesets = [x for x in zippedList if float(x[2]) <= threshold]

  ## In most cases, I will want to have the statistical data unzipped - this function will return the unzipped elements from the filteredGeneset
  if unzip == True: 
      def unzip(iterable):
          return zip(*iterable)
      geneList, nesList, fdrList = unzip(filteredGenesets)
      return geneList, nesList, fdrList

## In the event that I want to keep the elements zipped, unzip=False will return the native filteredGeneset object
  else:
      return filteredGenesets
 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
def fileAccessor(homeDirectory,jaccardFilter):
  ''' 
  This is the primary function that integrates all above functions. 
  fileAccessor accepts a home directory - specified by the user, in the form of a variable-assigned string - from which the OS module will walk through the available directories and files.
  Specifically, the function searches for the destination folders called << edb >> in which the raw data are kept; all other folders are skipped.
  Once within the << edb >> folder the function extracts the two raw data files (.gmt, .edb)
  The default output path for GSEA follows this same architecture, so the function should be generalisable.
  '''
  ## Survey directory trees to find the target folders << edb >>
  destinations = []
  for root, dirs, files in os.walk(homeDirectory, topdown=True):
      for folder in dirs:
          print(folder)
          if 'edb' not in folder:
              continue
          else:
              destination = str(f'{root}/{folder}')  #os.path.join(root,folder)
              print('destination:',destination)
              destinations.append(destination)
  
  ## Print out the identified file paths
  print(destinations)

  ## Build output containers 
  networkTemp = []
  metaTemp = []
  masterGeneLookup = []
  masterStats = []

  ## Iterate through the target file paths identified above
  for target in destinations:
    ## At this indentation level, we are within one of the target folder paths, so handling of gmt and edb data unique to one folder - i.e. unique to one pairing of conditions 
    print('target:',target)
    gmtData = []
    edbData = []

    for files in os.listdir(target):
      if '.gmt' in files:
        ## Open and extract data from the .gmt file - append to output container
        fileName = os.path.join(target,files)
        with open(fileName) as gmt_:
          rawGmt = gmt_.read()
          splitGmt = re.split(r'\n', rawGmt)
                  
        for i in range(len(splitGmt)):
          GMT = re.split(r'\t',splitGmt[i])
          gmtData.append(GMT)
                  
      elif '.edb' in files:
        ## Open and extract data from the .edb file - append to output container
        fileName = os.path.join(target,files)
        with open(fileName) as edb_:
          rawEdb = edb_.read()
          EDB = re.findall(r'<(.+?)>',rawEdb)
          edbData.append(EDB)
      else:
          continue
        
    ## First handle the .gmt data
    ## Build output containers for the full list of module names and associated genes from those modules (within .gmt data)
    moduleNames = []
    moduleGenes = []

    for i in range(len(gmtData)):
      ## Iterate through the .gmt data to extract the geneset/module name and the constituent genes of that module
      name = gmtData[i][0]
      genes = gmtData[i][2:]
      moduleNames.append(name)
      moduleGenes.append(genes)

    ## Zip together module names and module genes to avoid mismatching
    modulesZipped = list(zip(moduleNames,moduleGenes))
    
    ## Second handle the .edb data
    ## edbFileParser extracts the statistical data (NES and FDR) from the edb file, with the corresponding geneset name
    statsData = edbFileParser(edbData,label=target,mode='full')

    ## significantGenesets receives the zipped statistical data from edbFileParser and filters the full list with an FDR cutoff of <0.1 - it also has the functionality to unzip data
    statsDataThresholded = significantGenesets(statsData,0.1,unzip=False)

    ## Extract geneset names from the significant geneset data and use this to truncate the modulesZipped variable to include only the statistically significant genesets                
    sigGenesets = [x[0] for x in statsDataThresholded]
    sigModules = [x for x in modulesZipped if x[0] in sigGenesets]
    geneLookup = list(zip(moduleNames,moduleGenes,target))

    ## Calculate the Jaccard indices for the significantly enriched genesets
    jaccardIndices = Jaccard([x[0] for x in sigModules],[x[1] for x in sigModules])


    ## Produce output dataframe for data within this loop iteration
    ## Format of this dataframe is a node-node relationship table, which is packaged within the jaccardIndices zipped list
    networkData = pd.DataFrame({'node1':[x[0] for x in jaccardIndices],
                                'node2':[x[1] for x in jaccardIndices],
                                'jaccard_index':[x[2] for x in jaccardIndices],
                                'draw_edge':[1 if x[2] > jaccardFilter else 0 for x in jaccardIndices],
                                'condition':target})
    ## Format of this dataframe is a statistical information table, which is packaged within the statsDataThresholded zipped list
    metaData = pd.DataFrame({'node':[x[0] for x in statsDataThresholded],
                             'nes':[x[1] for x in statsDataThresholded],
                             'fdr':[x[2] for x in statsDataThresholded],
                             'upregulated?':[1 if x[1] > 0 else 0 for x in statsDataThresholded],
                             'condition':target})

    ## Append data from this loop to the output containers
    networkTemp.append(networkData)
    metaTemp.append(metaData)
    masterGeneLookup.append(geneLookup)
    masterStats.append(statsData)

  ## Concatenate output dataframes
  masterNetwork = pd.concat(networkTemp)
  masterMeta = pd.concat(metaTemp)
      
  
  return masterNetwork, masterMeta, masterStats, masterGeneLookup
 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
def iterativePlotting(network,meta,save=True):
  '''
  This function will iteratively plot the Network graphs for each condition within the data extracted from the previous function using the Python library NetworkX.
  '''
  ## Import specific libraries
  import sklearn  
  from sklearn.preprocessing import minmax_scale

  ## Option for automatic saving of plots to output path
  saves = [True,False]
  if save not in saves:
    raise ValueError('Save requires a boolean operator: True or False.')

  ## Iterate through the different conditions within the network analysis dataframe 
  for condition in network['condition'].unique():
    saveName = str(condition)

    ## Truncate to condition 
    shortNetwork = networkDf[networkDf['condition'] == condition]
    shortMeta = meta[meta['condition'] == condition]

    ## Generate edgelist from dataframe
    graph = nx.from_pandas_edgelist(shortNetwork, source='node1', target='node2', edge_attr='jaccard_index')

    ## Specify that no edges are created between nodes with a Jaccard index of < 0.5 - remove tham from the edge list 
    noEdge = list(filter(lambda e: e[2] < 0.5, (e for e in graph.edges.data('jaccard_index'))))
    noEdgePairs = list(e[:2] for e in noEdge)
    
    graph.remove_edges_from(noEdgePairs)

    ## Generate colour maps from the meta data provided - colour red if upregulated or blue if downregulated
    colorMap = []
    for node in graph:
      if node in list(shortMeta['node'][shortMeta['upregulated?'] == 1]):
        colorMap.append('#A43D40')
      else: 
        colorMap.append('#507D96')

    ## Extract the NES for each of the plotted nodes in the graph 
    nodeMap = []
    for node in graph:
      if node in list(shortMeta['node']):
        nes = shortMeta['nes'][shortMeta['node'] == node].values[0]
        nodeMap.append(nes)
      else: 
        continue     

    ## Normalise the scaling of the extracted NES values to use for node sizes
    normalisedNodeMap = minmax_scale(nodeMap,feature_range=(-2,2))

    ## Generate force-directed network graph using the NetworkX library
    pos = nx.spring_layout(graph, k=2/np.sqrt(len(graph.nodes())), iterations=90, seed=101, weight='jaccard_index',center=(0,0))
    
    ## Plotting
    fig, ax = plt.subplots(1,1,figsize=(12,12))
    fig = nx.draw(graph,pos=pos,node_color=colorMap,edgecolors='0.2',edge_color='0.2',linewidths=2,node_size=((normalisedNodeMap+3)*35),width=2)
    plt.tight_layout()

    ## Save if requested
    if save == True:
      plt.savefig(f'Images/{saveName}.png',dpi=150)
    else:
      continue

