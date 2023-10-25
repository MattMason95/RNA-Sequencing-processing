# AUTHOR: Matthew Mason
# EDITED: 10.2023
# DESC: R script for generating a single FeatureCounts matrix for submission to DeSeq2 

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> LIBRARY INSTANTIATION                      <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

library(magrittr)
library(dplyr)
library(purrr)

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> SCRIPT                                     <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

## Fetch files matching pattern
FILES<- list.files(pattern="*_featureCounts.txt$", recursive=True, full.names=True)

## Function for accessing the relevant columns of the individual files
countMerger<- function(f){
  counts<- read_tsv(f,col_names=True,comment=#)
  counts<- counts %>% dplyr::select(-c(Chr,Start,End,Strand,Length))
  return(counts)
}

rawCounts<- map(FILES, countMerger)
countMatrix<- purrr::reduce(rawCounts,full_join)

## FIN                  
