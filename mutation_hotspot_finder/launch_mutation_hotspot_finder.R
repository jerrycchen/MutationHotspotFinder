## This script is used to launch Mutation Hotspot Finder
## User needs to specifiy 3 parameter files in the directory "./input/" before running
##
## File I - "domain_coordinates.csv" specifies protein domain information
## ----
## User needs to specify domain names, the start and ending amino acid (AA) residue loci
## Each domain name and domain coordinates takes a new row. Please DO NOT change the header
## As an example, domain coordinates of KCNQ2 (Gene Bank: NP_742105.1) are provided
##
## File II - "mutation_positions.csv" specifies AA mutation loci
## ----
## User needs to specify the loci of mutations. Each row corresponds to one locus
## In case there are multiple mutations observed at one locus, record them in multiple rows. Please DO NOT change the header
## For example, there are two mutations observed at AA locus 114 in our example file
## 
## File III - "set_parameters.csv"
## ----
## User needs to specify the output file name, input domain_coordinates and mutation_position file names in this files
## In case user is only interested in finding mutation hotspot in a sub region of a protein, user may specify
## the start and ending position of such sub region, and the analysis will ignore domains and mutations out of this region
##
## In case there are multiple hypothesis to test, user may provide multiple File I's and File II's.
## Each hypotheis will take a new row in File III
## Please DO NOT change the header
## 

rm(list=ls())
source('src/functions.R')

parList <- read.csv(file='input/set_parameters.csv')

for (iter in 1:nrow(parList)) {
  output <- NA
  input.coordinates <- paste('input/', parList$Coord[iter], sep='')
  input.positions <- paste('input/', parList$Position[iter], sep='')
  simulation.range <- c(parList$Start[iter], parList$End[iter])
  outFile.path <- paste('output/', parList$Output[iter], sep='')
  
  output <- resamplingSimulation(input.coordinates, input.positions, simulation.range)
  
  write.csv(output, file=outFile.path, row.names=F)
  cat(iter, '\n')
}



