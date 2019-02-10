## This script is used to generate all results of the Zhang, Chen and Kim paper

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



