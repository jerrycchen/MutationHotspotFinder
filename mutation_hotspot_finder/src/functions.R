countMutations <- function(domain.info=NA, mutation.info=NA) {
  #' summarize the number of mutations located in each of the domains of interest
  #' return a numeric vector with the same rows of domains
  data.tmp <- domain.info
  out <- rep(0,nrow(data.tmp))
  for (tt in 1:nrow(data.tmp)) {
    loci <- seq(data.tmp[tt,2], data.tmp[tt,3])
    count.tmp <- sum(mutation.info %in% loci)
    out[tt] <- count.tmp
  }
  
  return(out)
  
}

resamplingSimulation <- function(input.coordinates='input/domain_coordinates.csv',
                                 input.positions='input/mutation_positions.csv',
                                 simulation.range=c(1,872),
                                 number.iteration=10000) {
  #' running resampling simulation
  #' return a data frame with resampling simulation statistics
                                   
  ## load original domain coordinates and mutation positions
  domain.info <- read.csv(file=input.coordinates, header=T, stringsAsFactors=F)
  mutation.info <- read.csv(file=input.positions, header=T)$POSITION
                                   
  ## filter out domains and mutations out of the simulation range
  domain.info <- domain.info[which((domain.info$START>=simulation.range[1]) & (domain.info$END<=simulation.range[2])), ]
  mutation.info <- mutation.info[mutation.info %in% (simulation.range[1]:simulation.range[2])]
  
  ## start simulations
  true.counts <- countMutations(domain.info, mutation.info)
  simulation.counts <- data.frame()
  simulation.counts <- rbind(true.counts, simulation.counts)
  colnames(simulation.counts) <- domain.info$DOMAIN
  
  resample.size <- length(mutation.info)
  resample.space <- simulation.range[1]:simulation.range[2]
  for (k0 in 2:number.iteration) {
    resample.tmp <- sample(resample.space, resample.size, replace=T)
    resample.count <- countMutations(domain.info, resample.tmp)
    simulation.counts[k0, ] <- resample.count
  }
  expected.counts <- as.numeric(round(apply(simulation.counts[(2:number.iteration), ], 2, FUN=mean),5))
  raw.epv <- as.numeric(round(apply(simulation.counts, 2, FUN=function(x){mean(x>=x[1])}),5))
  adj.epv <- raw.epv * length(domain.info$DOMAIN)
  adj.epv[adj.epv>1] <- 1
  out <- data.frame(DOMAIN=domain.info$DOMAIN, OBSERVED=true.counts, EXPECTED=expected.counts, RAW.P.VALUE=raw.epv, ADJ.P.VALUE=adj.epv)
  
  return(out)
                                   
}