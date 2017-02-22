##================  KCNQ3  ================##
## Perform analysis on the whole coding region

## Dataset preparation
Domain <- read.table("input/domain_coordinates_KCNQ3.txt", header=TRUE);
Position <- read.table("input/missense_position_KCNQ3.txt", header=TRUE)$AA_POSITION;
num.AA <- 872;
num.perm <- 10^5;
Range <- 230:381;  # Identified by previous test using the nmc method

## Design the mapping matrix to count # of mutations in each domain
MC <- matrix(0, nrow=nrow(Domain)+2, ncol=num.AA);
for (k in 1:nrow(Domain)) {
    c1 <- Domain[k,2];
    c2 <- Domain[k,3];
    MC[k,c1:c2] <- 1;
};

Group1 <- as.integer(colSums(MC[1:nrow(Domain),Range])==0);
MC[(nrow(Domain)+1),Range] <- Group1;

Group2 <- as.integer(colSums(MC[1:(nrow(Domain)+1), ])==0);
MC[(nrow(Domain)+2), ] <- Group2;

rownames(MC) <- c(as.character(Domain$DOMAIN), "Blank1", "Blank2");

## Summit the mutation frequencies at each position
True.Sum <- rep(0,num.AA);
for (k in 1:length(Position)) {
    True.Sum[Position[k]] <- True.Sum[Position[k]] + 1;
};

Output.Data <- as.data.frame(t(MC %*% True.Sum));
rownames(Output.Data)[1] <- "True.Value";

## Perform resampling and save the data
for (m in 1:num.perm) {
  Position.Sample <- sample(x=1:num.AA, size=length(Position), replace=TRUE);
  Perm.Sum <- rep(0,num.AA);
  
  for (k in 1:length(Position.Sample)) {
    Perm.Sum[Position.Sample[k]] <- Perm.Sum[Position.Sample[k]] + 1;
  };
  
  Output.Data[m+1, ] <- t(MC %*% Perm.Sum)[1, ];
};

## Compute the resampling extreme probability
EXT <- function(t) {
  t0 <- t[1];
  tp <- t[-1];
  percent1 <- mean(t0>=tp);
  percent2 <- mean(t0<=tp);
  return(min(c(percent1,percent2)));
};

p.Raw <- apply(Output.Data, 2, FUN = EXT);
p.Bon <- p.adjust(p.Raw, method="bonferroni");
p.Fdr <- p.adjust(p.Raw, method="fdr");

Obs <- as.vector(as.matrix(Output.Data[1, ]));
Exp <- round(apply(Output.Data[2:num.perm, ], 2, FUN=mean),1);
Output.Pval <- data.frame(Observed=Obs, Expected=Exp, Raw=p.Raw, Bonferroni=p.Bon, Fdr=p.Fdr);

## Write the simulated outcome into text
write.table(Output.Data, file="output/KCNQ3_resampling details.txt", sep="\t");
write.table(Output.Pval, file="output/KCNQ3_resampling summary.txt", sep="\t");




## Perform analysis within the region [226,389]
## No "S1", "S2", "S3", "Helix_B", "Helix_C", "Helix_D" or "Ankyrin_G_Binding"

## Dataset preparation
Domain <- read.table("input/domain_coordinates_KCNQ3.txt", header=TRUE);
Position <- read.table("input/missense_position_KCNQ3.txt", header=TRUE)$AA_POSITION;
Range <- 230:381;  # Identified by previous test using the nmc method

Range.New <- c(226,389);  #@__@ Specify the region to be analyzed;

num.AA <- (Range.New[2] - Range.New[1] + 1);
num.perm <- 10^5;  #-__-;

# Update the domain information
Domain <- Domain[(Domain$START>=Range.New[1] & Domain$END<=Range.New[2]), ];
Domain$START <- Domain$START-(Range.New[1]-1);
Domain$END <- Domain$END-(Range.New[1]-1);

# Update the position information
Position <- Position[(Position>=Range.New[1]&Position<=Range.New[2])];
Position <- Position - (Range.New[1] - 1);

## Design the mapping matrix to count the number of mutations in each domain
MC <- matrix(0, nrow=nrow(Domain)+1, ncol=num.AA);
for (k in 1:nrow(Domain)) {
  c1 <- Domain[k,2];
  c2 <- Domain[k,3];
  MC[k,c1:c2] <- 1;
};

Blank <- as.integer(colSums(MC[1:nrow(Domain), ])==0);
MC[(nrow(Domain)+1), ] <- Blank;

rownames(MC) <- c(as.character(Domain$DOMAIN), "Blank1");

## Summit the mutation frequencies at each position
True.Sum <- rep(0,num.AA);

for (k in 1:length(Position)) {
  True.Sum[Position[k]] <- True.Sum[Position[k]] + 1;
};

Output.Data <- as.data.frame(t(MC %*% True.Sum));
rownames(Output.Data)[1] <- "True.Value";

## Perform resampling and save the data
for (m in 1:num.perm) {
  Position.Sample <- sample(x=1:num.AA, size=length(Position), replace=TRUE);
  Perm.Sum <- rep(0,num.AA);
  
  for (k in 1:length(Position.Sample)) {
    Perm.Sum[Position.Sample[k]] <- Perm.Sum[Position.Sample[k]] + 1;
  };
  
  Output.Data[m+1, ] <- t(MC %*% Perm.Sum)[1, ];
};

## Compute the resampling extreme probability
EXT <- function(t) {
  t0 <- t[1];
  tp <- t[-1];
  percent1 <- mean(t0>=tp);
  percent2 <- mean(t0<=tp);
  return(min(c(percent1,percent2)));
};

p.Raw <- apply(Output.Data, 2, FUN = EXT);
p.Bon <- p.adjust(p.Raw, method="bonferroni");
p.Fdr <- p.adjust(p.Raw, method="fdr");

Obs <- as.vector(as.matrix(Output.Data[1, ]));
Exp <- round(apply(Output.Data[2:num.perm, ], 2, FUN=mean),1);
Output.Pval <- data.frame(Observed=Obs, Expected=Exp, Raw=p.Raw, Bonferroni=p.Bon, Fdr=p.Fdr);

## Write the simulated outcome into text
write.table(Output.Data, file="output/KCNQ3_resampling details_incluster.txt", sep="\t");
write.table(Output.Pval, file="output/KCNQ3_resampling summary_incluster.txt", sep="\t");