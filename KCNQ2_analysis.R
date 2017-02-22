##============Step 1: Data Preparation============##

##1-1 Input the number of amino acid residuals in the gene of interest
Protein_Length <- 872;  # Plug in the length of KCNQ2 protein sequence

##1-2 Fill the mutation position vector using the information from the database
Inputframe <- read.table("input/KCNQ2_missense position.txt", header=TRUE);  # Copy the mutation positions from this file

Mutant_Position <- Inputframe$AA_POSITION;  # There are 78 different positions for missense mutations

##1-3 Form the input matrix
Input_Data <- matrix(0, nrow=length(Mutant_Position), ncol=Protein_Length);  # Define the size and initial value of input

for (i in 1:length(Mutant_Position))  {
  
    j <- Mutant_Position[i];  
    Input_Data[i,j] <- 1;  
  
};

##============End of Step 1============##




##============Step 2: Run Code from NMC algorithm============##

source("nmc.R");

##============End of Step 2============##




##============Step 3: Analyze the input file============##

##3-1 Set for the desired methods
alpha <- 0.001;
multiple_test <- "Bonferroni";  # Alternatively, may set as "BH" or "None"

##3-2 Run test using these parameters from the last step
write.table(nmc(Input_Data, alpha=alpha, multtest=multiple_test), file="output/result_nmc method_KCNQ2.txt", sep="\t");  # Write the output into a separate file

##============End of Step 3============##



