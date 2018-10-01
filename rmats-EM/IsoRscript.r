#!/usr/bin/env Rscript

# Shihao Shen 05/22/18
# This is the script to input data into rMATS_ISO R function
args = commandArgs(trailingOnly=TRUE)

# The function "rMATS_Iso" takes the following inputs:

# IsoExon: String containing the name of the .IsoExon file to be read as input.
# IsoMatrix1: Vector of strings containing the names of the .IsoMatrix files for sample group 1 to be read as input.
# IsoMatrix2: Vector of strings containing the names of the .IsoMatrix files for sample group 2 to be read as input. If there is only 1 sample group, set this parameter to NULL. The default value is NULL.
# outfile: String specifying the name of the output file that will be generated with all of the statistical results.
# detect_clusters: TRUE or FALSE. Should the number of computing clusters be automatically determined? Default is FALSE.
# numClusters: Integer specifying how many computing clusters will be used. Default value is 4.
# seed: Random seed. Default value is 17495.
# S: Number of samples to draw for importance sampling. Default value is 500.
# pval_thresh: Below which module-wise p-value shall individual isoform differences be tested? Default value is 1.
# S_pairwise: Number of simulations to perform to compute pairwise isoform p-values. Default value is 1e+05.
# nIter: Maximum number of EM iterations to perform. Default value is 100.
# epsilon: Error tolerance for terminating the EM algorithm. Default value is 10^-2.

# test if there is at least one argument: if not, return an error
if (length(args)<4) {
  stop("At least four arguments must be supplied", call.=FALSE)
}else{
  IsoExon    <- args[1];
  IsoMatrix1 <- strsplit(args[2],',')[[1]];
  IsoMatrix2 <- strsplit(args[3],',')[[1]];
  outfile    <- args[4];
  
  #default
  detect_clusters <- FALSE;
  numClusters     <- 4;
  seed            <- 17495;
  S               <- 500;
  pval_thresh     <- 1;
  S_pairwise      <- 100000;
  nIter           <- 100;
  epsilon         <- 0.01;
  
  if (length(args)>4) {detect_clusters <- as.logical(args[5]);}
  if (length(args)>5) {numClusters <- as.integer(args[6]);}
  if (length(args)>6) {seed <- as.integer(args[7]);}
  if (length(args)>7) {S <- as.integer(args[8]);}
  if (length(args)>8) {pval_thresh <- as.numeric(args[9]);}
  if (length(args)>9) {S_pairwise <- as.integer(args[10]);}
  if (length(args)>10) {nIter <- as.integer(args[11]);}
  if (length(args)>11) {epsilon <- as.numeric(args[12]);}
  
  library(rMATSISO)
  results <- rMATS_Iso(IsoExon, IsoMatrix1, IsoMatrix2, outfile, detect_clusters, numClusters, seed, S, pval_thresh, S_pairwise, nIter, epsilon);

}
