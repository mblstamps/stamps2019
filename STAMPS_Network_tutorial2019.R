#---
#  title: "STAMPS network SPIEC-EASI/SPRING Tutorial 2019"
#  author: "Zachary D. Kurtz, Christian L. Mueller"
#  date: "7/26/2019"
#---
  

# Sparse InversE Covariance estimation for Ecological Association and Statistical Inference #
# 
# This package will be useful to anybody who wants to infer graphical models for all sorts of compositional data, 
# though primarily intended for microbiome relative abundance data (generated from 16S amplicon sequence data). 
# It also includes a generator for [overdispersed, zero inflated] multivariate, correlated count data. 
# Please see the paper published in [PLoS Comp Bio](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004226).

# One small point on notation: we refer to the method as "SPIEC-EASI" and reserve "SpiecEasi" for this package.

# Installation #
# 
# We  assume that all auxiliary packages are already installed - 
# for example pulsar, huge, MASS, etc. 
# If you get an unexpected error, you may need to download and install a missing dependency.

# From an interactive R session (already installed for the STAMPS course):
  
# library(devtools)
# install_github("zdk123/SpiecEasi")

library(SpiecEasi)

# News #
# The latest SpiecEasi (version 1.0.0 and higher) 
# now uses the [pulsar package](https://cran.r-project.org/package=pulsar) for stability-based model selection. 
# The methods are similar to previous releases, but contains some additional methods for 
# [speeding up computations](#pulsar-parallel-utilities-for-model-selection)

# The input arguments have changed slightly (but are backwards compatible) but the data structure 
# that is returned from `spiec.easi` has changed.
# The output to spiec.easi-fit models structure can be easily processed using new getter functions. 
# See `?getOptInd` for usage.
# You can revert to the previous release ([0.1.4](https://github.com/zdk123/SpiecEasi/releases/tag/v0.1.4)) 
# to avoid code-breaking changes.
 
################################  
# Basic Usage #
################################  

# Lets simulate some multivariate data under zero-inflated negative binomial model, 
# based on (high depth/count) round 1 of the American gut project, with a sparse network. 
# The basic steps are

  # 1. load the data and normalize counts to to common scale (min depth)
  # 2. fit count margins to the model
  # 3. generate a synthetic network
  # 4. generate some synthetic data
  # 5. clr transformation
  # 6. inverse covariance estimation along a lambda (sparsity) path
  # 7. stability selection using the StARS criterion
  # 8. evaluate performance
  
# If you analyze real data, skip 1-4.
  

# First, let's download data and R files necessary for the tutorial

## Set your working directory to your home directory
# setwd("~/")
#
# download.file("https://ndownloader.figshare.com/files/16638143", destfile = "STAMPS_Network_tutorial.tar.gz")
# untar("STAMPS_Network_tutorial.tar.gz")
#
# setwd("~/STAMPS_Network_tutorial")

# Filtered data of the American Gut project's first round
load("data/amgut1.filt.rda")

# Normalize data
depths <- rowSums(amgut1.filt)
amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))
  
# d: Number of OTUs
d <- ncol(amgut1.filt.cs)
# n: Number of samples
n <- nrow(amgut1.filt.cs)

################################  
# Warm-up: Generate synthetic data
################################  

set.seed(10010)
# e: Number of true edges in the synthetic graph
e <- d

# Generate graph with "cluster" topology
graph <- SpiecEasi::make_graph('cluster', d, e)

# Generate inverse covariance/precision matrix
Prec  <- graph2prec(graph)
# Derive correlation matrix
Cor   <- cov2cor(prec2cov(Prec))

# Generate synthetic microbiome data   
X <- synth_comm_from_counts(amgut1.filt.cs, mar=2, distr='zinegbin', Sigma=Cor, n=n)

# Let's run SPIEC-EASI on the data 
# The main SPIEC-EASI pipeline: 
# - Data transformation (CLR)
# - Sparse inverse covariance estimation (GLASSO or MB)
# - Stability-based model selection
 
se <- spiec.easi(X, method='mb', lambda.min.ratio=1e-2, nlambda=15)
  
# Examine ROC over lambda path 
huge::huge.roc(se$est$path, graph, verbose=FALSE)

# Examine precision-recall over the StARS edge probabilities for the selected graph
stars.pr(getOptMerge(se), graph, verbose=FALSE)
  
# StARS selected final network under: 
getRefit(se)

# The above example does not cover all possible options and parameters. 
# For example, other generative network models are available, 
# the lambda.min.ratio (the scaling factor that determines the minimum sparsity/lambda parameter) 
# shown here might not be right for your dataset, 
# and it's possible that you'll want more repetitions (number of subsamples) for StARS.
  
################################  
# Analysis of American Gut data 
################################  

# Now let's apply SpiecEasi directly to the American Gut data. 
# Don't forget that the normalization is performed internally in the `spiec.easi` function. 
# Also, we should use a larger number of StARS repetitions for real data. 
# We can pass in arguments to the inner StARS selection function as a list via the parameter `pulsar.params`. 
# If you have more than one processor available, you can also supply a number to `ncores`. 
# Also, let's compare results from the MB and glasso methods as well as SparCC (correlation).
  
# Meinshausen-Buehlmann Neighborhood selection
se.mb.amgut <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-2,
                            nlambda=20, pulsar.params=list(rep.num=50))

# "Graphical Lasso" method
se.gl.amgut <- spiec.easi(amgut1.filt, method='glasso', lambda.min.ratio=1e-2,
                            nlambda=20, pulsar.params=list(rep.num=50))
# "SparCC" method
sparcc.amgut <- sparcc(amgut1.filt)
## Define threshold 0.3 for SparCC correlation matrix for the graph
sparcc.graph <- abs(sparcc.amgut$Cor) >= 0.3
diag(sparcc.graph) <- 0

library(Matrix)
sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)

## Let's visualize the resulting graphs

## Create igraph objects
ig.mb     <- adj2igraph(getRefit(se.mb.amgut))
ig.gl     <- adj2igraph(getRefit(se.gl.amgut))
ig.sparcc <- adj2igraph(sparcc.graph)

library(igraph)
## Set size of vertex proportional to clr-mean
  
vsize    <- rowMeans(clr(amgut1.filt, 1))+6
am.coord <- layout.fruchterman.reingold(ig.mb)
  
par(mfrow=c(1,3))
plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
plot(ig.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")

# We can evaluate the weights on edges networks using the terms from the underlying model. 
# SparCC correlations can be used directly, while SpiecEasi networks need to be massaged a bit. 
# Note that since SPIEC-EASI is based on penalized estimators, 
# the edge weights are not directly comparable to SparCC (or Pearson/Spearman correlation coefficients)

library(Matrix)
  
secor  <- cov2cor(getOptCov(se.gl.amgut))

# Symmetrize the MB assymetric graph by taking the maximum value of the edge pairs e_ij / e_ji
sebeta <- symBeta(getOptBeta(se.mb.amgut), mode='maxabs')
elist.gl     <- summary(triu(secor*getRefit(se.gl.amgut), k=1))
elist.mb     <- summary(sebeta)
elist.sparcc <- summary(sparcc.graph*sparcc.amgut$Cor)
  
hist(elist.sparcc[,3], main='', xlab='edge weights')
hist(elist.mb[,3], add=TRUE, col='forestgreen')
hist(elist.gl[,3], add=TRUE, col='red')

## Let's look at the degree statistics from the networks inferred by each method.
  
dd.gl     <- degree.distribution(ig.gl)
dd.mb     <- degree.distribution(ig.mb)
dd.sparcc <- degree.distribution(ig.sparcc)
  
plot(0:(length(dd.sparcc)-1), dd.sparcc, ylim=c(0,.35), type='b',
       ylab="Frequency", xlab="Degree", main="Degree Distributions")
points(0:(length(dd.gl)-1), dd.gl, col="red" , type='b')
points(0:(length(dd.mb)-1), dd.mb, col="forestgreen", type='b')
legend("topright", c("MB", "glasso", "sparcc"), col=c("forestgreen", "red", "black"), pch=1, lty=1)
  
################################  
# Working with phyloseq #
################################  

# SpiecEasi includes some convenience wrappers to work directly with `phyloseq` objects.

library(phyloseq)
  
## Load round 2 of American gut project
data(amgut2.filt.phy)

se.mb.amgut2 <- spiec.easi(amgut2.filt.phy, method='mb', lambda.min.ratio=1e-2,
                             nlambda=20, pulsar.params=list(rep.num=50))
ig2.mb <- adj2igraph(getRefit(se.mb.amgut2),
                       vertex.attr=list(name=taxa_names(amgut2.filt.phy)))
plot_network(ig2.mb, amgut2.filt.phy, type='taxa', color="Rank3")

################################  
# SpiecEasi with multiple data sets - Cross-domain interactions
################################

# SpiecEasi now includes a convenience wrapper for dealing with multiple taxa sequenced on the same samples,
# such as 16S and ITS, 
# as seen in [Tipton, Müller, et. al. (2018)](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0393-0). It assumes that each taxa is in it's own data matrix and that all samples are in all data matrices in the same order.
  
# Here's an example run from the [HMP2 project](https://ibdmdb.org/tunnel/public/summary.html) 
# with 16S and Proteomics data.
  
## IMPORTANT: THIS DOES NOT WORK ON THE STAMPS servers!!!!

# library(phyloseq)
#data(hmp2)
#se.hmp2 <- spiec.easi(list(hmp216S, hmp2prot), method='mb', nlambda=40, lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05))
  
#dtype <- c(rep(1,ntaxa(hmp216S)), rep(2,ntaxa(hmp2prot)))
#plot(adj2igraph(getRefit(se.hmp2)), vertex.color=dtype+1, vertex.size=9)

################################
# pulsar: parallel utilities for model selection #
################################

# SpiecEasi is now using the [pulsar package](https://cran.r-project.org/package=pulsar) as the backend for performing model selection. In the default parameter setting, this uses the same [StARS](https://arxiv.org/abs/1006.3316) procedure as previous versions.
# As in the previous version of SpiecEasi, we can supply the `ncores` argument to the pulsar.params list to break up the subsampled computations into parallel tasks.

# In this example, we set the random seed to make consistent comparison across experiments.
  
## Default settings ##
pargs1 <- list(rep.num=50, seed=10010)
t1 <- system.time(
  se1 <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-3, nlambda=30,
                      sel.criterion='stars', pulsar.select=TRUE, pulsar.params=pargs1)
)
  
## Parallel multicore ##
pargs2 <- list(rep.num=50, seed=10010, ncores=4)
t2 <- system.time(
  se2 <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-3, nlambda=30,
                      sel.criterion='stars', pulsar.select=TRUE, pulsar.params=pargs2)
)

# We can further speed up StARS using the [bounded-StARS](https://arxiv.org/abs/1605.07072) ('bstars') method. The B-StARS approach computes network stability across the whole lambda path, but only for the first 2 subsamples. This is used to build an initial estimate of the summary statistic, which in turn gives us a lower/upper bound on the optimal lambda. The remaining subsamples are used to compute the stability over the restricted path. Since denser networks are more computational expensive to compute, this can significantly reduce computational time for datasets with many variables.
  
t3 <- system.time(
    se3 <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-3, nlambda=30,
                      sel.criterion='bstars', pulsar.select=TRUE, pulsar.params=pargs1)
)

t4 <- system.time(
    se4 <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-3, nlambda=30,
                      sel.criterion='bstars', pulsar.select=TRUE, pulsar.params=pargs2)
)

# We can see that in addition to the computational savings, the refit networks are identical.
  
## serial vs parallel
identical(getRefit(se1), getRefit(se2))

t1[3] > t2[3]
## stars vs bstars

identical(getRefit(se1), getRefit(se3))
t1[3] > t3[3]

identical(getRefit(se2), getRefit(se4))
t2[3] > t4[3]

## Batch Mode ##
# Pulsar gives us the option of running stability selection in batch mode, 
# using the [batchtools](https://mllg.github.io/batchtools/) package. 
# This will be useful to anyone with access to an hpc/distributing computing system. 
# Each subsample will be independently executed using a system-specific cluster function.
  
# This requires an external config file which will instruct the batchtools registry 
# how to construct the cluster function which will execute the individual jobs. 
# `batch.pulsar` has some built-in config files that are useful for 
# testing purposes (serial mode, "parallel", "snow", etc), 
# but it is advisable to create your own config file and pass in the absolute path. 
# See the [batchtools docs](https://mllg.github.io/batchtools/articles/batchtools.html#configuration-file) 
# for instructions on how to construct config file and template files (i.e. to interact with a queueing system 
# such as TORQUE or SGE).
                                                                                                                                                                                                                                                                                                                                                                                                         
bargs <- list(rep.num=50, seed=10010, conffile="parallel")
## See the config file stores:
pulsar::findConfFile('parallel')

## remove line below to turn on batchtools reporting
options(batchtools.verbose=FALSE)
se5 <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-3, nlambda=30,
                  sel.criterion='stars', pulsar.select='batch', pulsar.params=bargs)

################################
# SPIEC-EASI with latent variable graphical model (sparse+low-rank (slr) mode)
################################

# An alternative model of the precision matrix is a decomposition of sparse and low rank components 
# (stay tuned for a preprint manuscript in 2019).
# We can treat unobserved covariates as low rank and spread out, we can disentangle compositional, 
# technical artifacts and biological confounders from direct interactions, even if we can't observe them.

# This requires us to install an experimental branch of SpiecEasi:
library(devtools)
install_github("zdk123/SpiecEasi", ref='lowrank')

# Lets compare our SpiecEasi methods on the first two rounds of the American gut dataset.
args <- list(ncores=4)

se1 <- list(
  mb=spiec.easi(amgut1.filt, method='mb', pulsar.params=args,lambda.min.ratio=1e-2, sel.criterion='bstars'),
  glasso=spiec.easi(amgut1.filt, method='glasso', lambda.min.ratio=1e-2, pulsar.params=args, sel.criterion='bstars'),
  slr=spiec.easi(amgut1.filt, method='slr', r=2,pulsar.params=args, sel.criterion='bstars'))

se2 <- list(
  mb=spiec.easi(amgut2.filt.phy, method='mb', pulsar.params=args,
  lambda.min.ratio=1e-2, sel.criterion='bstars'),
  glasso=spiec.easi(amgut2.filt.phy, method='glasso', pulsar.params=args,
  lambda.min.ratio=1e-2, sel.criterion='bstars'),
  slr=spiec.easi(amgut2.filt.phy, method='slr', r=2, pulsar.params=args,
  lambda.min.ratio=1e-1, sel.criterion='bstars'))

# On the set of common taxa, we compare the recovery using the round 2 network as the 'true' network:

library(phyloseq)
common_tax <- intersect(taxa_names(amgut2.filt.phy), colnames(amgut1.filt))

ind1 <- na.exclude(match(common_tax, colnames(amgut1.filt)))
ind2 <- na.exclude(match(common_tax, taxa_names(amgut2.filt.phy)))

subset <- function(path) lapply(path, function(x) x[ind1,ind1])

pr.mb  <- stars.pr(subset(se1$mb$est$path), getRefit(se2$mb)[ind2,ind2], verbose=FALSE)
pr.gl  <- stars.pr(subset(se1$glasso$est$path), getRefit(se2$glasso)[ind2,ind2], verbose=FALSE)
pr.slr <- stars.pr(subset(se1$slr$est$path), getRefit(se2$slr)[ind2,ind2], verbose=FALSE)

## compare max F1 score
max(pr.slr$F1)
max(pr.mb$F1)
max(pr.gl$F1)

# Troubleshooting #

# A common issue that comes up with when running `spiec.easi` is coming up with an empty network after running StARS.
# For example:
pargs <- list(seed=10010)
se <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=5e-1,
nlambda=10, pulsar.params=pargs)
getOptInd(se)
sum(getRefit(se))/2


# As the warning indicates, the network stability could not be determined from the lambda path. 
# Looking at the stability along the lambda path, `se$select$stars$summary`, 
# we can see that the maximum value of the StARS summary statistic never crosses the default threshold (0.05).

# This problem we can fix by lowering `lambda.min.ratio` to explore denser networks.

se <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-1, nlambda=10, pulsar.params=pargs)

# We have now fit a network, but since we have only a rough, discrete sampling of networks along the lambda path, 
# we should check how far we are from the target stability threshold (0.05).

getStability(se)
sum(getRefit(se))/2

# To get closer to the mark, we should bump up `nlambda` to more finely sample of the lambda path, which gives a denser network.

se <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-1,
nlambda=100, pulsar.params=pargs)
getStability(se)
sum(getRefit(se))/2

# Save output if wanted
#save(amgut1.filt, depths, amgut1.filt.n, amgut1.filt.cs, d, n, e, graph, Prec, Cor, se, se.mb.amgut, se.gl.amgut, sparcc.amgut, sparcc.graph, ig.mb, ig.gl, 
#     ig.sparcc, file=".README.RData")



################################
# SPRING for quantitative microbiome data
################################
# Reproduce Figure 10 from the Yoon et al. 2019 paper
# Frontiers in Genetics, https://www.frontiersin.org/articles/10.3389/fgene.2019.00516/full

install.packages("remotes")
remotes::install_github("irinagain/mixedCCA")
library(mixedCCA) # For rank-based correlation

library(SpiecEasi) # For SpiecEasi method
library(pulsar) # For subsampling in parallel with StARS criterion
library(huge) # For graph estimation with mb method
library(phyloseq) # For otu_table, $Rank6, $Health_status, etc.

source("SPRINGhelpers.R") # contains SPRING, mclr, hugeKmb functions

load("data/qmphealthyrank6pruned.rdata")
# containing three data variable:
# X: copyadjusted count data
# QMP: quantitative count data (X -> RMP by dividing by total abundance -> QMP by multiplying by cell counts)
# qmphealthy6_only1filt: phyloseq class data.
# codes for the detail steps regarding pruning from QMPphy.RData to 106 healthy samples are in pruneQMP.R

n = dim(QMP)[1]
p = dim(QMP)[2]

rep.num = 50 # the repetition number of subsampling
nlam = 50 # the number of lambda sequence values
lambda.min.ratio = 1e-2 # the ratio of lambda.min over lambda.max
thresh = 0.1 # threshold for StARS criterion
subsample.ratio = 0.8 # subsample size ratio over the total samples
nc = 2 # numbere of cores for subsampling in parallel mode
seed = 10010 # seed for subsampling

# Commented out because of long run times
 
# SpiecEasi on copyadjusted count data
#ptm <- proc.time()
#fit.se <- spiec.easi(X, method = 'mb', sel.criterion = 'stars', lambda.min.ratio = lambda.min.ratio, nlambda = nlam, pulsar.params = list(thresh = thresh, subsample.ratio = subsample.ratio, seed = seed, rep.num = rep.num))
#proc.time() - ptm # 65 seconds

#save(fit.se, file = "data/QMP_SE.rdata")

# SPRING on quantitative count data
#ptm <- proc.time()
#fit.spring <- SPRING(QMP, lambda.min.ratio = lambda.min.ratio, nlambda = nlam, thresh = thresh, subsample.ratio = subsample.ratio, seed = seed, ncores = nc, rep.num = rep.num)
#proc.time() - ptm # 1353 seconds = 23 minutes

#save(fit.spring, file = "data/QMP_SPRING.rdata")

# SPRING on mclr-transformed data
#ptm <- proc.time()
#QMP_mclr <- mclr(QMP) # mclr(QMP) and mclr(X@.Data) are the same.
#fitmclr.spring <- SPRING(QMP_mclr, lambda.min.ratio = lambda.min.ratio, nlambda = nlam, thresh = thresh, subsample.ratio = subsample.ratio, seed = seed, ncores = nc, rep.num = rep.num)
#proc.time() - ptm # 1426 seconds = 24 minutes

#save(QMP_mclr, fitmclr.spring, file = "data/QMP_SPRING_onmclr.rdata")

# Similar to Figure 3 in the paper "Quantitative microbiome profiling..." by Vandeputte et al. in 2017.
# Create the tile figure (Figure 10) from Yoon et al. 2019 paper

library(phyloseq) # extract genera name from phyloseq format data
library(ggplot2) # draw tile figure
library(Matrix) # for triu and tril functions
library(gtable) # to add the box for ordered genera

# Figure 10 in our manuscript.
# SpiecEasi on compositional vs. SPRING on quantitative 

load("data/QMP_SE.rdata") # fit.se
load("data/QMP_SPRING.rdata") # fit.spring

# For StARS criterion
thresh = 0.1

# According to the threshold, find opt.index and save adjacency matrix and coefficient among the path.
opt.SE <- max(which(fit.se$select$stars$summary<thresh))
opt.K <- max(which(fit.spring$output$stars$summary<thresh))

adj.SE <- fit.se$est$path[[opt.SE]]; adj.SE <- as.matrix(adj.SE)
adj.K <- fit.spring$output$est$path[[opt.K]]; adj.K <- as.matrix(adj.K)

Coef.SE <- SpiecEasi::symBeta(fit.se$est$beta[[opt.SE]], mode='maxabs') # isSymmetric(Coef.SE)
Coef.SE <- as.matrix(Coef.SE)
Coef.K <- SpiecEasi::symBeta(fit.spring$output$est$beta[[opt.K]], mode='maxabs') # isSymmetric(Coef.K)
Coef.K <- as.matrix(Coef.K)

# For all genera, table for matching/mismatching signs between SE and SPRING
table(cbind.data.frame(SpiecEasi = as.numeric(sign(Coef.SE[upper.tri(Coef.SE)])), SPRING = as.numeric(sign(Coef.K[upper.tri(Coef.K)]))))

# Extract the genera names
taxainfo <- tax_table(qmphealthy6_only1filt)
genusnamesall <- substring(tax_table(qmphealthy6_only1filt)[, 6], 4)
### Change the genera name with parenthetis to with a number 2.
whereparenthesis <- which(substring(genusnamesall, 1, 1) == "[")
for ( ii in 1:length(whereparenthesis)){
  len <- length(unlist(strsplit(genusnamesall[whereparenthesis[ii]], "")))
  tmp <- substring(genusnamesall[whereparenthesis[ii]], 2, len-1)
  genusnamesall[whereparenthesis[ii]] <- paste(tmp, 2)
}
### Remove genera which do not have names or NAs.
p <- dim(Coef.SE)[1]
genus_NA <- which(is.na(tax_table(qmphealthy6_only1filt)[, 6]))
genus_noname <- which(substr(tax_table(qmphealthy6_only1filt)[, 6], start = 4, stop = 5) == "")
genusname <- (1:p)[-sort(union(genus_NA, genus_noname))]

GenusCoef.SE <- Coef.SE[genusname, genusname]
GenusCoef.K <- Coef.K[genusname, genusname]
QMP_genusname <- QMP[, genusname]
# Save the extracted genera name to their column names and row names
colnames(GenusCoef.SE) <- rownames(GenusCoef.SE) <- genusnamesall[genusname]
colnames(GenusCoef.K) <- rownames(GenusCoef.K) <- genusnamesall[genusname]

# Screen the genera which has at least one partial correlation larger than corthresh = 0.2.
corthresh = 0.2
isover_0.2 <- which(apply(GenusCoef.SE, 1, function(x) any( abs(x) > corthresh)) | apply(GenusCoef.K, 1, function(x) any(abs(x) > corthresh)))

# Order genera based on the total abundance
abundanceord <- order(colSums(QMP_genusname[, isover_0.2]), decreasing = TRUE)
# To check if ordering is correctly done. Should be zero
# sum(colSums(QMP_genusname[, isover_0.2[abundanceord]]) != sort(colSums(QMP_genusname[, isover_0.2[abundanceord]]), decreasing = TRUE))

# Coef matrices for only screened genera (has at least one partial correlation larger than corthresh = 0.2)
GenusCoef.SE_selected <- GenusCoef.SE[isover_0.2[abundanceord], isover_0.2[abundanceord]]
GenusCoef.K_selected <- GenusCoef.K[isover_0.2[abundanceord], isover_0.2[abundanceord]]


# Create a matrix to draw tile figure.
# SE is in lower matrix and SPRING in upper matrix.
TileMat <- matrix(0, ncol = length(isover_0.2), nrow = length(isover_0.2))
TileMat[lower.tri(TileMat)] <- as.matrix(GenusCoef.SE_selected)[lower.tri(GenusCoef.SE_selected)]
TileMat[upper.tri(TileMat)] <- as.matrix(GenusCoef.K_selected)[upper.tri(GenusCoef.K_selected)]

# create a coordinate vectors for geom_tile
coord <- cbind(rowv = rep(rownames(GenusCoef.K_selected), each = length(isover_0.2)), 
               colv = rep(rownames(GenusCoef.K_selected), length(isover_0.2)))
df_Tilemat <- cbind.data.frame(coord, partialcorr = c(TileMat))

# rows start from bottom to top and cols start from left to right.
# but we want start from the left top corner for both triangles. So make colv factor in reverse level.
df_Tilemat$rowv <- factor(df_Tilemat$rowv, levels = rownames(GenusCoef.K_selected))
df_Tilemat$colv <- factor(df_Tilemat$colv, levels = rev(rownames(GenusCoef.K_selected)))

### disagreement table
# For genera (which have at least one partial correlation larger than corthresh = 0.2), table for matching/mismatching signs between SE and SPRING
table(cbind.data.frame(SpiecEasi = as.numeric(sign(GenusCoef.SE_selected[upper.tri(GenusCoef.SE_selected)])),
                       SPRING = as.numeric(sign(GenusCoef.K_selected[upper.tri(GenusCoef.K_selected)]))))

#               SPRING
# SpiecEasi  -1   0   1
#         -1   1   5   0
#         0    2 440  16
#         1    0  13  19

library(xtable)
xtable(table(cbind.data.frame(SpiecEasi = as.numeric(sign(GenusCoef.SE_selected[upper.tri(GenusCoef.SE_selected)])),
                              SPRING = as.numeric(sign(GenusCoef.K_selected[upper.tri(GenusCoef.K_selected)]))))
)

p02 <- length(isover_0.2)
tile1 <- ggplot(df_Tilemat, aes(x = rowv, y = colv), fill = partialcorr) +        ## global aes
  geom_tile(aes(fill = partialcorr), color = "black", alpha = 1, width = 1, height = 1) +         ## to get the rect filled
  scale_fill_gradient2(high = "blue", mid = "white", low = "red", name = "Partial Correlation") +       ## color of the corresponding aes
  coord_cartesian(xlim = c(1, p02), ylim = c(1, p02), clip = 'off') + 
  theme_bw() +
  theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = c(-0.2, 1.25), legend.justification="left", legend.direction="horizontal", plot.margin = margin(1, 1, 1, 0, "cm")) + 
  annotate("text", x = (p02+2)/2, y = -0.7, label = "SpiecEasi") + 
  annotate("text", x = p02+1.5, y = (p02+2)/2, angle = 90, label = "SPRING on Quantitative counts") + 
  scale_x_discrete(position = "top") + 
  scale_alpha_continuous(guide=FALSE) +  scale_size_continuous(guide=FALSE) + 
  geom_segment(aes(x=0.5, y=(p02+0.5), xend = (p02+0.5), yend=0.5), color="grey", size = 0.5)

# Uncomment the following two lines when only you want to update the figure.
# pdf(file = "Tile_SE_vs_SPRING.pdf", width = 8, height = 8)
tile1
# dev.off()




