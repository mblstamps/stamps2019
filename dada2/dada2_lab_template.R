######################################################################
#####  GOAL: Choose and validate appropriate dada2 parameters    #####
#####         for processing the dataset in dada2/data/lab       #####
######################################################################
#####   DATA is taken from the following recent paper            #####
######################################################################
# Pyrethroid exposure alters internal and cuticle surface bacterial communities in Anopheles albimanus, ISMEJ, 2019.
# https://doi.org/10.1038/s41396-019-0445-5
# Sequencing data: First 18 samples from https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA512122
######################################################################

######################################################################
#####  IF DOING THIS ON THE CLOUD, change multi=TRUE to multi=2  #####
######################################################################

library(dada2); packageVersion("dada2") # 1.12 or later

path <- "XXX" # REPLACE XXX with the path to the directory the fastqs are in

# Read in forward and reverse fastq filenames
fnFs <- list.files(path, pattern="_1.fastq.gz", full.names=TRUE)
fnRs <- list.files(path, pattern="_2.fastq.gz", full.names=TRUE)

plotQualityProfile(fnFs)
plotQualityProfile(fnRs)
# Good quality/bad quality? What might make sense as truncation length parameters?

# Define filtered filenames
filtFs <- file.path(path, "filtered", basename(fnFs)) 
filtRs <- file.path(path, "filtered", basename(fnRs))

###################################################################
######  Are primers on these reads that need to be removed?  ######
######  How long is the sequenced amplicon? SEE PAPER        ######
###################################################################

track <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxEE=2, 
                       trimLeft=c(XXX, YYY), # REPLACE XXX/YYY with proper parameter choices
                       truncLen=c(XXX, YYY), # REPLACE XXX/YYY with proper parameter choices
                       multi=TRUE, verbose=TRUE)
# Were most reads retained during filtering? If not, would any parameter tweaks help?

errF <- learnErrors(filtFs, multi=TRUE)
errR <- learnErrors(filtRs, multi=TRUE)

plotErrors(errF)
plotErrors(errR)
# Do the fitted error models look reasonable?

ddF <- dada(filtFs, errF, pool="pseudo", multi=TRUE)
ddR <- dada(filtRs, errR, pool="pseudo", multi=TRUE)
# What pooling option makes sense? FALSE (default), "pseudo", or TRUE? See ?dada for more info
# For more pseudo-pooling detail: https://benjjneb.github.io/dada2/pseudo.html#pseudo-pooling

mm <- mergePairs(ddF, filtFs, ddR, filtRs, verbose=TRUE)
# Were most reads retained during merging. If not, why not?

sta <- makeSequenceTable(mm)
st <- removeBimeraDenovo(sta, multi=TRUE, verbose=TRUE)
sum(st)/sum(sta)
# Were most reads retained during chimera removal?

# Download rdp_train_set_14.fa.gz from here: https://benjjneb.github.io/dada2/training.html
tax <- assignTaxonomy(st, "~/tax/rdp_train_set_14.fa.gz", multi=TRUE)
unname(head(tax))
# Are reasonable taxonomic assignments being made to the abundant taxa?

library(phyloseq); library(ggplot2)
# Very simple data.frame, defining the 3 groups this data comes from (see paper for more info)
samdf <- data.frame(row.names=rownames(st), 
                    sampleID=rownames(st), 
                    treatment=rep(c("Unexposed", "Susceptible", "Resistant"), each=6))

ps <- phyloseq(sample_data(samdf), otu_table(st, taxa_are_rows=FALSE), tax_table(tax))
plot_ordination(ps, ordinate(ps, method="NMDS", distance="bray"), color="treatment") + 
  aes(size=4) + theme_bw()
# This is a small subset of the total data in this study...
