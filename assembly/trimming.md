
# Metagenome quality trimming
Authored by Jin Choi for EDAMAME

***
EDAMAME tutorials have a CC-BY [license](https://github.com/edamame-course/2015-tutorials/blob/master/LICENSE.md). _Share, adapt, and attribute please!_
***

## Overarching Goal
* This tutorial will contribute towards an understanding of **microbial metagenome analysis**

## Learning Objectives
* Assess the quality of "raw" metagenome data
* Trim raw reads to meet quality standards

***


## Background On Quality Trimming

When DNA sequencing takes place, errors are inevitable.  No sequencing method is perfect and some are drastically different than others in regard to sequence length and quality.  We're going to trim the poor quality tail end sections of our Illumina reads.  Although in general Illumina reads are very high quality, this degradation at the end of the sequencing run is typical of the Illumina sequencing platforms.

Some sequencing centers will remove library adapters (our sequencing center does), but you will have to check with your data provider to know what they give you and ALWAYS check for your self to verify what you have been told.

As always, you want to make sure you read the manual of any tool to be sure you know what the tool is doing to your data and if that is the right tool for the job.  Knowing which tool to use is very important -- you wouldn't use a saw to put a nail in a piece of wood, would you?

## Quality Trimming Your Sequence Data

The first step to starting an assembly or metagenomic read analysis is to remove bad quality sequences, a process we call quality trimming.  We'll start a cloud compute instance and practice trimming bad quality reads.

*Start your STAMPS instance.*

Install Trimmomatic - a program to quality trim reads
```
cd 
curl -O http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip
unzip Trimmomatic-0.38.zip
```

Download the data (You did on last session): The tutorial data is from [Sharon et al. 2013](http://www.ncbi.nlm.nih.gov/pubmed/22936250); it’s two data points from an infant gut sample. And it is a subsampled file so that we can do this tutorial in some reasonable amount of time, 100,000 sequences. Go to the home directory, make a directory named `metagenome`, change directores into the folder, download data, then unzip the data or decompress the data file. 
```
cd
mkdir metagenome
cd metagenome
wget https://s3.amazonaws.com/edamame/infant_gut.sub.tar.gz
tar -zxvf infant_gut.sub.tar.gz
```
This dataset contains paired end reads.  The paired end nature of datasets can be a huge advantage for assembly -- we know the approximate insert size between these reads.  However, when we trim bad quality sequences and/or reads, sometimes a previously paired read is left as an orphan (sad face).  However, for downstream analyses, we want to keep paired end reads that are together and separate these from single ended reads.  This is a bit of a process but important to understand.    

What's the quality of the dataset?  We are going to use the program [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) to both remove reads with a low quality threshold and remove sequencing adapters (non-biological sequences) prior to assembly.

For this dataset, we happen to know that our adapters are those in TruSeq3-PE-2.fa

Now, let's get to trimming.  To shorten this tutorial, let's trim only a few files to exemplify how we would do trimming.

```
java -jar ../Trimmomatic-0.38/trimmomatic-0.38.jar PE SRR492065_1.sub.fastq.gz SRR492065_2.sub.fastq.gz SRR492065.r1.pe SRR492065.r1.se SRR492065.r2.pe SRR492065.r2.se ILLUMINACLIP:../Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

If you look at the Trimmomatic manual, you'll see that this command performs the following:
This will perform the following:

Remove adapters
 
Remove leading low quality or N bases (below quality 3) (LEADING:3)

Remove trailing low quality or N bases (below quality 3) (TRAILING:3)

Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)

Drop reads below the 36 bases long (MINLEN:36)


A further breakdown of the command is here:

java: run java program, -jar: run jar program, /usr/local/bin/trimmomatic-0.36.jar: name of the program with path, PE: paired-end, SRR492066_1.sub.fastq.gz: first pared-end, SRR492066_2.sub.fastq.gz: second paired-end, s1_pe: output of first file(paired-end), s1_se: output of first file(single-end), s2_pe: output of second file(paired-end), s2_se: output of second file(single-end), ILLUMINACLIP: use illumina clip, /usr/local/share/adapters/TruSeq2-PE.fa: adapter file with path, 2:30:10 : <seed mismatches(specifies the maximum mismatch count which will still allow a full match to be performed)>:<palindrome clip threshold(specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment.)>:<simple clip threshold(specifies how accurate the match between any adapter etc. sequence must be against a read.)> [Here more detail](http://www.usadellab.org/cms/?page=trimmomatic)

```

Some handy quality and/or adapter trimming tools you might want to investigate are:   
   1. [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) - all purpose
   2. [cutadapt](https://code.google.com/p/cutadapt/) - adapter trimming
   3. [sickle](https://github.com/najoshi/sickle) - read quality trimming
   4. [scythe](https://github.com/vsbuffalo/scythe) - adapter contamination trimming


