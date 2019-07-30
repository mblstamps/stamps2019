Metagenome Workflow
=============

Today we will execute a *de novo* shotgun metagenome 
workflow. What makes a good metagenome workflow? As with 
everything we've covered so far, it depends :) 

But we have to start somewhere! We find that it's 
generally a good idea to do a "default" workflow. This allows you to start getting to know your data quickly, and to get **a result** quickly. After you have results from the default workflow, you can critically assess whether they make sense in light of your biological problem or whether you need to use different methods to answer your question. You can also measure attrition of information and determine if you need to do something about it. 

There are a lot of steps in a default *de novo* metagenome workflow. Below we have a diagram of each step one would genreally want to take. 

![](https://i.imgur.com/yGfZKpw.png)

There are 17 steps! Each of which requires a different tool. That's a lot of tools! 

Instead of orchestrating these 17 tools ourselves, we can put an automated workflow to work for us to do this. We chose to use the [Atlas](https://peerj.com/preprints/2843/) workflow because...it works. Which is actually a sort of high bar for stringing 17 tools together!

## Instance hygeine

Metagenomics workflows can take up quite a bit of space. We want to make sure that we have enough space and compute power on our instances to execute the whole thing. 

We need three things to run a workflow: harddrive space, ram, and cpus. 
+ **hard drive space** is phyiscal, persistent storage. We write our files to disk for long term storage
+ **RAM** is working memory. It temporarily stores all the information that your computer needs to work right that minute. Some tools need to load all of your data into RAM in order to work with it. Assembly is one of these steps -- to assemble something, you need at minimum the amount of RAM as you have in gigabytes of data.
+ **CPUs** are the electronic circuitry within a computer that carries out the instructions of a computer program by performing the basic arithmetic, logic, controlling, and input/output (I/O) operations specified by the instructions [article here](https://en.wikipedia.org/wiki/Central_processing_unit).

Log into your jetstream instance using `ssh` or the jupyter lab page.

First, let's check the amount of hard drive space we have. From your bash prompt, run:

```
df -h .
```

You should see an output like this:

```
Filesystem      Size  Used Avail Use% Mounted on
/dev/sda1        58G   26G   33G  44% /
```

You need at least 10G of space on your instance. If you don't have that much, put up a pink sticky and a helper will come around to help remove files. 

Let's also take a look at how much ram we have available:

```
free -mh
```

and you'll see something like this:
```
              total        used        free      shared  buff/cache   available
Mem:            15G        728M         12G         12M        2.0G         14G
Swap:            0B          0B          0B
```

And the amount of CPUs:

```
nproc
```

Which should return:

```
6
```

Ok! We now have an idea of how much compute is available to us. In general, it's really hard to know how much of each resource you need for each tool before you run it. You may be able to google around for this information, but you will also develop this intuition over time as you run more tools. We have enough compute to execute our tutorial because we are working with a subsampled dataset. 

The last bit of instance hygiene we need to take care of is making sure that we all have the right channels, or package managers, enabled for conda. We'll go over this in detail, but essentially if our computers don't know where to look to install software, or if they look at the wrong location first, things will go awry. 

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

## Getting started with ATLAS

We'll organize ourselves as if we're starting a brand new project. 

```
mkdir atlas_workflow
cd atlas_workflow/
```

We also need to install atlas to get this metagenome party started. We'll do that with this code, and we'll go over what it all means in our evening session.

```
conda create -y -n atlas
conda activate atlas
conda install -y -c bioconda metagenome-atlas=2.1.4
```

Ok! We now have atlas installed. Let's checked that our installation worked by running atlas.

```
atlas
```

You should see a help message that looks like this:

```
Usage: atlas [OPTIONS] COMMAND [ARGS]...

  ATLAS - workflows for assembly, annotation, and genomic binning of
  metagenomic and metatranscriptomic data.

  For updates and reporting issues, see: https://github.com/metagenome-
  atlas/atlas

Options:
  --version   Show the version and exit.
  -h, --help  Show this message and exit.

Commands:
  download  download reference files (need ~50GB)
  init      prepare configuration file and sample table for atlas run
  run       run atlas main workflow
```

## Retrieving the data

We will be working with a human stool shotgun metagenome dataset from the Human Microbiome Project. The accession number is SRS049959. This sample has about ~8 GB of data, and it would take too long to run atlas on the entire dataset (~2 days).

For this tutorial, we subsampled the SRS049959 to the first 1 million paired-end reads, and we will be working with this subsampled data. 

First, let's download our data.

```
mkdir data
wget -O data/SRS049959_R1.fq.gz https://osf.io/a5ekr/download
wget -O data/SRS049959_R2.fq.gz https://osf.io/daxze/download 
```

We can check that the download was successful by listing the files with `ls`:

```
ls -lh data
```

We see that we now have about 160 MB of data.

```
total 163M
-rw-rw-r-- 1 stamps19 stamps19 81M Jul 29 11:28 SRS049959_R1.fq.gz
-rw-rw-r-- 1 stamps19 stamps19 82M Jul 29 11:28 SRS049959_R2.fq.gz
```

## Configuring atlas!

We now need to tell atlas where our data is. We do this with the `atlas init` command. 

```
atlas init data/ --assembler megahit --data-type metagenome
```

In running this command, atlas has created a configuration file for itself using the path of our data. Let's take a look at this "config" file to see what decisions atlas is making under the hood. Note, here we are choosing to use megahit because it uses substantially less memory that (meta)spades.

```
less config.yaml
```

The first thing we see is that atlas would prefer that we have 32 GB of RAM for most tools, and 250 GB of RAM for assembly, and more than 6 threads for assembly. Sorry atlas, not gonna happen!!! We need to edit this file to set atlas straight about our resources. We will use nano to do this.

```
nano config.yaml
```

We need to edit these numbers so they reflect our actual resources. Make sure the top of your file looks like this:

```
########################
# Execution parameters
########################
# max cores per process
threads: 6
# Memory for most jobs especially from BBtools, which are memory demanding
java_mem: 16
# can be a subset of threads or altered if rule run_spades or run_megahit are being defined$
assembly_threads: 6
# in GB
assembly_memory: 16
```

To exit, type `ctrl X` type `y`, and press `enter`. This will write the changes we made to our file. 

Let's also look at the samples file and make sure that atlas got it right. 

```
cat samples.tsv
```

Everything looks good here! Atlas successfully realized that we are working with one sample with paired end reads.

## Running atlas!

Now that everything is set up, we could type `atlas run` and it would orchestrate everything for us. We're not going to do that though. We're going to run each step of the pipeline in pieces and look at the output. We're also going to skip the database download step, because it takes too much disk space and RAM to build the databases for annotation. 

### Quality Control

Let's start with quality control. We'll use the `-n` flag first, which tells atlas to run a `--dryrun`. This tells you everything that is about to be run but doesn't actually run it. (It also checks to be sure that your config files are formatted properly!)

```
atlas run qc -n
```

We see a bunch of green and gold output. This is good! Let's take a closer look at what this means:

```
Job counts:
	count	jobs
	1	apply_quality_filter
	1	build_decontamination_db
	1	build_qc_report
	1	calculate_insert_size
	1	combine_insert_stats
	1	combine_read_counts
	1	combine_read_length_stats
	1	deduplicate_reads
	1	finalize_sample_qc
	5	get_read_stats
	1	initialize_qc
	1	qc
	1	qcreads
	1	run_decontamination
	1	write_read_counts
	19
```

Let's run the QC command while we discuss what is happening (it takes about 10 minutes to run):

```
atlas run qc
```

Atlas is automated using a tool called [Snakemake](https://snakemake.readthedocs.io/en/v5.5.4/). We love snakemake in our lab :) - it does the thing, and orchestrates all the things for you. If this sounds appealing to you, here is a [tutorial](https://angus.readthedocs.io/en/2019/snakemake_for_automation.html) on how to use snakemake.

We see that 19 different steps are wrapped by the QC command. These steps control contamination, run quality control (erroneous k-mers, adapters, etc.), and generate reports about the different steps that are run. As each rule runs, we see green snakemake text alerting us to which step is running, and we see white text from the stdout ("standard output") from each tool that is run. Some of the longest steps are downloading and installing all of the software that's needed - it's using conda to do this, note!

When atlas is finished with QC, we'll see:

```
ATLAS finished
The last rule shows you the main output files
```

If we list the files with `ls -FC`, you'll see we have many more files and directories:

```
SRS049959    data       logs  reports      stats
config.yaml  databases  ref   samples.tsv
```

Our results are located in the `SRS049959` folder.

```
ls SRS049959/
```

And looking inside the folder in this folder:
```
ls SRS049959/sequence_quality_control
```
We see:
```
SRS049959_QC_R1.fastq.gz                       contaminants
SRS049959_QC_R2.fastq.gz                       finished_QC
SRS049959_QC_se.fastq.gz                       read_stats
```

These are our quality controlled reads. We also see a file that says `finished_QC`. This is how atlas alerts us that it finished running all quality control steps.

#### Challenge
Look inside the contaminants folder. What contaminant was removed by atlas?

### Assembly

QC. Done. Yay! Now that our reads are cleaned up and we have removed contaminants, let's assemble. Again, let's do a dry run first to see what steps atlas will run. 

```
atlas run assembly -n
```

```
Job counts:
        count   jobs
        1       align_reads_to_final_contigs
        1       align_reads_to_prefilter_contigs
        1       assembly
        1       assembly_one_sample
        1       bam_2_sam_contigs
        1       build_assembly_report
        2       calculate_contigs_stats
        1       convert_sam_to_bam
        1       error_correction
        1       filter_by_coverage
        1       finalize_contigs
        1       get_contigs_from_gene_names
        1       init_pre_assembly_processing
        1       merge_pairs
        1       merge_se_me_for_megahit
        1       passtrough_se_merged
        1       pileup
        1       pileup_prefilter
        1       predict_genes
        1       qc
        1       rename_contigs
        1       rename_megahit_output
        1       run_megahit
        24
```

Lots of things! Let's start the assembly and then we'll talk through these steps (takes ~12 min).

```
atlas run assembly
```

Megahit is the meat and potatoes of this pipeline. It orchestrates the assembly of our reads, generating long(er) contigs. (We ran it during the initial assembly tutorial, too.)

After the assembly is finished, the reads are mapped back to the assembly. This is good for quality control of our assembly (does our assembly fully capture the information in our reads?), and helps with downstream processes like binning.

The first step of annotation also occurs during the assembly phase. Here, the open reading frames (ORFs) are predicted. 

After these processes finish, we can take a look at the assembly and associated output files

```
ls SRS049959/
```

In this directory, we see the file `finished_assembly`, as well as new folders for annotation, assembly, and sequence alignment. 

We can check the size of our assembly

```
ls -lh SRS049959/assembly/
```

And take a look at the assembly itself

```
less SRS049959/assembly/SRS049959_final_contigs.fasta
```

Tons of contigs, most of them babies (#babycontigs).

**Challenge 1:** use a combination of `grep` and `wc` to count the number of contigs in our final assembly. 

**Challenge 2:** explore the other output files in the assembly directory and subfolders. Find the file with more assembly stats and view it with `less`.

### Binning


Next comes binning. Woohoo! 


```
atlas run binning -n
```

```

Job counts:
        count   jobs
        1       assembly
        1       bam_2_sam_binning
        1       binning
        1       build_bin_report
        1       combine_coverages
        1       convert_concoct_csv_to_tsv
        1       download_atlas_files
        1       get_bins
        1       get_contig_coverage_from_bb
        1       get_maxbin_cluster_attribution
        1       get_metabat_depth_file
        3       get_unique_bin_ids
        3       get_unique_cluster_attribution
        1       initialize_checkm
        1       maxbin
        1       metabat
        1       pileup_for_binning
        1       qc
        1       run_checkm_lineage_wf
        1       run_checkm_tree_qa
        1       run_concoct
        1       run_das_tool
        1       unpack_checkm_data
        27
```

Let's get this running as we break down what's actually happening.

```
atlas run binning
```

Atlas does both binning and bin evaluation. We can see from the dry run that atlas uses many binners -- maxbin, metabat, concoct, and DAS tool. DAS tool combines the output of the other three binning tools, and summarizes their overlaps.

For bin evaluation, atlas runs CheckM. CheckM also assigns taxonomy, so a little bit of annotation happens in this step as well.

After a bit of time passes, we will see this error on the screen:

```
*******************************************************************************
 [CheckM - tree] Placing bins in reference genome tree.
*******************************************************************************

  Identifying marker genes in 1 bins with 6 threads:
    Finished processing 1 of 1 (100.00%) bins.
  Saving HMM info to file.

  Calculating genome statistics for 1 bins with 6 threads:
    Finished processing 1 of 1 (100.00%) bins.

  Extracting marker genes to align.
  Parsing HMM hits to marker genes:
    Finished parsing hits for 1 of 1 (100.00%) bins.
  Extracting 43 HMMs with 6 threads:
    Finished extracting 43 of 43 (100.00%) HMMs.
  Aligning 43 marker genes with 6 threads:
    Finished aligning 43 of 43 (100.00%) marker genes.

  Reading marker alignment files.
  Concatenating alignments.
  Placing 1 bins into the genome tree with pplacer (be patient).
Killed
Uncaught exception: Sys_error("SRS049959/binning/DASTool/checkm/storage/tree/concatenated.pplacer.json: No such file or directory")
Fatal error: exception Sys_error("SRS049959/binning/DASTool/checkm/storage/tree/concatenated.pplacer.json: No such file or directory")
```

What happened? It looks as though `pplacer` failed. We can copy the relevant parts of the error message and google them to see if anyone has encountered these problems before. If you google:

```
Uncaught exception: Sys_error("SRS049959/binning/DASTool/checkm/storage/tree/concatenated.pplacer.json: No such file or directory")
Fatal error: exception Sys_error("SRS049959/binning/DASTool/checkm/storage/tree/concatenated.pplacer.json: No such file or directory")
```

You might come upon [this](https://github.com/Ecogenomics/CheckM/issues/37) github issue. 

The key turns out to be that `Killed` message in the output above, which generally tells you that the computer ran out of memory at this step. This is because we "only" have 16 GB of RAM on these computers, and pplacer requires more.

## Other

### Annotation

We will be skipping these sections because they require large databases that take ~100 GB of ram to build. However, on your samples, atlas would run gene annotation for each bin. 

## Other references

Don't like this pipeline? There are others that run in a very similar way! This [paper](https://www.biorxiv.org/content/biorxiv/early/2018/05/18/326363.full.pdf) does a great job of outlining different options. 

Want all the data? It's available from hmp dacc. 

