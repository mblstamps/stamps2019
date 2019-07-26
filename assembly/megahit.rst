Run the MEGAHIT assembler
=========================

`MEGAHIT <https://github.com/voutcn/megahit>`__ is a very fast, quite
good assembler designed for metagenomes.

First, install it::

   cd
   conda install -c bioconda megahit 

When asked to proceed, type 'y'.  MEGAHIT will then proceed to install and store its commands in such a way that you can access it.
Try typing::

   megahit


You will now see the version and command line options for MEGAHIT.  Note that MEGAHIT is looking for paired-end reads and a specified output directory.  

Now, download some data to try on an assembly::


----

Now, run the assembler! ::

   mkdir ~/assembly
   cd ~/assembly
   ln -fs ../data/*.subset.pe.fq.gz .

   ~/megahit/megahit --12 example.R1.fastq, example.R2.fastq

This will take about 5-10 minutes; at the end you should see output like
this::

   ... 4930 contigs, total 9642451 bp, min 201 bp, max 212315 bp, avg 1955 bp, N50 7771 bp
   ... ALL DONE. Time elapsed: 300 seconds

The output assembly will be in ``megahit_out/final.contigs.fa``.

After the assembly is finished
------------------------------

At this point we can do a bunch of things:

* annotate the assembly (:doc:`prokka_tutorial`);
* evaluate the assembly's inclusion of k-mers and reads;
* set up a BLAST database so that we can search it for genes of interest;
* quantify the abundance of the contigs or genes in the assembly, using the original read data set (:doc:`salmon_tutorial`);
* bin the contigs in the assembly into species bins;


