Run the MEGAHIT assembler
=========================
Tutorial by Titus Brown, Adina Howe, and Jihoon Yang

`MEGAHIT <https://github.com/voutcn/megahit>`__ is a very fast, quite
good assembler designed for metagenomes.

First, install it::

   cd
   conda install -c bioconda megahit 

When asked to proceed, type 'y'.  MEGAHIT will then proceed to install and store its commands in such a way that you can access it.

Try typing::

   megahit


You will now see the version and command line options for MEGAHIT.  Note that MEGAHIT is looking for paired-end reads and a specified output directory.  

Now, download some data to try on an assembly ::

   cd
   mkdir data
   #Mike's figshare get command here"

Now, run the assembler! ::

   cd
   mkdir assembly
   cd assembly
   mv ~/data/example*fastq .
   megahit --12 example.R1.fastq, example.R2.fastq

This will take about 5-10 minutes; at the end you should see output like
this::

   ... 4930 contigs, total 9642451 bp, min 201 bp, max 212315 bp, avg 1955 bp, N50 7771 bp
   ... ALL DONE. Time elapsed: 300 seconds

The output assembly will be in ``megahit_out/final.contigs.fa``.

After the assembly is finished
------------------------------

At this point we can do a bunch of things:

* annotate the assembly;
* evaluate the assembly's inclusion of k-mers and reads;
* set up a BLAST database so that we can search it for genes of interest;
* bin the contigs in the assembly into species bins;

Annotation with BLAST
---------------------

Install BLAST ::

   curl -O ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.9.0+-x64-linux.tar.gz
   tar -xvf ncbi-blast+-2.9.0-src.tar.gz
   sudo cp ncbi-blast-2.9.0+/bin/makeblastdb /usr/local/bin/

Let's look at two blast commands ::

   makeblastdb -h
   blastn -h

First, let's make are blast index database file from our reference database, ref.fa ::

   mkdir blast
   cd blast
   cp ~/data/ref.fa .
   cp ~/data/*fasta .
   makeblastdb -in ref.fa -dbtype nucl

Let's run a blast ::

   blastn -db ref.fa -query megahit_out/final.contigs.fa 

AHHHHH!  Press control-C to break that carnage!

Let's try this again, where now we specify a specific output type and filename. ::

   blastn -db ref.fa -query megahit_out/final.contigs.fa -outfmt 6 -out contigs.x.ref.blastnout

Let's take a look at this file together.  First, this is a nice `key <http://www.metagenomics.wiki/tools/blast/blastn-output-format-6>`_ for BLAST outputs in tabular format.

What can we see in this output?

**Note, use BLAST settings with `[caution] <https://academic.oup.com/bioinformatics/article/35/9/1613/5106166>`_

Checking against References
---------------------------

Let's say that you have some reference genomes to compare to and would like to know the rperesentation of the contigs in these genomes.  These three genomes are NC_003112.2, NC_004310.3, NC_009089.1 -- three well known pathogens.  Let's format each genome into a blast database and run three different blasts against our assembly.  

Wait...this seems like a lot of BLASTING!

Let's give this a try ::

   for x in NC*fasta; do echo $x; done

In this case, the echo command just "echos" the variable x within the for-loop.  We can change that command to "blastn".

   for x in NC*fasta; do makeblastdb -in $x -dbtype nucl; done
   
You should have 3 indexed databases for each reference file.  How could you do a blast of each database against the contigs?  Give it a try.

Here's the solution, with an extra parameter to speed things up :: 

   for x in NC*fasta; do blastn -db $x -query megahit_out/final.contigs.fa -outfmt 6 -out contigs.x.$x.blastnout -num_threads 6; done

This is a real pro-trick in doing bioinformatics swiftly!

   














