Polyidus --- Identifying viral integration sites from chimeric sequencing reads
==================================================================================


Introduction
------------

**The free Polyidus software identifies the exact genomic regions for integration of
a known virus.**


We developed Polyidus to identify viral integration sites with chimeric sequencing reads from any paired-end sequencing data.
First, Polyidus aligns reads to a viral genome.
It allows for partial mapping using local alignment, and removes any sequencing fragment where neither read maps to the virus.
Second, Polyidus aligns the selected reads to the host genome, permitting partial mapping.
Third, Polyidus identifies *chimeric reads*: those reads mapped partially to the host genome and partially to the virus genome.
Fourth, for each chimeric read, Polyidus reports the start and strand of integration in both the host and viral genomes.
Polyidus also reports the number of chimeric reads supporting each integration site.


Polyidus finds highly confident integration sites which contain chimeric sequencing reads.
Other methods perform the first two steps in reverse order, resulting in slower performance.
While some previous methods also align to the virus first, either the software no longer appears available where specified at publication, or they use BLAST instead of a faster short read aligner.
Unlike ViFi, Polyidus requires that the chimera match an existing viral genome reference.
Polyidus does not use non-chimeric fragments where one read maps entirely to host and one read maps entirely to virus genome.


Polyidus uses Bowtie2 (version 2.2.6) and vastly speeds up integration site finding.
Polyidus identified integration sites at an average of 8 core-hours on a 2.6 GHz Intel Xeon CPU E5-2650 v2 processor and 4 GB of RAM for whole genome sequencing data.
Previous methods require an average of 400 CPU core-hours.


Quick start
-----------

You can run the following commands to run Polyidus on a small dataset.
First clone the repository and download the example dataset (123 MB)::

    git clone https://github.com/hoffmangroup/polyidus.git
    cd polyidus
    wget https://www.pmgenomics.ca/hoffmanlab/proj/bellerophon/bellerophon-data-v1.tar.gz
    tar -xvf bellerophon-data-v1.tar.gz


Then, make sure you have bowtie2 (v2.2.6), samtools (v1.9), pysam (v0.8.4), and bedtools (v2.26.0) installed.
You need to run the following commands::

    mkdir data/polyidusOutput
    cd src
    python polyidus.py ../data/hg38/hg38_bwt2_index ../data/hpv16/hpv16_bowt_ind --fastq ../data/fastqfiles/SiHa_R1.fastq.gz ../data/fastqfiles/SiHa_R2.fastq.gz --outdir ../data/polyidusOutput


If you have any issues with using the software, please post it to GitHub issues and do not contact us directly.


If you prefer not to post on issues, you can contact the maintainer with email: mehran.karimzadehreghbati@mail.utoronto.ca
