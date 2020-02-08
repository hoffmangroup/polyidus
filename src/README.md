# Bellerophon test run

You can download the example data for
running Bellerophon from https://www.pmgenomics.ca/hoffmanlab/proj/bellerophon/bellerophon-data-v1.tar.gz .



```
usage: Identify viral integration sites in FASTQ files.
       Requires a host reference genome index,
       a viral reference genome index, and the FASTQ file(s)
       of the experiment to identify the exact
       integration sites.
       Make sure you have generated genome indices for
       both the host and virus genomes.
       You can use bwa instead of bowtie2 as well.
       Make sure the aligner (bowtie2 or bwa), samtools,
       and bedtools also exist in $PATH.
       [-h] [--fastq [FASTQ [FASTQ ...]]] [--outdir OUTDIR]
       [--aligner {bwa,bowtie2}]
       hostindex viralindex

positional arguments:
  hostindex             Path to host index prefix (either bwa or bowtie2)
  viralindex            Path to viral index prefix (either bwa or bowtie2)

optional arguments:
  -h, --help            show this help message and exit
  --fastq [FASTQ [FASTQ ...]]
                        Path to fastq files (if more than 1, assumes second is
                        the second pair.
  --outdir OUTDIR       Path to output folder
  --aligner {bwa,bowtie2}
                        Choose from bwa or bowtie2 (default)

Citation: Karimzadeh M., Arlidge C., Rostami A., Lupien M., Bratman S., and
Hoffman M. M. Integration of human papillomavirus into the human genome
promotes cancer by modifying local chromatin and transcription.


VIRALINDEX=../data/hpv16/hpv16_bowt_ind
HOSTINDEX=../data/hg38/hg38_bwt2_index
FASTQFILES=(../data/fastqfiles/SiHa_R1.fastq.gz ../data/fastqfiles/SiHa_R2.fastq.gz)
OUTDIR=../data/bellerophonOutput
mkdir -p $OUTDIR
echo -e '#!/bin/sh' > testRun.sh
echo "python bellerophon.py $HOSTINDEX $VIRALINDEX --fastq ${FASTQFILES[@]} --outdir $OUTDIR" >> testRun.sh
sbatch -c 1 -p hoffmangroup --mem=16G -t 36:00:00 -e bellerophon.ERR -o bellerophon.LOG testRun.sh
```

# Interpreting the results

The output file exactHpvIntegrations.tsv contains the exact information of where the integration site starts
in the virus as well as the host.
The fourth column IntegrationSite shows the exact position in the host
that the integration occurs.
Second and third columns show a 10,000 bp window upstream of the integration host in the host genome
that is not affected by that integration.


The eigth column ViralIntegrationSite shows the exact nucleotide in the virus with integration.
Sixth and 7th columns indicate a 10,000 bp window in the virus genome that is not affected by that integration.
Since the 6th and 7th column are calculated by adding or subtracting the value of 10,000, we assume the virus has a circular genome.
Therefore if for a virus of length 7,904 bp an integration occurs at position 879, the start site of -9121 indicates that positions
upstream of 879 are not affected by that integration.
A start site of -2,975 (7,904 - 10,879), however, will indicate that the positions downstream of 879 are not affected by that integration.


