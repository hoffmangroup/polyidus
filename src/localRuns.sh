## ATAC-seq run -- SiHa
```
usage: Bellerophon identifies viral integration sites.
       The script requires a host reference genome index,
       a viral reference genome index, and the fastq file(s)
       of the experiment to identify the exact integration sites.
       Bellerophon uses bowtie2.Make sure you have generated genome
       indices using bowtie2-build and bowtie2 exists in $PATH.
       You can use bwa instead of bowtie2 as well.

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


```
