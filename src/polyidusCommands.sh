bowtie2 -p 1 --local -x ../data/hpv16/hpv16_bowt_ind -1 ../data/fastqfiles/SiHa_R1.fastq.gz -2 ../data/fastqfiles/SiHa_R2.fastq.gz | samtools sort -o ../data/bellerophonOutput/viral/virusAligned-temp.bam -
samtools view -bS -f 4 -F 264 ../data/bellerophonOutput/viral/virusAligned-temp.bam -o ../data/bellerophonOutput/viral/pair1mapped.bam
samtools view -bS -f 8 -F 260 ../data/bellerophonOutput/viral/virusAligned-temp.bam -o ../data/bellerophonOutput/viral/pair2mapped.bam
samtools view -bS -f 1 -F 12 ../data/bellerophonOutput/viral/virusAligned-temp.bam -o ../data/bellerophonOutput/viral/bothPairsMapped.bam
samtools merge - ../data/bellerophonOutput/viral/pair1mapped.bam ../data/bellerophonOutput/viral/pair2mapped.bam ../data/bellerophonOutput/viral/bothPairsMapped.bam | samtools sort -n - -o ../data/bellerophonOutput/viral/virusAligned-final.bam
bedtools bamtofastq -i ../data/bellerophonOutput/viral/virusAligned-final.bam -fq ../data/bellerophonOutput/viral/ViralAligned_1.fastq -fq2 ../data/bellerophonOutput/viral/ViralAligned_2.fastq
bowtie2 -p 1 --local -x ../data/hg38/hg38_bwt2_index -1 ../data/bellerophonOutput/viral/ViralAligned_1.fastq -2 ../data/bellerophonOutput/viral/ViralAligned_2.fastq | samtools sort -o ../data/bellerophonOutput/host/hostAligned.bam -
samtools sort ../data/bellerophonOutput/host/hostAligned.bam -o ../data/bellerophonOutput/host/hostAligned_sorted.bam
samtools index ../data/bellerophonOutput/host/hostAligned_sorted.bam
