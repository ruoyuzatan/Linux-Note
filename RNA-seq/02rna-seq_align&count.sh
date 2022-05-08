#！/bin/bash
#用于RNA-seq数据比对和计算reads
index=/home/hanlei/hanlei/quinoa_genome
REFannotation_PATH=/home/hanlei/hanlei/quinoa_genome/quinoa_genome_annotation.gtf
DATA_PATH=/home/hanlei/hanlei/quinoa_genome/data/quinoa_salt_300mM

mkdir $DATA_PATH/04align
mkdir $DATA_PATH/05counts

for ((i=11921127;i<=11921149;i++));
do
hisat2 -p 10 -x $index/quinoa -1 $DATA_PATH/01clean/SRR${i}_1.paired.fq -2 $DATA_PATH/01clean/SRR${i}_2.paired.fq -S $DATA_PATH/04align/SRR${i}.sam
samtools view -S $DATA_PATH/04align/SRR${i}.sam -b > $DATA_PATH/04align/SRR${i}.bam
samtools sort -n $DATA_PATH/04align/SRR${i}.bam -o $DATA_PATH/04align/SRR${i}_nsorted.bam
featureCounts -T 10 -p -a $REFannotation_PATH  -o  $DATA_PATH/05counts/SRR${i}.count $DATA_PATH/04align/SRR${i}_nsorted.bam
done