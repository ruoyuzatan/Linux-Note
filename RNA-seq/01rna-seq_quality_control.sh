#！/bin/bash
#用于对RNA-seq数据清洗和质量控制
adapter=/home/hanlei/app/Trimmomatic-0.38/adapters/TruSeq3-PE.fa
DATA_PATH=/home/hanlei/hanlei/quinoa_genome/data/quinoa_salt_300mM #数据目录的上一级
mkdir $DATA_PATH/01clea $DATA_PATH/02unpair $DATA_PATH/03QC
mkdir $DATA_PATH/04align
mkdir $DATA_PATH/05counts

#质量控制
for ((i=11921132;i<=11921144;i++));
do
trimmomatic PE -threads 5 -phred33 $DATA_PATH/raw/SRR${i}_1.fastq.gz $DATA_PATH/raw/SRR${i}_2.fastq.gz  $DATA_PATH/01clean/SRR${i}_1.paired.fq.gz  $DATA_PATH/02unpair/SRR${i}_1.unpaired.fq.gz  $DATA_PATH/01clean/SRR${i}_2.paired.fq.gz  $DATA_PATH/02unpair/SRR${i}_2.unpaired.fq.gz ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:10 MINLEN:36;
pigz -d -k -p 5 $DATA_PATH/01clean/SRR${i}_*.fq.gz #解压cleanreads
fastqc -f fastq -t 5 $DATA_PATH/01clean/SRR${i}_*.fq -o $DATA_PATH/03QC #生成质量控制报告哦
done

multiqc $DATA_PATH/03QC #合并质量控制报告
