https://zhuanlan.zhihu.com/p/144830963
PRE=/home/hanlei/app/sra/sratoolkit.2.11.1-ubuntu64/bin
DATA=/home/hanlei/hanlei/quinoa_genome/data/quinoa_salt_200mM
nohup $PRE/prefetch -f yes PRJNA605324 -O $DATA >> $DATA/downdata.log 2>&1 &
SRR2989985 PRJNA305752
$PRE/prefetch SRR10914949 SRR10914948 SRR10914947 SRR10914946 SRR10914945 SRR10914944 -O $DATA

for ((i=529;i<=532;i++)); do $PRE/fastq-dump --split-files SRR7755${i}/SRR7755${i}.sra;done

#解压 
$PRE/fastq-dump --split-files SRR10216549.sra

