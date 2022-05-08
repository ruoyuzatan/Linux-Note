#一、数据质量控制
#Trimmomatic

#search
conda search trimmomatic
#安装
conda install trimmomatic

#Trimmomatic 质控用法
adapter=/home/hanlei/app/Trimmomatic-0.38/adapters/TruSeq3-PE.fa
# SE 模式
#SE模式下，只有一个输入文件和一个质控后的输出文件，运行命令如下
trimmomatic SE –threads <线程数> <input> <output> <step1> <step2> …<step1><step2>… 表示每一步的质控参数

# PE模式
#PE 模式下，有两个输入文件（正向测序reads和反向测序reads）和四个质控后的输出文件（双端序列都保留的paired序列文件和只保留一端序列的unpaired序列文件），运行命令如下：
trimmomatic PE -threads 12 -phred33 C4_1.fq.gz C4_2.fq.gz C4_1.paired.fq.gz C4_1.unpaired.fq.gz C4_2.paired.fq.gz C4_2.unpaired.fq.gz ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:10 MINLEN:36
#其中R1.fq.gz以及 R2.fq.gz 为两个输入文件
#R1.paired.fq.gz 、R1.unpaired.fq.gz、 R2.paired.fq.gz 、R2.unpaired.fq.gz 为四个对应的输出文件,如果在接下来要进行序列比对的话用的文件只需要使用到两个paired文件。
#Phred33 设置碱基的质量格式，默认的是-phred64。
#ILLUMINACLIP:$adapter.fa:2:30:10 adapter.fa为接头文件，2表示最大mismatch数，30表示palindrome模式下碱基的匹配阈值，10表示simple模式下碱基的匹配阈值。
#LEADING: 3 表示切除reads 5’端碱基质量低于3的碱基。
#TRAILING:3 表示切除3’ 端碱基质量低于3的碱基。
#SLIDINGWINDOW:4:15 表示以4个碱基为窗口进行滑动，切除窗口内碱基平均质量小于15的。
#MINLEN:36 丢弃以上步骤处理后，序列长度小于36的reads。
#HEADCROP:切掉从起始开始的10个碱基
#批量进行质量控制
mkdir 01cleandata
mkdir 02unpair
adapter=/home/hanlei/app/Trimmomatic-0.38/adapters/TruSeq3-PE.fa
DATA_PATH=/home/hanlei/hanlei/quinoa_genome/data/dro_qui
PAIROUT_PATH=/home/hanlei/hanlei/quinoa_genome/data/dro_qui/01cleandata
UNPAIROUT_PATH=/home/hanlei/hanlei/quinoa_genome/data/dro_qui/02unpair
FILE=/home/hanlei/hanlei/Transcriptomes/Lihong_raw/filelist.txt #文件的换行符需要从CRLF设置成LF，设置方法为在Ultraedit中将^r^n替换为^n
cd $DATA_PATH
cat $FILE | while read line
do
nohup trimmomatic PE -threads 12 -phred33 ${line}_1.fq.gz ${line}_2.fq.gz $OUT_PATH/${line}_1.paired.fq.gz ${line}_1.unpaired.fq.gz $OUT_PATH/${line}_2.paired.fq.gz ${line}_2.unpaired.fq.gz ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:15 MINLEN:36 >$OUT_PATH/${line}.log 2>&1 & 
done

for ((i=984;i<=985;i++));
do
trimmomatic PE -threads 12 -phred33 $DATA_PATH/SRR2989${i}_1.fastq.gz $DATA_PATH/SRR2989${i}_2.fastq.gz  $PAIROUT_PATH/SRR2989${i}_1.paired.fq.gz  $UNPAIROUT_PATH/SRR2989${i}_1.unpaired.fq.gz  $PAIROUT_PATH/SRR2989${i}_2.paired.fq.gz  $UNPAIROUT_PATH/SRR2989${i}_2.unpaired.fq.gz ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:10 MINLEN:36; 
done

trimmomatic PE -threads 12 -phred33 SRR2989985_1.fastq.gz SRR2989985_2.fastq.gz  $PAIROUT_PATH/SRR2989985_1.paired.fq.gz  $PAIROUT_PATH/SRR2989985_1.unpaired.fq.gz  $UNPAIROUT_PATH/SRR2989985_2.paired.fq.gz $UNPAIROUT_PATH/SRR2989985_2.unpaired.fq.gz ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:10 MINLEN:36; 

#批量解压paired文件
pigz -d *.gz -k -p 12

#fastq质量汇报
# 安装fastqc工具包
$ sudo apt-get install fastqc
# 执行fastqc
#raw_data
fastqc -f fastq -t 30 AS1_1.fq 
#clean_data
fastqc -f fastq -not 30 -o /home/shengli1/Rabbit/llh_raw_data_2year/results C2_1.paired.fq, C2_2.paired.fq
#multiqc质量报告（需要在conda中启动，multiqc可以对几个fastqc报告文件进行总结并汇总到一个报告文件中，以更直观到防止展示。）
multiqc <analysis directory>
#批量进行质量控制

FASTA_RAWDATA_PATH=/home/hanlei/hanlei/quinoa_genome/data/dro_qui/01cleandata
FASTQC_RAWOUTPATH=/home/hanlei/hanlei/quinoa_genome/data/dro_qui/03QC
cd $FASTA_RAWDATA_PATH
for id in $(ls *.fastq); do nohup fastqc -f fastq -t 5 $id -o $FASTQC_RAWOUTPATH >$FASTQC_RAWOUTPATH/${id}.log 2>&1 &  done

FASTA_DATA_PATH=/home/hanlei/hanlei/quinoa_genome/data/dro_qui/01cleandata
FASTQC_OUTPATH=/home/hanlei/hanlei/quinoa_genome/data/dro_qui/03QC
cd $FASTA_DATA_PATH
for id in $(ls *.fq); do nohup fastqc -f fastq -t 5 $id -o $FASTQC_OUTPATH >$FASTQC_OUTPATH/${id}.log 2>&1 &  done

二、序列比对
#STAR进行序列比对
#先要进入conda环境，见脚本conda
#建立索引
STAR --runThreadN 6 

--runMode genomeGenerate #工作模式（建立索引）

--genomeDi output_folder #索引文件存储位置

--genomeSAindexNbases 12 #报错，14太大，需要用12

--genomeFastaFiles 00ref/TAIR10_Chr.all.fasta #参考基因组

--sjdbOverhang 149 #可变剪切的预测

--sjdbGTFfile 00ref/Araport11_GFF3_genes_transposons.201606.gtf#注释文件
STAR --runThreadN 6 --runMode genomeGenerate --genomeDi output_folder --genomeSAindexNbases 12 --genomeFastaFiles 00ref/TAIR10_Chr.all.fasta --sjdbGTFfile 00ref/Araport11_GFF3_genes_transposons.201606.gtf#注释文件

#进行比对
STAR --runThreadN 5 \

--genomeDir output_folder \

--readFilesCommand zcat \#解压缩文件

--readFilesIn 02clean_data/output_forward_paired.fq.gz \#输入文件位置

 02clean_data/output_reverse_paired.fq.gz \

--outFileNamePrefix 03align_out/sample2_ \ #结果文件

--outSAMtype BAM SortedByCoordinate \ #输出BAM文件并排序

--outBAMsortingThreadN 5 \

--quantMode TranscriptomeSAM GeneCounts #定量分析每个基因上的reads数，为使用RSEM进行定量分析做准备

#实例 STAR --runThreadN 20--genomeDir zhangdir --readFilesIn C1_1.paired.fq C1_2.paired.fq --outFileNamePrefix align_out/C1 --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 10

#bowtie2进行序列比对
#(1)安装
sudo apt-get install bowtie2
#(2)建立索引（或下载官方索引）
REFGenome_PATH=~/hanlei/Transcriptomes/reference/genome/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta #基因组序列
REFannotation_PATH=~/hanlei/Transcriptomes/reference/gtf #注释文件目录
index=~/hanlei/Transcriptomes/reference/index #索引文件目录
bowtie2-build --large-index $REFGenome_PATH $index/iwgsc1.0 
#文件格式默认为fasta
#前面为基因组文件，后面为输出索引文件目录及文件名
#--large-index建立大基因组索引


#(3)进行序列比对
bowtie2 -p 6 -3 5 --local -x $REF_PATH/Triticum_aestivum_IWGSC -1 $DATA_PATH/CR_1.paired.fq -2 $DATA_PATH/CR_2.paired.fq -S SRR3208744
#$REF_PATH 索引文件
#(4)批量进行序列比对
DATA_PATH=/home/hanlei/hanlei/Transcriptomes/plosone/01cleandata
index=~/hanlei/Transcriptomes/reference/index
OUT_PATH=~/hanlei/Transcriptomes/plosone/03align
FILE=/home/shengli1/Rabbit/HL_wheat/HL_wheat/raw/filelist.txt #RAS12_1 RAS12_2 RAS48_1 RAS48_2 RCK12_1 RCK12_2 RCK48_1 RCK48_2 
#文件的换行符需要从CRLF设置成LF，设置方法为在Ultraedit中将^r^n替换为^n
cd $DATA_PATH
nohup cat $FILE | while read line
do
bowtie2 -p 16 -3 5 --local -x $REF_PATH/Triticum_aestivum_IWGSC -1 $DATA_PATH/${line}_1.paired.fq -2 $DATA_PATH/${line}_2.paired.fq -S $OUT_PATH/${line}.sam
done
> rundata.log 2>&1 &

for ((i=529;i<=532;i++));
do
bowtie2 -p 16 -3 5 --local -x $index/Ta_IWGSC1.0_index -1  $DATA_PATH/SRR7755${i}_1.paired.fq -2  $DATA_PATH/SRR7755${i}_2.paired.fq -S $OUT_PATH/SRR7755${i}.sam
done


#hisat2进行序列比对
#（1）安装hisat2
sudo apt-get install hisat2

#(2)将GTF注释文件转换成hisat2能用的格式
gffread **.gff3 -T -o **.gtf   #gff3格式需转换成gtf格式才能用
REFGenome_PATH=/home/hanlei/hanlei/quinoa_genome/quinoa_genome.fa
REFannotation_PATH=/home/hanlei/hanlei/quinoa_genome
index=/home/hanlei/hanlei/quinoa_genome/index
hisat2_extract_exons.py $REFannotation_PATH/quinoa_genome_annotation.gtf > quinoa.exon #外显子/转录组
hisat2_extract_splice_sites.py $REFannotation_PATH/quinoa_genome_annotation.gtf > quinoa.ss #可变剪切
hisat2_extract_snps.py snp.txt > genome.snp #SNP序列
#(3)创建索引文件（即将序列打断成片段进行多种组合，以便提高比对效率）,必须要在基因组文件存放的文件夹建立索引
hisat2-build quinoa_genome.fa quinoa –p 4#建立基因组索引
nohup hisat2-build -p 30 --large-index $REFGenome_PATH  --ss $REFannotation_PATH/Triticum_aestivum_IWGSC.ss --exon $REFannotation_PATH/Triticum_aestivum_IWGSC.exon $index/Triticum_aestivum_IWGSC_index  > index.log 2>&1 & #建立基因组+转录组
#-p 30 p为线程, 一般多少核去运行，这个看自己电脑的内存，
#--snp genome.snp 如果要建立SNP索引，则将这一命令放到上面一行命令中

nohup hisat2  -t -p 16 -x  $REF_PATH/RNA-seq_hisat2_index  -1  RAS12_1_1.paired.fq   -2  RAS12_1_2.paired.fq -S test.sam > program_2.log 2>&1 &
hisat2 -t -p 30 -x  $REF_PATH/RNA-seq_hisat2_index  -1  RAS12_1_1.paired.fq   -2  RAS12_1_2.paired.fq -s test.sam | samtools view -bS 1>$OUT_PATH/test.bam
#利用hisat2进行批量序列比对
DATA_PATH=/home/hanlei/hanlei/quinoa_genome/data/dro_qui/01cleandata
index=/home/hanlei/hanlei/quinoa_genome
OUT_PATH=/home/hanlei/hanlei/quinoa_genome/data/dro_qui/04align
FILE=/home/shengli1/Rabbit/xbb_At_transcriptome_raw_data/filelist.txt #RAS12_1 RAS12_2 RAS48_1 RAS48_2 RCK12_1 RCK12_2 RCK48_1 RCK48_2 
#文件的换行符需要从CRLF设置成LF，设置方法为在Ultraedit中将^r^n替换为^n
cd $DATA_PATH
cat $FILE | while read line
do
hisat2 -p 16 -x $REF_PATH/RNA-seq_hisat2_index -1 $DATA_PATH/${line}_1.paired.fq -2 $DATA_PATH/${line}_2.paired.fq |samtools view -bS 1>$OUT_PATH/${line}.bam
done
#输出sam
for ((i=984;i<=985;i++));
do
hisat2 -p 30 -x $index/quinoa -1 $DATA_PATH/SRR2989${i}_1.paired.fq -2 $DATA_PATH/SRR2989${i}_2.paired.fq -S $DATA_PATH/SRR2989${i}.sam
done

for ((i=984;i<=985;i++));
do
hisat2 -p 30 -x $index/quinoa -1 $DATA_PATH/SRR2989${i}_1.paired.fq -2 $DATA_PATH/SRR2989${i}_2.paired.fq |samtools view -bS 1>$OUT_PATH/SRR2989${i}.bam
done

#使用samtools实现sam向BAM文件转换，并排序、添加index
#将HISAT2处理的结果输出到samtools转化为bam格式
#此处使用6核，约使用6.4G内存，平均每文件处理需30min
DATA_PATH=/home/hanlei/hanlei/Transcriptomes/Lihong_raw/mapping
OUT_PATH=/home/hanlei/hanlei/Transcriptomes/Lihong_raw/lihong_clean/sorted
FILE=~/hanlei/Transcriptomes/Lihong_raw/filelist.txt #RAS12_1 RAS12_2 RAS48_1 RAS48_2 RCK12_1 RCK12_2 RCK48_1 RCK48_2 
#使用samtools实现sam向BAM文件转换，并排序、添加index
cat $FILE | while read line
do 
nohup samtools view -S $DATA_PATH/${line}.sam -b > $OUT_PATH/${line}.bam $OUT_PATH/logs/${line}_bam.log 2>&1 &
done

for ((i=984;i<=985;i++));do samtools view -S SRR2989${i}.sam -b > SRR2989${i}.bam;done
# 将所有的bam文件按默认的染色体位置进行排序
for ((i=984;i<=985;i++));do samtools sort SRR2989${i}.bam -o SRR2989${i}_sorted.bam;done
# 将所有的排序文件建立索引，索引文件.bai后缀
for ((i=984;i<=985;i++));do samtools index SRR2989${i}_sorted.bam;done

#一步运行
for i in `seq 56 62`
do
    samtools view -S SRR35899${i}.sam -b > SRR35899${i}.bam
    samtools sort SRR35899${i}.bam -o SRR35899${i}_sorted.bam
    samtools index SRR35899${i}_sorted.bam
done

#-n 按read name 排序 ，如果不指定则按染色体位置排序
for ((i=984;i<=985;i++));do samtools sort -n SRR2989${i}.bam -o SRR2989${i}_nsorted.bam;done 

DATA_PATH=~/hanlei/Transcriptomes/plosone/03align
OUT_PATH=~/hanlei/Transcriptomes/plosone/04sorted
FILE=~/hanlei/Transcriptomes/Lihong_raw/filelist.txt
cat $FILE | while read line
do 
nohup samtools sort -n $DATA_PATH/${line}.bam -o $OUT_PATH/sorted_${line}.bam > $OUT_PATH/logs/${line}.log 2>&1 &
done
# -n 按read name 排序 ，如果不指定则按染色体位置(pos)排序
# -@ 线程
# -m 每个线程运行时的内存大小，默认500M，可以根据自己的情况调大一点
for ((i=529;i<=532;i++));
do
samtools sort -n $DATA_PATH/SRR7755${i}.bam -o $OUT_PATH/SRR7755${i}.bam
done

index=/home/hanlei/hanlei/Transcriptomes/Lihong_raw/mapping/sorted/index

cd $OUT_PATH
for i in $(ls sorted_*.bam)
do 
nohup samtools index $i > $OUT_PATH/logs/${i}_sorted.log 2>&1 & 
done

for i in $(ls *.bam)
do 
samtools index $i
done


# 对已经排序的bam文件进行简单质量控制
对BAM文件进行统计分析
#启动python2.7环境
source activate python2
#安装numpy
conda install numpy
sudo apt-get install python2.7-dev
#安装gcc
sudo apt  install gcc
#查看是否安装
gcc --version
# 用pip命令安装
pip install RSeQC
# 对bam文件进行质控，其余都同样的进行
bam_stat.py  -i sorted_*.bam

#将gtf转换成bed
/home/hanlei/app/gtf2bed gencode.v19.long_noncoding_RNAs.gtf >gencode.v19.long_noncoding_RNAs.bed

#查看基因组覆盖率

REFannotation_PATH=~/hanlei/Transcriptomes/reference/gtf
$REFannotation_PATH/xxc_genome.gtf
for i in $(ls sorted*.bam)
do
read_distribution.py -i $i -r $REFannotation_PATH/xxc_genome.bed >  $i.count
done

三、对reads进行计数并添加注释

（一）featurecounts(计数前需要按照名称进行排序)
#1.使用
featureCounts -T 10 -a $gtf -o read.count -p  *.bam
#2.参数
 #input file	输入的bam/sam文件，支持多个文件输入
 #-a < string >	参考gtf文件名，支持Gzipped文件格式
 #-F	参考文件的格式，一般为GTF/SAF，C语言版本默认的格式为GTF格式
 #-A	提供一个逗号分割为两列的文件，一列为gtf中的染色体名，另一列为read中对应的染色体名，用于将gtf和read中的名称进行统一匹配，注意该文件提交时不需要列名
 #-J	对可变剪切进行计数
 #-G < string >	当-J设置的时候，通过-G提供一个比对的时候使用的参考基因组文件，辅助寻找可变剪切
 #-M	如果设置-M，多重map的read将会被统计到
 #-O	允许多重比对，即当一个read比对到多个feature或多个metafeature的时候，这条read会被统计多次
 #-o < string >	输出文件的名字，输出文件的内容为read 的统计数目
 #-T	线程数目，1~32
 #下面是有关featrue/metafeature选择的参数
 #-p	只能用在paired-end的情况中，会统计fragment而不统计read
 #-B	在-p选择的条件下，只有两端read都比对上的fragment才会被统计
 #-C	如果-C被设置，那融合的fragment（比对到不同染色体上的fragment）就不会被计数，这个只有在-p被设置的条件下使用
 #-d < int >	最短的fragment，默认是50
 #-D < int >	最长的fragmen，默认是600
 #-f	如果-f被设置，那将会统计feature层面的数据，如exon-level，否则会统计meta-feature层面的数据，如gene-levels
 #-g < string >	当参考的gtf提供的时候，我们需要提供一个id identifier 来将feature水平的统计汇总为meta-feature水平的统计，默认为gene_id，注意！选择gtf中提供的id identifier！！！
 #-t < string >	设置feature-type，-t指定的必须是gtf中有的feature，同时read只有落到这些feature上才会被统计到，默认是“exon”
 #参考网址：https://www.jianshu.com/p/9cc4e8657d62

#3.实例
FILE_PATH=/home/hanlei/hanlei/quinoa_genome/data/dro_qui/04align
REFannotation_PATH=/home/hanlei/hanlei/quinoa_genome/quinoa_genome_annotation.gtf #注释文件目录
OUT_PATH=/home/hanlei/hanlei/quinoa_genome/data/dro_qui/05counts
REFGE_PATH=/home/hanlei/hanlei/Transcriptomes/reference/genome/iwgsc_refseqv1.0_all_chromosomes #基因组文件目录

cd $FILE_PATH
for ((i=984;i<=985;i++));
do
featureCounts -T 10 -p -a $REFannotation_PATH  -o $OUT_PATH/SRR2989${i}.count  $FILE_PATH/SRR2989${i}_nsorted.bam
done
#> $OUT_PATH/${i}.count.log 2>&1 &

（二）#使用HTSEQ得到COUNT值（入门）
#安装
conda install htseq
# 利用htseq-count对sort之后的bam文件进行reads计数
htseq-count -s yes -r name -f bam test.bam /home/shengli1/Rabbit/sequence/A.thaliana/genome.gtf > count.txt
#f bam/sam： 指定输入文件格式，默认SAM
#-r name/pos: 你需要利用samtool sort对数据根据read name或者位置进行排序，默认是name
#-s yes/no/reverse: 数据是否来自于strand-specific assay。DNA是双链的，所以需要判断到底来自于哪条链。如果选择了no， 那么每一条read都会跟正义链和反义链进行比较。默认的yes对于双端测序表示第一个read都在同一个链上，第二个read则在另一条链上。
#-a 最低质量， 剔除低于阈值的read
#-m 模式 union（默认）, intersection-strict and intersection-nonempty。一般而言就用默认的，作者也是这样认为的。
#-i id attribute: 在GTF文件的最后一栏里，会有这个基因的多个命名方式（如下）， RNA-Seq数据分析常用的是gene_id， 当然你可以写一个脚本替换成其他命名方式。

FILE_PATH=/home/hanlei/hanlei/quinoa_genome/data/dro_qui/04align
REFannotation_PATH=/home/hanlei/hanlei/quinoa_genome/quinoa_genome_annotation.gtf #注释文件目录
OUT_PATH=/home/hanlei/hanlei/quinoa_genome/data/dro_qui/05counts
cd $FILE_PATH

for ((i=984;i<=985;i++));
do htseq-count -r name -f bam $FILE_PATH/SRR2989${i}_nsorted.bam $REFannotation_PATH > $OUT_PATH/SRR2989${i}.count;
done 

for i in $(ls sorted_*.bam)
do
nohup htseq-count -r name -f bam $FILE_PATH/${i} $REFannotation_PATH/iwgsc_refseqv1.0_HighConf_2017Mar13.gtf  >$OUT_PATH/${i}.count \
2 > $OUT_PATH/${i}.log 2>&1 &
done 
#--idattr=gene 可以用这个改程序要寻找的id名字

#在R中合并文件，并将ID转换为Gene_symbol
#运行R设置工作路径并查看当前目录下的文件
setwd("H:/Plant Physiology/生物信息学/运行数据/xbb_Athaliana")
list.files()
# 防止R自动把字符串string的列辨认成factor
options(stringAsFactor =  FALSE)  
# 从19R576_paired.count文件中读取数据，并添加列名，之后生成新文件命名为control_W55
fun1 <- function(x){
read.table(x,,sep = "\t",col.names=c("gene_id",x))
}
filename <-  list.files(pattern=".txt")
xx <- lapply(filename,fun1)
raw_count1<-merge(merge(xx[[1]],xx[[2]], by = "gene_id"),merge(xx[[3]],xx[[4]], by = "gene_id"),by ="gene_id")
raw_count2<-merge(merge(xx[[5]],xx[[6]], by = "gene_id"),merge(xx[[7]],xx[[8]], by = "gene_id"),by ="gene_id")
raw_count<-merge(raw_count1,raw_count2, by = "gene_id")
#删除显示未匹配数据的前5行，删之前确定一下
raw_count_filt <- raw_count[-1:-5,]
#因为我们无法在EBI数据库上直接搜索找到ENSMUSG00000024045.5这样的基因，只能是ENSMUSG00000024045的整数，没有小数点，所以需要进一步替换为整数的形式。
# 第一步通过gsub工具将匹配到的.以及后面的数字连续匹配并替换为空，并赋值给ENSEMBL
ENSEMBL <- gsub("\\.\\d*", "", raw_count_filt$gene_id) 
# 将ENSEMBL重新添加到raw_count_filt1矩阵
row.names(raw_count_filt) <- ENSEMBL
# 查看>head(raw_count_filt)

#对基因进行注释-获取gene_symbol
# 首先检查BioManager包是否安装
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# 如果BioManager包已经安装，则通过BioManager包安装biomaRt包
BiocManager::install("biomaRt")
# 因为我的文件是放在桌面上的，所以先设置R的工作目录为桌面路径
> setwd("C://Users/My/Desktop/")
# 加载放在桌面上的.Rdata文件
> load("raw_count.Rdata")
# 加载所需要的R包
> library('biomaRt')
> library("curl")
# 用bioMart对差异表达基因进行注释
#载入biomaRt并显示一下能连接的数据库
>library("biomaRt")
>listMarts()
#用useMart函数选定数据库
listMarts(host="plants.ensembl.org")
a<-useMart(biomart="plants_mart",host="plants.ensembl.org")
datasets<-listDatasets(a) #显示当前数据库所含的基因组注释
mart <- useDataset("athaliana_eg_gene",useMart(biomart="plants_mart",host="plants.ensembl.org"))
my_ensembl_gene_id<-row.names(raw_count_filt)
mms_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
#合并数据，先修改行名
c<-colnames(raw_count_filt)
c[1]="ensembl_gene_id"
colnames(raw_count_filt)<-c
diff_name<-merge(raw_count_filt,mms_symbols,by="ensembl_gene_id")

四、DEseq2筛选差异表达基因
> library(tidyverse)
> library(DESeq2)
> #import data
> setwd("F:/rna_seq/data/matrix")
> mycounts<-read.csv("readcount.csv")
> head(mycounts)
                   X control1 control2 treat1 treat2
1 ENSMUSG00000000001     1648     2306   2941   2780
2 ENSMUSG00000000003        0        0      0      0
3 ENSMUSG00000000028      835      950   1366   1051
4 ENSMUSG00000000031       65       83     52     53
5 ENSMUSG00000000037       70       53     94     66
6 ENSMUSG00000000049        0        3      4      5
#这里有个x，需要去除，先把第一列当作行名来处理
> rownames(mycounts)<-mycounts[,1]
#把带X的列删除
> mycounts<-mycounts[,-1]
> head(mycounts)
                   control1 control2 treat1 treat2
ENSMUSG00000000001     1648     2306   2941   2780
ENSMUSG00000000003        0        0      0      0
ENSMUSG00000000028      835      950   1366   1051
ENSMUSG00000000031       65       83     52     53
ENSMUSG00000000037       70       53     94     66
ENSMUSG00000000049        0        3      4      5
# 这一步很关键，要明白condition这里是因子，不是样本名称；小鼠数据有对照组和处理组，各两个重复
> condition <- factor(c(rep("control",2),rep("treat",2)), levels = c("control","treat"))
> condition
[1] control control treat   treat  
Levels: control treat
#colData也可以自己在excel做好另存为.csv格式，再导入即可
> colData <- data.frame(row.names=colnames(mycounts), condition)
> colData
         condition
control1   control
control2   control
treat1       treat
treat2       treat