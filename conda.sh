RNAseq
#conda的安装与使用
wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod 777 Miniconda3-latest-Linux-x86_64.sh #给执行权限
bash Miniconda3-latest-Linux-x86_64.sh #运行
#一路yes，安装成功
#找到你刚才安装的miniconda，如果没有更改过安装位置的话应该是在/home下面，cd到miniconda3的bin目录下面，能看到有一个activate
#添加权限
chmod 777 activate 
#启动conda
. ./activate
. /home/shengli1/miniconda3/bin/activate
#这两个点不是连在一起的
#当命令行前面出现(base)的时候说明现在已经在conda的环境中了。
conda list #查看命令
#利用conda安装生物信息软件
conda install samtools #安装命令
#可以用这个命令进行搜索看想安装的软件存不存在
conda search samtools
#添加生物信息分析常用的channel
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/

#对channel常用的操作

#显示安装的频道
 conda config --set show_channel_urls yes 
#查看已经添加的channels
conda config --get channels
#已添加的channel在哪里查看
vim ~/.condarc
#利用conda安装生物信息软件
conda install samtools
#安装完成后，可以用“which 软件名”来查看该软件安装的位置：
which samtools
#搜索目前软件包有哪几个版本
conda search samtools
#如需要安装特定的版本:
conda install 软件名=版本号
conda samtools=1.10
#查看已安装软件:
conda list
#更新指定软件:
conda update samtools
#卸载指定软件:
conda remove gatk
#退出conda环境
. ./deactivate
