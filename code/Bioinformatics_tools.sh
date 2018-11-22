#! /usr/bin/ bash

# This file starts by listing all the programs you will need to preprocess your raw data before metabarcoding.
# We first create the directory where we will add all the tools and software programs needed. We will also add the programs to the computer PATH.
# This is done so you can call the tool you need regardless of where in the computer you are, just like the usual bash commands.


# Lets BEGIN the journey — your life is about to get better. A good thing to remember is that you always need to let the computer know “where things are”

# Open terminal. You can find it in your applications folder from finder. Might be handy to add tour dock.
mkdir ~/bio_tools
echo “PATH=$PWD/bio_tools:$PATH” >> ~/.bash_profile

cd ~/bio_tools

# Install miniconda
# This package is python based and may take a little while to finish. once it is finished it will make your life easy with regards to dependencies and bioinformatics.
# in the files you've downloaded from cloudstor there is a script that will download the python 3.7 version of miniconda into your computer. You can always go to the
# website https://conda.io/miniconda.html and download the version or type you prefer. for now i have added it to the Code folder. This version is for Mac OS X and accessed Oct 2018.

bash ~/WAHMCC_metabarcoding_data_analysis/code/Miniconda3-latest-MacOSX-x86_64.sh

# Adding vsearch 
## This is the free alternative to usearch (we didn't download Usearch because R. Edgar has not yet made his tool open source and applies a restriction for anything above 32 bit) - Vsearch can do most functions (not Unnoise3 algorithm for denoising reads) of Usearch but might be slower and lesser accuracy !! I Don’t know…

curl -LO https://github.com/torognes/vsearch/releases/download/v2.5.1/vsearch-2.5.1-macos-x86_64.tar.gz
tar -xzvf vsearch-2.5.1-macos-x86_64.tar.gz
rm vsearch-*.tar.gz
chmod +x vsearch/*
echo “PATH=$PWD/vsearch:$PATH” >> ~/.bash_profile

# Adding fastQC 
# This is an important sequence QC program that can be run both CLI and GUI. It tells you what your raw Fastq files from illumina look like. you can basically set all your parameters using the information learned from it's output.

wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
unzip fastqc_v0.11.7.zip
chmod +x FastQC/*
sudo ln -s $PWD/FastQC/fastqc /usr/local/bin/fastqc

# Adding fastp
# This tool works the same way as fastqc but in my opinion it is better and has functionality to do some filtering.

# note: the fastp version in bioconda may be not the latest (bug fixes)
conda install -c bioconda fastp

# Alternatively, you can download fastp using git and compiling from source, as follows; remove the #s and run.
#git clone https://github.com/OpenGene/fastp.git
#cd fastp
#make
#make install
#cd ..

# Adding Trimmomatic
# This is a java tool - it is like the standard in the field with regards to applying your parameters from fastQC output to your raw data.

curl -LO http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
unzip Trimmomatic-0.36.zip
rm Trimmomatic-0.36.zip
echo “PATH=$PWD/Trimmomatic_0.36/trimmomatic-0.36.jar:$PATH” >> ~/.bash_profile

# Adding BBtools 
# This tool is very handy for QC and other manipulations - would be great to know how to use these tools well because I think they will easily converted in Pawsey )

curl -L https://sourceforge.net/projects/bbmap/files/latest/download -o bbtools.tar.gz
tar -xzvf bbtools.tar.gz
rm -rf bbtools.tar.gz 
chmod +x bbmap/*
echo “PATH=$PWD/bbmap:$PATH” >> ~/.bash_profile


# Adding a shi7
# This tool will download a few tools one of them is "gotta_split" which is a paired-end demultiplexer..very very fast and accurate, 

wget https://github.com/knights-lab/shi7/releases/download/v0.9.9/shi7_0.9.9_mac_release.zip
unzip  shi7_0.9.9_mac_release.zip
rm shi7_0.9.9_mac_release.zip
chmod +x shi7_0.9.9_mac_release/*
echo “PATH=$PWD/shi7_0.9.9_mac_release:$PATH” >> ~/.bash_profile

#Adding cutadapt
# This tool is essential when you need to trim and QC your raw fastq files. It is awesome at removing adapters and primers
conda install -c bioconda cutadapt
# Alternatively, pip install --user --upgrade cutadapt

#Adding seqkit & seqtk
# Both of these tools are extremely fast that can deal with a variety of raw sequence data like fastq/fasta. They are also used for subsampling if you want to do that.
conda install -c bioconda seqkit
git clone https://github.com/lh3/seqtk.git;
cd seqtk; make
echo “PATH=$PWD/seqtk:$PATH” >> ~/.bash_profile

source ~/.bash_profile

#END
