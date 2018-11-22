# WAHMCC_metabarcoding_data_analysis
This repository contains code and information on how the TrEnD_Lab conducts exploratory metabarcoding data analysis. R and shell are used for code. 

**A summary of what and how we do metabarcoding**

The metabarcoding bioinformatics pipeline that we have developed in our laboratory has been designed to address the biases and inconsistencies observed in metabarcoding research. We identify these biases in the early stages of the bioinformatics process by removing all non-biological sequences such as adapters, primers and low quality bases as well as during the clustering and taxonomy allocation stages. 

Our stringent protocols are in harmony with maintaining a low level of contamination as well as sufficient recovery rate of target sequences due to our unique dual sequence barcoding that we associated with each sample for the purposes of multiplexing and contamination avoidance. Finally we use DADA2 to pick our amplicon sequence variants not the usual 97% clustering approach operational taxonomic units (OTUs). 

DADA2 is a sequence denoising algorithm that applies sophisticated statistical and machine learning approaches to learn the error profiles of each of the target sequences and uses that information to make the best inference about the true bases of each short read sequence which are then binned together at 100% similarity as exact sequence variants (ASVs). Then using a trained naïve bayes classifier the ASVs are assigned taxonomy to Genus and sometimes to Species level (dependent of the length of and region of the marker gene sequenced) against a curated database of microbial reference sequences such as the RDP, SILVA, GTDB and NCBI's 16SMicrobial. 

The tools we use include FastQC for read quality visualization, Shi7 and GHAPv2 for demultiplexing and preprocessing of samples, cutadapt for removal of all non-biological sequences. Finally we use the R programming language to invoke DADA2 for quality filtering, error correction, ASV picking and the package Phyloseq for preliminary statistical data exploration.


## Bioinformatics tools installation
* There is a script in the code/ folder that will do this for you.*

Please be aware you will need to have wget, curl, python2/3 (***you may need to read up about what it will do to your computer to have both installed unless you can figure out a hack that can help you do that without dramas you can use stackoverflow it is a good place for solutions***), and conda installed (*conda we will install together*) inorder to continue with this workshop.

###The following text will starts by listing all the programs you will need to preprocess your raw target gene sequencing data.

We first create the directory where we will add all the tools and software programs needed. We will also add the programs to the computer PATH. This is done so you can call the tool you need regardless of where in the computer you are, just like the usual bash commands.


Lets BEGIN the journey — your life is about to get better. A good thing to remember is that you always need to let the computer know “where things are”

*   1.   Open terminal. *You can find it in your applications folder from finder. Might be handy to add tour dock.*
mkdir ~/bio_tools
echo “PATH=$PWD/bio_tools:$PATH” >> ~/.bash_profile

cd ~/bio_tools

*   2.  Install miniconda
This package is python based and may take a little while to finish. once it is finished it will make your life easy with regards to dependencies and bioinformatics. You can always go to the website https://conda.io/miniconda.html and download the version or type you prefer. For now i have added the bash script for miniconda to the Code/ folder. This version is for Mac OS X and accessed Oct 2018. **Copy and paste the following and it will download miniconda for you.**

bash ~/WAHMCC_metabarcoding_data_analysis/code/Miniconda3-latest-MacOSX-x86_64.sh

*   3.  **Adding vsearch** 
This is the free alternative to usearch (we didn't download Usearch because R. Edgar has not yet made his tool open source and applies a restriction for anything above 32 bit) - Vsearch can do most functions (not Unnoise3 algorithm for denoising reads) of Usearch but might be slower and lesser accuracy !! I Don’t know…

cd ~/bio_tools/
curl -LO https://github.com/torognes/vsearch/releases/download/v2.5.1/vsearch-2.5.1-macos-x86_64.tar.gz
tar -xzvf vsearch-2.5.1-macos-x86_64.tar.gz
rm vsearch-*.tar.gz
chmod +x vsearch/*
echo “PATH=$PWD/vsearch:$PATH” >> ~/.bash_profile

*   4.  **Adding fastQC** 
This is an important sequence QC program that can be run both CLI and GUI. It tells you what your raw Fastq files from illumina look like. you can basically set all your parameters using the information learned from it's output.

cd ~/bio_tools/
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
unzip fastqc_v0.11.7.zip
chmod +x FastQC/*
sudo ln -s $PWD/FastQC/fastqc /usr/local/bin/fastqc

*   5.  **Adding fastp**
This tool works the same way as fastqc but in my opinion it is better and has functionality to do some filtering.
*note: the fastp version in bioconda may be not the latest (bug fixes)*
conda install -c bioconda fastp
*Alternatively, you can download fastp using git clone and compile locally, as followss. make sure you're in the bio_tools/ dir.*

cd ~/bio_tools/
git clone https://github.com/OpenGene/fastp.git
cd fastp
make
make install
cd ..

*   5.  **Adding Trimmomatic**
This is a java tool - it is considered the trimming standard in the field of bioinformatics with regards to applying your parameters from fastQC output to your raw data.

cd ~/bio_tools/
curl -LO http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
unzip Trimmomatic-0.36.zip
rm Trimmomatic-0.36.zip
echo “PATH=$PWD/Trimmomatic_0.36/trimmomatic-0.36.jar:$PATH” >> ~/.bash_profile

*   6. **Adding BBtools** 
This tool is very handy for QC and other manipulations - would be great to know how to use these tools well because I think they can be easily run on HPC like Pawsey)

cd ~/bio_tools/
curl -L https://sourceforge.net/projects/bbmap/files/latest/download -o bbtools.tar.gz
tar -xzvf bbtools.tar.gz
rm -rf bbtools.tar.gz 
chmod +x bbmap/*
echo “PATH=$PWD/bbmap:$PATH” >> ~/.bash_profile


*   7.  **Adding a shi7**
This tool is my favourit it has excelent documentation, easy to use, speedy. When you install it will download a few tools that will work together nicely for data preprocessing. A tool in shi7 i like to use is shi7_learning.py, given the directory of the raw data it will take a good subset of your reads and learn everything you need to know to proparly clean your reads, such as adapters you know of and maybe some you don't, trimming parameters, merging and will give you a bash one liner using their shi7.py tool to get it all done with a results being a nicely cleaned and combined fasta output file. Another of their tools is "gotta_split", it is program that takes a MOTHUR formated dualbarcode .txt file and will demultiplex your paired-end file..very very fast and accurate.. I encourge you to vist the Knoghtslab github page, your life will get better.

cd ~/bio_tools/
wget https://github.com/knights-lab/shi7/releases/download/v0.9.9/shi7_0.9.9_mac_release.zip
unzip  shi7_0.9.9_mac_release.zip
rm shi7_0.9.9_mac_release.zip
chmod +x shi7_0.9.9_mac_release/*
echo “PATH=$PWD/shi7_0.9.9_mac_release:$PATH” >> ~/.bash_profile

*   8.  **Adding cutadapt**
This tool is essential when you need to trim and QC your raw fastq files. It is awesome at removing adapters and primers in a multithreaded fasion, yay. The catch is you'll need to install Python3 to use it's multithread flag(-j). 

cd ~/bio_tools/
conda install -c bioconda cutadapt

*Alternatively you can use* 

pip install --user --upgrade cutadapt

*   9.    **Adding seqkit & seqtk**
Both of these tools are extremely fast and can deal with a variety of raw sequence data like fastq/fasta in many different ways. They are also used for subsampling if you want to do that. recommend seqtk...

cd ~/bio_tools/
conda install -c bioconda seqkit
git clone https://github.com/lh3/seqtk.git;
cd seqtk; make
echo “PATH=$PWD/seqtk:$PATH” >> ~/.bash_profile

*   10.  **Adding Magic-Blast**
This tool is made specifically to blast short read sequences against the NCBI database, or *alternatively a localling installed refernce database formatted in blast format.* If you know how to use Blastn you can use Magic blast too.

cd ~/bio_tools/
curl -LO ftp://ftp.ncbi.nlm.nih.gov/blast/executables/magicblast/LATEST/ncbi-magicblast-1.3.0-x64-macosx.tar.gz
tar -xzvf ncbi-magicblast-1.3.0-x64-macosx.tar.gz
rm ncbi-magicblast-1.3.0-x64-macosx.tar.gz
cp ncbi-magicblast-1.3.0/bin/* . # moving the executables to our ~/happy_bin so they are in our PATH

magicblast -version # for me, magicblast: 1.3.0
source ~/.bash_profile

*I hope you will find it usefull and be aware that these tools may not be applicable to your particular project and by no means they are the best or only tools out there. these are just ones that have worked for me. **ENJOY***

### *Most of the code is recycled from other resources in the web. very small amount is what i have written my self.* 

#   END
