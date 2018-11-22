#!/bin/sh

# BaseSpace sequencing run downloader

# This shell script will download a paired-end fastq.gz files from a particular sequencing run conducted in the TrEnD_Lab. To successfully download your raw data you will need the r2 fastq file ID number that you can obtain from the BaseSpace website. Ideally you will click on the run tab on the top of the BaseSpace page (after logining in) and then select your particular run and once you're in the run page click on the sample id with your library name. Then a page will open, please scroll down and you will see your files in the bottom right of teh page. click on the R2 file and then copy the id number on the URL that appears above. Once obtained please paste that id number in the part of this script that says #r2 id=.

echo "start at:`date`"
cd ~/pipeline_folder/raw_fq #Change me to where you want the raw files to go in your computer.
runname=Lig25 # change me to the name that appears at the beginning of your fastq_files or library.
token=XXXXXXXXXXXXXXXXXXXXXX # this special number is unique, so if you want to use this code to access a BaseSpace account please change to the token generated for that account.

#r2
id=11737184383

for readtype in R2 R1
do
if [ $readtype == "R1" ]
then
id="$(($id + 1))"
fi

wget -O ${runname}.${readtype}.fastq.gz 'https://api.basespace.illumina.com/v1pre3/files/'$id'/content?access_token='$token
gunzip ${runname}.${readtype}.fastq.gz
done

echo "End at:`date`"


