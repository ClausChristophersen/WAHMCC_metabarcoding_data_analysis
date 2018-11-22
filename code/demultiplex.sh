#!/bin/sh

#Paired-end fastq file demultiplex script. This script depends
#on the barcodes.txt file being correctly formatted in the form of the Mothur barcode file.

#Demultiplexing using Shi7
gotta_split ./raw_fq/Lig25.R1.fastq ./raw_fq/Lig25.R2.fastq ../mapping_file/LigationLib_mapping_files.txt ../demux_fq/ --deep
#the --deep will tell the demultiplexer to also look into the reads. you can also add --errors with the desired mismatch number allowed. I recommend you don't use this flag but if you must because you're not demultiplexing enough reads/sample. Then make sure that the hamming distances (sequence difference) between individual barcodes is sufficently far. So as the uniqueness of your barcodes is maintained even if you allow for a mismatch. For combination tagging approach, such that the combination is unique not each single barcode alone then gotta_split may not be for you (you can try anyway) and you should consider "SeekDeep or AXE" or write your own (make it multithread or you'll suffer). If you use Python3 for the job then import the "multiprocessing library".


#End
