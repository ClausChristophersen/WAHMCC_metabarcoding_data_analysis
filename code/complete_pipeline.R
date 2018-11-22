############################################################################################################################Using TrEnD lab workflow and DADA2######################################

title: "V4_vaginal_analysis"
author: "Alishum_Ali"
date: "18/11/2018"
output:
  html_document: default
  pdf_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown


**NOTE**
When you download and open RStudio it will automatically navigate to your $HOME directory. If you'd rather another spot, then you can do that manually in the bottom right pane. There click on Files and navigate the directory to your prefered folder.

# WE BEGIN

The first thing i would usually do is set up the directory tree for the project.  
```{r setting up project folder, warning=FALSE, message=FALSE, tidy=TRUE, results='hold'}

dir.create("~/pipeline_folder") 

subfolder_names <- c("raw_fq", "QC_data", "GHAP_demux", "mapping_file", "code", "primer-trim_fq", "adapt-trim_fq", "DADA_output") 

for (j in 1:length(subfolder_names)){
  folder<-dir.create(paste0("~/pipeline_folder/", subfolder_names[j]))
}

setwd("~/pipeline_folder/")
path <- "~/pipeline_folder/GHAP_demux"
list.files(path)
```

The mapping files & code files can be downloaded from cloudstor use the link i sent in the email. In the folder downloaded i have provided the codes needed, the trend lab's barcode sequence list, primers/adapters sequence files and the read mapping file to demultiplex the fastq files that we'll use for this workshop.  
If you have downloaded the folder, just adjust the file path in the code below.

```{r copying files to project folder,  warning=FALSE, message=FALSE, tidy=TRUE, results='hold'}
file.copy(c("~/Project_folder/mapping_file", "~/Project_folder/code"), "~/pipeline_folder")
```
# Adding all the tools you'll need for metabarcoding in your computer.

# First of all we need to set your computer up so it has all the tools that we need for all the analysis. We do that as such. Remove the # to run code.

```{r setting up computer, warning=FALSE, message=FALSE, tidy=TRUE, results='hold'}
bash <- "/usr/bin/bash"
#system2(bash, "./code/Bioinformatics_tools.sh")
```

# Downloading fastq.gz files from BaseSpace into "raw_fq" folder. This will take about 10 min

```{r fastq download, warning=FALSE, message=FALSE, tidy=TRUE, results='hold'}
system2(bash, "./code/raw_fq_download.sh")
```

# Pre-Demultiplexing the fastq files we have to do a couple of crucial preprocessing steps.

We need to unzip our fastq.gz files, this may take a bit of space but it will make it easy for all tools to use as input.

```{r fastq gunzip, warning=FALSE, message=FALSE, tidy=TRUE, results='hold'}
gunzip <- "/usr/bin/gunzip"
system2(gunzip, args = "./raw_fq/*fastq")
```

Run fastqc on the raw sequences and assess what you are dealing with. We can do that in the following way.

```{r quality check pre trimming, warning=FALSE, message=FALSE, tidy=TRUE, results='hold'}
fastqc <- "/usr/local/bin/fastqc"

system2(fastqc, args = c("./raw_fq/Lig25*.fastq", "-f", "fastq", "-t", "14", "-q", "-o", "./QC_data")) # all you need is the html format in your QC_data folder
```

To maximise recovery of demultiplexed files we will remove any adapters at the begining of our reads using the cutadapt tool.

#Adapter removal using cutadapt
```{r removing adapters, warning=FALSE, message=FALSE, tidy=TRUE, results='hold'}

cutadapt <- "/Users/18746334/miniconda3/bin/cutadapt"
adapters <- c("-a", "file:/Users/18746334/pipeline_folder/mapping_files/adapters.fa", "-A", "file:/Users/18746334/pipeline_folder/mapping_files/adapters.fa", "-g", "file:/Users/18746334/pipeline_folder/mapping_files/adapters.fa", "-G", "file:/Users/18746334/pipeline_folder/mapping_files/adapters.fa")

system2(cutadapt, args = c(adapters, "--quiet", "-e", "0.05", "-O", "8", "--trim-n", "-j", "14", "-o", "./raw_fq/cLig25_R1.fastq", "-p", "./raw_fq/cLig25_R2.fastq", "./raw_fq/Lig25.R1.fastq", "./raw_fq/Lig25.R1.fastq"))

```

# Now we confirm that we have removed the first of the non-biological sequences in our fastq data.
```{r quality check post trimming, warning=FALSE, message=FALSE, tidy=TRUE, results='hold'}
system2(fastqc, args = c("./raw_fq/cLig25*.fastq", "-f", "fastq", "-t", "10", "-q", "-o", "./QC_data")) # all you need is the html format in your QC_data folder
```

# Now we demultiplex using the GHAP pipeline module
```{r demultiplexing, needed, warning=FALSE, message=FALSE, tidy=TRUE, results='hold'}
demultiplexer <- "/Users/18746334/GHAPv2/DemultiplexCustomBC"
system2(demultiplexer, args = c("-m", "./mapping_files/rLig25_mapping_file.txt", "-name", "16S", "-merge", "-o", "./primer-trim_fq", "-id", "0", "-fbc", "1", "-fpc", "2", "-rbc", "3", "-rpc", "4", "./raw_fq/Lig25*.fastq"))
```

# Loading in the packages we need for this workshop
```{r packages, warning=FALSE, message=FALSE, tidy=TRUE, results='hold'}
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("dada2", version = "3.8")

source("http://bioconductor.org/biocLite.R")
biocLite(suppressUpdates = FALSE)
biocLite("ShortRead", suppressUpdates = FALSE)
biocLite("Biostrings", suppressUpdates = FALSE)
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

list.files(path)
```
Here we specify the full path where our raw data is situated

```{r path and raw files}
path <- "/Users/18746334/pipeline_folder/GHAP_demux"
list.files(path)
fnFs <- sort(list.files(path, pattern = "dm_16S_Lig25.R1_", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "dm_16S_Lig25.R2_", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs),"dm_16S_Lig25.R1_"),'[',2)
sample.names ## Lists the names of all the sequence files
```
# Starting primer removal steps

```{r primers}
FWD <- "CWACGCGARGAACCTTACC"
REV <- "ACRACACGAGCTGACGAC"
```

## Primer orientation and removal function

```{r orient&remove, warning=FALSE, message=FALSE, tidy=TRUE}
allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
    RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}
```
### reading an example of our fastq files to count the presence of primers

```{r primer counts, warning=FALSE, message=FALSE, tidy=TRUE}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]),
FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]),
REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]),
REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))
```

### The actual removal step
in the next code chunk we do a coupleof things, we specify teh path to the tool that will do the trimming and then give it the flags specific to our files. All that will be put together in a loop to get all files cleaned and stored in a folder called "cutadapt".

```{r cutadapt tool, warning=FALSE, message=FALSE, tidy=TRUE, comment="NA"}

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, paste0(sample.names,"_R1.fastq"))
fnRs.cut <- file.path(path.cut, paste0(sample.names,"_R2.fastq"))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)
# Run Cutadapt
for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, "-j", "14", "-m", "100",
    "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
    fnFs[i], fnRs[i])) # input files
}
```

```{r removal confirmation, warning=FALSE, message=FALSE, tidy=TRUE}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
```

# DADA2 pipeline
Now that we have cleaned our fastq files from primer contamination we can happily use the dada2 pipeline without fear.

```{r Import files, warning=FALSE, message=FALSE, tidy=TRUE}
#Forward and reverse fastq filenames have the format: 
fnFs <- sort(list.files(path.cut, pattern = "R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path.cut, pattern = "R2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: 

sample.names <- sapply(strsplit(basename(fnFs.cut),"."),'[',1)
sample.names
```

##Inspect read quality profiles

We start by visualizing the quality profiles of the forward reads:

```{r Quality Profile forward, warning=FALSE, message=FALSE, tidy=TRUE}
plotQualityProfile(fnFs.cut[1:2])
```

We visualize the quality profile of the reverse reads:

```{r Quality Profile reverse, warning=FALSE, message=FALSE, tidy=TRUE}
plotQualityProfile(fnRs.cut[1:2])
```

##Filter and trim

Assigning the filenames for the output of the filtered reads to be stored as fastq.gz files. 

```{r filt_output, warning=FALSE, message=FALSE, tidy=TRUE}
filtFs <- file.path(path.cut, "filtered", paste0(sample.names, "_R1_filt.fastq.gz"))
filtRs <- file.path(path.cut, "filtered", paste0(sample.names, "_R2_filt.fastq.gz"))
```

```{r Filter and Trim, warning=FALSE, message=FALSE, tidy=TRUE}
out<-filterAndTrim(fnFs.cut, filtFs, fnRs.cut, filtRs, maxN = 0, maxEE = c(2,2), truncQ =2, truncLen = c(175, 220), rm.phix = TRUE, compress = TRUE, multithread = TRUE) 
(out)
```

## Learning the errors in our reads

```{r Error rate, warning=FALSE, message=FALSE, tidy=TRUE}
errF <- learnErrors(filtFs,multithread = TRUE)
```

```{r Error rate reverse, warning=FALSE, message=FALSE, tidy=TRUE}
errR <- learnErrors(filtRs,multithread = TRUE)
```

As a sanity check, it is worth visualizing the estimated error rates:

```{r Plot Quality Profile, warning=FALSE, message=FALSE, tidy=TRUE}
plotErrors(errF, nominalQ=TRUE)
```

```{r Plot Quality Profile, warning=FALSE, message=FALSE, tidy=TRUE}
plotErrors(errF, nominalQ=TRUE)
```

##Dereplication

Dereplication combines all identical sequencing reads into their “unique sequences” with a corresponding “abundance” equal to the number of reads with that unique sequence. Dereplication substantially reduces computation time by eliminating redundant comparisons.

Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: DADA2 *retains a summary of the quality information associated with each unique sequence*. The consensus quality profile of a unique sequence is the average of the positional qualities from the dereplicated reads. These quality profiles inform the error model of the subsequent sample inference step, significantly increasing DADA2’s accuracy.


```{r Dereplication, warning=FALSE, message=FALSE, tidy=TRUE}
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
```

##Sample Inference

At this step, [the core sample inference algorithm](https://www.nature.com/articles/nmeth.3869#methods) is applied to the dereplicated data. 

```{r dada2, warning=FALSE, message=FALSE, tidy=TRUE}
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```

Inspecting the returned `data-class` object:

```{r dada2 object, warning=FALSE, message=FALSE, tidy=TRUE}
dadaFs[1]
```

The DADA2 algorithm inferred 76 true sequence variants from the 61773 unique sequences in the first sample. There is much more to the `dada-class` return object than this (see help(`"dada-class"`) for some info), including multiple diagnostics about the quality of each denoised sequence variant, but that is beyond the scope of an introductory tutorial.

##Merge paired reads

After applying the core sample inference algorithm, the forward aand reverse reads are merged to obtain the full denoise sequences. Merging is performed by aligning the denoised forward reads with the reverse-complement of the corresponding denoised reverse reads, and then constructing the merged "contig" sequences. By default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical to each other in the overlap region. 


```{r Mergers, warning=FALSE, message=FALSE}
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```

The `mergers` object is a list of `data.frames` from each sample. Each `data.frame` contains the merged `$sequence`, its `$abundance`, and the indices of the `$forward` and `$reverse` sequence variants that were merged. Paired reads that did not exactly overlap were removed by `mergePairs`, further reducing spurious output.

##Construct Sequence Table

We can now construct an amplicon sequence variant table (ASV) table, a higher-resolution version of the OTU table produced by traditional methods.

```{r Seqtab, warning=FALSE, message=FALSE, tidy=TRUE}
seqtab <- makeSequenceTable(mergers)
```


```{r seqtab dimension, warning=FALSE,message=FALSE, tidy=TRUE}
dim(seqtab)
```


```{r Distribution of lengths, warning=FALSE,message=FALSE, tidy=TRUE}
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
```


The sequence table is a `matrix` with rows corresponding to (and named by) the samples, and columns corresponding to (and named by) the sequence variants. This table contains 41 ASVs. 

##Remove Chimeras

The core `dada` method corrects substitution and indel errors, but chimeras remain. Fortunately, the accuracy of the sequence variants after denoising makes identifying chimeras simpler than it is when dealing with fuzzy OTUs. Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant “parent” sequences.

```{r seqtab chimera removal, warning=FALSE, message=FALSE}
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
```

```{r seqtab sum, warning=FALSE, message=FALSE, tidy=TRUE}
sum(seqtab.nochim)/sum(seqtab)
```

The frequency of chimeric sequences varies substantially from dataset to dataset, and depends on on factors including experimental procedures and sample complexity. Here chimeras make up about 32% of the merged sequence variants, but when we account for the abundances of those variants we see they account for less than 1% of the merged sequence reads.


##Track reads through the pipeline

As a final check of our progress, we will look at the number of reads that made it through each step in the pipeline: 

```{r Track reads, warning=FALSE, message=FALSE, tidy=TRUE}
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

##Assigning taxonomy with the newer more accurate DECIPHER package using it's IdTaxa module


```{r loading iDTAXA}
library(DECIPHER); packageVersion("DECIPHER")
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("~/Fasta_databases/DADA2_db/GTDB_r86-mod_September2018.rdata") 
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
taxa <- taxid
```

Inspecting the taxonomic assignments:

```{r Taxonomy inspection, warning=FALSE, message=FALSE, tidy=TRUE}
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

# Writing out files we will need for BLAST/MEGAN and External statistics work
```{r Constructing nicer count and taxa tables, warning=FALSE, message=FALSE, tidy=TRUE}
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(t(seqtab.nochim))[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "DADA_output/ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "DADA_output/ASVs_counts.txt", sep="\t", quote=F)

# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "DADA_output/ASVs_taxonomy.txt", sep="\t", quote=F)
```


# Exploratory analysis
```{r packages for analysis, message=FALSE, warning=FALSE, tidy=TRUE }
library("knitr")
library("BiocStyle")
.cran_packages <- c("ggplot2", "gridExtra", "devtools")
.bioc_packages <- c("phyloseq", "DECIPHER", "phangorn")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
    install.packages(.cran_packages[!.inst])
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite(.bioc_packages[!.inst], ask = F)
}
# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
library(devtools)
devtools::install_github("benjjneb/decontam")
library(ggplot2)
library(phyloseq)
library(gridExtra)
library(phangorn)
library(decontam)

theme_set(theme_bw())
```
## Creating phyloseq object

```{r Phyloseq object, message=FALSE, warning=FALSE, tidy=TRUE}
samples.out <- rownames(seqtab.nochim)
samdf <- read.table(file= "v4SamplesInfo.txt", header = TRUE)
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

(ps)
```
Removing ATCC standard from further analysis

```{r removing standards, message=FALSE, warning=FALSE, tidy=TRUE}
ps1 <- prune_samples(sample_names(ps) != "ATCC_v4", ps) 
(ps1)
```

#decontamination

```{r using decontam algorithm, message=FALSE, warning=FALSE, tidy=TRUE}
df <- as.data.frame(sample_data(ps1)) # putting our sample metadata into a ggplot form.
df$LibrarySize <- sample_sums(ps1)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=qubit)) + geom_point()

summary(df)
## Identifying contaminating species using the "either" method of decontam.

sample_data(ps1)$is.neg <- sample_data(ps1)$SampleID == c("SCExt_v4", "MBExt_v4")
contamdf.freq <- isContaminant(ps1, method="either", neg="is.neg", conc="rfu", threshold=0.1)
head(contamdf.freq)
table(contamdf.freq$contaminant)
head(which(contamdf.freq$contaminant))

plot_frequency(ps1, taxa_names(ps1)[which(contamdf.freq$contaminant)], conc="rfu") + 
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")

set.seed(100)

plot_frequency(ps1, taxa_names(ps1)[sample(which(contamdf.freq$contaminant))], 
               conc="rfu") + xlab("DNA Concentration (PicoGreen fluorescent intensity)")

# Removing contaminants from the phyloseq object

ps.noncontam <- prune_taxa(!contamdf.freq$contaminant, ps1)
ps.noncontam # now you can take this object and do some analysis using phyloseq package

ps.noncontam1 <- prune_samples(sample_names(ps.noncontam) != "SCExt_v4", ps.noncontam)
ps.noncontam1 <- prune_samples(sample_names(ps.noncontam1) != "MBExt_v4", ps.noncontam1)
ps.noncontam1
```

#Plot alpha & beta diversity

```{r Plotting diversity, message=FALSE, warning=FALSE, tidy=TRUE}
# plotting all samples
plot_richness(ps, x="samples", measures=c("Observed", "Shannon", "Simpson"), shape = "centrifuge", color="extraction", title = "DNA extraction methodology & vaginal microbiome alpha diversity using 16S rRNA V4 target region")
```
# plotting decontaminated and samples only
```{r alpha diversity}
plot_richness(ps.noncontam1, x="samples", measures=c("Observed", "Shannon", "Simpson"), shape = "centrifuge", color="extraction", title = "DNA extraction methodology & vaginal microbiome alpha diversity using 16S rRNA V4 target region (samplesOnly)")
```
# Plot beta diversity
```{r plotting Beta diversity, }
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="extraction", label = "SampleID", title="Beta diversity Bray-NMDS all samples")

ps.prop <- transform_sample_counts(ps.noncontam, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="extraction", label = "SampleID", title="Beta diversity Bray NMDS (no Mock)")

ps.prop <- transform_sample_counts(ps.noncontam1, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="extraction", label = "SampleID", title="V4 Vaginal swab Beta diversity Bray NMDS (samplesOnly)")
```

```{r other plots,  message=FALSE, warning=FALSE, tidy=TRUE}
# Bar chart with all samples
top100 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:100]
ps.top100 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top100 <- prune_taxa(top100, ps.top100)
plot_bar(ps.top100, x="SampleID", fill="Genus", title = "Relative abundance at species level V4 AllSamples") + facet_wrap(~extraction, scales="free_x")


# Bar chart with no mock & blanks Genus
top100 <- names(sort(taxa_sums(ps.noncontam1), decreasing=TRUE))[1:100]
ps.top100 <- transform_sample_counts(ps.noncontam1, function(OTU) OTU/sum(OTU))

ps.top100 <- prune_taxa(top100, ps.top100)
plot_bar(ps.top100, x="SampleID", fill="Genus", title = "Relative abundance at species level V4 SamplesOnly") + facet_wrap(~extraction, scales="free_x")

#Plot Network graph everything
g <- make_network(ps, type="samples", distance="bray", max.dist = 0.9, 
                  keep.isolates=FALSE)
plot_network(g, ps, type="samples", 
             color="centrifuge", shape="extraction", point_size=4, alpha=1,
             label="SampleID", hjust = 1.35 , 
             line_weight=0.5, line_alpha=0.4,
             title= "MiSeq V4 Vaginal swab extraction protocol comparison")

#Plot Network graph no blanks or Mock
g <- make_network(ps.noncontam1, type="samples", distance="bray", max.dist = 0.4, 
                  keep.isolates=FALSE)
plot_network(g, ps.noncontam1, type="samples", 
             color="centrifuge", shape="extraction", point_size=4, alpha=1,
             label="SampleID", hjust = 1.35 , 
             line_weight=0.5, line_alpha=0.4,
             title= "MiSeq V4 Vaginal swab extraction protocol comparison")


# plot scree to explain axis difference
plot_ordination(ps.noncontam1, ordinate(ps.noncontam1, "DCA"), type="scree")

# Heatmap with no cleaning what it would look like.
plot_heatmap(ps, sample.label="SampleID", species.label="Genus") 

# Heatmap with cleaning what it would look like.
plot_heatmap(ps.noncontam1, sample.label="SampleID", species.label="Genus")


```

##**Finished, now you can go ahead and do some real statistics**
