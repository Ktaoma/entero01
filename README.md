# Installation

First, two environments must be created to avoid any conflicts among softwares.

``` 
#1. first environment for genetal tasks
conda create --name General_env python=3.7.11 bioconda::nanofilt bioconda::diamond bioconda::rasusa bioconda::seqkit bioconda::minimap2 bioconda::samtools -y

#2. second enviroment for genome assembly task
conda create --name SPAdes_env python=3.8.20 -y

#3. install SPAdes software for second enviroment
wget https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz
tar -xf SPAdes-4.0.0-Linux.tar.gz
```

Also, please install the following R packages before running the main script.
```
install.packages(c("dplyr","gtools","data.table","ggplot2","reshape2","stringr","BiocManager"))
BiocManager::install("msa")
BiocManager::install("Biostrings")
```

# Fastq data
Raw sequencing reads of Coxsackievirus A16 are available under BioProject accession number [PRJNA1214984](https://www.ncbi.nlm.nih.gov/sra/?term=SRR32105493).

# Usage
After installing all required packages, you can run the script with following command

```
#1. making diamond database
gzip -d db/*
diamond makedb --in db/ref_AA.fasta -d db/ref_AA

#2. running assembly from main script
./01_assembly.sh /path/to/fastq_folder output_name db /path/to/SPAdes
```
