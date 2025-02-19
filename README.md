# Installation & Requirements.

First, two enviroments must be created to avoid any conflicts from packages.

``` 
#1. first environment for genetal tasks
conda create --name env bioconda::nanofilt bioconda::diamond bioconda::rasusa bioconda::seqkit bioconda::minimap2 bioconda::samtools python=3.7.11 -y

#2. second enviroment for genome assembly
conda create -n SPades_env python=3.8.20 -y

#3. Then, download SPades from source
wget https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz
tar -xf SPAdes-4.0.0-Linux.tar.gz
```

Please install the following R packages before running the script.
1. msa (https://github.com/UBod/msa)
2. Biostrings (https://github.com/Bioconductor/Biostrings)
3. gtools (https://github.com/r-gregmisc/gtools)
4. data.table (https://github.com/Rdatatable/data.table)
5. dplyr (https://github.com/tidyverse/dplyr)
6. reshape2 (https://github.com/cran/reshape2)
7. ggplot2 (https://github.com/tidyverse/ggplot2)
8. stringr (https://github.com/tidyverse/stringr)


# Fastq data
Raw sequencing reads were deposited in the Sequence Read Archive (SRA) under BioProject accession number [PRJNA1214984](https://www.ncbi.nlm.nih.gov/sra/?term=SRR32105493).

# Usage
After installing all required packages, you can run the script with following command

```
#1. making diamond database
diamond makedb --in ref_AA.fasta -d db/ref_AA

#2. running assembly from main script
./01_assembly.sh /path/to/fastq_folder output_name db /path/to/SPades
```
