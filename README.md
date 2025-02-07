# Installation.

``` 
conda create --name env bioconda::nanofilt bioconda::diamond bioconda::rasusa bioconda::seqkit bioconda::spades bioconda::minimap2 bioconda::samtools -y
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
Raw sequencing reads were deposited in the Sequence Read Archive (SRA) under BioProject accession number PRJNA1214984.

# Usage
After install all required packages, you can run the script with following command

```
./01_assembly.sh /path/to/fastq_folder output_name db
```
