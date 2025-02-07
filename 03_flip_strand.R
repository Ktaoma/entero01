library(dplyr)
library(msa)
library(Biostrings)
library(reshape2)
library(stringr)

#https://stackoverflow.com/questions/48218332/how-to-output-results-of-msa-package-in-r-to-fasta
alignment2Fasta <- function(alignment, filename) {
  sink(filename)
 
  n <- length(rownames(alignment))
  for(i in seq(1, n)) {
    cat(paste0('>', rownames(alignment)[i]))
    cat('\n')
    the.sequence <- toString(unmasked(alignment)[[i]])
    cat(the.sequence)
    cat('\n')  
  }
 
  sink(NULL)
}

args <- commandArgs(trailingOnly = TRUE)

# Flip all contigs to plus strand 
df <- readAAStringSet(paste0(args[1],"/QC_contigs.fa"))
df_name <- do.call(rbind,names(df) %>% strsplit("_","")) %>% as.data.frame() 
df_DNA <- readDNAStringSet(paste0(args[1],"/filter_contigs.fa"))

all_DNA <- c()
for (i in 1:nrow(df_name)) {
  name_ <- df_name[i,1]
  frame_ <- df_name[i,2] %>% gsub("frame=","",.) %>% as.numeric()
  
  if (frame_ < 0) {
    DNA <- df_DNA[name_] %>% reverseComplement()
    all_DNA <- c(all_DNA,DNA)
  }else{
    all_DNA <- c(all_DNA,df_DNA[name_])
  }
}

do.call(c,all_DNA) %>% writeXStringSet(.,paste0(args[1],"/optimal_read_DNA_plus.fa"))

x_aln <- do.call(c,all_DNA) %>% msa(method="Muscle",verbose=T)
alignment2Fasta(x_aln,paste0(args[1],"/optimal_read_DNA_plus_align.fa"))

# Filter UTR region by self alignment
xcon <- consensusMatrix(x_aln,as.prob = T)
df_con <- xcon %>% as.data.frame()
idx <- row.names(df_con)
vec_consensus <- c()
for (x in 1:ncol(df_con)) {
  dft <- df_con[,x] %>% as.data.frame()
  dft$r <- rownames(dft)
  dft$code <- idx
  cnt <- dft %>% filter(. != 0)  
  
  max_ <- cnt %>% filter(. == max(.))
  
  if (max_$. %>% unique() >= 0.75) {
    
    print("haha")
    con <- cnt %>% filter(. == max(.))
    vec_consensus <- c(vec_consensus,con$code)  
    
  } else {
    
    print("hoho")
    vec_consensus <- c(vec_consensus,",")
    
  }
}

merge_vec <- paste(vec_consensus,collapse = "") %>% gsub("-","",.) %>% gsub(",","N",.)
mb <- merge_vec%>% DNAStringSet()
min_start = unlist(lapply(gregexpr(pattern = 'A|T|C|G', merge_vec), min))
max_end = unlist(lapply(gregexpr(pattern = 'A|T|C|G', merge_vec), max))
consensus_final <- mb %>% subseq(min_start,max_end)
names(consensus_final) <- args[2]
consensus_final %>% writeXStringSet(.,paste0(args[1],"/final_consensus_DNA_tmp.fa"))
  

# Filter UTR region by reference 
ref_ <- readDNAStringSet(paste0(args[1],"/ref_AA_best.fa"))
ref_name <- names(ref_) %>% strsplit(" ") %>% unlist()
db_cross <- read.table(paste0(args[3],"/cross_id.tsv"),header = T) %>% 
  filter(Accession == ref_name[1]) 
ref_DNA <- readDNAStringSet(paste0(args[3],"/ref_DNA.fasta"))
ref_db <- names(ref_DNA) %>% strsplit(" ") %>% do.call(rbind,.) %>% as.data.frame() 
ref_db$idx <- rownames(ref_db)
ref_idx <- ref_db %>% filter(str_detect(V1,db_cross$Nucleotide))
aln_both <- c(consensus_final,ref_DNA[as.numeric(ref_idx$idx)]) %>% msa(method="Muscle",verbose=T)
aln_df <- aln_both %>% as.matrix() %>% t() %>% as.data.frame()
colnames(aln_df) <- c("sample","ref")
aln_df %>% filter(ref != '-') %>% filter(sample != '-') %>% select(sample) %>% unlist() %>% as.vector() %>% paste(.,collapse = "")  -> curate_read
obj <- DNAStringSet(curate_read)
names(obj) <- args[2]

# Export final consensus
obj %>% writeXStringSet(.,paste0(args[1],"/final_consensus_DNA.fa"))

