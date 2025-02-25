library(dplyr)
library(msa)
library(Biostrings)
library(reshape2)
library(ggplot2)

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


#Import AA from the complete contigs
args <- commandArgs(trailingOnly = TRUE)
AA_seq <- readAAStringSet(paste0(args[1],"/contigs_full_length_AA.fa")) %>% msa()
alignment2Fasta(AA_seq,paste0(args[1],"/candidate_plot.fa"))
x <- readAAStringSet(paste0(args[1],"/candidate_plot.fa")) %>% as.matrix()  %>% as.data.frame() 

#Assign the type of AA based on frequency 
for (u in 1:ncol(x)) {
  
  x2 <- table(x[,u]) %>% as.data.frame()
  
  x2$type <- ifelse(x2$Freq == max(x2$Freq),"REF",
                    ifelse(x2$Var1 == '-',"DEL",
                           ifelse(x2$Var1 == 'X','X',"ALT")))
  
  x3 <- x[,u] %>% as.data.frame()
  colnames(x3) <- "Var1"
  
  x4 <- inner_join(x3,x2)
  x[,u] <- x4$type
}

#Curate result
x$sample <- rownames(x)
x_plot <- x %>% melt(id.vars = "sample")
x_vec <- x_plot$variable %>% unlist() %>% as.vector()
x_plot$dis <- gsub("V","",x_vec) %>% as.numeric()

xt <- table(x_plot$sample,x_plot$value) %>% as.data.frame() %>% dcast(Var1~Var2)
xt_ <- xt %>% select(-REF)
if (ncol(xt_) == 1) {
  xt$sum <- 0
}else{
  xt$sum <- rowSums(xt_[2:ncol(xt_)])  
}

#Export the AA figure
reor <- xt %>% arrange(desc(sum),ascending=F) %>% select(Var1) %>% unlist() %>% as.vector()
x_plot$sample <- factor(x_plot$sample,levels=reor)
x_plot$value <- factor(x_plot$value,levels=c("REF","ALT","DEL","X"))

x_plot  %>%
  ggplot(aes(dis, sample)) +
  geom_tile(aes(fill = value),width=3, height=0.9)+
  scale_fill_manual(values = c('DEL' = 'brown1', 'REF' = 'grey80',"ALT" = 'royalblue',"X"='black'))+
  theme_minimal()+
  xlab("AA position")+
  ylab("")+
  labs(fill="AA types")

ggsave(paste0(args[1],"/AA_plot.pdf"))
