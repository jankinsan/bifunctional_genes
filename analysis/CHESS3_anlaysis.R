library(ggplot2)
library(ggrepel)

#CHESS match to NRs and expression?
chess_dat<- read.delim("C:/Users/DELL/Downloads/chess3.1.3.GRCh38.gtf.gz", header=FALSE)
colnames(chess_dat) <-c("chr", "source", "type", "start", "end", "misc", "strand", "misc2", "details")
  
sapply(chess_dat[1:5,9], function(x){strsplit(x, split=";")})

#get only transcript data!
chess_transcripts <- chess_dat[chess_dat$type=="transcript",]
strings <- data.frame()
idxs <- data.frame()
for (i in 1:15){
  x=chess_transcripts[i,9] 
  string <- strsplit(x, split="; ")
  strings <- rbind(strings, c(string[[1]]))
  idxs <- rbind(idxs, c(if(length(grep("gene_id", string[[1]]))>0) grep("gene_id", string[[1]]) else 0, #get gene id or send 0
                        if(length(grep("gene_name", string[[1]]))>0) grep("gene_name", string[[1]]) else 0, #get gene name or send 0
                        if(length(grep("db_xref", string[[1]]))>0) grep("db_xref", string[[1]]) else 0, #get accession from source sequence or 0 idx
                        if(length(grep("db_xref", string[[1]]))>0) grep("db_xref", string[[1]]) else 0))#get gene type or 0
}

chess_transcipts 

#copilot-generated code
# Define a function to parse the attributes column
parse_attributes <- function(attributes) {
  attr_list <- strsplit(attributes, ";")[[1]]
  attr_list <- trimws(attr_list)
  attr_list <- attr_list[attr_list != ""]
  attr_df <- data.frame(matrix(ncol = 2, nrow = length(attr_list)))
  colnames(attr_df) <- c("key", "value")
  for (i in 1:length(attr_list)) {
    key_value <- strsplit(attr_list[i], " ")[[1]]
    key <- key_value[1]
    value <- paste(key_value[-1], collapse = " ")
    attr_df$key[i] <- key
    attr_df$value[i] <- gsub('"', '', value)
  }
  return(attr_df)
}

# Apply the function to the last column of the data frame
attributes_list <- lapply(chess_transcripts[,9], parse_attributes)

# Get all unique keys
all_keys <- unique(unlist(lapply(attributes_list, function(df) df$key)))

# Create a data frame to hold the parsed attributes
attributes_df <- data.frame(matrix(ncol = length(all_keys), nrow = dim(chess_transcripts)[1]))
colnames(attributes_df) <- all_keys

# Fill the attributes data frame
for (i in 1:length(attributes_list)) {
  if (i%%1000 ==0){print(i)}
  for (j in 1:nrow(attributes_list[[i]])) {
    key <- attributes_list[[i]]$key[j]
    value <- attributes_list[[i]]$value[j]
    attributes_df[i, key] <- value
  }
}

# Combine the original data frame with the parsed attributes data frame
chess_transcripts <- cbind(chess_transcripts, attributes_df)

# Remove the original attributes column
chess_transcripts <- chess_transcripts[, -9]

# Display the resulting data frame
print(head(chess_transcripts))

#get only genes with NR/XR and NM/XM transcripts!
nr_xr_genes <- na.omit(unique(c(chess_transcripts$gene_name[grep("NR_", chess_transcripts$db_xref)],chess_transcripts$gene_name[grep("XR_", chess_transcripts$db_xref)])))
nm_xm_genes <- na.omit(unique(c(chess_transcripts$gene_name[grep("NM_", chess_transcripts$db_xref)],chess_transcripts$gene_name[grep("XM_", chess_transcripts$db_xref)])))

bifunc_genes <- intersect(nm_xm_genes, nr_xr_genes)

imp_dat<- matrix()
imp_dat <- do.call(rbind, lapply(bifunc_genes, function(gene){chess_transcripts[which(chess_transcripts$gene_name==gene),]}))

#add refseq and gencode or other transcript id's here
imp_dat$acc <- sapply(imp_dat$db_xref, function(x){strsplit(strsplit(x, split=",")[[1]][1], split="RefSeq:")[[1]][2]})
imp_dat$acc[is.na(imp_dat$db_xref)] <- imp_dat$transcript_id[is.na(imp_dat$db_xref)] # for chess only transcripts
imp_dat$acc[is.na(imp_dat$acc)] <- sapply(imp_dat$db_xref[is.na(imp_dat$acc)], function(x){strsplit(strsplit(x, split=",")[[1]][1], split="GENCODE:")[[1]][2]})


imp_dat<- imp_dat[, -20] # removing the tag column: not needed for now!
#keep only data columns without NAs (NA values: I wouldn't be able to plot)!
imp_dat_plot <- na.omit(imp_dat)

#wd for outputs
setwd("F:/Janki/bifunctional_genes")
dir.create("CHESS3")

#plot gene-wise to gt all important information in one place!
gene_plot <- function(gene_list){
  for (gene in gene_list){
  print(paste0("Running for ", gene))
  gene_dat <- imp_dat_plot[imp_dat_plot$gene_name==gene,c(11,20,18,19)]

  #if gene doesn't exist!
  if (dim(gene_dat)[1]==0){
    print("gene not in dataset!")
    next
  } 
  
  #add status for non-coding
  gene_dat$status <- "coding"
  gene_dat$status[grep("NR_", gene_dat$acc)] <- "noncoding"
  gene_dat$status[grep("XR_", gene_dat$acc)] <- "noncoding"
  
  gene_dat[,3:4]<- apply(gene_dat[,3:4], 2, as.numeric)
  print("saving plots!")
  plot_genes<- reshape2::melt(gene_dat[,-1], id=c("acc", "status"))
  pdf(paste0("./CHESS3/", gene, ".pdf"))
  
  if(gene %in% c("PLEKHA5", "TIA1")){
    barplot <- ggplot(plot_genes, aes(x=acc, y=value, fill=value))+ geom_bar(stat='identity') +
      facet_wrap(~variable,scales = "free_y", ncol = 1) + 
      labs(x="Transcripts", y="Expression", title=gene, fill="")+
      theme_bw()+theme(axis.text.x = element_text(size =10, angle = 90, vjust = 0.5, hjust=1, face="bold"),
                       axis.text.y = element_text(size =18, face="bold"),
                       strip.text.x = element_text(size = 18, face="bold"),
                       axis.title=element_text(size = 22, face="bold"),
                       plot.title=element_text(size = 22, face="bold", hjust=0.5))
    barplot2 <- ggplot(plot_genes, aes(x=acc, y=value, fill=status))+ geom_bar(stat='identity') +
      facet_wrap(~variable,scales = "free_y", ncol = 1) + 
      labs(x="Transcripts", y="Expression", title=gene, fill="Transcript Type")+
      theme_bw()+theme(axis.text.x = element_text(size =10, angle = 90, vjust = 0.5, hjust=1, face="bold"),
                       axis.text.y = element_text(size =18, face="bold"),
                       strip.text.x = element_text(size = 18, face="bold"),
                       legend.title = element_text(size = 18, face="bold"),
                       legend.text = element_text(size = 18, face="bold"),
                       axis.title=element_text(size = 22, face="bold"),
                       plot.title=element_text(size = 22, face="bold", hjust=0.5))+
      scale_fill_manual(labels = c("coding", "noncoding"),
                        values = c("#00bfc4", "#f8766d"))
    
  } else{
    barplot <- ggplot(plot_genes, aes(x=acc, y=value, fill=value))+ geom_bar(stat='identity') +
      facet_wrap(~variable,scales = "free_y", ncol = 1) + 
      labs(x="Transcripts", y="Expression", title=gene, fill="")+
      theme_bw()+theme(axis.text.x = element_text(size =18, angle = 90, vjust = 0.5, hjust=1, face="bold"),
                       axis.text.y = element_text(size =18, face="bold"),
                       strip.text.x = element_text(size = 18, face="bold"),
                       axis.title=element_text(size = 22, face="bold"),
                       plot.title=element_text(size = 22, face="bold", hjust=0.5))
    
    barplot2 <- ggplot(plot_genes, aes(x=acc, y=value, fill=status))+ geom_bar(stat='identity') +
      facet_wrap(~variable,scales = "free_y", ncol = 1) + 
      labs(x="Transcripts", y="Expression", title=gene, fill="Transcript Type")+
      theme_bw()+theme(axis.text.x = element_text(size =18, angle = 90, vjust = 0.5, hjust=1, face="bold"),
                       axis.text.y = element_text(size =18, face="bold"),
                       strip.text.x = element_text(size = 18, face="bold"),
                       legend.title = element_text(size = 18, face="bold"),
                       legend.text = element_text(size = 18, face="bold"),
                       axis.title=element_text(size = 22, face="bold"),
                       plot.title=element_text(size = 22, face="bold", hjust=0.5))+
      scale_fill_manual(labels = c("coding", "noncoding"),
                        values = c("#00bfc4", "#f8766d"))
  }
 
  plot(barplot)
  plot(barplot2)

  dev.off()
  
  ggsave(paste0("./CHESS3/", gene, "_exp.svg"), barplot, width = 7, height = 9)
  ggsave(paste0("./CHESS3/color_", gene, "_exp.svg"), barplot2, width = 7, height = 9)
}}
 
gene_plot(gene_list=c("CLASRP", "HNRNPDL", "KHDRBS1", "KHDRBS2", "NOP56", "PLEKHA5",
                      "RPS9", "SNX25", "SRSF1", "SRSF2", "TASP1", "UBP1", "VCPKMT", "SNX16", 
                      "XIAP", "TIA1", "USP15", "FOXP1", "TGFBR1", "DDX3Y", "FUBP1", "ACVR1B" ))

write.csv(imp_dat, file="./CHESS3/20012026_CHESS3_bifunc_genes_transcripts.csv")
imp_dat <- read.csv("./CHESS3/20012026_CHESS3_bifunc_genes_transcripts.csv")[,-1]
