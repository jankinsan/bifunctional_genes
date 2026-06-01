#cd /scratch/bioschool/phd/blz227562/nm_nr
#module load compiler/R/4.1.2/intel2020
#R

#parallel work
library(parallel)

#get genomic data and load it!
hg38_dat <- read.delim("./data/GCF_000001405.40_GRCh38.p14_genomic.gtf", skip=5, header=FALSE)
head(hg38_dat)

#keep only genes
hg38_dat$gene <- sapply(hg38_dat[,9], function(x){
  str_trans <- strsplit(x, split=";")[[1]][1]
  split = unlist(strsplit(str_trans, split="gene_id "))
  return(split[length(split)])
})

hg38_dat$gene_biotype <- sapply(hg38_dat[,9], function(x){
  str_trans <- strsplit(x, split=";")[[1]]
  split = unlist(strsplit(str_trans[grep(" gene_biotype ", str_trans)], split="gene_biotype "))
  return(split[length(split)])
})

hg38_genes =  hg38_dat[which(hg38_dat[,3]=="gene"),]
colnames(hg38_genes) <- c("chr", "annotation", "type", "start", "end", "misc", "strand", "misc2", "details", "gene", "gene_biotype")

#getting simplified chromosome names for matching!
unique(hg38_genes$chr)
chr_dat <- read.delim("GCF_000001405.40_GRCh38.p14_assembly_report.txt", skip=63)
hg38_genes$chr <- chr_dat$UCSC.style.name[match(hg38_genes$chr, chr_dat$RefSeq.Accn)]

#keeping only bifunc genes
hg38_counts <- read.csv("23122023_GRCh38.p14_bifunc_gene_counts_status.csv")
gene_idx <- na.omit(unlist(sapply(hg38_counts$X, function(x){
  which(hg38_genes$gene==x)
})))
hg38_genes_bifunc <- hg38_genes[gene_idx,]
dim(hg38_genes_bifunc)

#match with circular RNAs
circRNA_dat <- read.delim("./data/human_bed_v3.0.txt")
head(circRNA_dat)

#matching and writing to file 
fileName="circRNA_bifunc"
circRNA_dat_matched <- data.frame()
for(chromosome in unique(hg38_genes_bifunc[,1])){
  print(paste0("chr: ", chromosome))
  chr_gene_pos <- hg38_genes_bifunc[intersect(which(hg38_genes_bifunc[,"chr"]==chromosome), which(hg38_genes_bifunc[,"strand"]=="+")),]
  chr_gene_neg <- hg38_genes_bifunc[intersect(which(hg38_genes_bifunc[,"chr"]==chromosome), which(hg38_genes_bifunc[,"strand"]=="-")),]
  chr_circ_pos <- circRNA_dat[intersect(which(circRNA_dat$Chro==chromosome), which(circRNA_dat$Strand=="+")),]
  chr_circ_neg <- circRNA_dat[intersect(which(circRNA_dat$Chro==chromosome), which(circRNA_dat$Strand=="-")),]
  
  #select miRNAs matching each NM/NR/XM/XR
  res_pos <- mclapply(1:dim(chr_gene_pos)[1], function(gene_idx){
    circ_idx_match <- intersect(which(chr_circ_pos$Start>= chr_gene_pos$start[gene_idx]),
                               which(chr_circ_pos$End<=chr_gene_pos$end[gene_idx]))
    if(length(circ_idx_match)>0){
      dat_matched <- data.frame()
      for(circ_idx in circ_idx_match){
        dat_matched <- rbind(dat_matched, data.frame(circRNA=chr_circ_pos$circAltas_ID[circ_idx],chr=chr_circ_pos$Chro[circ_idx],
                                                     strand=chr_circ_pos$Strand[circ_idx], circRNA_start=chr_circ_pos$Start[circ_idx],
                                                     circRNA_end=chr_circ_pos$End[circ_idx], gene=chr_gene_pos$gene[gene_idx], 
                                                     biotype_NCBI = unname(unlist(chr_gene_pos$gene_biotype[gene_idx])), gene_start=chr_gene_pos$start[gene_idx], 
                                                     gene_end=chr_gene_pos$end[gene_idx]))
      }
      return(dat_matched)
    } else{
      return(NA)
    }
  }, mc.cores=4)
  circRNA_dat_matched_pos <- do.call(rbind.data.frame, res_pos)
  
  res_neg<- mclapply(1:dim(chr_gene_neg)[1], function(gene_idx){
    circ_idx_match <- intersect(which(chr_circ_pos$Start>= chr_gene_neg$start[gene_idx]),
                                which(chr_circ_pos$End<=chr_gene_neg$end[gene_idx]))
    if(length(circ_idx_match)>0){
      dat_matched <- data.frame()
      for (circ_idx in circ_idx_match){
        dat_matched <- rbind(dat_matched, data.frame(circRNA=chr_circ_neg$circAltas_ID[circ_idx],chr=chr_circ_neg$Chro[circ_idx],
                                                     strand=chr_circ_neg$Strand[circ_idx], circRNA_start=chr_circ_neg$Start[circ_idx],
                                                     circRNA_end=chr_circ_neg$End[circ_idx], gene=chr_gene_neg$gene[gene_idx], 
                                                     biotype_NCBI = unname(unlist(chr_gene_neg$gene_biotype[gene_idx])), gene_start=chr_gene_neg$start[gene_idx], 
                                                     gene_end=chr_gene_neg$end[gene_idx]))
      }
      return(dat_matched)
    } else{
      return(NA)
    }}, mc.cores=4)
  circRNA_dat_matched_neg <- do.call(rbind.data.frame, res_neg)
  
  circRNA_dat_matched <- rbind.data.frame(circRNA_dat_matched,
                                      na.omit(circRNA_dat_matched_pos),
                                      na.omit(circRNA_dat_matched_neg))
}
print(paste0("saving file for", fileName))
write.table(circRNA_dat_matched, file=paste0(Sys.Date(), "_", fileName, "_matched_all_genes.csv"), row.names=FALSE, 
            sep=",")



#to run this after all processing! converts the data to atomic vectors and make the '$' operator problematic!
hg38_genes_bifunc <- apply(hg38_genes_bifunc,2,as.character)
write.table(hg38_genes_bifunc, file="23042024_bifunc_genes_coord_chrom_biotype.txt", row.names=FALSE)

hg38_genes <- apply(hg38_genes,2,as.character)
write.table(hg38_genes, file="23042024_genes_coord_chrom_biotype.txt", row.names=FALSE)


#circRNA check genes
circRNA_dat_matched <- read.csv("./circRNA_matched/2024-04-24_circRNA_bifunc_matched_all_genes.csv")
circRNA_count <- data.frame()
for (gene in unique(circRNA_dat_matched$gene)){
  circRNA_count<- rbind.data.frame(circRNA_count, 
                                   data.frame(gene, length(which(circRNA_dat_matched$gene==gene))))
}
