#cd /scratch/bioschool/phd/blz227562/nm_nr
#module load compiler/R/4.1.2/intel2020
#R

#parallel work
library(parallel)

#get genomic data and load it!
hg38_dat <- read.delim("GCF_000001405.40_GRCh38.p14_genomic.gtf", skip=5, header=FALSE)
hg38_dat$acc <- sapply(hg38_dat[,9], function(x){
  str_trans <- strsplit(x, split=";")[[1]][2]
  split = unlist(strsplit(str_trans, split="transcript_id "))
  return(split[length(split)])
})
hg38_dat$gene <- sapply(hg38_dat[,9], function(x){
  str_trans <- strsplit(x, split=";")[[1]][1]
  split = unlist(strsplit(str_trans, split="gene_id "))
  return(split[length(split)])
})

print(head(hg38_dat))
print(dim(hg38_dat))
#removing the last row since it's useless and coming from the gtf format!
hg38_dat <- hg38_dat[-4693812,]

#separating only exon data since that's what we'll be matching TargetScan and microCLIP targets to!
hg38_exons<- hg38_dat[which(hg38_dat[,3]=="exon"),]
colnames(hg38_exons) <- c("chr", "annotation", "type", "start", "end", "misc", "strand", "misc2", "details", "acc", "gene")

#get accessions from my counts document
#match them with the type of gene!
acc<- read.csv("18122023_GRCh38.p14_gene_accessions.csv")
hg38_counts <- read.csv("18122023_GRCh38.p14_genewise_counts_status.csv")
acc$ids <- sapply(acc$accession, function(x){
  substr(x, 2, nchar(x))
})
acc$gene_type <- hg38_counts$status[match(acc$gene, hg38_counts[,1])]
print(head(acc))
hg38_exons$gene_type <- acc$gene_type[match(hg38_exons$acc, acc$ids)]
print(head(hg38_exons))
#replaces NA for gene_type with "not_assigned"
hg38_exons$gene_type[is.na(hg38_exons$gene_type)] = "not_assigned"
length(unique(hg38_exons$acc[-which(hg38_exons$gene_type=="not_assigned")]))

#getting simplified chromosome names for matching!
unique(hg38_exons$chr)
chr_dat <- read.delim("GCF_000001405.40_GRCh38.p14_assembly_report.txt", skip=63)
hg38_exons$chr <- chr_dat$UCSC.style.name[match(hg38_exons$chr, chr_dat$RefSeq.Accn)]


#run the matching function!
mirMatchComplete <- function(mir_dat, fileName){
  mir_dat_matched <- data.frame()
  for(chromosome in unique(mir_dat$chr)){
    print(paste0("chr: ", chromosome))
    chr_acc_pos <- hg38_exons[intersect(which(hg38_exons$chr==chromosome), which(hg38_exons$strand=="+")),]
    chr_acc_neg <- hg38_exons[intersect(which(hg38_exons$chr==chromosome), which(hg38_exons$strand=="-")),]
    chr_mir_pos <- mir_dat[intersect(which(mir_dat$chr==chromosome), which(mir_dat$strand=="+")),]
    chr_mir_neg <- mir_dat[intersect(which(mir_dat$chr==chromosome), which(mir_dat$strand=="-")),]
    
    #select miRNAs matching each NM/NR/XM/XR
    res_pos <- mclapply(1:dim(chr_mir_pos)[1], function(mir_idx){
      mir_idx_match <- intersect(which(chr_acc_pos$start<chr_mir_pos$start[mir_idx]),
                                 which(chr_acc_pos$end>chr_mir_pos$end[mir_idx]))
      if(length(mir_idx_match)>0){
        dat_matched <- data.frame()
        for(x in mir_idx_match){
          dat_matched <- rbind(dat_matched, data.frame(microRNA=chr_mir_pos$name[mir_idx],chr=chr_mir_pos$chr[mir_idx],
                     strand=chr_mir_pos$strand[mir_idx], microRNA_start=chr_mir_pos$start[mir_idx],
                     microRNA_end=chr_mir_pos$end[mir_idx],accession=chr_acc_pos$acc[x],
                     gene=chr_acc_pos$gene[x], gene_type = chr_acc_pos$gene_type[x],
                     transcript_start=chr_acc_pos$start[x], transcript_end=chr_acc_pos$end[x],
                     score=chr_mir_pos$score[mir_idx]))
        }
        return(dat_matched)
      } else{
        return(NA)
      }
    }, mc.cores=6)
    mir_dat_matched_pos <- do.call(rbind.data.frame, res_pos)

    res_neg<- mclapply(1:dim(chr_mir_neg)[1], function(mir_idx){
      mir_idx_match <- intersect(which(chr_acc_neg$start<chr_mir_neg$start[mir_idx]),
                                 which(chr_acc_neg$end>chr_mir_neg$end[mir_idx]))
      if(length(mir_idx_match)>0){
        dat_matched <- data.frame()
        for (x in mir_idx_match){
          dat_matched <- rbind(dat_matched, data.frame(microRNA=chr_mir_neg$name[mir_idx], chr=chr_mir_neg$chr[mir_idx],
                            strand=chr_mir_neg$strand[mir_idx], microRNA_start=chr_mir_neg$start[mir_idx],
                            microRNA_end=chr_mir_neg$end[mir_idx], accession=chr_acc_neg$acc[x],
                            gene=chr_acc_neg$gene[x], gene_type = chr_acc_pos$gene_type[x],
                            transcript_start=chr_acc_pos$start[x], transcript_end=chr_acc_pos$end[x],
                            score = chr_mir_neg$score[mir_idx]))
        }
        return(dat_matched)
      } else{
        return(NA)
      }}, mc.cores=6)
    mir_dat_matched_neg <- do.call(rbind.data.frame, res_neg)
    
    mir_dat_matched <- rbind.data.frame(mir_dat_matched,
                                        na.omit(mir_dat_matched_pos),
                                        na.omit(mir_dat_matched_neg))
  }
  print(paste0("saving file for", fileName))
  write.table(mir_dat_matched, file=paste0(Sys.Date(), "_", fileName, "_matched_all_genes.csv"), row.names=FALSE, 
              sep=",")
  return(mir_dat_matched)
}

#loading and running!
microCLIP_pre <- read.delim("16122023_microCLIPPredicted_liftover_hg38.bed", header=FALSE)
colnames(microCLIP_pre)<- c("chr", "start", "end", "name", "score", "strand", "method")
microCLIP_trial<- mirMatchComplete(microCLIP_pre[1:20000,], file="microCLIP_trial")

TS_default_pred <- read.delim("26022024_TargetScan_predicted_default_hg38.bed", sep="\t", header=FALSE)
colnames(TS_default_pred)<- c("chr", "start", "end", "name", "score", "strand", 
                              "thickStart", "thickEnd", "itemRgb", "blockCount", 
                              "blockSizes", "blockStarts")
TS_default_matched <- mirMatchComplete(TS_default_pred, fileName= "TS_default_pred")
microCLIP_matched <- mirMatchComplete(microCLIP_pre, fileName= "microCLIP_pred")


