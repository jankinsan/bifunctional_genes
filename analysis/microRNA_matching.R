#read RefSeq data to get the exon coordinates!
hg38 <- read.delim("hg38.000001405.40-RS_2023_03.ncbiRefSeq.gtf", header = FALSE)
#separate transcripts and exons to different tables; exons could be later used for the miRNA matching!
hg38_exons <- hg38[which(hg38[,3]=="exon"),]

#processing for mapping to miRNAs! 
colnames(hg38_exons) <- c("chr", "date", "type", "start", "end", "misc", "strand", "misc2", "details")
hg38_exons_break<-sapply(hg38_exons$details, function(x){print(strsplit(x, split=";"))})
hg38_exons$gene <- sapply(1:length(hg38_exons_break), function(x){
  strsplit(hg38_exons_break[[x]][1], split="gene_id ")[[1]][2]
})
hg38_exons$acc <- sapply(1:length(hg38_exons_break), function(x){
  strsplit(hg38_exons_break[[x]][2], split="transcript_id ")[[1]][2]
})
hg38_counts <- read.csv("/lustre/sonam.dhamija/janki_nr/23122023_GRCh38.p14_bifunc_gene_counts_status.csv")
gene_idx <- na.omit(unlist(sapply(hg38_counts$X, function(x){
  which(hg38_exons$gene==x)
})))
hg38_exons_bifunc <- hg38_exons[gene_idx,]

#read microCLIP experimental dataset
#microRNA predictions and experimentally validated
#load
microCLIP_exp <- read.delim("./hg38_liftOver/16122023_microCLIPExperimental_liftover_hg38.bed", header=FALSE)
colnames(microCLIP_exp)<- c("chr", "start", "end", "microRNA", "score", "strand")
mir_exp_acc <- data.frame()


#some predicted sequences have uncommon labels for chromosomes at places, I am not usin them for the time being!
#subsetting datasets chromosome-wise, makes the finding more efficient!
for (chromosome in unique(microCLIP_exp$chr)){
  print(paste0("chr: ", chromosome))
  #positive strand
  chr_acc_pos <- hg38_exons_bifunc[intersect(which(hg38_exons_bifunc$chr==chromosome), which(hg38_exons_bifunc$strand=="+")),]
  chr_mir_pos <- microCLIP_exp[intersect(which(microCLIP_exp$chr==chromosome), which(microCLIP_exp$strand=="+")),]
  
  #negative strand
  chr_acc_neg <- hg38_exons_bifunc[intersect(which(hg38_exons_bifunc$chr==chromosome), which(hg38_exons_bifunc$strand=="-")),]
  chr_mir_neg <- microCLIP_exp[intersect(which(microCLIP_exp$chr==chromosome), which(microCLIP_exp$strand=="-")),]
  
  #select miRNAs matching each NM/NR/XM/XR
  for(mir_idx in 1:dim(chr_mir_pos)[1]){
    mir_idx_match <- intersect(which(chr_acc_pos$start<chr_mir_pos$start[mir_idx]),
                               which(chr_acc_pos$end>chr_mir_pos$end[mir_idx]))
    if(length(mir_idx_match)>0){
      for (x in mir_idx_match){
        mir_exp_acc <- rbind.data.frame(mir_exp_acc, data.frame(microRNA=chr_mir_pos$microRNA[mir_idx],
                                                                chr=chr_mir_pos$chr[mir_idx],
                                                                strand=chr_mir_pos$strand[mir_idx],
                                                                microRNA_start=chr_mir_pos$start[mir_idx],
                                                                microRNA_end=chr_mir_pos$end[mir_idx], 
                                                                accession=chr_acc_pos$acc[x],
                                                                gene=chr_acc_pos$gene[x],
                                                                transcript_start=chr_acc_pos$start[x],
                                                                transcript_end=chr_acc_pos$end[x]))
      }
    }
  }
  for(mir_idx in 1:dim(chr_mir_neg)[1]){
    mir_idx_match <- intersect(which(chr_acc_neg$start<chr_mir_neg$start[mir_idx]),
                               which(chr_acc_neg$end>chr_mir_neg$end[mir_idx]))
    if(length(mir_idx_match)>0){
      for (x in mir_idx_match){
        mir_exp_acc <- rbind.data.frame(mir_exp_acc, data.frame(microRNA=chr_mir_neg$microRNA[mir_idx],
                                                                chr=chr_mir_neg$chr[mir_idx],
                                                                strand=chr_mir_neg$strand[mir_idx],
                                                                microRNA_start=chr_mir_neg$start[mir_idx],
                                                                microRNA_end=chr_mir_neg$end[mir_idx],
                                                                accession=chr_acc_neg$acc[x],
                                                                gene=chr_acc_neg$gene[x],
                                                                transcript_start=chr_acc_pos$start[x],
                                                                transcript_end=chr_acc_pos$end[x]))
      }
    }
  }
}
write.table(mir_exp_acc, file="01012024_microCLIP_exp_matched_bifunc.csv", row.names=FALSE, 
            sep=",")
#microCLIP Predicted
microCLIP_pre <- read.delim("./hg38_liftOver/16122023_microCLIPPredicted_liftover_hg38.bed", header=FALSE)
colnames(microCLIP_pre)<- c("chr", "start", "end", "microRNA", "score", "strand")
mir_pre_acc <- data.frame()
#predicted
for (chromosome in unique(microCLIP_pre$chr)){
  print(paste0("chr: ", chromosome))
  #positive strand
  chr_acc_pos <- hg38_exons_bifunc[intersect(which(hg38_exons_bifunc$chr==chromosome), which(hg38_exons_bifunc$strand=="+")),]
  chr_mir_pos <- microCLIP_pre[intersect(which(microCLIP_pre$chr==chromosome), which(microCLIP_pre$strand=="+")),]
  
  #negative strand
  chr_acc_neg <- hg38_exons_bifunc[intersect(which(hg38_exons_bifunc$chr==chromosome), which(hg38_exons_bifunc$strand=="-")),]
  chr_mir_neg <- microCLIP_pre[intersect(which(microCLIP_pre$chr==chromosome), which(microCLIP_pre$strand=="-")),]
  
  #select miRNAs matching each NM/NR/XM/XR
  for(mir_idx in 1:dim(chr_mir_pos)[1]){
    mir_idx_match <- intersect(which(chr_acc_pos$start<chr_mir_pos$start[mir_idx]),
                               which(chr_acc_pos$end>chr_mir_pos$end[mir_idx]))
    if(length(mir_idx_match)>0){
      for (x in mir_idx_match){
        mir_pre_acc <- rbind.data.frame(mir_pre_acc, data.frame(microRNA=chr_mir_pos$microRNA[mir_idx],
                                                                chr=chr_mir_pos$chr[mir_idx],
                                                                strand=chr_mir_pos$strand[mir_idx],
                                                                microRNA_start=chr_mir_pos$start[mir_idx],
                                                                microRNA_end=chr_mir_pos$end[mir_idx],
                                                                score= chr_mir_pos$score[mir_idx],
                                                                accession=chr_acc_pos$acc[x],
                                                                gene=chr_acc_pos$gene[x],
                                                                transcript_start=chr_acc_pos$start[x],
                                                                transcript_end=chr_acc_pos$end[x]))
      }
    }
  }
  for(mir_idx in 1:dim(chr_mir_neg)[1]){
    mir_idx_match <- intersect(which(chr_acc_neg$start<chr_mir_neg$start[mir_idx]),
                               which(chr_acc_neg$end>chr_mir_neg$end[mir_idx]))
    if(length(mir_idx_match)>0){
      for (x in mir_idx_match){
        mir_pre_acc <- rbind.data.frame(mir_pre_acc, data.frame(microRNA=chr_mir_neg$microRNA[mir_idx],
                                                                chr=chr_mir_neg$chr[mir_idx],
                                                                strand=chr_mir_neg$strand[mir_idx],
                                                                microRNA_start=chr_mir_neg$start[mir_idx],
                                                                microRNA_end=chr_mir_neg$end[mir_idx],
                                                                score= chr_mir_neg$score[mir_idx],
                                                                accession=chr_acc_neg$acc[x],
                                                                gene=chr_acc_neg$gene[x],
                                                                transcript_start=chr_acc_pos$start[x],
                                                                transcript_end=chr_acc_pos$end[x]))
      }
    }
  }
}
write.table(mir_pre_acc, file="01012024_microCLIP_pre_matched_bifunc.csv", row.names=FALSE, 
            sep=",")
#TargetScan Conserved Site, Broadly conserved families
TS_broad_cons <- read.delim("./hg38_liftOver/16122023_Targets_CS_pctiles_hg38_broadConsFam_consSite.bed", header=FALSE)
colnames(TS_broad_cons)<- c("chr", "start", "end", "name", "score", "strand", 
                            "thickStart", "thickEnd", "itemRgb", "blockCount", 
                            "blockSizes", "blockStarts")
TS_broad_non <- read.delim("./hg38_liftOver/16122023_Targets_CS_pctiles_hg38_broadConsFam_nonConsSite.bed", header=FALSE)
colnames(TS_broad_non)<- c("chr", "start", "end", "name", "score", "strand", 
                            "thickStart", "thickEnd", "itemRgb", "blockCount", 
                            "blockSizes", "blockStarts")
#matching with genome tracks!
mirMatch <- function(mir_dat, fileName){
mir_dat_matched <- data.frame()
for (chromosome in unique(mir_dat$chr)){
  print(paste0("Running for", fileName))
  print(paste0("chr: ", chromosome))
  #positive strand
  chr_acc_pos <- hg38_exons_bifunc[intersect(which(hg38_exons_bifunc$chr==chromosome), which(hg38_exons_bifunc$strand=="+")),]
  chr_mir_pos <- mir_dat[intersect(which(mir_dat$chr==chromosome), which(mir_dat$strand=="+")),]
  
  #negative strand
  chr_acc_neg <- hg38_exons_bifunc[intersect(which(hg38_exons_bifunc$chr==chromosome), which(hg38_exons_bifunc$strand=="-")),]
  chr_mir_neg <- mir_dat[intersect(which(mir_dat$chr==chromosome), which(mir_dat$strand=="-")),]
  
  #select miRNAs matching each NM/NR/XM/XR
  for(mir_idx in 1:dim(chr_mir_pos)[1]){
    mir_idx_match <- intersect(which(chr_acc_pos$start<chr_mir_pos$start[mir_idx]),
                               which(chr_acc_pos$end>chr_mir_pos$end[mir_idx]))
    if(length(mir_idx_match)>0){
      for (x in mir_idx_match){
        mir_dat_matched <- rbind.data.frame(mir_dat_matched, data.frame(microRNA=chr_mir_pos$name[mir_idx],
                                                                        chr=chr_mir_pos$chr[mir_idx],
                                                                        strand=chr_mir_pos$strand[mir_idx],
                                                                        microRNA_start=chr_mir_pos$start[mir_idx],
                                                                        microRNA_end=chr_mir_pos$end[mir_idx],
                                                                        accession=chr_acc_pos$acc[x],
                                                                        gene=chr_acc_pos$gene[x],
                                                                        transcript_start=chr_acc_pos$start[x],
                                                                        transcript_end=chr_acc_pos$end[x],
                                                                        TargetScan_score=chr_mir_pos$score[mir_idx]))
      }
    }
  }
  for(mir_idx in 1:dim(chr_mir_neg)[1]){
    mir_idx_match <- intersect(which(chr_acc_neg$start<chr_mir_neg$start[mir_idx]),
                               which(chr_acc_neg$end>chr_mir_neg$end[mir_idx]))
    if(length(mir_idx_match)>0){
      for (x in mir_idx_match){
        mir_dat_matched <- rbind.data.frame(mir_dat_matched, data.frame(microRNA=chr_mir_neg$name[mir_idx],
                                                                                    chr=chr_mir_neg$chr[mir_idx],
                                                                                    strand=chr_mir_neg$strand[mir_idx],
                                                                                    microRNA_start=chr_mir_neg$start[mir_idx],
                                                                                    microRNA_end=chr_mir_neg$end[mir_idx],
                                                                                    accession=chr_acc_neg$acc[x],
                                                                                    gene=chr_acc_neg$gene[x],
                                                                                    transcript_start=chr_acc_pos$start[x],
                                                                                    transcript_end=chr_acc_pos$end[x],
                                                                                    TargetScan_score = chr_mir_neg$score[mir_idx]))
      }
    }
  }
}
write.table(mir_dat_matched, file=paste0(Sys.Date(), fileName, "_matched_bifunc.csv"), row.names=FALSE, 
              sep=",")
return(mir_dat_matched)
}

TS_broad_cons_matched <- mirMatch(TS_broad_cons, fileName = "TS_broad_cons")
TS_broad_non_matched <- mirMatch(TS_broad_non, fileName= "TS_broad_non")

