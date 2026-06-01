#hpc!
#module load compiler/R/4.1.2/intel2020
#R

library(parallel)
setwd("/scratch/bioschool/phd/blz227562/nm_nr")
#reading microCLIP predictions matched to hybrid genes
mir_pre_acc <- read.csv("2024-03-22_microCLIP_pred_matched_all_genes.csv")
#Target Scan default predictions
TS_default_matched <- read.csv("2024-03-22_TS_default_pred_matched_all_genes.csv")
TS_default_matched$microRNA<- unlist(sapply(TS_default_matched$microRNA, function(x){strsplit(x, split=":")[[1]][2]}))
#melt to add all family members; first taking out only human microRNAs!
mirna_info <- read.delim("miR_Family_Info.txt")
hsa_mirna_idx <- unlist(sapply(mirna_info$MiRBase.ID, function(x){
  length(unlist(strsplit(x, split="hsa-")))
}))
hsa_mirna_info<-mirna_info[which(hsa_mirna_idx>1),]
hsa_mirna_info$miR.family <- unlist(sapply(hsa_mirna_info$miR.family, 
                                           function(x){strsplit(x, split=' ')[[1]][1]}))

TS_melt_list <- mclapply(unique(TS_default_matched$microRNA), function(mir){
  print(paste0("miR Family: ", mir))
  #getting positions!
  mir_idx <- which(hsa_mirna_info$miR.family==mir)
  TS_idx <- which(TS_default_matched$microRNA==mir)
  
  #check if the mir Family in the dataframe
  if(length(TS_idx)==0){
    return(NA)
  } else if(length(TS_idx)>0){
    TS_melt_dat <- data.frame()
    TS_melt_dat <-TS_default_matched[rep(TS_idx,length(mir_idx)),]
    TS_melt_dat[, 1] <- unlist(lapply(hsa_mirna_info$MiRBase.ID[mir_idx], function(x){
    rep(x,length(TS_idx))}))
    return(TS_melt_dat)
  }
}, mc.cores=15)
TS_default_melt <- na.omit(do.call(rbind.data.frame, TS_melt_list))
write.csv(TS_default_melt, file=paste0(Sys.Date(), "_TS_default_pred_all_mirna_members.csv"))

#Add the predicted kd and other scores for each match!
#read the score file!
TS_scores <- read.delim("Predicted_Targets_Context_Scores.default_predictions.txt")
hsa_TS_scores <- unlist(sapply(TS_scores$miRNA, function(x){
  length(unlist(strsplit(x, split="hsa-")))
}))
hsa_TS_scores<-TS_scores[which(hsa_TS_scores>1),]
rm(TS_scores)
#combine gene and mirna to make unique strings to combine all scores from TargetScan in one place!
hsa_TS_scores$mir_gene = paste0(hsa_TS_scores$miRNA, "_", hsa_TS_scores$Gene.Symbol)
TS_default_melt$mir_gene = paste0(TS_default_melt$microRNA, "_", TS_default_melt$gene)
TS_added_scores <- data.frame()
TS_to_add <- data.frame()
idx_loop=0
for (str_mir in unique(hsa_TS_scores$mir_gene)){
  if(idx_loop%%1000==0){
    print(idx_loop)
  }
  gen_mir_idx <- which(hsa_TS_scores$mir_gene == str_mir)
  melt_idx <- which(TS_default_melt$mir_gene == str_mir)
  if (length(melt_idx)==0){
    #print(paste0("Gene-miRNA pair: ", str_mir,  " is not in the target scan matched dataset.."))
  } else if(length(gen_mir_idx)==1){
    TS_to_add <- TS_default_melt[melt_idx, ]
    TS_to_add[,13:18] <- hsa_TS_scores[gen_mir_idx, c(6,9:13)]
  } else {
    #get the order of the miRNA target sites sorted..
    TS_spliced <- TS_default_melt[melt_idx,]
    info_spliced <- hsa_TS_scores[gen_mir_idx,] 
    
    #get the different starts; from both the dataframes that are to be added together!
    mir_start<- sort(unique(TS_spliced$microRNA_start))
    mir_start_utr <- sort(info_spliced$UTR_start)
    #check the number of binding sites are the same! 
    if(length(mir_start)==length(mir_start_utr)){
      TS_to_add<- TS_spliced
      #to add other score parameters along with the coordinates
      idxs<- sapply(1:length(mir_start_utr), function(x){
        idx_site<- which(TS_spliced$microRNA_start==mir_start[x])})
      for (i in 1:length(mir_start_utr)){
        TS_to_add[idxs[[i]],13:18] <- info_spliced[i, c(6,9:13)] 
      }
    } else if (length(mir_start_utr)>length (mir_start)){ 
      #get the score of the microRNA start site closest to the UTR start position in the scores dataframe!
      diff <- mir_start_utr - (TS_spliced$microRNA_start[1]- min(TS_spliced$transcript_start))
      mir_keep<- match(min(diff), diff)
      TS_to_add<- TS_spliced
      TS_to_add[idxs[[i]],13:18] <- info_spliced[mir_keep, c(6,9:13)] 
    } else{
      print(paste0("You have a case still in question: ", str_mir))
    }
  }
  #combining it all in a data frame
  TS_added_scores <- rbind.data.frame(TS_added_scores, TS_to_add)
  TS_to_add <- data.frame()
  idx_loop = idx_loop+1
}
write.csv(TS_added_scores, file=paste0(Sys.Date(), "_TS_default_melt_scores_complete.csv"))

#27-08-2024 Rerunning from here!
#read the file!
library(parallel)
TS_added_scores<- read.csv("2024-03-28_TS_default_melt_scores_complete.csv")
TS_added_scores<- TS_added_scores[,-1]
#match with microCLIP predictions!
mir_pre_acc$mir_gene_acc = paste0(mir_pre_acc$microRNA, "__", mir_pre_acc$gene, "__", mir_pre_acc$accession)
TS_added_scores$mir_gene_acc = paste0(TS_added_scores$microRNA, "__", TS_added_scores$gene, "__", TS_added_scores$accession)
#take all possible gene-miRNA-transcript combinations and calculate the number of binding sites!
all_combinations <- unique(c(mir_pre_acc$mir_gene_acc, TS_added_scores$mir_gene_acc))
write.table(all_combinations, file=paste0(Sys.Date(), "_all_combinations_gene-mir-transcript.txt"))
common_targets <- mclapply(1:length(all_combinations), function(i){
  if(i%%1000==0){print(i)}
  x= all_combinations[i]
  return(c(x, length(which(mir_pre_acc$mir_gene_acc==x)), 
           length(which(TS_added_scores$mir_gene_acc==x))))
}, mc.cores=15)
#convert list to data frame
common_dat <- as.data.frame(do.call(rbind, common_targets))

common_dat[,4:6 ]<- t(sapply(common_dat[,1], function(x){
  c(strsplit(x, split="__")[[1]][1], strsplit(x, split="__")[[1]][2], 
    strsplit(x, split="__")[[1]][3])
}))
colnames(common_dat)<- c("mir_gene_acc", "microCLIP_sites", "TargetScan_Sites", "microRNA", "gene", "accession")
write.csv(common_dat, file=paste0(Sys.Date(), "_all_targets_TS_microCLIP_combined.csv"))

common_sites <- common_dat[intersect(which(common_dat$microCLIP_sites>0),which(common_dat$TargetScan_Sites>0)),]
write.csv(common_sites, file=paste0(Sys.Date(), "_common_sites_TS_microCLIP_combined.csv"))

#my laptop
#microRNA downstream analysis!
setwd("E:/MSR BLY/THESIS WORK/NP NR Project/15122023")
#read dataset from humans!
gene_status <- read.csv("./counts/18122023_GRCh38.p14_genewise_counts_status.csv")

mir_combined_all<-read.csv("./microRNA/2024-08-27_all_targets_TS_microCLIP_combined.csv")
mir_combined_all$gene_status <- gene_status$status[match(mir_combined_all$gene, gene_status$X)]

#summarize microRNA statistics!
mir_summary <- data.frame(details=c("sites", "microRNAs", "genes", "transcripts"))
mir_microCLIP <- which(mir_combined_all$microCLIP_sites>=1)
mir_TargetScan<- which(mir_combined_all$TargetScan_Sites>=1)

mir_summary$microCLIP <- c(length(mir_microCLIP),
                           length(unique(mir_combined_all$microRNA[mir_microCLIP])),
                           length(unique(mir_combined_all$gene[mir_microCLIP])),
                           length(unique(mir_combined_all$accession[mir_microCLIP])))

mir_summary$TargetScan <- c(length(mir_TargetScan),
                            length(unique(mir_combined_all$microRNA[mir_TargetScan])),
                            length(unique(mir_combined_all$gene[mir_TargetScan])),
                            length(unique(mir_combined_all$accession[mir_TargetScan])))
#sites in both predictions
mir_common <- intersect(mir_TargetScan, mir_microCLIP)
mir_summary$common<- c(length(intersect(mir_TargetScan, mir_microCLIP)),
                       length(unique(mir_combined_all$microRNA[intersect(mir_TargetScan, mir_microCLIP)])),
                       length(unique(mir_combined_all$gene[intersect(mir_TargetScan, mir_microCLIP)])),
                       length(unique(mir_combined_all$accession[intersect(mir_TargetScan, mir_microCLIP)]))) 
#all sites which have a
mir_summary$combined <- c(length(mir_combined_all[,1]),
                       length(unique(mir_combined_all$microRNA)),
                       length(unique(mir_combined_all$gene)),
                       length(unique(mir_combined_all$accession)))

for (gene_type in c('hybrid', "only_ncRNA", "only_mRNA")){
  print(paste0("Running for ", gene_type, "..."))
  mir_select_type <- which(mir_combined_all$gene_status==gene_type)
  select_dat_mc <- c(length(intersect(mir_microCLIP, mir_select_type)),
                     length(unique(mir_combined_all$microRNA[intersect(mir_microCLIP, mir_select_type)])),
                     length(unique(mir_combined_all$gene[intersect(mir_microCLIP, mir_select_type)])),
                     length(unique(mir_combined_all$accession[intersect(mir_microCLIP, mir_select_type)])))
  select_dat_ts <- c(length(intersect(mir_TargetScan, mir_select_type)),
                     length(unique(mir_combined_all$microRNA[intersect(mir_TargetScan, mir_select_type)])),
                     length(unique(mir_combined_all$gene[intersect(mir_TargetScan, mir_select_type)])),
                     length(unique(mir_combined_all$accession[intersect(mir_TargetScan, mir_select_type)])))
  select_dat <- c(length(intersect(mir_common, mir_select_type)),
                  length(unique(mir_combined_all$microRNA[intersect(mir_common, mir_select_type)])),
                  length(unique(mir_combined_all$gene[intersect(mir_common, mir_select_type)])),
                  length(unique(mir_combined_all$accession[intersect(mir_common, mir_select_type)])))
  select_dat_all <- c(length(mir_select_type),
                  length(unique(mir_combined_all$microRNA[mir_select_type])),
                  length(unique(mir_combined_all$gene[mir_select_type])),
                  length(unique(mir_combined_all$accession[mir_select_type])))
  
  colnames_bind <- sapply(c("sites", "microRNAs", "genes", "transcripts"), function(x){paste0(gene_type, "_", x)})
  all_bind_dat<- cbind.data.frame(colnames_bind, select_dat_mc, select_dat_ts, select_dat, select_dat_all)
  colnames(all_bind_dat)<- colnames(mir_summary)
  mir_summary <- rbind.data.frame(mir_summary, all_bind_dat)
}
row.names(mir_summary)<- 1:16
write.csv(mir_summary, file=paste0("./microRNA/", Sys.Date(), "_microRNA_summary.csv"), row.names = FALSE)


