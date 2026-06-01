#Match unique ORFs for ORFs from mRNAs and also the KSS
matchORF <- function(ORFfile, KSSfile, mRNAfile, outFile){
  #read files
  ORF_dat <- read.csv(ORFfile)
  KSS_dat <- read.csv(KSSfile)
  mRNA_ORFs_dat <- read.csv(mRNAfile)
  print("Files read!")
  
  #removing useless first two columns
  ORF_dat <- ORF_dat[,-c(1,2)]
  ORF_dat$exact_mRNA <- "No"
  #adding an ORF length column to match it!
  KSS_dat$ORF_length <- KSS_dat$end_position-KSS_dat$start_position
  #match the mRNA portion first!
  for (idx in 1:dim(ORF_dat)[1]){
    if(idx%%1000 == 0) {
      print(idx)
    } 
    mat <- match(ORF_dat$x[idx], mRNA_ORFs_dat$orf_sequence)
    if(!is.na(mat)){
      ORF_dat$exact_mRNA[idx] <- "Yes" 
    }
  }
  print("If exact matches found in mRNA, then matched!")
  
  #KSS matching next!
  combined_dat<- data.frame()
  #to keep track of how many times each accession is coming!
  acc_dat <- data.frame()
  for (idx in 1:dim(ORF_dat)[1]){
    if(idx%%1000 == 0) {
      print(idx)
     }
    if(ORF_dat$num_transcripts[idx]==1){
      acc <- paste0(ORF_dat$accs[idx], "_", ORF_dat$len_strip[idx])
      acc_dat <- rbind.data.frame(acc_dat, acc)
      KSS_dat_sliced <- KSS_dat[intersect(which(KSS_dat$gene==ORF_dat$genes[idx]), which(KSS_dat$sequence_id==ORF_dat$accs[idx])),]
      
      #match size and take the first one!
      all_same_len<- which(KSS_dat_sliced$ORF_length==ORF_dat$len_strip[idx])
      if(length(all_same_len)>1){
        to_take_idx <- all_same_len[length(grep(acc, acc_dat[,1]))]
      } else {
        to_take_idx<- all_same_len
      }
      combined_dat <- rbind.data.frame(combined_dat,
                                       c(ORF_dat[idx,], KSS_dat_sliced[to_take_idx, 6:10]))
     } else{
      #matching similarity score with the first transcript for the time-being!
      #need to modify this later!
      acc <- strsplit(ORF_dat$accs[idx], split=",")[[1]][1]
      acc2 <- paste0(acc, "_", ORF_dat$len_strip[idx])
      acc_dat <- rbind.data.frame(acc_dat, paste0(acc, "_", ORF_dat$len_strip[idx]))
      
      #taking the first gene too!
      gene <- strsplit(ORF_dat$genes[idx], split=",")[[1]][1]
      KSS_dat_sliced <- KSS_dat[intersect(which(KSS_dat$gene==gene), which(KSS_dat$sequence_id==acc)),]
      
      #match size and take the first one!
      all_same_len<- which(KSS_dat_sliced$ORF_length==ORF_dat$len_strip[idx])
        if(length(all_same_len)>1){
      to_take_idx <- all_same_len[length(grep(acc2, acc_dat[,1]))]
      } else {
        to_take_idx<- all_same_len
      }
      combined_dat <- rbind.data.frame(combined_dat,
                                       c(ORF_dat[idx,], KSS_dat_sliced[to_take_idx, 6:10]))
    } 
  }
  print("KSS matching finished, writing to file!")
  write.csv(combined_dat, file=outFile, row.names = FALSE)
}

matchORF(ORFfile="E:/MSR BLY/THESIS WORK/NP NR Project/15122023/ORFs/uniqueORFs_nmdMatch/2024-05-21_unique0RFs_bifunc_ATG_nmd_blast.csv",
KSSfile="E:/MSR BLY/THESIS WORK/NP NR Project/15122023/kss/14062024_kss-len_formatted_pep_ATG_bifunc_genes.csv",
mRNAfile="E:/MSR BLY/THESIS WORK/NP NR Project/15122023/coding_transcripts/orfs/15072024_formatted_bifuncGenes_codingPEP_ATGorfs.csv",
outFile="E:/MSR BLY/THESIS WORK/NP NR Project/15122023/ORFs/uniqueORFs_nmdMatch/2024-07-17_unique0RFs_bifunc_ATG_nmd_blast_mRNAs_KSS.csv")


print("Running for CTG")
matchORF(ORFfile="E:/MSR BLY/THESIS WORK/NP NR Project/15122023/ORFs/uniqueORFs_nmdMatch/2024-05-21_unique0RFs_bifunc_CTG_nmd_blast.csv",
         KSSfile="E:/MSR BLY/THESIS WORK/NP NR Project/15122023/kss/14062024_kss-len_formatted_pep_CTG_bifunc_genes.csv",
         mRNAfile="E:/MSR BLY/THESIS WORK/NP NR Project/15122023/coding_transcripts/orfs/15072024_formatted_bifuncGenes_codingPEP_CTGorfs.csv",
         outFile="E:/MSR BLY/THESIS WORK/NP NR Project/15122023/ORFs/uniqueORFs_nmdMatch/2024-07-17_unique0RFs_bifunc_CTG_nmd_blast_mRNAs_KSS.csv")

print("Running for TTG")
matchORF(ORFfile="E:/MSR BLY/THESIS WORK/NP NR Project/15122023/ORFs/uniqueORFs_nmdMatch/2024-05-21_unique0RFs_bifunc_TTG_nmd_blast.csv",
         KSSfile="E:/MSR BLY/THESIS WORK/NP NR Project/15122023/kss/14062024_kss-len_formatted_pep_TTG_bifunc_genes.csv",
         mRNAfile="E:/MSR BLY/THESIS WORK/NP NR Project/15122023/coding_transcripts/orfs/15072024_formatted_bifuncGenes_codingPEP_TTGorfs.csv",
         outFile="E:/MSR BLY/THESIS WORK/NP NR Project/15122023/ORFs/uniqueORFs_nmdMatch/2024-07-17_unique0RFs_bifunc_TTG_nmd_blast_mRNAs_KSS.csv")

print("Running for GTG")
matchORF(ORFfile="E:/MSR BLY/THESIS WORK/NP NR Project/15122023/ORFs/uniqueORFs_nmdMatch/2024-05-21_unique0RFs_bifunc_GTG_nmd_blast.csv",
         KSSfile="E:/MSR BLY/THESIS WORK/NP NR Project/15122023/kss/14062024_kss-len_formatted_pep_GTG_bifunc_genes.csv",
         mRNAfile="E:/MSR BLY/THESIS WORK/NP NR Project/15122023/coding_transcripts/orfs/15072024_formatted_bifuncGenes_codingPEP_GTGorfs.csv",
         outFile="E:/MSR BLY/THESIS WORK/NP NR Project/15122023/ORFs/uniqueORFs_nmdMatch/2024-07-17_unique0RFs_bifunc_GTG_nmd_blast_mRNAs_KSS.csv")
