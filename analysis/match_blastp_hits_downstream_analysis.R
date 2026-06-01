match_blastp_analysis <- function (ORF_file, blastp_file, outFile, wd){
  print("Setting directory and reading files..")
  setwd(wd)
  seq_orf = read.delim(file=ORF_file, sep="\t")
  blastp_bifunc<- read.delim(blastp_file, header=FALSE, sep=",")
  colnames(blastp_bifunc)<- unlist(strsplit("qseqid sseqid stitle pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore qcovs qcovhsp qcovus", 
                                            split=" "))
  #get the numbers!
  print(length(unique(blastp_bifunc$qseqid)))
  matORF_idx <- match(unique(blastp_bifunc$qseqid), seq_orf$acc)
  blastedORF_genes <- unlist(sapply(seq_orf$genes[matORF_idx], function(x){strsplit(x, split=",")}))
  length(unique(blastedORF_genes))
  
  matID <- match(blastp_bifunc$qseqid, seq_orf$acc)
  pred_blast_comb<- cbind.data.frame(qseqid = blastp_bifunc[,1] , qgenes = seq_orf$genes[matID],
                                     qaccs = seq_orf$accs[matID], qseq = seq_orf$x[matID],
                                     qstart_codon= seq_orf$start[matID],
                                     qstop_codon= seq_orf$stop[matID], blastp_bifunc[,-1])
  pred_blast_comb$C_term_status<- "FALSE"
  C_term_idx <- which(pred_blast_comb$send==pred_blast_comb$slen)
  pred_blast_comb$C_term_status[C_term_idx] <- "TRUE"
  length(unique(pred_blast_comb$qseq))
  
  write.table(pred_blast_comb, outFile, sep="\t", row.names=FALSE,
              quote=FALSE)
  print("Saving the comnined file...")
  
  
  #C-terminus hitting queries count
  print("Blasting ORFs.. ")
  print(length(unique(pred_blast_comb$qseqid[which(pred_blast_comb$C_term_status=="FALSE")])))
  print("Non-blasiting ORFs..")
  print(length(unique(pred_blast_comb$qseqid[-which(pred_blast_comb$C_term_status=="FALSE")])))
  
  
  #Are any ORFs encoded by more than one gene?
  gene_number<- data.frame("number_genes"=c(1:12), "number_ORFs" = rep(0,12), "number_smORFs"=rep(0,12),
                           "number_non_blastedORFs" = rep(0,12), "number_Cterm_blastedORFs"=rep(0,12), 
                           "number_nonCterm_blastedORFs" = rep(0,12), 
                           "number_NMD_transcript_ORFs"=rep(0,12), 
                           "number_nonNMD_transcript_ORFs"=rep(0,12))
  
  
}

match_blastp(ORF_file = "unique_ORFs/2024-05-05_uniqueORFs_bifunc_ATG.tsv",
             blastp_file = "./blastp/05052024_csv/05052024_blastp-fast_hybridGenes_uniqueORFs_ATG.tsv", 
             outFile="unique_ORFs/2024-05-21_blastpHybridGenes_ATG_uniqueORFsMatched.tsv",
             wd= "F:/Janki/bifunctional_genes/data/ORFs")
match_blastp(ORF_file = "unique_ORFs/2024-05-05_uniqueORFs_bifunc_CTG.tsv",
             blastp_file = "./blastp/05052024_csv/05052024_blastp-fast_hybridGenes_uniqueORFs_CTG.tsv", 
             outFile="unique_ORFs/2024-05-21_blastpHybridGenes_CTG_uniqueORFsMatched.tsv",
             wd= "F:/Janki/bifunctional_genes/data/ORFs")
match_blastp(ORF_file = "unique_ORFs/2024-05-05_uniqueORFs_bifunc_GTG.tsv",
             blastp_file = "./blastp/05052024_csv/05052024_blastp-fast_hybridGenes_uniqueORFs_GTG.tsv",
             outFile="unique_ORFs/2024-05-21_blastpHybridGenes_GTG_uniqueORFsMatched.tsv",
             wd= "F:/Janki/bifunctional_genes/data/ORFs")
match_blastp(ORF_file = "unique_ORFs/2024-05-05_uniqueORFs_bifunc_TTG.tsv",
             blastp_file = "./blastp/05052024_csv/05052024_blastp-fast_hybridGenes_uniqueORFs_TTG.tsv", 
             outFile="unique_ORFs/2024-05-21_blastpHybridGenes_TTG_uniqueORFsMatched.tsv",
             wd= "F:/Janki/bifunctional_genes/data/ORFs")


#add nmd transcript info to the unique ORFs!
#matching predicted ORFs to NMD transcripts
nmd_info=read.csv("F:/Janki/NM NR Project/15122023/nmd_transcripts/13052024_bifunc_gene_accessions_nmd.csv")
check_nmd_transcripts<- function(wd, ORF_file, blastp_file, out_string, out_string_blastp){
  
  print("Setting directory and reading files..")
  setwd(wd)
  seq_orf = read.delim(file=ORF_file)
  print(paste0("matching predicted ORFs to NMD transcripts for ", out_string))
  seq_orf$nmd_transcripts = 0
  seq_orf$nmd_status <- "nonNMD"
  for (idx in 1:dim(seq_orf)[1]){
    if (idx%%1000==0){
      print(idx)
    }
    #if the ORF is encoded by only one transcript
    if(seq_orf$num_transcripts[idx]==1){
      mat <- match(seq_orf$accs[idx], nmd_info$accession)
      if(nmd_info$status[mat]==""){
        seq_orf$nmd_transcripts[idx]=0
        if(isTRUE(grep("XR_", seq_orf$accs[idx]))){
          seq_orf$nmd_status[idx] = "only_XR"
        }
      } else if(nmd_info$status[mat]=="nmd_transcript"){
        seq_orf$nmd_transcripts[idx]=1
        seq_orf$nmd_status[idx]="NMD"
      } else if(nmd_info$status[mat]=="changed_to_NM"){
        seq_orf$nmd_transcripts[idx]=0 
        seq_orf$nmd_status[idx]="changed_to_NM"
      }
    } else{
      nmd_status <- sapply(strsplit(seq_orf$accs[idx], split=",")[[1]], function(x){
        mat <- match(x, nmd_info$accession)
        if(nmd_info$status[mat]==""){
          if(length(grep("XR_",x))>0){
            return (c(0, "only_XR"))
          } else return(c(0, "nonNMD"))
        }else if(nmd_info$status[mat]=="nmd_transcript"){
          return(c(1, "NMD"))
        } else if(nmd_info$status[mat]=="changed_to_NM"){
          return(c(0, "changed_to_NM"))
        }
        })
      seq_orf$nmd_transcripts[idx] = sum(as.numeric(nmd_status[1,]))
      
      #check the number of nmd status and match!
      nonnmd_sum <- length(which(nmd_status[2,]=="nonNMD"))
      xr_sum <- length(which(nmd_status[2,]=="only_XR"))
      changed_sum <- length(which(nmd_status[2,]=="changed_to_NM"))
      
      if(seq_orf$num_transcripts[idx]== seq_orf$nmd_transcripts[idx]){
        seq_orf$nmd_status[idx]="NMD"
      } else if(seq_orf$num_transcripts[idx]!= seq_orf$nmd_transcripts[idx]){
        #some NR non-NMD
        if(nonnmd_sum>0){
          seq_orf$nmd_status[idx]="nonNMD"
        }
        #only XR non-NMD
        if(seq_orf$num_transcripts[idx]- seq_orf$nmd_transcripts[idx]==xr_sum){
          seq_orf$nmd_status[idx]="onlyXR_nonNMD"
        }
        
        if(seq_orf$num_transcripts[idx]==changed_sum){
          seq_orf$nmd_status[idx]="changed_to_NM"
        }
      } 
    }
  }
  
 
  dir.create("uniqueORFs_nmdMatch")
  write.csv(seq_orf, file = paste0("./uniqueORFs_nmdMatch/", Sys.Date(), "_", out_string, "_", "nmd.csv"))
  
  #read blastp file 
  blastp_bifunc <- read.delim(blastp_file)
  
  seq_orf$blast_status<- "nonblast"
  seq_orf$blast_status[which(seq_orf$acc %in% blastp_bifunc$qseqid[blastp_bifunc$C_term_status==TRUE])] <- "blast_Cterm"
  seq_orf$blast_status[which(seq_orf$acc %in% blastp_bifunc$qseqid[blastp_bifunc$C_term_status!=TRUE])] <- "blast_nonCterm"
  write.csv(seq_orf, file = paste0("./uniqueORFs_nmdMatch/", Sys.Date(), "_", out_string, ".csv"))
  }

check_nmd_transcripts(ORF_file = "unique_ORFs/2024-05-05_uniqueORFs_bifunc_ATG.tsv",
                      wd= "F:/Janki/bifunctional_genes/data/ORFs",
                      blastp_file="unique_ORFs/2024-05-21_blastpHybridGenes_ATG_uniqueORFsMatched.tsv",
                      out_string = "bifunc_ATG_uniqueORFs", 
                      out_string_blastp="unique0RFs_bifunc_ATG_nmd_blast")
check_nmd_transcripts(ORF_file = "unique_ORFs/2024-05-05_uniqueORFs_bifunc_CTG.tsv",
                      wd= "F:/Janki/bifunctional_genes/data/ORFs",
                      blastp_file="unique_ORFs/2024-05-21_blastpHybridGenes_CTG_uniqueORFsMatched.tsv",
                      out_string = "bifunc_CTG_uniqueORFs", 
                      out_string_blastp="unique0RFs_bifunc_CTG_nmd_blast")
check_nmd_transcripts(ORF_file = "unique_ORFs/2024-05-05_uniqueORFs_bifunc_TTG.tsv",
                      wd= "F:/Janki/bifunctional_genes/data/ORFs",
                      blastp_file="unique_ORFs/2024-05-21_blastpHybridGenes_TTG_uniqueORFsMatched.tsv",
                      out_string = "bifunc_TTG_uniqueORFs", 
                      out_string_blastp="unique0RFs_bifunc_TTG_nmd_blast")
check_nmd_transcripts(ORF_file = "unique_ORFs/2024-05-05_uniqueORFs_bifunc_GTG.tsv",
                      wd= "F:/Janki/bifunctional_genes/data/ORFs",
                      blastp_file="unique_ORFs/2024-05-21_blastpHybridGenes_GTG_uniqueORFsMatched.tsv",
                      out_string = "bifunc_GTG_uniqueORFs", 
                      out_string_blastp="unique0RFs_bifunc_GTG_nmd_blast")

#match nmd status with blastp status! 

atg_nonNMD_orfs <- read.csv("F:/Janki/bifunctional_genes/data/ORFs/uniqueORFs_nmdMatch/2024-05-21_unique0RFs_bifunc_ATG_nmd_blast.csv")
kss_atg <- read.csv("F:/Janki/bifunctional_genes/data/kss/14062024_kss-len_formatted_pep_ATG_bifunc_genes.csv")
 

