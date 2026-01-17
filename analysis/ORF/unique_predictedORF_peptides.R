library(ggplot2)

setwd("E:/MSR BLY/THESIS WORK/NP NR Project/15122023/ORFs")
#read orf prediction formatted data!
ORF_DAT <- read.csv("predicted/23122023_formatted_allpep_ATG-TTG-GTG-CTC-bifunc_genes.csv")
#match gene status with our annotation: no need here since they are bifunc genes!
ORF_DAT$X <- paste0(ORF_DAT$acc_orf, "_", ORF_DAT$orf_number)

out_list <- lapply(unique(ORF_DAT$orf_sequence), function(sequ){
  x = which(ORF_DAT$orf_sequence==sequ)
  ORF_DAT[x,]
})
names(out_list)<- unique(ORF_DAT$orf_sequence)

#get this to an output file
dir.create("unique_ORFs")
lapply(out_list, write.table, "unique_ORFs/26122023_bifunc_ORFs_sorted_by_sequences.txt", append=T)

seq_orf_dat <- apply(rbind(names(out_list)), 2, function(x){
  num_transcripts <- dim(out_list[[x]])[1]
  accs <- unlist(strsplit(c(out_list[[x]]$acc_orf), split=">"))[c(FALSE, TRUE)]
  genes <- c(unique(out_list[[x]]$gene))
  return(cbind(x, num_transcripts, accs=paste(accs, collapse = ","), genes=paste(genes, collapse = ","), out_list[[x]][1,c(7,9,10)]))
})
seq_orf <- data.frame()

#getting the numbers!
length(seq_orf$x)
length(unique(seq_orf$x))

for (i in 1:length(seq_orf_dat)){
  seq_orf <- rbind.data.frame(seq_orf, seq_orf_dat[[i]])
}
length(which(seq_orf$start=="start:ATG"))
length(unique(unlist(sapply(seq_orf$gene[which(seq_orf$start=="start:ATG")], function(x){strsplit(x, split=",")}))))

#separate the microproteins (50aa-120aa)
ORF_DAT$len_strip <- as.numeric(sapply(1:dim(ORF_DAT)[1], function(i){
  strsplit(ORF_DAT$length[i], split="length:")[[1]][2]}))
seq_orf$len_strip <- as.numeric(sapply(1:dim(seq_orf)[1], function(i){
  strsplit(seq_orf$length[i], split="length:")[[1]][2]}))

microprot_dat <- seq_orf[intersect(which(seq_orf$len_strip <=300), which(seq_orf$len_strip>=150)),]
length(unique(microprot_dat$genes))
microprot_dat2 <- ORF_DAT[intersect(which(ORF_DAT$len_strip <=300), which(ORF_DAT$len_strip>=150)),]

#unique number of genes from microproteins
length(unique(unlist(sapply(microprot_dat$gene, function(x){strsplit(x, split=",")}))))

#check the start:ATG from microproteins
length(which(microprot_dat$start=="start:ATG"))
length(unique(unlist(sapply(microprot_dat$gene[which(microprot_dat$start=="start:ATG")], function(x){strsplit(x, split=",")}))))
#write the files above for later use!
write.csv(seq_orf, file="unique_ORFs/26122023_ORF_sequences_unique_bifunc_allSTART.csv")
write.csv(microprot_dat, file="unique_ORFs/26122023_50-120aa_ORF_sequences_unique_bifunc_allSTART.csv")
write.csv(microprot_dat2, file="unique_ORFs/26122023_50-120aa_formatted_allpep_ATG-TTG-GTG-CTG-bifunc_genes.csv")

#giving a unique accession to all predicted ORFs
seq_orf$acc <- sapply(1:dim(seq_orf)[1], function(x){paste0("predictedORF_", x)})
#writing to a fasta file for blastp
sapply(1:dim(seq_orf)[1], function(x){
  write(paste(paste0(">", seq_orf[x,9]), paste(seq_orf[x,c(5, 6, 7)], collapse = " "), 
              paste0("genes:", seq_orf$genes[x]),
              paste0("accs:", seq_orf$accs[x]),
              paste0("num_transcripts:", seq_orf[x,2])), 
        file="26122023_fasta_uniqueORFsbySeq.txt", append=TRUE)
  write(paste0(seq_orf[x,c(1)]), file="26122023_fasta_uniqueORFsbySeq.txt", append=TRUE)
})

#10-01-2024
#Reading the files back!
seq_orf<- read.csv(file="unique_ORFs/26122023_ORF_sequences_unique_bifunc_allSTART.csv")
seq_orf$acc <- sapply(1:dim(seq_orf)[1], function(x){paste0("predictedORF_", x)})
microprot_dat <- seq_orf[intersect(which(seq_orf$len_strip <=300), which(seq_orf$len_strip>=150)),]
microprot_dat2<- read.csv(file="unique_ORFs/26122023_50-120aa_formatted_allpep_ATG-TTG-GTG-CTG-bifunc_genes.csv")
length(unique(microprot_dat$genes))
length(unique(microprot_dat$X))

#check the number of ORFs mapping to each category of transcript number
#check the number of ORFs mapping to each category of transcript number

#Number of transcripts encoding the same ORF!
num_unique_ORFs <- data.frame(transcripts_encoding_ORF=1:max(microprot_dat$num_transcripts))
num_unique_ORFs$num_ORFs<- sapply(num_unique_ORFs[,1],function(x){
  length(which(seq_orf$num_transcripts==x))
})
num_unique_ORFs$num_smORFs<- sapply(num_unique_ORFs[,1],function(x){
  length(which(microprot_dat$num_transcripts==x))
})
write.csv(num_unique_ORFs, file="unique_ORFs/11012024_unique_ORFs_numbers.csv", row.names = FALSE)

size_unique_ORFs <- cbind("ORF_size" = c("50-100", 
                                              "101-500", 
                                              "501-1000", 
                                              "1001-2000", 
                                              ">2000"), 
                          number = c(length(which(seq_orf$len_strip<301)),
                                     length(intersect(which(seq_orf$len_strip>300),which(seq_orf$len_strip<1501))),
                                     length(intersect(which(seq_orf$len_strip>1500),which(seq_orf$len_strip<3000))),
                                     length(intersect(which(seq_orf$len_strip>3000),which(seq_orf$len_strip<6001))),
                                     length(which(seq_orf$len_strip>6000))))
write.csv(num_unique_ORFs, file="unique_ORFs/11012024_unique_ORFs_numbers.csv", row.names = FALSE)

#all start codons
orf_nums <- data.frame("number_ORF"=length(unique(seq_orf$X)), 
                       "number_genes"=length(unique(seq_orf$gene)), 
                       "number_transcripts"=length(unique(unlist(sapply(seq_orf$accs, function(x){
                         strsplit(x, split=",")
                       })))))
orf_nums["ATG",] <- c(length(unique(seq_orf$X[which(seq_orf$start=="start:ATG")])),
                      length(unique(seq_orf$gene[which(seq_orf$start=="start:ATG")])),
                      length(unique(unlist(sapply(seq_orf$accs[which(seq_orf$start=="start:ATG")],function(x){
                        strsplit(x, split=",")
                      })))))
orf_nums["GTG",] <- c(c(length(unique(seq_orf$X[which(seq_orf$start=="start:GTG")])),
                        length(unique(seq_orf$gene[which(seq_orf$start=="start:GTG")])),
                        length(unique(unlist(sapply(seq_orf$accs[which(seq_orf$start=="start:GTG")], function(x){
                          strsplit(x, split=",")
                        }))))))
orf_nums["TTG",] <- c(length(unique(seq_orf$X[which(seq_orf$start=="start:TTG")])),
                      length(unique(seq_orf$gene[which(seq_orf$start=="start:TTG")])),
                      length(unique(unlist(sapply(seq_orf$accs[which(seq_orf$start=="start:TTG")],function(x){
                        strsplit(x, split=",")
                      })))))
orf_nums["CTG",] <- c(length(unique(seq_orf$X[which(seq_orf$start=="start:CTG")])),
  length(unique(seq_orf$gene[which(seq_orf$start=="start:CTG")])),
  length(unique(unlist(sapply(seq_orf$accs[which(seq_orf$start=="start:CTG")], function(x){
    strsplit(x, split=",")
  })))))

#sort by stop codons
orf_nums$num_TAA <- c(length(unique(seq_orf$X[which(seq_orf$stop=="stop:TAA")])), 
                      length(unique(seq_orf$X[intersect(which(seq_orf$start=="start:ATG"),
                                                        which(seq_orf$stop=="stop:TAA"))])),
                      length(unique(seq_orf$X[intersect(which(seq_orf$start=="start:GTG"),
                                                        which(seq_orf$stop=="stop:TAA"))])),
                      length(unique(seq_orf$X[intersect(which(seq_orf$start=="start:TTG"),
                                                        which(seq_orf$stop=="stop:TAA"))])),
                      length(unique(seq_orf$X[intersect(which(seq_orf$start=="start:CTG"),
                                                        which(seq_orf$stop=="stop:TAA"))])))
orf_nums$num_TGA <- c(length(unique(seq_orf$X[which(seq_orf$stop=="stop:TGA")])), 
                      length(unique(seq_orf$X[intersect(which(seq_orf$start=="start:ATG"),
                                                        which(seq_orf$stop=="stop:TGA"))])),
                      length(unique(seq_orf$X[intersect(which(seq_orf$start=="start:GTG"),
                                                        which(seq_orf$stop=="stop:TGA"))])),
                      length(unique(seq_orf$X[intersect(which(seq_orf$start=="start:TTG"),
                                                        which(seq_orf$stop=="stop:TGA"))])),
                      length(unique(seq_orf$X[intersect(which(seq_orf$start=="start:CTG"),
                                                        which(seq_orf$stop=="stop:TGA"))])))
orf_nums$num_TAG <- c(length(unique(seq_orf$X[which(seq_orf$stop=="stop:TAG")])), 
                      length(unique(seq_orf$X[intersect(which(seq_orf$start=="start:ATG"),
                                                        which(seq_orf$stop=="stop:TAG"))])),
                      length(unique(seq_orf$X[intersect(which(seq_orf$start=="start:GTG"),
                                                        which(seq_orf$stop=="stop:TAG"))])),
                      length(unique(seq_orf$X[intersect(which(seq_orf$start=="start:TTG"),
                                                        which(seq_orf$stop=="stop:TAG"))])),
                      length(unique(seq_orf$X[intersect(which(seq_orf$start=="start:CTG"),
                                                        which(seq_orf$stop=="stop:TAG"))])))
#plot ORF Numbers
colnames(orf_nums)[1]<- "Total"

#04-01-2024
#Read the blastp results and continue analysis!
blastp_bifunc<- read.delim("./blastp/29122023_blastp-fast_hybridGenes_uniqueORFs.tsv", header=FALSE)
colnames(blastp_bifunc)<- unlist(strsplit("qseqid sseqid stitle pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore qcovs qcovhsp qcovus", 
                     split=" "))

#check the number of predicted ORFs that are blasting!
length(unique(blastp_bifunc$qseqid))
