## set command line arguments ----
args <- commandArgs(trailingOnly = TRUE)

#stop the script if no command line argument
if(length(args)==0){
  print("Please add arguments as needed!")
  stop("Requires command line argument.")
}

#write function
unique_ORFs_get <- function(start, inFile, outFileString){
  
  ORF_DAT <- read.csv(inFile)

  seq_orf_dat <- lapply(unique(ORF_DAT$orf_sequence), function(x){
    unique_orf_dat <- ORF_DAT[which(ORF_DAT$orf_sequence==x),]
    num_transcripts <- dim(unique_orf_dat)[1]
    accs <- unlist(strsplit(c(unique_orf_dat$acc_orf), split=">"))[c(FALSE, TRUE)]
    genes <- c(unique(unique_orf_dat$gene))
    return(cbind(x, num_transcripts, accs=paste(accs, collapse = ","), genes=paste(genes, collapse = ","), unique_orf_dat[1,c(7,9,10)]))
  })
  seq_orf <- do.call(rbind, seq_orf_dat)
  seq_orf <- cbind.data.frame(acc=sapply(1:dim(seq_orf)[1], function(x){paste0("predictedORF_",start, "_",x)}), seq_orf)
  #getting the numbers!
  print(paste0("Original number of ORFs with ", start, ": ", length(ORF_DAT$X)))
  print(paste0("Unique number of ORFs with ", start, ": ", length(unique(seq_orf$x))))
  print(paste0("Unique ORFs are from ",
               length(unique(unlist(sapply(seq_orf$gene, function(x){strsplit(x, split=",")})))),
               " genes!"))

  #separate the smORFs (50aa-100aa)
  seq_orf$len_strip <- as.numeric(sapply(1:dim(seq_orf)[1], function(i){
    strsplit(seq_orf$length[i], split="length:")[[1]][2]}))
  
  smORF_dat <- seq_orf[intersect(which(seq_orf$len_strip <=300), which(seq_orf$len_strip>=150)),]
  print(paste0("Unique number of smORFs (<=300 nt) with ", start, ": ", length(unique(smORF_dat$x))))
  print(paste0("Unique smORFs (<=300 nt) are from ",
               length(unique(unlist(sapply(smORF_dat$gene, function(x){strsplit(x, split=",")})))),
               " genes!"))

   #write the files above for later use!
  write.table(seq_orf, file=paste0("./unique_ORFs/", Sys.Date(), "_uniqueORFs_", outFileString, "_", start, ".tsv"), sep="\t")
  write.csv(smORF_dat, file=paste0("./unique_ORFs/", Sys.Date(), "_unique_smORFs_", outFileString, "_", start, ".tsv"), sep="\t")

  #writing to a fasta file for blastp
  sapply(1:dim(seq_orf)[1], function(x){
    write(paste(paste0(">", seq_orf[x,1]), paste(seq_orf[x,c(6, 7, 8)], collapse = " "), 
                paste0("genes:", seq_orf$genes[x]),
                paste0("accs:", seq_orf$accs[x]),
                paste0("num_transcripts:", seq_orf[x,3])), 
          file=paste0("./unique_ORFs/", Sys.Date(), "_fasta_uniqueORFs_", outFileString, "_", start, ".txt"), append=TRUE)
    write(paste0(seq_orf[x,2]), paste0("./unique_ORFs/", Sys.Date(), "_fasta_uniqueORFs_", outFileString, "_", start, ".txt"), append=TRUE)
  })
  
}
#run function
unique_ORFs_get(start=args[1], inFile=args[2], outFileString=args[3])
