library(recount3)

#Match and separate data for specific junction types!
#the RefSeq junctions have been annotated separately previously and they are just being loaded
setwd("/scratch/bioschool/phd/blz227562/")

dir.create("recount3_dat")
shared_juncs <- read.csv("bifunc_genes_shared_junctions.csv")
nc_only_juncs <- read.csv("bifunc_genes_noncoding_junctions.csv")
coding_only_juncs <- read.csv("bifunc_genes_coding_junctions.csv")

coding_genes_juncs <- read.csv("coding_genes_junctions.csv")
noncoding_genes_juncs <- read.csv("noncoding_genes_junctions.csv")

#add a status label: Bifunc_coding, Bifunc_noncoding, Bifunc_shared
coding_only_juncs$status <- "Bifunc_coding"
nc_only_juncs$status <- "Bifunc_noncoding"
shared_juncs$status <- "Bifunc_shared"
coding_genes_juncs$status <- "Coding_gene"
noncoding_genes_juncs$status <- "Noncoding_gene"

#make NA's to 0 for transcript frequencies!
coding_only_juncs$coding_freq [is.na(coding_only_juncs$coding_freq)]<- 0
nc_only_juncs$coding_freq[is.na(nc_only_juncs$coding_freq)] <- 0
shared_juncs$coding_freq[is.na(shared_juncs$coding_freq)] <- 0
coding_genes_juncs$coding_freq[is.na(coding_genes_juncs$coding_freq)] <- 0
noncoding_genes_juncs$coding_freq[is.na(noncoding_genes_juncs$code_freq)] <- 0

#noncoding frequency 
coding_only_juncs$nc_freq [is.na(coding_only_juncs$nc_freq)]<- 0
nc_only_juncs$nc_freq[is.na(nc_only_juncs$nc_freq)] <- 0
shared_juncs$nc_freq[is.na(shared_juncs$nc_freq)] <- 0
coding_genes_juncs$nc_freq[is.na(coding_genes_juncs$nc_freq)] <- 0
noncoding_genes_juncs$nc_freq[is.na(noncoding_genes_juncs$code_freq)] <- 0

#label 'none' where there might be no transcript matching to that junction
coding_only_juncs$nc_transcript[is.na(coding_only_juncs$nc_transcript)]<- "none"
nc_only_juncs$coding_transcript[is.na(nc_only_juncs$coding_transcript)] <- "none"
shared_juncs$coding_transcript[is.na(shared_juncs$coding_transcript)] <- "none"
shared_juncs$nc_transcript[is.na(shared_juncs$nc_transcript)] <- "none"
coding_genes_juncs$nc_transcript[is.na(coding_genes_juncs$nc_transcript)] <- "none"
noncoding_genes_juncs$coding_transcript[is.na(noncoding_genes_juncs$coding_transcript)] <- "none"

#Combining junction information for all bifunctional genes
combined_juncs<- rbind.data.frame(shared_juncs,nc_only_juncs, coding_only_juncs)

#Start data retrieval from recount3 and map the annotated junctions!
#Get a list of all projects from recount and subset only gtex and tcga data!
human_projects <- available_projects()
gtex_tcga_projects <- human_projects[c(which(human_projects$file_source=="tcga"), which(human_projects$file_source=="gtex")),]

for (idx in 1: dim(gtex_tcga_projects)[1]){
  print(paste0("Downloading ", gtex_tcga_projects$project[idx], " data from ", gtex_tcga_projects$file_source[idx], "..."))
  
  #get code from TCGA/GTEx study to download data!
  project_code <- gtex_tcga_projects$project[idx]
  project_info <- subset(human_projects, project == project_code & project_type == "data_sources")
  
  #retrieving junction data from recount
  print("Retrieving junction data.. ")
  junc_data <- create_rse(project_info, type = "jxn")
  #getting coordinates and counts
  junc_coords <- rowRanges(junc_data)
  junc_counts <- as.matrix(junc_data@assays@data@listData[["counts"]])
  
  #mapping the junctions!
  print("Matching the annotated junctions with recount3 data..")
  #for bifunctional genes
  combined_junc_counts <- rbind.data.frame(junc_counts[match(as.character(shared_juncs$junction), rownames(junc_counts)),], 
                                           junc_counts[match(as.character(nc_only_juncs$junction), rownames(junc_counts)),],
                                           junc_counts[match(as.character(coding_only_juncs$junction), rownames(junc_counts)),])
  combined_junc_counts[is.na(combined_junc_counts)] <- 0
  combined_junc_counts <- cbind.data.frame(combined_juncs, combined_junc_counts)
  
  #for coding and non-coding genes
  coding_junc_counts <- as.data.frame(junc_counts[match(as.character(coding_genes_juncs$junction), rownames(junc_counts)),])
  coding_junc_counts[is.na(coding_junc_counts)] <- 0
  coding_junc_counts <- cbind.data.frame(coding_genes_juncs, coding_junc_counts)
  noncoding_junc_counts <- as.data.frame(junc_counts[match(as.character(noncoding_genes_juncs$junction), rownames(junc_counts)),])
  noncoding_junc_counts[is.na(noncoding_junc_counts)] <- 0  
  noncoding_junc_counts <- cbind.data.frame(noncoding_genes_juncs, noncoding_junc_counts)

  print("Writing data to files..")
  write.table(combined_junc_counts, file=paste0("/recount3_dat/", gtex_tcga_projects$file_source[idx], "_", gtex_tcga_projects$project[idx], "_bifunc_combined_junc_counts.csv"), sep=",")
  write.table(coding_junc_counts, file=paste0("/recount3_dat/", gtex_tcga_projects$file_source[idx], "_", gtex_tcga_projects$project[idx], "_coding_genes_junc_counts.csv"), sep=",")
  write.table(noncoding_junc_counts, file=paste0("/recount3_dat/", gtex_tcga_projects$file_source[idx], "_", gtex_tcga_projects$project[idx], "_noncoding_genes_junc_counts.csv"), sep=",")

}

