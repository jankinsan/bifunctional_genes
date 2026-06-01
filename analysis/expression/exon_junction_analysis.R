#conda activate exp_env
# load required libraries
library(dplyr)
library(tidyr)

#set working directory!
setwd("/scratch/bioschool/phd/blz227562/exp_noncoding/recount3_dat")
junc_files <- files <- dir()[intersect(grep("gtex_", dir()), grep("bifunc_combined", dir()))]

dir.create("tissue_medians")
dir.create("stats")

tissue_medians <- list()

for (junc_file in junc_files){
  print(paste0("Reading ", junc_file, "..."))
  junc_dat <- read.csv(junc_file)
  tissue_details <- strsplit(junc_file, split="_bifunc_combined")[[1]][1]
  
  print("Computing key statistics..")
  print("Calculating noncoding weights for shared junctions..")
  #Rlative Noncoding Share
  #calculate the noncoding weight for each junction
  #number of noncoding transcripts in that junction divided by total number of transcripts having that annotated junction
  #assign zero weight to codung junctions!
  noncoding_weight <- junc_dat$nc_freq/(junc_dat$nc_freq + junc_dat$coding_freq)
  junc_dat <- cbind.data.frame(junc_dat[,1:7], noncoding_weight, junc_dat[,8:dim(junc_dat)[2]])
  
  #pivot longer so each sample has a row!
  junc_dat_long <- junc_dat %>%
    pivot_longer(
      cols = 9:last_col(),, 
      names_to = "sample_id", 
      values_to = "reads"
    )
  
  #A sample should have minimum 20 reads mapping to splice junctions of a gene!
  
  MIN_READS = 20
  #calculate relative noncoding isoform share!
  #rel_nc = (nc_reads + (shared_reads* nc_weight))/ sum of all reads
  print("calculate relative noncoding isoform share!")
  rel_nc_results <- junc_dat_long %>%
    group_by(sample_id, gene) %>%
    # 1. Calculate the total reads for the gene/sample group first
    mutate(total_spliced_reads = sum(reads, na.rm = TRUE)) %>%
    
    # 2. Drop the group if it doesn't meet your threshold
    filter(total_spliced_reads >= MIN_READS) %>%
    
    # 3. Summarize the remaining high-confidence groups
    summarize(
      total_spliced_reads = unique(total_spliced_reads),
      weighted_nc_reads = sum(reads * noncoding_weight, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    
    mutate(
      rel_nc = ifelse(total_spliced_reads > 0, weighted_nc_reads / total_spliced_reads, 0)
    )
  
  #Calculate unique noncoding read share!
  # Set your minimum unique read threshold per gene per sample
  print("Calculating unqiue noncoding fraction/share..")
  MIN_UNIQUE_READS <- 10 
  unc_results <- junc_dat_long %>%
    # Filter out shared junctions immediately
    filter(status != "Bifunc_shared") %>%
    
    group_by(sample_id, gene) %>%
    mutate(
      # Denominator: Sum of unique coding + unique noncoding reads
      total_unique_reads = sum(reads, na.rm = TRUE)
    ) %>%
    
    # Filter out samples that don't have enough unique raw evidence
    filter(total_unique_reads >= MIN_UNIQUE_READS) %>%
    
    summarize(
      total_unique_reads = unique(total_unique_reads),
      # Numerator: Sum of reads where weight is exactly 1 (purely non-coding)
      unique_nc_reads = sum(reads[noncoding_weight == 1], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      uncf = unique_nc_reads / total_unique_reads
    )
  
  #compute unique coding fraction!
  print("Computing unqiue coding fraction..")
  uc_results <- junc_dat_long %>%
    # Filter out shared junctions immediately
    filter(status != "Bifunc_shared") %>%
    
    group_by(sample_id, gene) %>%
    mutate(
      # Denominator: Sum of unique coding + unique noncoding reads
      total_unique_reads = sum(reads, na.rm = TRUE)
    ) %>%
    
    # Filter out samples lacking sufficient unique raw evidence
    filter(total_unique_reads >= MIN_UNIQUE_READS) %>%
    
    summarize(
      total_unique_reads = unique(total_unique_reads),
      # Numerator: Sum of reads where weight is exactly 0 (purely coding)
      unique_coding_reads = sum(reads[noncoding_weight == 0], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      ucf = unique_coding_reads / total_unique_reads
    )
  
  #Calculate for each gene & add!
  df_rel_nc_core <- rel_nc_results %>% select(sample_id, gene, rel_nc, total_spliced_reads, weighted_nc_reads)
  df_unc_core  <- unc_results %>% select(sample_id, gene, uncf, total_unique_reads, unique_nc_reads)
  df_uc_core   <- uc_results %>% select(sample_id, gene, ucf, unique_coding_reads)
  
  # Merge everything into a single master long-format dataframe
  print("Combining datasets..")
  combined_scores <- df_rel_nc_core %>%
    # Join the unique non-coding fraction
    full_join(df_unc_core, by = c("sample_id", "gene")) %>%
    # Join the unique coding fraction
    full_join(df_uc_core, by = c("sample_id", "gene"))
  
  write.csv(combined_scores , file=paste0("./stats/", tissue_details, "_combined_stats.csv"))
  
  # Preview the combined dataset..
  print("Combined dataset..")
  print(head(combined_scores))
  
  sig_bifunc <- combined_scores  %>%
    # 1. Group by gene to evaluate across its entire sample population
    group_by(gene) %>%
    
    # 2. Filter: keep genes where the proportion of active samples is > 30%
    # (mean() of a logical condition gives the proportion of TRUE values)
    filter(mean(uncf > 0.05 & !is.na(uncf)) > 0.30) %>%
    
    # 3. Remove grouping layer so downstream steps behave normally
    ungroup()
  
  print(paste0(length(unique(sig_bifunc$gene)), " genes have >5% splicing output in at least 30% of the samples.."))
  write.csv(sig_bifunc, file=paste0("./stats/", tissue_details,"_sig_bifunc_tissue.csv"))
  
  #Get the median tissue exp values using the unique noncoding share.. 
  uncf_tissue_matrix <- combined_scores %>%
    # 1. Group by gene to evaluate across the tissue's sample population
    group_by(gene) %>%
    
    # 2. Dynamically filter: Keep genes where >30% of the samples in this tissue have uncf > 0.05
    # Over 30% of ALL samples (including NAs) must be > 0.05
    filter(mean(uncf > 0.05 & !is.na(uncf)) > 0.30) %>%
    
    # 3. Keep only the active samples to calculate a clean representative value
    filter(uncf > 0.05) %>%
    
    # 4. Calculate the median for the surviving high-confidence genes
    summarize(
      median_uncf = median(uncf, na.rm = TRUE),
      .groups = "drop"
    )
  
  write.csv(uncf_tissue_matrix, file=paste0("./tissue_medians/", tissue_details, "_median_unc.csv"))
  
  #add details to list for further manipulation and plotting!
  tissue_medians[[tissue_details]]<- uncf_tissue_matrix
}

#tissue_medians <- lapply(tissue_medians, function(tissue_medians) tissue_medians[, 1:2])
#Now, combine the tissue-wise medians to use for heatmaps
tissue_medians_long <- data.frame()
combine_medians <- data.frame()

for (tissue in names(tissue_medians)){
  print(tissue)
  tissue_medians[[tissue]] <- cbind.data.frame(tissue_medians[[tissue]], tissue)
  tissue_medians_long <- rbind.data.frame(tissue_medians_long, tissue_medians[[tissue]])
} 
write.csv(tissue_medians_long, file="./stats/GTEX_COMBINED_TISSUE_MEDIANS_LONG_UNCF.csv")

#combine in a way to make a table for heatmap
combined_medians <- tissue_medians_long |> pivot_wider (names_from = tissue, 
                                                             values_from = median_uncf)
#write to file
write.csv(combined_medians, file="./stats/GTEX_COMBINED_TISSUE_MEDIANS_UNCF.csv")

#calculate tau specificity for each gene!