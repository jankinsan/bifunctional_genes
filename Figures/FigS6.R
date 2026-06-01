#longread data analysis 
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
setwd("D:/Janki/bifunctional genes/longread/")

#read expression counts and tpms!
counts <- read.delim("OUT.transcript_grouped_counts.tsv")
tpm <- read.delim("OUT.transcript_grouped_tpm.tsv")

#get CPC2 labels for the transcripts
#NOT USED IN THE ANALYSIS AS novel/discovered transcripts not analyzed in analysis!
chunk_files <- dir("./cpc")
cpc_labels <- data.frame()
for (chunk in chunk_files){
  print(paste0("Reading ", chunk))
  chunk_dat <- read.delim(paste0("./cpc/",chunk))
  cpc_labels <- rbind.data.frame(cpc_labels, chunk_dat)
}

#match cpc labels and move to a different data.frame!
matched_tpm <- cbind.data.frame(transcript=tpm[,1], label=cpc_labels$label[match(tpm[,1], cpc_labels$X.ID)], tpm[,-1])

#map this to bifunctional genes first!
acc <- read.csv("18122023_GRCh38.p14_gene_accessions.csv")
acc$accession <- sapply(acc$accession, function(x){
  strsplit(x, split=">")[[1]][2]})
gene_status <- read.csv("20062024__GRCh38.p14_genewise_status.csv")
gene_status$status[gene_status$status=="only_mRNA"] <- "Coding"
gene_status$status[gene_status$status=="only_ncRNA"] <- "Noncoding"
gene_status$status[gene_status$status=="hybrid"] <- "Bifunctional"

matched_tpm <- cbind.data.frame(transcript=matched_tpm[,1], 
                                mat_acc = gsub("_[0-9]+$", "", matched_tpm$transcript),
                                matched_tpm[,-1])

matched_tpm <- cbind.data.frame(matched_tpm[,c(1,2)], 
                                gene=acc$gene[match(matched_tpm$mat_acc, acc$accession)],
                                matched_tpm[,-c(1,2)])
#remove the unassigned RefSeq transcripts!
gene_matched_tpm <- matched_tpm[-which(is.na(matched_tpm$gene)),]
gene_matched_tpm2 <- cbind.data.frame(gene_matched_tpm[,c(1:3)], 
                                     gene_type = gene_status$status[match(gene_matched_tpm$gene, gene_status$X)],
                                     gene_matched_tpm[,-c(1:3)])

#formatting sample names for better analysis!
colnames(gene_matched_tpm2)[6:32]<- c("A549_rep1_run3", "A549_rep2_run1", "A549_rep3_run2",
                                "A549_rep5_run3","H9_rep2_run3", "H9_rep3_run3", "H9_rep4_run3", 
                                "HEYA8_rep1_run3", "HEYA8_rep2_run3", "HEYA8_rep3_run2",
                                "HCT116_rep1_run4", "HCT116_rep3_run2", "HCT116_rep4_run1", 
                                "HCT116_rep5_run1", "HEPG2_rep1_run1", "HEPG2_rep1_run2",
                                "HEPG2_rep2_run1", "HEPG2_rep4_run1", "HEPG2_rep4_run2", 
                                "HEPG2_rep5_run3", "K562_rep1_run2", "K562_rep2_run1", "K562_rep3_run1",
                                "K562_rep4_run2", "MCF7_rep1_run2", "MCF7_rep3_run3", 
                                "MCF7_rep4_run2")

#get the unique number of genes that are expressed in at least one cell line!
length(unique(gene_matched_tpm2$gene))

#keep only rows where rowSums TPM >=1
subset_gene_matched_tpm <- gene_matched_tpm2[which(rowSums(gene_matched_tpm2[,6:32])>=1),]
length(unique(subset_gene_matched_tpm$gene))

#get tpm fractions for type of genes!
bifunc_idx <- which(subset_gene_matched_tpm$gene_type=="Bifunctional")
coding_idx <- which(subset_gene_matched_tpm$gene_type=="Coding")
noncoding_idx <- which(subset_gene_matched_tpm$gene_type=="Noncoding")

#sum of TPMs for RefSeq transcripts
tpm_sums <- apply(subset_gene_matched_tpm[,-c(1:5)], 2, function(x){
  c(Bifunctional=sum(x[bifunc_idx]),
    Coding=sum(x[coding_idx]),
    Noncoding = sum(x[noncoding_idx]))
})
tpm_sums_plot <- reshape2::melt(tpm_sums)
tpm_sums_plot$cell_line <- as.factor(sapply(tpm_sums_plot$Var2, function(x){
  strsplit(as.character(x), split="_")[[1]][1]
}))

write.csv(tpm_sums_plot, file="tpm_sums_plot.csv")
plot_stacked_annot <- ggplot(data = tpm_sums_plot, 
                             aes(x =Var2, y=log2(value), 
                                 fill=Var1, width=0.65)) +
  geom_bar(position="stack", stat="identity")+ #facet_wrap(~cell_line)+
  labs(y="Log2(TPM)", x="Replicates", fill="Gene Annotation") +
  theme(axis.text.x = element_text(size=13, colour="black",angle = 90, vjust = 0.5, hjust=1), 
        axis.text.y = element_text(size=13, colour="black"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="bottom",
        legend.text = element_text(size=13, colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("Coding"="#e8e288", "Bifunctional" = "#82c0cc", 
                               "Noncoding" = "#ffa5a5"))
plot(plot_stacked_annot)  

ggplot(data = tpm_sums_plot, 
       aes(x =Var2, y=value, 
           fill=Var1, width=0.65)) +
  geom_bar(position="stack", stat="identity")+ #facet_wrap(~cell_line)+
  labs(y="Log2(TPM)", x="Replicates", fill="Gene Annotation") +
  theme(axis.text.x = element_text(size=13, colour="black",angle = 90, vjust = 0.5, hjust=1), 
        axis.text.y = element_text(size=13, colour="black"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="bottom",
        legend.text = element_text(size=13, colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("Coding"="#e8e288", "Bifunctional" = "#82c0cc", 
                               "Noncoding" = "#ffa5a5"))
  

bifunc_genes <- unique(gene_matched_tpm2$gene[which(gene_matched_tpm2$gene_type=="Bifunctional")])

#get number of transcripts and kind of transcripts reads are mapping to for bifunctional genes!
gene_status_longread<- t(sapply(bifunc_genes, function(gene){
  matched_dat <- gene_matched_tpm2[gene_matched_tpm2$gene==gene,]
  if(length(grep("^[NX]M_", matched_dat$transcript))>=1){
    if(length(grep("^[NX]R_", matched_dat$transcript))>=1){
      return(c(gene, "Bifunctional both", length(grep("^[NX]M_", matched_dat$transcript)), length(grep("^[NX]R_", matched_dat$transcript)),
               length(grep("^NXR_", matched_dat$transcript))))
    } else {
      return(c(gene, "Bifunctional coding", length(grep("^[NX]M_", matched_dat$transcript)), 0,0))
    }
  } else {
    return(c(gene, "Bifunctional noncoding", 0, length(grep("^[NX]R_", matched_dat$transcript)),length(grep("^XR_", matched_dat$transcript))))
  }
  }))
colnames(gene_status_longread) <- c("gene", "transcript_evidence", "num_coding_trans", "num_nc_trans", "num_xr_trans")

gene_status_longread[,2] <- gsub("Bifunctional ", "", gene_status_longread[,2])
#make a pie chart here!
gene_prop_plot <- data.frame(table(gene_status_longread[,2]))
ggplot(gene_prop_plot, aes(x=" ", y=Freq, fill=Var1)) +
  geom_bar(width =1, stat = "identity", color="white", linewidth=1, alpha=0.8) +
  coord_polar("y", start=0) + labs(title="Bifunctional genes with detectable TPMs",
                                   fill = "Transcript type expressed")+
  geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), size=6.75)+
  theme(
    plot.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(size=15, colour = "black", face="bold", hjust=0.5),
    axis.line.x = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position="bottom",
    strip.text.y = element_text(size=13, colour = "black", face="bold"),
    legend.title = element_text(size=13, color="black", face="bold"),
    legend.text = element_text(size=12,colour = "black"))+ 
  scale_fill_manual(values = c("both"="#04C1C5", "coding" = "#F8766D", 
                               "noncoding" = "darkgrey"))

#calculate median across replicates and then plot for the cell lines!
#get genes with coding and noncoding expression in a gene!
bifunc_both <- gene_status_longread[which(gene_status_longread[,2]=="both"), 1]
plot_bifunctional_both <- subset_gene_matched_tpm %>%
  filter(gene %in%  bifunc_both)

#pivot longer to get tpm values separately!
plot_bifunc <- pivot_longer(plot_bifunctional_both, 
                            cols = 6:last_col(),  
                            names_to = "Cell_Line_Rep",
                            values_to = "TPM")


print("Computing unqiue coding fraction..")
unc_results <- plot_bifunc %>%
  group_by(gene, Cell_Line_Rep) %>%
  summarize(
    total_tpm = sum(TPM, na.rm = TRUE),
    unique_nc_tpm = sum(TPM[grepl("^[NX]R_", transcript)], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    uncf = if_else(total_tpm > 0, unique_nc_tpm / total_tpm, 0)
  )

#get cell lines to get median noncoding fraction!
unc_results$Cell_Line <- gsub("_rep[0-9]+_run[0-9]+", "", unc_results$Cell_Line_Rep)

#get cell-line medians for each gene!
unc_median_tpm_ratio <- unc_results %>%
  group_by(gene, Cell_Line) %>%
  summarize(median_nc_tpm_ratio=median(uncf, na.rm=TRUE), .groups="drop")

#make a compact dataframe!
heatmap_unc <- pivot_wider(unc_median_tpm_ratio, names_from = Cell_Line, 
                           values_from=median_nc_tpm_ratio)
#keep only genes which have non-zero median noncoding TPM ratios
heatmap_unc<- data.frame(heatmap_unc[rowSums(heatmap_unc[,-1])>0,])
rownames(heatmap_unc )<- heatmap_unc[,1]

#save file for later use!
write.csv(heatmap_unc, file="nc_tpm_ratio_bifunctional_genes.csv")

#add a heatmap here!
library("ComplexHeatmap")
col_fun2=ComplexHeatmapcol_fun2= colorRamp2(c(0, 0.10, 1), c("white", "#F8766D", "#04C1C5"))
#svg("trial_heatmap_top_genes.svg", width=700, height=550)
Heatmap(t(heatmap_unc[,-1]), name="Median NC TPM", col=col_fun2, na_col = "white",
        show_column_names = FALSE)


#get dominant transcript for each gene in each cell line based on median tpm across replicates!
plot_bifunc$Cell_Line <- gsub("_rep[0-9]+_run[0-9]+", "", plot_bifunc$Cell_Line_Rep)
transcript_medians <- plot_bifunc %>%
  group_by(gene, transcript, Cell_Line) %>%
  summarize(Median_TPM = median(TPM, na.rm = TRUE), .groups = "drop") 

dominant_transcript <- transcript_medians %>% #select the most expressed transcript
  group_by(gene, Cell_Line) %>%
  filter(Median_TPM == max(Median_TPM) & Median_TPM > 0)

dominant_transcript$transcript_type <- "coding"
dominant_transcript$transcript_type[grep("^[NX]R_", dominant_transcript$transcript)] <- "noncoding"
#get number of genes for which you have the dominant transcript
length(unique(dominant_transcript$gene))

length(unique(dominant_transcript$gene[dominant_transcript$transcript_type=="noncoding"]))

#save the data for resuse!
write.csv(dominant_transcript, file="dominant_transcripts_median_TPM_longread.csv")

#add labels to top 5 genes for plotting!
lab_plot <- dominant_transcript %>% group_by(Cell_Line, transcript_type) %>%
  mutate(label = ifelse(rank(dplyr::desc(Median_TPM)) <= 5, as.character(gene), "")) %>%
  ungroup()


ggplot(lab_plot, aes(x=Cell_Line, y=log2(Median_TPM), color=transcript_type))+
  geom_point(size=2.5, alpha=0.5,
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),)+ 
  geom_text_repel(aes(label = label, color = transcript_type),
                  # This MUST match the position_jitterdodge parameters used above
                  position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
                  size = 3.5,                  # Font size of labels
                  fontface = "bold",
                  max.overlaps = 20,         # Prevents clutter by hiding overlapping labels
                  box.padding = 0.3,         # Space around the text box
                  point.padding = 0.2,       # Space between point and text
                  segment.color = "grey50",  # Color of the line connecting text to dot
                  segment.alpha = 0.6,
                  show.legend = FALSE) +
  labs(y="Log2(Median_TPM)", x="Cell Line", color="Dominant Transcript Type") +
  theme(axis.text.x = element_text(size=13, colour="black"), 
        axis.text.y = element_text(size=13, colour="black"),
        axis.title = element_text(size=13, face="bold"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="bottom",
        legend.title =  element_text(size=13, colour = "black", face="bold"),
        legend.text = element_text(size=13, colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())#+
  #scale_color_manual(values = c("noncoding" = "#414487FF", "coding"="#E3e418ff"))
