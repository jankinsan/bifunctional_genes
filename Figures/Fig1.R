## RUNNING THIS FROM RStudio
#loading libraries
#install.packages("PupillometryR")
library(ggplot2)
library(viridis)
library(PupillometryR)


setwd("E:/Janki/bifunctional_genes/data")

#get genes, match status and plot!
hg38_genes <- read.csv("GCF_000001405.40_GRCh38.p14_genes.csv")
hg38_genes$length <- hg38_genes$end - hg38_genes$start
#keep only genes transcribed!
hg38_genes2 <- hg38_genes[hg38_genes$gene_type!="not_assigned",]

#read the transcript data!
hg38_transcript <- read.csv("GCF_000001405.40_GRCh38.p14_transcripts.csv")
hg38_transcript$length <- hg38_transcript$end - hg38_transcript$start

gene_status <- read.csv("./counts/18122023_GRCh38.p14_genewise_counts_status.csv")
#update the symbols so that the gene symbols match!
changed_idx <- which(is.na(match(gene_status$X, hg38_genes$gene)))
gene_status$final_gene <- gene_status$X
gene_status$final_gene[changed_idx] <- c("HSALR1", "GAR1-DT")
gene_status$status[gene_status$status=="only_mRNA"] <- "Coding"
gene_status$status[gene_status$status=="only_ncRNA"] <- "Noncoding"
gene_status$status[gene_status$status=="hybrid"] <- "Bifunctional"

write.csv(gene_status, file="./counts/23112024_GRCh38.p14_genewise_status.csv", row.names=FALSE)
hg38_genes$gene_status<-gene_status$status[match(hg38_genes$gene, gene_status$final_gene)]

#Getting gene numbers to plot!
#total summary including genes that have no transcripts
hg38_genes$gene_biotype <- sapply(1:dim(hg38_genes)[1], function(x){
  strsplit(strsplit(hg38_genes[x,9], split="gene_biotype ")[[1]][2], split=";")[[1]][1]
})
gene_data_summary <- data.frame(gene_type=unique(hg38_genes$gene_biotype))
gene_data_summary$number <- sapply(gene_data_summary$gene_type, function(x){
  length(which(hg38_genes$gene_biotype==x))
})
#getting the biotypes of transcribed genes
hg38_genes2<-hg38_genes[-which(is.na(hg38_genes$gene_status)),]
gene_data_summary$number_trans <- sapply(gene_data_summary$gene_type, function(x){
  length(which(hg38_genes2$gene_biotype==x))
})
gene_data_summary$percent_num <- gene_data_summary$number/sum(gene_data_summary$number)
gene_data_summary$percent_num_trans <- gene_data_summary$number_trans/sum(gene_data_summary$number_trans)

#from transcript data, remove unassigned transcripts!
hg38_transcript2 <- hg38_transcript[-grep("unassigned_transcript", hg38_transcript$acc),]
#this still leaves some transcripts!
acc<- read.csv("./counts/18122023_GRCh38.p14_gene_accessions.csv")
#remove zeroes
acc<- acc[-which(acc$accession==0),]
#match only the IDs from this file!
acc$ids <- sapply(acc$accession, function(x){
  substr(x, 2, nchar(x))
})
#some transcripts are on multiple loci, so they have multiple coordinates.
#but for this analysis, we are using the first match, the exact match. 
#will take care of these in the later stages and during expression quantification!
hg38_transcripts_matched<- hg38_transcript[match(acc$ids, hg38_transcript$acc),]
unmat <- hg38_transcript[-match(acc$ids, hg38_transcript$acc),]
hg38_transcripts_matched$gene_type <- gene_status$status[match(acc$gene[match(hg38_transcripts_matched$acc, acc$ids)], gene_status$X)]

#Figure 1A
#STACKED BAR PLOT OF GENE TYPE---------------------------
summ_plot <- reshape2::melt(gene_data_summary[,c(1,5,4)])
write.csv(gene_data_summary, file="counts/02072024_NCBI_Annotation_numbers_human.csv")

#this is for horizontal, not running this now!
plot_stacked_type <- ggplot(data = summ_plot, 
                            aes(y = variable, x=value, 
                                fill=gene_type, width=0.65)) +
  geom_bar(position="stack", stat="identity", colour="black")+
  labs(x="Fraction of Genes", y="", fill="Gene Type") +
  theme(axis.text.x = element_text(size=13, colour="black"),
        axis.title.x=element_text(size=13, color="black", face="bold"), 
        axis.text.y=element_text(size=13, color="black", face="bold"), 
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="bottom",
        legend.title=element_text(size=13, color="black", face="bold"),
        legend.text = element_text(size=12, colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())+
  #scale_fill_viridis(discrete=TRUE, direction=-1)+
  scale_y_discrete(labels=c("Genes with transcripts","All Genes"))

#pdf("F:/Janki/bifunctional_genes/figures/Fig1A.pdf")
plot(plot_stacked_type)
#dev.off()

#Figure 1B
#alluvial plot to represent the ncbi biotype and our annotation
library("ggalluvial")
library("ggrepel")
library("data.table")
library("dplyr")

#making combinations
#putting gene categories with less than 1000 genes as others
hg38_genes2$gene_status <- recode(hg38_genes2$gene_status,
                                  "ncRNA encoding" = "Noncoding",
                                  "mRNA encoding" = "Coding",
                                  "bifunctional" = "Bifunctional")

hg38_genes2$gene_biotype2 <-hg38_genes2$gene_biotype
idx <- unique(c(which(hg38_genes2$gene_biotype=="protein_coding"), 
                which(hg38_genes2$gene_biotype=="lncRNA"),
                which(hg38_genes2$gene_biotype=="miRNA"),
                which(hg38_genes2$gene_biotype=="transcribed_pseudogene"),
                which(hg38_genes2$gene_biotype=="snoRNA")))
hg38_genes2$gene_biotype2[-idx] <- "others"
combs <- CJ(hg38_genes2$gene_biotype2, hg38_genes2$gene_status, unique=TRUE)
allu_dat <- data.frame(combs)

allu_dat$freq <- sapply(1:dim(combs)[1], function(x){
  length(intersect(which(hg38_genes2$gene_biotype==allu_dat[x,1]),
                   which(hg38_genes2$gene_status==allu_dat[x,2])))
}) 
colnames(allu_dat) <- c("ncbi", "new_annot", "freq")

#alluvial plot
allu_plot <- ggplot(data = allu_dat,
       aes(axis1 = ncbi, axis2 = new_annot, y = freq)) +
  scale_x_discrete(limits = c("NCBI Biotype", "Our Annotation"), expand = c(.2, .05)) +
  geom_flow(aes(fill = new_annot), width=0.3) + 
  geom_stratum(aes(fill = new_annot), width=0.3) + 
  labs(y="Number of Genes", fill="Gene Type") +
  geom_label_repel(stat = "stratum", aes(label = after_stat(stratum),size=11)) + guides (size="none")+
  scale_fill_manual(values = c("Coding"="#e8e288", "Bifunctional" = "#82c0cc", "Noncoding" = "#ffa5a5"))+
  theme_bw() + theme(legend.title=element_text(size=13,face="bold"),
                    legend.text=element_text(size=13),
                    legend.position="bottom",
                     axis.title = element_text(size=13, face="bold"),
                     axis.text.x = element_text(size=13,face="bold", color="black"),
                    axis.text.y = element_text(size=13, color="black"),
                    legend.key.size = unit(1, "lines"))
plot(allu_plot)

#pdf("F:/Janki/bifunctional_genes/figures/Figure1B.pdf")
plot(allu_plot)
#dev.off()

##Figure 1G 
###CHROMOSOMAL DISTRIBUTION OF ALL KINDS OF GENES
#remove the genes with both X and Y chromosomes 
#adding these numbers to either of these chromosomes
hg38_genes2$chr_simplified <- sapply(hg38_genes2$chr, function(x){
 strsplit(x, split="_")[[1]][1] 
})
chr_dist <- data.frame(chr=unique(hg38_genes2$chr_simplified))
chr_dist$all_genes <- sapply(chr_dist$chr, function(x){
  length(which(hg38_genes2$chr_simplified==x))
})
chr_dist$Coding <- sapply(chr_dist$chr, function(x){
  idx <- which(hg38_genes2$gene_status=="Coding")
  length(which(hg38_genes2$chr_simplified[idx]==x))
})
chr_dist$Bifunctional <- sapply(chr_dist$chr, function(x){
  idx <- which(hg38_genes2$gene_status=="Bifunctional")
  length(which(hg38_genes2$chr_simplified[idx]==x))
})
chr_dist$Noncoding <- sapply(chr_dist$chr, function(x){
  idx <- which(hg38_genes2$gene_status=="Noncoding")
  length(which(hg38_genes2$chr_simplified[idx]==x))
})

chr_dist$percent_mrna <- chr_dist$Coding/chr_dist$all_genes
chr_dist$percent_ncrna <- chr_dist$Noncoding/chr_dist$all_genes
chr_dist$percent_bifunctional <- chr_dist$Bifunctional/chr_dist$all_genes


#FIGURE 1H
##strandedness of bifunctional genes chromosome-wise!
hg38_bifunc <- hg38_genes2[hg38_genes2$gene_status=="Bifunctional",]
chr_dist$bifunc_pos <- sapply(chr_dist$chr, function(x){
  length(intersect(which(hg38_bifunc$chr_simplified==x), which(hg38_bifunc$strand=="+")))
})
chr_dist$bifunc_neg <- sapply(chr_dist$chr, function(x){
  length(intersect(which(hg38_bifunc$chr_simplified==x), which(hg38_bifunc$strand=="-")))
})

chr_dist$percent_bifunc_pos <- chr_dist$bifunc_pos/(chr_dist$bifunc_pos+chr_dist$bifunc_neg)
chr_dist$percent_bifunc_neg <- chr_dist$bifunc_neg/(chr_dist$bifunc_pos+chr_dist$bifunc_neg)
#mosaic and stacked bar plots for chromosomal distribution
#ordering for the plot to be chromosome-wise!
plot_chr_dist<- reshape2::melt(chr_dist[,c(1,3:5)])
plot_chr_dist$variable<- as.factor(plot_chr_dist$variable)
plot_chr_dist$chr<- factor(plot_chr_dist$chr, levels=chr_dist$chr)
chr_plot<- ggplot(data = plot_chr_dist,aes(x = chr, y=value, 
                                           fill=variable, width=0.75)) +
  geom_bar(position="stack", stat="identity", colour="black")+
  labs(y="Number of Genes", x="", fill="Gene Annotation") +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size=13, face="bold", color="black"),
        axis.text.y = element_text(size=13, color="black"),
        axis.ticks.x=element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="none",
        legend.text = element_text(size=13, colour = "black", face="bold"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("Coding"="#e8e288", "Bifunctional" = "#82c0cc", "Noncoding" = "#ffa5a5"))

plot_chr_dist2<- reshape2::melt(chr_dist[,c(1,6:8)])
plot_chr_dist2$variable <- sapply(plot_chr_dist2$variable, function(x){
  if (x=="percent_mrna"){
    return ("Coding")
  } else if (x == "percent_ncrna") {
    return ("Noncoding")
  } else { return ("Bifunctional")}
})
plot_chr_dist2$variable<- factor(plot_chr_dist2$variable, levels=c("Coding", "Bifunctional", "Noncoding"))
plot_chr_dist2$chr<- factor(plot_chr_dist2$chr, levels=chr_dist$chr)
chr_plot2<- ggplot(data = plot_chr_dist2,aes(x = chr, y=value, 
                                           fill=variable, width=0.75)) +
  geom_bar(position="stack", stat="identity", colour="black")+
  labs(y="Fraction of Genes", x = "Chromosome", fill="Gene Annotation") +
  theme(axis.text.x = element_text(size=13, colour="black", angle = 90),
        axis.text.y = element_text(size=13, colour="black"),
        axis.title = element_text(size=13, face="bold", color="black"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="bottom",
        legend.text = element_text(size=13),
        legend.title = element_text(size=13, colour = "black", face="bold"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("Coding"="#e8e288", "Bifunctional" = "#82c0cc", "Noncoding" = "#ffa5a5"))

#strandedness plot 
plot_chr_strand<- reshape2::melt(chr_dist[,c(1,11,12)])
plot_chr_strand$variable <- sapply(plot_chr_strand$variable, function(x){
  if (x=="percent_bifunc_pos"){
    return ("+")
  } else { return ("-")}
})
plot_chr_strand$variable<- factor(plot_chr_strand$variable, levels=c("+", "-"))
plot_chr_strand$chr<- factor(plot_chr_strand$chr, levels=chr_dist$chr)
chr_plot_strand<- ggplot(data = plot_chr_strand,aes(x = chr, y=value, 
                                             fill=variable, width=0.75)) +
  geom_bar(position="stack", stat="identity", colour="black")+
  labs(y="Bifunctional Gene Fraction", x= "Chromosome", fill="Strand") +
  theme(axis.text.x = element_text(size=13, colour="black", angle = 90),
        axis.text.y = element_text(size=13, colour="black"),
        axis.title = element_text(size=13, color="black", face="bold"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="bottom",
        legend.title = element_text(size=13, face="bold"),
        legend.text = element_text(size=13, colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank()) #+
 # scale_fill_viridis(direction = -1, discrete = TRUE)

#Numbers of stranded genes
plot_chr_strand2<- reshape2::melt(chr_dist[,c(1,9,10)])
plot_chr_strand2$variable <- sapply(plot_chr_strand2$variable, function(x){
  if (x=="bifunc_pos"){
    return ("+")
  } else { return ("-")}
})
plot_chr_strand2$variable<- factor(plot_chr_strand2$variable, levels=c("+", "-"))
plot_chr_strand2$chr<- factor(plot_chr_strand2$chr, levels=chr_dist$chr)
chr_plot_strand2<- ggplot(data = plot_chr_strand2,aes(x = chr, y=value, 
                                                    fill=variable, width=0.75)) +
  geom_bar(position="stack", stat="identity", colour="black")+
  labs(y="Number of Bifunctional Genes", fill="Strand") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=13, color="black"),
        axis.title.y = element_text(size=13, face="bold"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="none",
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())

#pdf("F:/Janki/bifunctional_genes/figures/Fig1G_1H.pdf")
plot(chr_plot)
plot(chr_plot2)
plot(chr_plot_strand)
plot(chr_plot_strand2)
#dev.off()

#save data
write.csv(chr_dist, file="./23112024_New_annotation_chr_dist.csv", row.names = FALSE)

##Figure 1D
#Transcript Numbers per gene
#p-value to calculate between groups
library(ggpubr)
library(rstatix)
my_comparisons <- list( c("Bifunctional", "Coding"), c("Coding", "Noncoding"),
                        c("Bifunctional", "Noncoding") )

#plot
gene_status$status <- recode(gene_status$status,
                              "ncRNA encoding" = "Noncoding",
                              "mRNA encoding" = "Coding",
                              "bifunctional" = "Bifunctional")
trans_num_plot <- ggplot(data = gene_status,
                  aes(x = status,y=log2(total), fill=status)) +
  geom_violin(width=1.6)+
  geom_boxplot(width=0.1)+
  labs(y="log2 (Total number of transcripts)", fill="Gene Annotation", x="") +
  ylim(0,10)+
  theme(axis.text = element_text(size=13, colour="black"),
        axis.title.y = element_text(size=13, face="bold"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="none",
        legend.title = element_text(size=12, colour = "black", face="bold"),
        legend.text = element_text(size=11, colour = "black", face="bold"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())+ 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "Bifunctional", size=6, label.y = 9.5)+
  scale_fill_manual(values = c("Coding"="#e8e288", "Bifunctional" = "#82c0cc", "Noncoding" = "#ffa5a5"))

plot(trans_num_plot)

#calculate effect sizes for number of transcripts
##bifunctional versus other groups
bifunc_cod_dat <- gene_status[c(which(gene_status$status=="Bifunctional"),
                                which(gene_status$status=="Coding")), ]
bifunc_cod_dat$log2_total <- log2(bifunc_cod_dat$total)
cohens_d(log2_total ~ status, data = bifunc_cod_dat)

bifunc_noncod_dat <- gene_status[c(which(gene_status$status=="Bifunctional"),
                                which(gene_status$status=="Noncoding")), ]
bifunc_noncod_dat$log2_total <- log2(bifunc_noncod_dat$total)
cohens_d(log2_total ~ status, data = bifunc_noncod_dat)

cod_noncod_dat <- gene_status[c(which(gene_status$status=="Coding"),
                                   which(gene_status$status=="Noncoding")), ]
cod_noncod_dat$log2_total <- log2(cod_noncod_dat$total)
cohens_d(log2_total ~ status, data = cod_noncod_dat)


##Figure 1E
#Gene lengths
gene_len_plot <- ggplot(data = hg38_genes2,
                        aes(x = gene_status, y=log2(length), fill=gene_status)) +
  geom_violin(width=0.8)+ ylim(5, 23)+
  geom_boxplot(width=0.1)+ labs(y="log2 (Gene length)", fill="Gene Annotation", x="") +
  theme(axis.text = element_text(size=13, colour="black"),
        axis.title.y = element_text(size=13, face="bold"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="none",
        legend.title = element_text(size=12, colour = "black", face="bold"),
        legend.text = element_text(size=11, colour = "black", face="bold"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())+ 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "Bifunctional", size=6, label.y = 21.5)+
  scale_fill_manual(values = c("Coding"="#e8e288", "Bifunctional" = "#82c0cc", "Noncoding" = "#ffa5a5"))
plot(gene_len_plot)

#calculate effect sizes for gene lengths 
##bifunctional versus other groups
bifunc_cod_len <- hg38_genes2[c(which(hg38_genes2$gene_status=="Bifunctional"),
                                which(hg38_genes2$gene_status=="Coding")), ]
bifunc_cod_len$log2_length <- log2(bifunc_cod_len$length)
cohens_d(log2_length ~ gene_status, data = bifunc_cod_len)

bifunc_noncod_len <- hg38_genes2[c(which(hg38_genes2$gene_status=="Bifunctional"),
                                   which(hg38_genes2$gene_status=="Noncoding")), ]
bifunc_noncod_len$log2_length <- log2(bifunc_noncod_len$length)
cohens_d(log2_length ~ gene_status, data = bifunc_noncod_len)

cod_noncod_len <- hg38_genes2[c(which(hg38_genes2$gene_status=="Coding"),
                                which(hg38_genes2$gene_status=="Noncoding")), ]
cod_noncod_len$log2_length <- log2(cod_noncod_len$length)
cohens_d(log2_length ~ gene_status, data = cod_noncod_len)

#RUN EXON SPECIFIC ANALYSIS ON HPC!
#back to pc for plotting!
library(dplyr)
library(tidyverse)
number_exons<- read.csv("10072024_number_exons_per_transcript.csv")
exons_keep <- number_exons[-which(number_exons$gene_type=="not_assigned"),]

#update annotation
exons_keep$gene_type[exons_keep$gene_type=="only_ncRNA"]<- "Noncoding"
exons_keep$gene_type[exons_keep$gene_type=="only_mRNA"]<- "Coding"
exons_keep$gene_type[exons_keep$gene_type=="hybrid"]<- "Bifunctional"

#keeping only the transcript with the highest number of exons for each gene
exons_keep_plot <- na.omit(exons_keep %>% group_by(gene) %>%
  slice_max(num, n = 1, with_ties = FALSE)  %>% ungroup())

hg38_genes2$num_exons <- exons_keep_plot$num[match(hg38_genes2$gene, exons_keep_plot$gene)]
plot_exons <- ggplot(data = exons_keep_plot,
                     aes(x = as.factor(gene_type), y=log2(num), fill=gene_type)) +
  geom_violin(width=0.8)+ ylim(0,9)+
  geom_boxplot(width=0.1)+ labs(y="log2 (Number of exons)", fill="Gene Annotation", x="") +
  theme(axis.text = element_text(size=13, colour="black"),
        axis.title.y = element_text(size=13, face="bold"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="right",
        legend.title = element_text(size=13, colour = "black", face="bold"),
        legend.text = element_text(size=12, colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())+ 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "Bifunctional", size=6, label.y = 8.75)+
  scale_fill_manual(values = c("Coding"="#e8e288", "Bifunctional" = "#82c0cc", "Noncoding" = "#ffa5a5"))
plot(plot_exons)

#calculate effect sizes
##bifunctional versus other groups
bifunc_cod_num_exon <- exons_keep_plot[c(which(exons_keep_plot$gene_type=="Bifunctional"),
                                         which(exons_keep_plot$gene_type=="Coding")), ]
bifunc_cod_num_exon$log2_num_exon <- log2(bifunc_cod_num_exon$num)
cohens_d(log2_num_exon ~ gene_type, data = bifunc_cod_num_exon)

bifunc_noncod_num_exon <- exons_keep_plot[c(which(exons_keep_plot$gene_type=="Bifunctional"),
                                            which(exons_keep_plot$gene_type=="Noncoding")), ]
bifunc_noncod_num_exon$log2_num_exon <- log2(bifunc_noncod_num_exon$num)
cohens_d(log2_num_exon ~ gene_type, data = bifunc_noncod_num_exon)

cod_noncod_num_exon <- exons_keep_plot[c(which(exons_keep_plot$gene_type=="Coding"),
                                         which(exons_keep_plot$gene_type=="Noncoding")), ]
cod_noncod_num_exon$log2_num_exon <- log2(cod_noncod_num_exon$num)
cohens_d(log2_num_exon ~ gene_type, data = cod_noncod_num_exon)

exons_keep$transcript_type <- ""
exons_keep$transcript_type <- sapply(exons_keep$acc, function(x){
  if (length(grep("NR_", x))==1){
    return("noncoding_val")
  } else if (length(grep("NM_", x))==1){
    return("coding_val")
  } else if (length(grep("XR_", x))==1){
    return("noncoding_pred")
  } else if (length(grep("XM_", x))==1){
    return("coding_pred")
  }
})
plot_exons_trans<- ggplot(data = exons_keep,
                      aes(x = as.factor(transcript_type), y=num, fill=gene_type)) +
  geom_boxplot()+ labs(y="Number of exons", fill="Transcript Type", x="") +
  theme(axis.text = element_text(size=13, colour="black"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="right",
        legend.text = element_text(size=13, colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("Coding"="#e8e288", "Bifunctional" = "#82c0cc", "Noncoding" = "#ffa5a5"))
plot(plot_exons_trans)
plot_exons_trans2<- ggplot(data = exons_keep,
                          aes(x = as.factor(transcript_type), y=log2(num), fill=gene_type)) +
  geom_boxplot()+ labs(y="Log2(Number of exons)", fill="Transcript Type", x="") +
  theme(axis.text = element_text(size=13, colour="black"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="right",
        legend.text = element_text(size=13, colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("Coding"="#e8e288", "Bifunctional" = "#82c0cc", "Noncoding" = "#ffa5a5"))
plot(plot_exons_trans2)


print(sessionInfo())
