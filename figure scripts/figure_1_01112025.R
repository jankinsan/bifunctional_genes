##27-06-2024
#IITD HPC--------------------
#get genomic data and load it!
hg38_dat <- read.delim("GCF_000001405.40_GRCh38.p14_genomic.gtf", skip=5, header=FALSE)
hg38_dat$acc <- sapply(hg38_dat[,9], function(x){
  str_trans <- strsplit(x, split=";")[[1]][2]
  split = unlist(strsplit(str_trans, split="transcript_id "))
  return(split[length(split)])
})
hg38_dat$gene <- sapply(hg38_dat[,9], function(x){
  str_trans <- strsplit(x, split=";")[[1]][1]
  split = unlist(strsplit(str_trans, split="gene_id "))
  return(split[length(split)])
})

print(head(hg38_dat))
print(dim(hg38_dat))
#removing the last row since it's useless and coming from the gtf format!
hg38_dat <- hg38_dat[-4693812,]

#separating gene, transcript and exon data to plot!
hg38_genes<- hg38_dat[which(hg38_dat[,3]=="gene"),]
hg38_exons<- hg38_dat[which(hg38_dat[,3]=="exon"),]
hg38_transcripts<- hg38_dat[which(hg38_dat[,3]=="transcript"),]
colnames(hg38_genes) <- c("chr", "annotation", "type", "start", "end", "misc", "strand", "misc2", "details", "acc", "gene")
colnames(hg38_exons) <- c("chr", "annotation", "type", "start", "end", "misc", "strand", "misc2", "details", "acc", "gene")
colnames(hg38_transcripts) <- c("chr", "annotation", "type", "start", "end", "misc", "strand", "misc2", "details", "acc", "gene")

#get accessions from my counts document
#match them with the type of gene!
acc<- read.csv("18122023_GRCh38.p14_gene_accessions.csv")
hg38_counts <- read.csv("20062024__GRCh38.p14_genewise_status.csv")
acc$ids <- sapply(acc$accession, function(x){
  substr(x, 2, nchar(x))
})
acc$gene_type <- hg38_counts$status[match(acc$gene, hg38_counts[,1])]
print(head(acc))

hg38_genes$gene_type <- hg38_counts$status[match(hg38_genes$gene, hg38_counts$final_gene)]
hg38_exons$gene_type <- acc$gene_type[match(hg38_exons$acc, acc$ids)]
hg38_transcripts$gene_type <- acc$gene_type[match(hg38_transcripts$acc, acc$ids)]

print("Exons...")
print(head(hg38_exons))
print("Genes...")
print(head(hg38_genes))
print("Transcripts...")
print(head(hg38_transcripts))

#replaces NA for gene_type with "not_assigned"
hg38_genes$gene_type[is.na(hg38_genes$gene_type)] = "not_assigned"
hg38_exons$gene_type[is.na(hg38_exons$gene_type)] = "not_assigned"
hg38_transcripts$gene_type[is.na(hg38_transcripts$gene_type)] = "not_assigned"
print(paste0("Exons assigned: ", length(hg38_exons$acc[-which(hg38_exons$gene_type=="not_assigned")])))
print(paste0("Genes assigned: ", length(unique(hg38_genes$gene[-which(hg38_genes$gene_type=="not_assigned")]))))
print(paste0("Transcripts assigned: ", length(unique(hg38_transcripts$acc[-which(hg38_transcripts$gene_type=="not_assigned")]))))


#getting simplified chromosome names for matching!
unique(hg38_exons$chr)
chr_dat <- read.delim("GCF_000001405.40_GRCh38.p14_assembly_report.txt", skip=63)
hg38_genes$chr <- chr_dat$UCSC.style.name[match(hg38_genes$chr, chr_dat$RefSeq.Accn)]
hg38_exons$chr <- chr_dat$UCSC.style.name[match(hg38_exons$chr, chr_dat$RefSeq.Accn)]
hg38_transcripts$chr <- chr_dat$UCSC.style.name[match(hg38_transcripts$chr, chr_dat$RefSeq.Accn)]

write.csv(hg38_genes, file="GCF_000001405.40_GRCh38.p14_genes.csv", row.names = FALSE)
write.csv(hg38_exons, file="GCF_000001405.40_GRCh38.p14_exons.csv", row.names = FALSE)
write.csv(hg38_transcripts, file="GCF_000001405.40_GRCh38.p14_transcripts.csv", row.names = FALSE)

##MOVED THE OUTPUT FILES TO MY LAPTOP!---------------
## RUNNING THIS FROM RStudio
#loading libraries
#install.packages("PupillometryR")
library(ggplot2)
library(viridis)
library(PupillometryR)

#get genes, match status and plot!
setwd("E:/MSR BLY/THESIS WORK/NP NR Project/15122023")
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

#Figure 2A
#STACKED BAR PLOT OF GENE TYPE---------------------------
summ_plot <- reshape2::melt(gene_data_summary[,c(1,5,4)])
write.csv(gene_data_summary, file="figure_drafts/02072024_NCBI_Annotation_numbers_human.csv")

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

plot_stacked_type2 <- ggplot(data = summ_plot, 
                             aes(x = variable, y=value, 
                                 fill=gene_type, width=0.65)) +
  geom_bar(position="stack", stat="identity", colour="black")+
  labs(y="Fraction of Genes", x="", fill="Gene Type") +
  theme(axis.text = element_text(size=13, colour="black"),
        axis.title = element_text(size=13, colour="black", face = "bold"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="bottom",
        legend.text = element_text(size=13, colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())+
  #scale_fill_brewer(palette = "Spectral")+
  scale_x_discrete(labels=c("Genes with transcripts", "All Genes"))

dir.create("figure_drafts")
pdf("figure_drafts/23112024_figure2a_human.pdf")
plot(plot_stacked_type2)
plot(plot_stacked_type)
dev.off()

#Figure 2B
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

#Protein-coding versus non-protein coding genes
#STACKED/MOSAIC FOR NUMBERS OF THREE CATEGORIES:hybrid, mRNA and ncRNA
#colours: #e8e288, #82c0cc, #ffa5a5
stacked_plot<- data.frame(annot=c("Bifunctional", "Noncoding", "Coding"))
stacked_plot$annot = factor(stacked_plot$annot, levels = c("Bifunctional", "Noncoding", "Coding"))
stacked_plot$number<- c(length(which(hg38_genes2$gene_status=="Bifunctional")),
                        length(which(hg38_genes2$gene_status=="Noncoding")),
                        length(which(hg38_genes2$gene_status=="Coding")))
stacked_plot$labels_genes <-c(length(which(hg38_genes2$gene_status=="Bifunctional"))/length(hg38_genes2$gene_status),
                              length(which(hg38_genes2$gene_status=="Noncoding"))/length(hg38_genes2$gene_status),
                              length(which(hg38_genes2$gene_status=="Coding"))/length(hg38_genes2$gene_status))

#genes that are said to be protein-coding
prot_idx <- which(hg38_genes2$gene_biotype=="protein_coding")
stacked_plot$prot_number<- c(length(which(hg38_genes2$gene_status[prot_idx]=="Bifunctional")),
                             length(which(hg38_genes2$gene_status[prot_idx]=="Noncoding")),
                             length(which(hg38_genes2$gene_status[prot_idx]=="Coding")))
stacked_plot$prot_coding_labels <- c(length(which(hg38_genes2$gene_status[prot_idx]=="Bifunctional"))/length(hg38_genes2$gene_status[prot_idx]),
                                     length(which(hg38_genes2$gene_status[prot_idx]=="Noncoding"))/length(hg38_genes2$gene_status[prot_idx]),
                                     length(which(hg38_genes2$gene_status[prot_idx]=="Coding"))/length(hg38_genes2$gene_status[prot_idx]))

#all other genes
stacked_plot$non_prot_number<- c(length(which(hg38_genes2$gene_status[-prot_idx]=="Bifunctional")),
                                 length(which(hg38_genes2$gene_status[-prot_idx]=="Noncoding")),
                                 length(which(hg38_genes2$gene_status[-prot_idx]=="Coding")))
stacked_plot$non_prot_coding_labels <- c(length(which(hg38_genes2$gene_status[-prot_idx]=="Bifunctional"))/length(hg38_genes2$gene_status[-prot_idx]),
                                         length(which(hg38_genes2$gene_status[-prot_idx]=="Noncoding"))/length(hg38_genes2$gene_status[-prot_idx]),
                                         length(which(hg38_genes2$gene_status[-prot_idx]=="Coding"))/length(hg38_genes2$gene_status[-prot_idx]))
write.csv(stacked_plot, file="figure_drafts/23112024_New_Annotations_numbers_human.csv")

#plot
summ_plot2 <- reshape2::melt(stacked_plot[,c(1,3, 5, 7)])

#ggplot
plot_stacked_annot <- ggplot(data = summ_plot2, 
                             aes(x = variable, y=value, 
                                 fill=annot, width=0.65)) +
  geom_bar(position="stack", stat="identity")+
  labs(y="Fraction of Genes", x="NCBI Gene Biotype", fill="") +
  theme(axis.text.x = element_text(size=13, colour="black",angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=13, colour="black"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="bottom",
        legend.text = element_text(size=13, colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("Coding"="#e8e288", "Bifunctional" = "#82c0cc", "Noncoding" = "#ffa5a5"))+
  scale_x_discrete(labels=c("All Genes", "Protein-Coding Genes", "Non-Protein Coding Genes"))

plot_stacked_annot2 <- ggplot(data = summ_plot2, 
                              aes(y = variable, x=value, 
                                  fill=annot, width=0.85)) +
  geom_bar(position="stack", stat="identity")+
  labs(x="Fraction of Genes", y="NCBI Gene Biotype", fill="") +
  theme(axis.text = element_text(size=13, colour="black"),
        axis.title = element_text(size=13, colour="black", face="bold"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="right",
        legend.text = element_text(size=13, colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("Coding"="#e8e288", "Bifunctional" = "#82c0cc", "Noncoding" = "#ffa5a5"))+
  scale_y_discrete(labels=c("All Genes", "Protein-Coding Genes", "Non-Protein Coding Genes"))

pdf("figure_drafts/23112024_figure2b.pdf")
plot(allu_plot)
plot(plot_stacked_annot)
plot(plot_stacked_annot2)
dev.off()

##Figure 2C 
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

pdf("./figure_drafts/23112024_human_chr_dist.pdf")
plot(chr_plot)
plot(chr_plot2)
plot(chr_plot_strand)
plot(chr_plot_strand2)
dev.off()

#save data
write.csv(chr_dist, file="./figure_drafts/23112024_New_annotation_chr_dist.csv", row.names = FALSE)

##Figure 2E
#Transcript Numbers per gene
#unusable dataframe
trans_numbers <- data.frame(status=c("Coding", "Bifunctional", "Noncoding"),
                            num_total = c(sum(gene_status$total[gene_status$status=="mRNA encoding"]),
                                          sum(gene_status$total[gene_status$status=="bifunctional"]),
                                          sum(gene_status$total[gene_status$status=="ncRNA enncoding"])),
                            num_coding = c(sum(gene_status$NM_XM[gene_status$status=="mRNA encoding"]),
                                           sum(gene_status$NM_XM[gene_status$status=="bifunctional"]),
                                           sum(gene_status$NM_XM[gene_status$status=="ncRNA encoding"])),
                            num_noncoding = c(sum(gene_status$NR_XR[gene_status$status=="mRNA encoding"]),
                                           sum(gene_status$NR_XR[gene_status$status=="bifunctional"]),
                                           sum(gene_status$NR_XR[gene_status$status=="ncRNA encoding"])))

#p-value to calculate between groups
library(ggpubr)
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

##Figure 2F
#Gene lengths
gene_len_plot <- ggplot(data = hg38_genes2,
                         aes(x = gene_status, y=length, fill=gene_status)) +
  geom_boxplot()+
  labs(y="Gene Lengths (nucleotides)", fill="Gene Annotation", x="") +
  theme(axis.text = element_text(size=13, colour="black", face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="right",
        legend.title = element_text(size=12, colour = "black", face="bold"),
        legend.text = element_text(size=11, colour = "black", face="bold"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())+
  stat_compare_means(comparisons = my_comparisons)+
  scale_fill_manual(values = c("Coding"="#e8e288", "Bifunctional" = "#82c0cc", "Noncoding" = "#ffa5a5"))+
  scale_color_manual(values = c("Coding"="#e8e288", "Bifunctional" = "#82c0cc", "Noncoding" = "#ffa5a5"))
plot(gene_len_plot)

#Gene lengths
gene_len_plot2 <- ggplot(data = hg38_genes2,
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
plot(gene_len_plot2)

##Figure 2G
#Transcript lengths
hg38_transcripts_matched$gene_type <- recode(hg38_transcripts_matched$gene_type,
                             "ncRNA encoding" = "Noncoding",
                             "mRNA encoding" = "Coding",
                             "bifunctional" = "Bifunctional")
trans_len_plot <- ggplot(data = hg38_transcripts_matched,
                        aes(x = gene_type, y=length, fill=gene_type)) +
  geom_boxplot()+
  labs(y="Transcript span (nucleotides)", fill="Gene Annotation", x="") +
  theme(axis.text = element_text(size=13, colour="black", face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="right",
        legend.title = element_text(size=12, colour = "black", face="bold"),
        legend.text = element_text(size=11, colour = "black", face="bold"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("Coding"="#e8e288", "Bifunctional" = "#82c0cc", "Noncoding" = "#ffa5a5"))+
  scale_color_manual(values = c("Coding"="#e8e288", "Bifunctional" = "#82c0cc", "Noncoding" = "#ffa5a5"))
  plot(trans_len_plot)

#Transcript lengths
trans_len_plot2 <- ggplot(data = hg38_transcripts_matched,
                           aes(x = gene_type, y=log2(length), fill=gene_type)) +
    geom_boxplot()+ labs(y="Log2(Transcript span)", fill="Gene Annotation", x="") +
  theme(axis.text = element_text(size=13, colour="black"),
        axis.title.y = element_text(size=13, face="bold"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="right",
        legend.title = element_text(size=12, colour = "black", face="bold"),
        legend.text = element_text(size=11, colour = "black", face="bold"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())+ 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "Bifunctional", size=6, label.y = 21.5)+
  scale_fill_manual(values = c("Coding"="#e8e288", "Bifunctional" = "#82c0cc", "Noncoding" = "#ffa5a5"))
plot(trans_len_plot2)

#plots
ggplot(hg38_transcripts_matched, aes(x=length, y=gene_type, colour=gene_type))+
  geom_jitter()+theme_bw()+
  scale_fill_manual(values = c("Coding"="#e8e288", "Bifunctional" = "#82c0cc", "Noncoding" = "#ffa5a5"))
ggplot(hg38_transcripts_matched, aes(x=length, y=gene_type, colour=gene_type))+
  geom_boxplot()+theme_bw()

#################################################################################
#Number of exons, average number of exons per transcript
#IGIB HPC!: transferred files from laptop to lustre folder

hg38_exons <- read.csv("./15122023_hg38/GCF_000001405.40_GRCh38.p14_exons.csv")
head(hg38_exons)

number_exons <- data.frame(acc=unique(hg38_exons$acc))
number_exons$num <- sapply(1:dim(number_exons)[1], function(idx){
  print(idx)
  x = number_exons$acc[idx]
  length(which(hg38_exons$acc==x))
})
number_exons$gene <- hg38_exons$gene[match(number_exons$acc, hg38_exons$acc)]
number_exons$gene_type<- hg38_exons$gene_type[match(number_exons$acc, hg38_exons$acc)]

gene_status <- read.csv("20062024__GRCh38.p14_genewise_status.csv")
number_exons$gene_status <- gene_status$status[match(number_exons$gene, gene_status$final_gene)]
head(number_exons)
write.csv(number_exons, file="10072024_number_exons_per_transcript.csv", row.names=FALSE)

###############################################################################################
#LAPTOP:plotting
number_exons<- read.csv("figure_drafts/10072024_number_exons_per_transcript.csv")
exons_keep <- number_exons[-which(number_exons$gene_type=="not_assigned"),]
#update annotation
exons_keep$gene_type[exons_keep$gene_type=="only_ncRNA"]<- "Noncoding"
exons_keep$gene_type[exons_keep$gene_type=="only_mRNA"]<- "Coding"
exons_keep$gene_type[exons_keep$gene_type=="hybrid"]<- "Bifunctional"

plot_exons <- ggplot(data = exons_keep,
                     aes(x = as.factor(gene_type), y=num, fill=gene_type)) +
  geom_boxplot()+ labs(y="Number of exons", fill="Gene Annotation", x="") +
  theme(axis.text = element_text(size=13, colour="black"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="right",
        legend.text = element_text(size=13, colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())+
  stat_compare_means(comparisons = my_comparisons)+
  scale_fill_manual(values = c("Coding"="#e8e288", "Bifunctional" = "#82c0cc", "Noncoding" = "#ffa5a5"))
plot(plot_exons)
plot_exons2 <- ggplot(data = exons_keep,
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
plot(plot_exons2)

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

##############################################################################################################
#NEW system: add a supplementary figure for number and sitribution of noncoding transcripts 
library(ggplot2)
library(UpSetR)
library(cowplot)

#read transcripts counts and status as per refseq
counts_hg38 <- read.csv("F:/Janki/bifunctional_genes/data/counts/23112024_GRCh38.p14_genewise_status.csv")
counts_hg38$status<- sub("mRNA encoding", "coding", counts_hg38$status)
counts_hg38$status<- sub("ncRNA encoding", "noncoding", counts_hg38$status)

bifunc <- counts_hg38[counts_hg38$status=="bifunctional", ]
bifunc$percent_noncoding_transcripts <- bifunc$NR_XR/bifunc$total


ggplot(bifunc, aes(x=percent_noncoding_transcripts))+
  #geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, color="blue", fill="blue")+
  geom_vline(aes(xintercept=0.5),
             color="red", linetype="dashed", size=1.25)+ 
  labs(x="Ratio of Noncoding Transcripts", y= "Density of bifunctional genes")+
  theme(axis.text = element_text(size=13, colour="black"),
        axis.title = element_text(size=13, face="bold"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="right",
        legend.title = element_text(size=13, colour = "black", face="bold"),
        legend.text = element_text(size=12, colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())

#upset plots for hg38 vs hs1
counts_hs1<- read.csv("F:/Janki/bifunctional_genes/data/counts/18122023_T2T-CHM13v2.0_genewise_counts_status.csv")
counts_hs1$status<- sub("hybrid", "bifunctional", counts_hs1$status)
counts_hs1$status<- sub("only_mRNA", "coding", counts_hs1$status)
counts_hs1$status<- sub("only_ncRNA", "noncoding", counts_hs1$status)

#get lists to compare!
plot_all <- list(GRCh38.p14= counts_hg38$X, 
                    T2T_CHM13v2.0= counts_hs1$X)

plot_bifunc <- list(GRCh38.p14= counts_hg38$X[counts_hg38$status=="bifunctional"], 
                    T2T_CHM13v2.0= counts_hs1$X[counts_hs1$status=="bifunctional"])


plot_coding <- list(GRCh38.p14= counts_hg38$X[counts_hg38$status=="coding"], 
                    T2T_CHM13v2.0= counts_hs1$X[counts_hs1$status=="coding"])

plot_noncoding <- list(GRCh38.p14= counts_hg38$X[counts_hg38$status=="noncoding"], 
                    T2T_CHM13v2.0= counts_hs1$X[counts_hs1$status=="noncoding"])

#make upset plots
upset(fromList(plot_all),
      mainbar.y.label="Intersecting genes (all genes)",
      point.size = 3.5, line.size = 1.5, 
      sets.bar.color = "darkblue", main.bar.color = "darkred",
      text.scale = c(1.75, 1.75, 1.5, 1.25,2,2.5),
      mb.ratio=c(0.6, 0.4), set_size.angles = 90)

upset(fromList(plot_bifunc),
      mainbar.y.label="Intersecting bifunctional genes",
      point.size = 3.5, line.size = 1.5, 
      sets.bar.color = "darkblue", main.bar.color = "darkred",
      text.scale = c(1.75, 1.75, 1.5, 1.2,2,2.5),
      mb.ratio=c(0.6, 0.4), set_size.angles = 90)

upset(fromList(plot_coding),
      mainbar.y.label="Intersecting coding genes",
      point.size = 3.5, line.size = 1.5, 
      sets.bar.color = "darkblue", main.bar.color = "darkred",
      text.scale = c(1.75, 1.75, 1.5, 1.2,2,2.5),
      mb.ratio=c(0.6, 0.4), set_size.angles = 90)

upset(fromList(plot_noncoding),
      mainbar.y.label="Intersecting noncoding genes",
      point.size = 3.5, line.size = 1.5, 
      sets.bar.color = "darkblue", main.bar.color = "darkred",
      text.scale = c(1.75, 1.75, 1.5, 1.2,2,2.5),
      mb.ratio=c(0.6, 0.4), set_size.angles = 90)

#TRASH STUFF!
na.omit(hg38_genes[,12:13])
plot_dat$gene_status <- as.factor(plot_dat$gene_status)
ggplot(plot_dat, aes(x=gene_status, y=log10(length), fill=gene_status))+geom_boxplot()

length(na.omit(hg38_genes$gene_status))

hg38_transcripts <- read.csv("GCF_000001405.40_GRCh38.p14_transcripts.csv")

#remove unassigned transcripts!
hg38_trans <- hg38_transcripts[-grep("unassigned", hg38_transcripts$acc),]

