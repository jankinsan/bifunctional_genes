library("ggplot2")
library("ggalluvial")
library("ggrepel")
library("data.table")
library("dplyr")
library("rstatix")

setwd("F:/Janki/bifunctional_genes/data")

#read counts_updated
gene_status<- read.csv("./counts/23112024_GRCh38.p14_genewise_status.csv")
#add validated gene status for only genes with validated transcripts!
gene_status$val_status<- "No validated transcript"
gene_status$val_status[gene_status$NR>0]<- "noncoding"
gene_status$val_status[gene_status$NM>0]<- "coding"
gene_status$val_status[intersect(which(gene_status$NM>0),
                                     which(gene_status$NR>0))]<- "bifunctional"

#reding previously separated gene, transcript and exon data to plot!
hg38_genes<- read.csv("GCF_000001405.40_GRCh38.p14_genes.csv")
hg38_exons<- read.csv("GCF_000001405.40_GRCh38.p14_exons.csv")
hg38_transcripts<- read.csv("GCF_000001405.40_GRCh38.p14_transcripts.csv")

#keep only genes and exons 
hg38_genes$gene_type <- gene_status$status[match(hg38_genes$gene, gene_status$final_gene)]
hg38_exons$gene_type <- gene_status$status[match(hg38_exons$gene, gene_status$final_gene)]
hg38_transcripts$gene_type <- gene_status$status[match(hg38_transcripts$gene, gene_status$final_gene)]
#add status based on validated transcripts!
hg38_genes$val_gene_type <- gene_status$val_status[match(hg38_genes$gene, gene_status$final_gene)]
hg38_exons$val_gene_type <- gene_status$val_status[match(hg38_exons$gene, gene_status$final_gene)]
hg38_transcripts$val_gene_type <- gene_status$val_status[match(hg38_transcripts$gene, gene_status$final_gene)]

#replaces NA for gene_type with "not_assigned"
hg38_genes$gene_type[is.na(hg38_genes$gene_type)] = "not_assigned"
hg38_exons$gene_type[is.na(hg38_exons$gene_type)] = "not_assigned"
hg38_transcripts$gene_type[is.na(hg38_transcripts$gene_type)] = "not_assigned"

hg38_genes$val_gene_type[is.na(hg38_genes$val_gene_type)] = "not_assigned"
hg38_exons$val_gene_type[is.na(hg38_exons$val_gene_type)] = "not_assigned"
hg38_transcripts$val_gene_type[is.na(hg38_transcripts$val_gene_type)] = "not_assigned"

print(paste0("Exons assigned: ", length(hg38_exons$acc[-which(hg38_exons$gene_type=="not_assigned")])))
print(paste0("Genes assigned: ", length(unique(hg38_genes$gene[-which(hg38_genes$gene_type=="not_assigned")]))))
print(paste0("Transcripts assigned: ", length(unique(hg38_transcripts$acc[-which(hg38_transcripts$gene_type=="not_assigned")]))))

#getting the NCBI biotype
hg38_genes$gene_biotype <- sapply(1:dim(hg38_genes)[1], function(x){
  strsplit(strsplit(hg38_genes[x,9], split="gene_biotype ")[[1]][2], split=";")[[1]][1]
})

#keep only genes with a transcript!
hg38_genes2<- hg38_genes[hg38_genes$gene_type!="not_assigned",]
hg38_genes2$gene_biotype2 <-hg38_genes2$gene_biotype
idx <- unique(c(which(hg38_genes2$gene_biotype=="protein_coding"), 
                which(hg38_genes2$gene_biotype=="lncRNA"),
                which(hg38_genes2$gene_biotype=="miRNA"),
                which(hg38_genes2$gene_biotype=="transcribed_pseudogene"),
                which(hg38_genes2$gene_biotype=="snoRNA")))
hg38_genes2$gene_biotype2[-idx] <- "others"


combs <- CJ(hg38_genes2$gene_biotype2, hg38_genes2$gene_typ, hg38_genes2$val_gene_type, unique=TRUE)
allu_dat <- data.frame(combs)

allu_dat$freq <- sapply(1:dim(combs)[1], function(x){
  length(intersect(which(hg38_genes2$val_gene_type==allu_dat[x,3]),
    intersect(which(hg38_genes2$gene_biotype2==allu_dat[x,1]),
                   which(hg38_genes2$gene_type==allu_dat[x,2]))))
}) 
colnames(allu_dat) <- c("ncbi", "new_annot", "val_annot", "freq")

#alluvial plot to look at how different genes are annotated according to NCBI and our annotation!
allu_plot <- ggplot(data = allu_dat,
                    aes(axis1 = ncbi, axis2 = new_annot, axis3 = val_annot, y = freq)) +
  scale_x_discrete(limits = c("NCBI Biotype", "Our Annotation", "Validated only"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = val_annot), width=1/8, curve_type = "sine") + 
  geom_stratum(aes(fill = val_annot),alpha = 1, width = 1/8, reverse = TRUE) + 
  labs(y="Number of Genes", fill="Gene Type") +
  geom_label_repel(stat = "stratum", aes(label = after_stat(stratum),size=11)) + guides (size="none")+
  scale_fill_manual(values = c("coding"="#e8e288", "bifunctional" = "#82c0cc", "noncoding" = "#ffa5a5", "No validated transcript"="darkolivegreen3"))+
  theme_bw() + theme(legend.title=element_text(size=11,face="bold"),
                     legend.text=element_text(size=11),
                     legend.position="bottom",
                     axis.title = element_text(size=11, face="bold"),
                     axis.text.x = element_text(size=11,face="bold", color="black"),
                     axis.text.y = element_text(size=11, color="black"),
                     legend.key.size = unit(1, "lines"))
plot(allu_plot)


#make a subset with only genes with at least one validated transcript
validated_genes <- gene_status[(gene_status$NR + gene_status$NM)>0,]
validated_genes$val_total <- validated_genes$NM + validated_genes$NR
trans_numbers <- data.frame(val_status=c("coding", "bifunctional", "noncoding"),
                            num_total = c(sum(validated_genes$val_total[validated_genes$val_status=="coding"]),
                                          sum(validated_genes$val_total[validated_genes$val_status=="bifunctional"]),
                                          sum(validated_genes$val_total[validated_genes$val_status=="noncoding"])),
                            num_coding = c(sum(validated_genes$NM[validated_genes$val_status=="coding"]),
                                           sum(validated_genes$NM[validated_genes$val_status=="bifunctional"]),
                                           sum(validated_genes$NM[validated_genes$val_status=="noncoding"])),
                            num_noncoding = c(sum(validated_genes$NR[validated_genes$val_status=="coding"]),
                                              sum(validated_genes$NR[validated_genes$val_status=="bifunctional"]),
                                              sum(validated_genes$NR[validated_genes$val_status=="noncoding"])))

#p-value to calculate between groups
library("ggpubr")
my_comparisons <- list( c("bifunctional", "coding"), c("coding", "noncoding"),
                        c("bifunctional", "noncoding") )

#plot
trans_num_plot <- ggplot(data = validated_genes,
                         aes(x = val_status,y=log2(val_total), fill=val_status)) +
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
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "bifunctional", size=6, label.y = 9.5)+
  scale_fill_manual(values = c("coding"="#e8e288", "bifunctional" = "#82c0cc", "noncoding" = "#ffa5a5"))

plot(trans_num_plot)

#calculate effect sizes for number of transcripts
##bifunctional versus other groups
bifunc_cod_dat <- validated_genes[c(which(validated_genes$val_status=="bifunctional"),
                                which(validated_genes$val_status=="coding")), ]
bifunc_cod_dat$log2_total <- log2(bifunc_cod_dat$total)
cohens_d(log2_total ~ val_status, data = bifunc_cod_dat)

bifunc_noncod_dat <- validated_genes[c(which(validated_genes$val_status=="bifunctional"),
                                   which(validated_genes$val_status=="noncoding")), ]
bifunc_noncod_dat$log2_total <- log2(bifunc_noncod_dat$total)
cohens_d(log2_total ~ val_status, data = bifunc_noncod_dat)

cod_noncod_dat <- validated_genes[c(which(validated_genes$val_status=="coding"),
                                which(validated_genes$val_status=="noncoding")), ]
cod_noncod_dat$log2_total <- log2(cod_noncod_dat$total)
cohens_d(log2_total ~ val_status, data = cod_noncod_dat)

##Figure 2F
#Gene lengths
hg38_genes3 <- hg38_genes2[!hg38_genes2$val_gene_type=="No validated transcript",]
hg38_genes3$length <- hg38_genes3$end - hg38_genes3$start

gene_len_plot <- ggplot(data = hg38_genes3,
                        aes(x = val_gene_type, y=log2(length), fill=val_gene_type)) +
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
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "bifunctional", size=6, label.y = 21.5)+
  scale_fill_manual(values = c("coding"="#e8e288", "bifunctional" = "#82c0cc", "noncoding" = "#ffa5a5"))
plot(gene_len_plot)

#calculate effect sizes for gene lengths 
##bifunctional versus other groups
bifunc_cod_len <- hg38_genes3[c(which(hg38_genes3$val_gene_type=="bifunctional"),
                                which(hg38_genes3$val_gene_type=="coding")), ]
bifunc_cod_len$log2_length <- log2(bifunc_cod_len$length)
cohens_d(log2_length ~ val_gene_type, data = bifunc_cod_len)

bifunc_noncod_len <- hg38_genes3[c(which(hg38_genes3$val_gene_type=="bifunctional"),
                                   which(hg38_genes3$val_gene_type=="noncoding")), ]
bifunc_noncod_len$log2_length <- log2(bifunc_noncod_len$length)
cohens_d(log2_length ~ val_gene_type, data = bifunc_noncod_len)

cod_noncod_len <- hg38_genes3[c(which(hg38_genes3$val_gene_type=="coding"),
                                which(hg38_genes3$val_gene_type=="noncoding")), ]
cod_noncod_len$log2_length <- log2(cod_noncod_len$length)
cohens_d(log2_length ~ val_gene_type, data = cod_noncod_len)

#Exon numbers
number_exons <- data.frame(acc=unique(hg38_exons$acc))
number_exons$num <- sapply(1:dim(number_exons)[1], function(idx){
  print(idx)
  x = number_exons$acc[idx]
  length(which(hg38_exons$acc==x))
})
number_exons$gene <- hg38_exons$gene[match(number_exons$acc, hg38_exons$acc)]
number_exons$val_gene_type<- hg38_exons$val_gene_type[match(number_exons$acc, hg38_exons$acc)]
write.csv(number_exons, file="F:/Janki/bifunctional_genes/revision/22052026_exon_numbers_validated_only.csv")

#read data on later date!
number_exons <- read.csv(file="E:/Janki/bifunctional_genes/revision/22052026_exon_numbers_validated_only.csv")
exons_keep <- number_exons[-c(which(number_exons$val_gene_type=="No validated transcript"),
                              which(number_exons$val_gene_type=="not_assigned")),]
#keeping only the transcript with the highest number of exons for each gene
exons_keep_plot <- na.omit(exons_keep %>% group_by(gene) %>%
                             slice_max(num, n = 1, with_ties = FALSE)  %>% ungroup())


plot_exons <- ggplot(data = exons_keep_plot,
                      aes(x = as.factor(val_gene_type), y=log2(num), fill=val_gene_type)) +
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
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "bifunctional", size=6, label.y = 8.75)+
  scale_fill_manual(values = c("coding"="#e8e288", "bifunctional" = "#82c0cc", "noncoding" = "#ffa5a5"))
plot(plot_exons)

#calculate effect sizes
##bifunctional versus other groups
bifunc_cod_num_exon <- exons_keep_plot[c(which(exons_keep_plot$val_gene_type=="bifunctional"),
                                         which(exons_keep_plot$val_gene_type=="coding")), ]
bifunc_cod_num_exon$log2_num_exon <- log2(bifunc_cod_num_exon$num)
cohens_d(log2_num_exon ~ val_gene_type, data = bifunc_cod_num_exon)

bifunc_noncod_num_exon <- exons_keep_plot[c(which(exons_keep_plot$val_gene_type=="bifunctional"),
                                            which(exons_keep_plot$val_gene_type=="noncoding")), ]
bifunc_noncod_num_exon$log2_num_exon <- log2(bifunc_noncod_num_exon$num)
cohens_d(log2_num_exon ~ val_gene_type, data = bifunc_noncod_num_exon)

cod_noncod_num_exon <- exons_keep_plot[c(which(exons_keep_plot$val_gene_type=="coding"),
                                         which(exons_keep_plot$val_gene_type=="noncoding")), ]
cod_noncod_num_exon$log2_num_exon <- log2(cod_noncod_num_exon$num)
cohens_d(log2_num_exon ~ val_gene_type, data = cod_noncod_num_exon)

#FIGURE S1E & s1F
###CHROMOSOMAL DISTRIBUTION OF ALL KINDS OF GENES
#remove the genes with both X and Y chromosomes 
#adding these numbers to either of these chromosomes
##strandedness of bifunctional genes chromosome-wise!
hg38_genes3$chr_simplified <- sapply(hg38_genes3$chr, function(x){
  strsplit(x, split="_")[[1]][1] 
})
chr_dist <- data.frame(chr=unique(hg38_genes3$chr_simplified))
chr_dist$all_genes <- sapply(chr_dist$chr, function(x){
  length(which(hg38_genes3$chr_simplified==x))
})
chr_dist$Coding <- sapply(chr_dist$chr, function(x){
  idx <- which(hg38_genes3$val_gene_type=="coding")
  length(which(hg38_genes3$chr_simplified[idx]==x))
})
chr_dist$Bifunctional <- sapply(chr_dist$chr, function(x){
  idx <- which(hg38_genes3$val_gene_type=="bifunctional")
  length(which(hg38_genes3$chr_simplified[idx]==x))
})
chr_dist$Noncoding <- sapply(chr_dist$chr, function(x){
  idx <- which(hg38_genes3$val_gene_type=="noncoding")
  length(which(hg38_genes3$chr_simplified[idx]==x))
})

chr_dist$percent_mrna <- chr_dist$Coding/chr_dist$all_genes
chr_dist$percent_ncrna <- chr_dist$Noncoding/chr_dist$all_genes
chr_dist$percent_bifunctional <- chr_dist$Bifunctional/chr_dist$all_gene
hg38_bifunc <- hg38_genes3[hg38_genes3$val_gene_type=="bifunctional",]
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


plot(chr_plot)
plot(chr_plot2)
plot(chr_plot_strand)
plot(chr_plot_strand2)

print(sessionInfo())