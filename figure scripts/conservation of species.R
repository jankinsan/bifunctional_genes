#conservation of bifunctional genes
BiocManager::install("orthogene")
library("orthogene")
setwd("C:/Users/DELL/OneDrive - IIT Delhi/phd work/NM NR Project/15122023")
#read dataset from humans!
gene_status <- read.csv("./counts/18122023_GRCh38.p14_genewise_counts_status.csv")

all_genes <- gene_status
rownames(all_genes) <- all_genes$X
#convert to mouse orthologs
gene_df <- orthogene::convert_orthologs(gene_df = all_genes,
                                        gene_input = "rownames", 
                                        gene_output = "dict", 
                                        input_species = "human",
                                        output_species = "mouse",
                                        non121_strategy = "kbs",
                                        method = "gprofiler") 
all_genes$mouse_orthologs <-NA
all_genes$mouse_orthologs <- gene_df[match(all_genes$X, names(gene_df))]

#convert to zebrafish orthologs
gene_df2 <- orthogene::convert_orthologs(gene_df = all_genes,
                                        gene_input = "rownames", 
                                        gene_output = "dict", 
                                        input_species = "human",
                                        output_species = "zebrafish",
                                        non121_strategy = "kbs",
                                        method = "gprofiler") 
all_genes$zfish_orthologs <-NA
all_genes$zfish_orthologs <- gene_df2[match(all_genes$X, names(gene_df2))]

#convert to chimpanzee orthologs
gene_df3 <- orthogene::convert_orthologs(gene_df = all_genes,
                                         gene_input = "rownames", 
                                         gene_output = "dict", 
                                         input_species = "human",
                                         output_species = "chimpanzee",
                                         non121_strategy = "kbs",
                                         method = "gprofiler") 
all_genes$chimp_orthologs <-NA
all_genes$chimp_orthologs <- gene_df3[match(all_genes$X, names(gene_df3))]

#convert to c. elegans ortholgs
gene_df4 <- orthogene::convert_orthologs(gene_df = all_genes,
                                         gene_input = "rownames", 
                                         gene_output = "dict", 
                                         input_species = "human",
                                         output_species = "celegans",
                                         non121_strategy = "kbs",
                                         method = "gprofiler") 
all_genes$celegans_orthologs <-NA
all_genes$celegans_orthologs <- gene_df4[match(all_genes$X, names(gene_df4))]

#read and match status
mouse<- read.csv("./counts/other species/2024-03-05_counts_mouse_genewise_status.csv")
zebrafish<- read.csv("./counts/other species/2024-03-05_counts_zebrafish_genewise_status.csv")
chimp<- read.csv("./counts/other species/2024-03-05_counts_chimpanzee_genewise_status.csv")
celegans <- read.csv("./counts/other species/2024-03-05_counts_zebrafish_genewise_status.csv")

#Match genes and add status 
all_genes$mouse_status <- NA
all_genes$mouse_status <-sapply(1:dim(all_genes)[1], function(x){
  gene <- all_genes$mouse_orthologs[x]
  if(!is.na(gene)){
   mat<-match(gene, mouse$X)
   if(!is.na(mat)){
   return(mouse$status[mat])
   } else { return ("unmatched")}
   } else { return(NA)}
})
all_genes$zfish_status <- NA
all_genes$zfish_status <-sapply(1:dim(all_genes)[1], function(x){
  gene <- all_genes$zfish_orthologs[x]
  if(!is.na(gene)){
    mat<-match(gene, zebrafish$X)
    if(!is.na(mat)){
      return(zebrafish$status[mat])
    } else { return ("unmatched")}
  } else { return(NA)}
})
all_genes$chimp_status <- NA
all_genes$chimp_status <-sapply(1:dim(all_genes)[1], function(x){
  gene <- all_genes$chimp_orthologs[x]
  if(!is.na(gene)){
    mat<-match(gene, chimp$X)
    if(!is.na(mat)){
      return(chimp$status[mat])
    } else { return ("unmatched")}
  } else { return(NA)}
})
all_genes$celegans_status <- NA
all_genes$celegans_status <-sapply(1:dim(all_genes)[1], function(x){
  gene <- all_genes$celegans_orthologs[x]
  if(!is.na(gene)){
    mat<-match(gene, celegans$X)
    if(!is.na(mat)){
      return(chimp$status[mat])
    } else { return ("unmatched")}
  } else { return(NA)}
})


#write to file
write.csv(all_genes, file=paste0("./counts/", Sys.Date(), "_counts_status_Across_species.csv"), row.names = FALSE)

#02-11-2025
#reading back the data to plot it!
all_genes <- read.csv("./counts/2025-04-29_counts_status_Across_species.csv")
#get genes which are hybrid in all genes!
all_hybrid <- all_genes$X[intersect(intersect(which(all_genes$status=="hybrid"), which(all_genes$mouse_status=="hybrid")), 
                        intersect(which(all_genes$zfish_status=="hybrid"), which(all_genes$chimp_status=="hybrid")))]
human_hybrid <- all_genes$X[intersect(intersect(which(all_genes$status=="hybrid"), which(all_genes$mouse_status!="hybrid")), 
                                    intersect(which(all_genes$zfish_status!="hybrid"), which(all_genes$chimp_status!="hybrid")))]
write.csv(all_hybrid, file="conserved_hybrid_genes.csv", col.names=NULL, row.names = FALSE)
write.csv(human_hybrid, file="only_human_hybrid_genes.csv", col.names=NULL, row.names = FALSE)

#upset diagram to represent this information!
library("UpSetR")
library("ComplexHeatmap")
#get the hybrid genes for each species!
hybrid_all <- unique(c(which(all_genes$status=="hybrid"),
         which(all_genes$mouse_status=="hybrid"),
         which(all_genes$zfish_status=="hybrid"),
         which(all_genes$chimp_status=="hybrid")))
bifunc_genes <- na.omit(all_genes) # omits genes which are NA's in any of the species!

plot_genes <- list(Human=bifunc_genes$X[bifunc_genes$status=="hybrid"], 
                  Mouse=bifunc_genes$X[bifunc_genes$mouse_status=="hybrid"],
                  Chimpanzee=bifunc_genes$X[bifunc_genes$chimp_status=="hybrid"],
                  Zebrafish=bifunc_genes$X[bifunc_genes$zfish_status=="hybrid"])
upset(fromList(plot_genes), mainbar.y.label="Number of bifunctional genes",
      point.size = 3.5, text.scale=1.5)

#not removing NA's
bifunc_genes2 <- all_genes[hybrid_all,] # does not omit genes which are NA's in any of the species!

plot_genes2 <- list(Human=bifunc_genes2$X[bifunc_genes2$status=="hybrid"], 
                   Mouse=bifunc_genes2$X[bifunc_genes2$mouse_status=="hybrid"],
                   Chimpanzee=bifunc_genes2$X[bifunc_genes2$chimp_status=="hybrid"],
                   Zebrafish=bifunc_genes2$X[bifunc_genes2$zfish_status=="hybrid"])
upset(fromList(plot_genes2),
      mainbar.y.label="Intersecting bifunctional genes",
      point.size = 5, line.size = 2.5, 
      sets.bar.color = "#92D050", main.bar.color = "coral",
      text.scale = c(2, 2, 1.5, 1.75,2,2.5),
      mb.ratio=c(0.6, 0.4), set_size.angles = 90)


#get conserved orthologs and whether or not they are bifunctional!
all_genes$conserv_type <- 0
all_genes$conserv_type <- sapply(1:dim(all_genes)[1], function(x){
  length(which(!is.na(all_genes[x, match(c("zfish_orthologs","chimp_orthologs","mouse_orthologs"),
  colnames(all_genes))])))+1
})
#conserv_type gives the number of species where the genes and orthologs are present!
#plot this!
library(ggplot2)
plot_dat <- data.frame(table(all_genes[,c(9,18)]))
plot_dat$gene_fraction <- 0
plot_dat$gene_fraction[which(plot_dat$status=="hybrid")]= plot_dat$Freq[which(plot_dat$status=="hybrid")]/sum(plot_dat$Freq[which(plot_dat$status=="hybrid")])
plot_dat$gene_fraction[which(plot_dat$status=="only_ncRNA")]= plot_dat$Freq[which(plot_dat$status=="only_ncRNA")]/sum(plot_dat$Freq[which(plot_dat$status=="only_ncRNA")])
plot_dat$gene_fraction[which(plot_dat$status=="only_mRNA")]= plot_dat$Freq[which(plot_dat$status=="only_mRNA")]/sum(plot_dat$Freq[which(plot_dat$status=="only_mRNA")])

#replace hybrid with bifunctional 
plot_dat$status <- gsub("hybrid", "Bifunctional genes", plot_dat$status)
plot_dat$status <- gsub("only_mRNA", "Coding genes", plot_dat$status)
plot_dat$status <- gsub("only_ncRNA", "Noncoding genes", plot_dat$status)

#PLOT FOR CONSERVATION OF BIFUNCTIONAL GENES!
ggplot(plot_dat, aes(x=as.factor(conserv_type), y=gene_fraction, fill=status))+ 
  geom_bar(stat = "identity", color="black", width=0.5)+ 
  labs(x="Gene Orthologs (in species)", y="Gene Fraction")+
  theme()+facet_wrap(~status, ncol=1)+
  scale_x_discrete(labels = c("Only humans", "Human + 1", "Human + 2", "Human + 3"))+
  theme(axis.text.x=element_text(angle = 90),
        plot.background = element_blank(),
        axis.text= element_text(size=13, color = 'black'),
        axis.title = element_text(size=13, color = 'black', face = "bold"),
        panel.border = element_rect( fill=NA), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "lightgrey"),
        panel.grid.minor = element_blank(),
        legend.position="none",
        legend.title=element_blank(),
        strip.text = element_text(size = 13))+
  scale_fill_manual(values = c("Coding genes"="#e8e288", "Bifunctional genes" = "#82c0cc", "Noncoding genes" = "#ffa5a5"))

#Get circular RNAs
circRNA_matched <- read.csv("F:/Janki/NM NR Project/15122023/circRNA_matched/2024-04-24_circRNA_bifunc_matched_all_genes.csv")

#match genes for each category
#first get bifunc genes only for human separately for this 
bifunc_dat <- all_genes[all_genes$status=="hybrid",]
bifunc_dat$circRNA_status <- 'no'
bifunc_dat$circRNA_num <- 0
bifunc_dat$circRNA_status[match(unique(circRNA_matched$gene), bifunc_dat$X)] <- 'yes'
bifunc_dat$circRNA_num<- sapply(bifunc_dat$X, function(gene){
  length(which(circRNA_matched$gene==gene))
})

#get summary of genes and plot!
circ_summ <- as.data.frame(table(bifunc_dat[,c(18,19)]))
circ_summ$percent <- c(circ_summ[1,3]/sum(circ_summ[c(1,5), 3]),
                       circ_summ[2,3]/sum(circ_summ[c(2,6), 3]),
                       circ_summ[3,3]/sum(circ_summ[c(3,7), 3]),
                       circ_summ[4,3]/sum(circ_summ[c(4,8), 3]),
                       circ_summ[5,3]/sum(circ_summ[c(1,5), 3]),
                       circ_summ[6,3]/sum(circ_summ[c(2,6), 3]),
                       circ_summ[7,3]/sum(circ_summ[c(3,7), 3]),
                       circ_summ[8,3]/sum(circ_summ[c(4,8), 3]))
circ_plot <- reshape2::melt(circ_summ[,c(1,2,4)])

#plot for the genes in each category!
library(dplyr)
circ_plot$conserv_type <- recode(circ_plot$conserv_type,
                    "1" = "Only humans",
                    "2" ="Human + 1",
                    "3"= "Human + 2",
                    "4" ="Human + 3")
ggplot(data=circ_plot, aes(x=" ", y=value, group=circRNA_status, fill=circRNA_status)) +
  geom_bar(width = 1, stat = "identity", color="black", alpha=0.8) +
  coord_polar("y", start=0) + labs (fill = "circRNA")+
  facet_wrap(.~ conserv_type, ncol=1) +
  theme(
    plot.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(color = 'black', face = "bold", hjust=0.5),
    panel.background = element_rect(fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 13, colour = "black", margin = margin(t = 5, b = 5, l = 2, r = 2)),
    legend.position="bottom",
    legend.title= element_text(size = 13, colour = "black"),
    legend.text = element_text(size=12, colour = "black"))

set_names = comb_name(m)
genes<- vector()
for (set in set_names){
  genes_comb <- extract_comb(m, set)
  genes<- c(genes, genes_comb)
}
