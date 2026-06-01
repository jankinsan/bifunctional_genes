library(ggplot2)
library(ComplexHeatmap)
library(magick)
library(circlize)
library(grid)
library(dplyr)
library(tidyr)
library(UpSetR)

combined_medians<- read.csv("E:/Janki/bifunctional_genes/data/GTEx exon junctions/GTEX_COMBINED_TISSUE_MEDIANS_UNCF.csv")

#remove the row numbers
plot_combined <- as.matrix(combined_medians[,-c(1:2)])
row.names(plot_combined) <- as.character(combined_medians[,2])

#replace NAs with zero 
plot_combined[is.na(plot_combined)] <- 0

#change rownames to remove the gtex label!
colnames(plot_combined)<- sapply(colnames(plot_combined), function(x){
  strsplit(x, split="gtex_")[[1]][2]
})

setwd("E:/Janki/bifunctional_genes/data/GTEx exon junctions/")

#make heatmap of 2609 bifunctional genes with non-zero median noncoding read fraction
col_fun= colorRamp2(c(0, 0.05, 0.5, 1), c("white", "#414487FF", "#22A884FF", "#E3e418ff"))
png("figure_exp_heatmap.png", width=450, height=1250, units="px")
Heatmap(plot_combined[,-28], name="Median NC", na_col = "white", use_raster = TRUE, 
        raster_quality = 5, col=col_fun, show_row_names = FALSE)
dev.off()

Heatmap(plot_combined[,-28], name="UNCF Median", na_col = "white", use_raster = TRUE, 
        raster_quality = 5, show_row_names = FALSE)

#get genes with noncoding fraction> 0.5 in at least one tissue!
high_nc <- plot_combined[rowSums(plot_combined>0.5)>= 1, -28]


#get genes with noncoding fraction > 0.5 in at least 30 tissues!
plot_top <- plot_combined[rowSums(plot_combined>0.5)>= 30, -28]
col_fun2= colorRamp2(c(0, 0.5, 0.75, 1), c("white", "#414487FF", "#22A884FF", "#E3e418ff"))
#svg("heatmap_top_genes.svg", width=700, height=550)
Heatmap(plot_top, name="Median NC", col=col_fun2, na_col = "white", column_names_gp = gpar(fontsize = 4),
        row_names_gp = gpar(fontsize = 4), legend())
#dev.off()

#plot for conserved genes and genes probed in our study!
genes_plot <- match(c("XIAP", "SNX16", "TIA1", "FOXP1", "TGFBR1", "USP15", "FUBP1", "ACVR1B", "DDX3Y"), rownames(plot_combined))
genes_plot <- match(c("CLASRP", "HNRNPDL", "KHDRBS1", "KHDRBS2", "NOP56", "PLEKHA5", "RPS9", "SNX25", 
                      "SRSF1", "SRSF2", "TASP1", "UBP1", "VCPKMT"), rownames(plot_combined))
Heatmap(na.omit(plot_combined[genes_plot,-28]), name="Median NC", col=col_fun2)

#get gene-wise expression! How many tissues is noncoding fraction expressed in?
num_tissues <- data.frame(genes=unique(rownames(plot_combined)))
num_tissues$five_percent <- sapply(num_tissues$genes, function(x){
  length(which(plot_combined[x,]>0.05))
})
num_tissues$twentyfive_percent <- sapply(num_tissues$genes, function(x){
  length(which(plot_combined[x,]>0.25))
})
rank_genes <- data.frame(freq=c("1-5", "6-10", "11-15", "16-20", "21-25", "25+"),
                         num_genes_five=0)
rank_genes$num_genes_five <- c(length(which(num_tissues$five_percent<=5)), 
                               length(intersect(which(num_tissues$five_percent>5),
                                                which(num_tissues$five_percent<=10))), 
                               length(intersect(which(num_tissues$five_percent>10),
                                                which(num_tissues$five_percent<=15))),
                               length(intersect(which(num_tissues$five_percent>15),
                                                which(num_tissues$five_percent<=20))),
                               length(intersect(which(num_tissues$five_percent>20),
                                                which(num_tissues$five_percent<=25))),
                               length(which(num_tissues$five_percent>25)))
rank_genes$num_genes_twentyfive <- c(length(which(num_tissues$twentyfive_percent<=5)), 
                                     length(intersect(which(num_tissues$twentyfive_percent>5),
                                                      which(num_tissues$twentyfive_percent<=10))), 
                                     length(intersect(which(num_tissues$twentyfive_percent>10),
                                                      which(num_tissues$twentyfive_percent<=15))),
                                     length(intersect(which(num_tissues$twentyfive_percent>15),
                                                      which(num_tissues$twentyfive_percent<=20))),
                                     length(intersect(which(num_tissues$twentyfive_percent>20),
                                                      which(num_tissues$twentyfive_percent<=25))),
                                     length(which(num_tissues$twentyfive_percent>25)))
plot_num <- reshape2::melt(rank_genes)

#change variables
plot_num$variable <- gsub("num_genes_five", ">5% non-coding reads", plot_num$variable)
plot_num$variable <- gsub("num_genes_twentyfive", ">25% non-coding reads", plot_num$variable)

# barplot showcasing genes with median noncoding fraction > 5% or >25% and how many tissues these are expressed in!
plot_num$freq <- factor(plot_num$freq, 
                        levels = unique(plot_num$freq))
ggplot(plot_num, aes(x=freq, y=value, fill=variable))+
  geom_bar(position="dodge", stat="identity", alpha=0.75)+ ylim(0,2000)+
  labs(x= "Number of tissues noncoding reads were detected in", y="Numbers of genes")+
  theme(plot.background = element_blank(),
        axis.text= element_text(size=13, color = 'black'),
        axis.title = element_text(size=13, color = 'black', face = "bold"),
        panel.border = element_rect( fill=NA), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "lightgrey"),
        panel.grid.minor = element_blank(),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=13),
        strip.text = element_text(size = 13))+
  geom_text(aes(label = value), 
            position=position_dodge(width=1.0), 
            vjust = -0.5)+
  scale_fill_manual(values=c(
    ">5% non-coding reads"= "#414487FF", 
      ">25% non-coding reads"="#E3e418ff"))

#for all genes, calculating the tau specificity score for noncoding splice variants
tau_genes <- data.frame(genes=combined_medians$gene, 
                        tau=0)
tau_genes$tau <- sapply(combined_medians$gene, function(x){
  # get all median values for a gene 
  row_data <- combined_medians[combined_medians$gene == x, -c(1,2)]
  val <- as.numeric(unlist(row_data))
  
  # Convert NAs to 0 
  val[is.na(val)] <- 0
  
  # give zeroes to genes that have 0 noncoding fraction across tissues!
  max_val <- max(val)
  if (max_val == 0) {
    return(0) 
  }
  
  # Compute tau
  tau_score <- sum(1 - (val / max_val)) / 30
  return(tau_score)
})

#get the number of tissues this gene is expressed in >25% noncoding fraction
tau_genes$num_tissues <- num_tissues$twentyfive_percent[match(tau_genes$genes,num_tissues$genes)]

#calculate global medians for each tissue for each gene!
tau_genes$global_median_nc <- sapply(combined_medians$gene, function(x){
  median(na.omit(as.numeric(combined_medians[which(combined_medians$gene==x),-c(1,2)])))
})

tau_genes$global_mean_nc <- sapply(combined_medians$gene, function(x){
  mean(na.omit(as.numeric(combined_medians[which(combined_medians$gene==x),-c(1,2)])))
})

tau_genes$min_tissue_nc <- sapply(combined_medians$gene, function(x){
  min(na.omit(as.numeric(combined_medians[which(combined_medians$gene==x),-c(1,2)])))
})


#plot tau versus global median noncoding fraction
ggplot(tau_genes, aes(y=tau, x=global_median_nc, color=global_median_nc))+
  geom_point(size=2.5,alpha=0.175)+
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed", size = 1)+
  geom_hline(yintercept = 0.25, color="red", linetype = "dashed", size = 1)+
  geom_hline(yintercept=0.75, color="red", linetype="dashed", size=1)+
  labs(x="Global Median Noncoding Fraction", y="Tau score")+
  theme(plot.background = element_blank(),
        axis.text= element_text(size=13, color = 'black'),
        axis.title = element_text(size=13, color = 'black', face = "bold"),
        panel.border = element_rect( fill=NA), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "lightgrey"),
        panel.grid.minor = element_blank(),
        legend.position="right",
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        strip.text = element_text(size = 13))+
  scale_color_viridis_c()

#combine the combined medians with the tau scores as well as global medians!
rownames(tau_genes) <- tau_genes$genes

combined_dat <- merge(tau_genes, plot_combined[,-28], by="row.names")
write.csv(combined_dat[,-1], file="combined_nc_medians_tau.csv")


library(corrplot)
library(ggcorrplot)
library(viridisLite)
#correlation plot
corrplot_dat <- cor(plot_combined[,-28])
ggcorrplot(corrplot_dat, hc.order = TRUE, outline.col = "white",
           ggtheme=theme(plot.background = element_blank(),
                        axis.text= element_text(size=2.5, color = 'black'),
                        panel.border = element_rect(color="darkgrey", fill=NA), 
                        panel.background = element_rect(fill="white"),
                        panel.grid.minor = element_blank(),
                        legend.position="right",
                        legend.title=element_blank(),
                        legend.text=element_text(size=10))) +
  scale_fill_gradientn(colors = viridis(256, option = 'D'))
#save data!
write.csv(corrplot_dat, file="correlation_median_nc_GTEX_tissues.csv")


corrplot_dat2 <- cor(t(plot_combined[,-28]))
png("correlation_plot_genes.png")
ggcorrplot(corrplot_dat2, outline.col = "white",
           ggtheme=theme(plot.background = element_blank(),
                         axis.text= element_text(size=2.5, color = 'black'),
                         panel.border = element_rect(color="darkgrey", fill=NA), 
                         panel.background = element_rect(fill="white"),
                         panel.grid.minor = element_blank(),
                         legend.position="right",
                         legend.title=element_blank(),
                         legend.text=element_text(size=10))) +
  scale_fill_gradientn(colors = viridis(256, option = 'D'))
dev.off()
