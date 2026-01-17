#microRNA Analysed data visualization!
#load the libraries needed!
library(ggplot2)


#read data
mir_gene_summary <- read.csv("F:/Janki/bifunctional_genes/data/microRNA/2024-08-28_counts_common_microRNA_targets.csv")

#arrange data into gene_type targets!
mir_gene_plot <- cbind.data.frame(mir=unique(mir_gene_summary$microRNA), 
                                  coding_gene_targets=0, 
                                  noncoding_gene_targets=0,
                                  bifunctional_gene_targets=0,
                                  total_gene_targets=0,
                                  percent_coding_targets=0,
                                  percent_noncoding_targets=0,
                                  percent_bifunc_targets=0)

#get number for gene targets and other things!
mir_gene_plot$total_gene_targets <- sapply(mir_gene_plot$mir, function(mir){
return(length(which(mir_gene_summary$microRNA==mir)))
})
mir_gene_plot$coding_gene_targets <- sapply(mir_gene_plot$mir, function(mir){
  idx <- which(mir_gene_summary$microRNA==mir)
  return(length(which(mir_gene_summary$gene_type[idx]=="only_mRNA")))
})
mir_gene_plot$noncoding_gene_targets <- sapply(mir_gene_plot$mir, function(mir){
  idx <- which(mir_gene_summary$microRNA==mir)
  return(length(which(mir_gene_summary$gene_type[idx]=="only_ncRNA")))
})
mir_gene_plot$bifunctional_gene_targets <- sapply(mir_gene_plot$mir, function(mir){
  idx <- which(mir_gene_summary$microRNA==mir)
  return(length(which(mir_gene_summary$gene_type[idx]=="hybrid")))
})

#get percentage of each gene type targets!
mir_gene_plot$percent_coding_targets <- mir_gene_plot$coding_gene_targets/mir_gene_plot$total_gene_targets
mir_gene_plot$percent_noncoding_targets <- mir_gene_plot$noncoding_gene_targets/mir_gene_plot$total_gene_targets
mir_gene_plot$percent_bifunc_targets <- mir_gene_plot$bifunctional_gene_targets/mir_gene_plot$total_gene_targets

plot_dat <- reshape2::melt(mir_gene_plot[,1:4])
#factorize mir name and type of gene targets!
plot_dat$mir <- as.factor(plot_dat$mir)
plot_dat$variable <- as.factor(plot_dat$variable)

#plot
ggplot(data=plot_dat, aes(x=" ", y=value, group=mir, fill=variable)) +
  geom_bar(width = 1, stat = "identity", color="black") +
  coord_polar("y", start=0) + labs(title="Number of target genes", fill="Target Gene Type")+
  facet_wrap(.~ mir, ncol=18) +theme(
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
        strip.text = element_text(size = 7.5, face = "bold", colour = "black", margin = margin(t = 5, b = 5, l = 2, r = 2)),
        legend.position="bottom",
        legend.text = element_text(size=10, colour = "black"))+
  scale_fill_manual(values = c("coding_gene_targets"="#e8e288", "bifunctional_gene_targets" = "#82c0cc", "noncoding_gene_targets" = "#ffa5a5"))


#plots for percentage 
plot_percent <- reshape2::melt(mir_gene_plot[,c(1,6:8)])
#factorize mir name and type of gene targets!
plot_percent$mir <- as.factor(plot_dat$mir)
plot_percent$variable <- as.factor(plot_percent$variable)

#plot
ggplot(data=plot_percent, aes(x=" ", y=value, group=mir, fill=variable)) +
  geom_bar(width =0.5, stat = "identity", color="black") +
  coord_polar("y", start=0) + labs(title="Ratio of target genes", fill="Target Gene Type")+
  facet_wrap(.~ mir, nrow=15) +theme(
    plot.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(color = 'black', face = "bold", hjust=0.5, size=10),
    panel.background = element_rect(fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    #strip.text = element_text(size = 7.5, face = "bold", colour = "black"), 
                              #margin = margin(t = 5, b = 5, l = 2, r = 2),
    legend.position="right",
    legend.title = element_text(size=8.5, colour = "black", face="bold"),
    legend.text = element_text(size=7.5, colour = "black"))+
  scale_fill_manual(values = c("percent_coding_targets"="#e8e288", "percent_bifunc_targets" = "#82c0cc", "percent_noncoding_targets" = "#ffa5a5"))


ggplot(data=plot_percent, aes(x=" ", y=value, group=mir, fill=variable)) +
  geom_bar(width =1, stat = "identity", color="black", linewidth=0.65) +
  coord_polar("y", start=0) + labs(title="Ratio of target genes", fill="Target Gene Type")+
  facet_wrap(.~ mir, ncol=8) +theme(
    plot.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(color = 'black', face = "bold", hjust=0.5, size=10),
    panel.background = element_rect(fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position="",
    legend.title = element_text(size=8, colour = "black", face="bold"),
    legend.text = element_text(size=8,colour = "black"))+
  scale_fill_manual(values = c("percent_coding_targets"="#e8e288", "percent_bifunc_targets" = "#82c0cc", "percent_noncoding_targets" = "#ffa5a5"))


microRNA_dat <- read.csv("2024-08-28_counts_common_microRNA_targets.csv")
#mir-21
mir21 <- microRNA_dat[which(microRNA_dat$microRNA=="hsa-miR-21-5p"),]

mir21_plot <- reshape2::melt(mir21[,c(2,10,11)])
mir21_plot$variable <- factor(mir21_plot$variable, levels=c("NR_XR","NM_XM"))
ggplot(data=mir21_plot, aes(x=gene, y=value, fill=variable, width=0.65))+
  geom_bar(position="stack", stat="identity", colour="black")+
  theme_bw()+theme(axis.text.x=element_text(angle=90, face = "bold", colour="black"))+
  labs(x="Target Genes", y="Number of transcripts", fill="Transcript Type")+
  scale_fill_discrete(labels=c("Noncoding", "Coding"))

#mir-129
mir129 <- microRNA_dat[which(microRNA_dat$microRNA=="hsa-miR-129-5p"),]
mir129_plot <- reshape2::melt(mir129[,c(2,10,11)])
mir129_plot$variable <- factor(mir129_plot$variable, levels=c("NR_XR","NM_XM"))
ggplot(data=mir129_plot, aes(x=gene, y=value, fill=variable, width=0.65))+
  geom_bar(position="stack", stat="identity", colour="black")+
  theme_bw()+theme(axis.text.x=element_text(angle=90, face = "bold", colour="black"))+
  labs(x="Target Genes", y="Number of transcripts", fill="Transcript Type")+
  scale_fill_discrete(labels=c("Noncoding", "Coding")) 

#mir-320a
mir320a <- microRNA_dat[which(microRNA_dat$microRNA=="hsa-miR-320a"),]
mir320_plot <- reshape2::melt(mir320a[,c(2,10,11)])
mir320_plot$variable <- factor(mir320_plot$variable, levels=c("NR_XR","NM_XM"))
ggplot(data=mir320_plot, aes(x=gene, y=value, fill=variable, width=0.65))+
  geom_bar(position="stack", stat="identity", colour="black")+
  theme_bw()+theme(axis.text.x=element_text(angle=90, face = "bold", colour="black"))+
  labs(x="Target Genes", y="Number of transcripts", fill="Transcript Type")+
  scale_fill_discrete(labels=c("Noncoding", "Coding")) 



print(sessionInfo())