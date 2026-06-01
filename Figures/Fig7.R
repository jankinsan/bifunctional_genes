#figure for ORFs!
library("ggplot2")
library("RColorBrewer")


#numbers checked from predicted files to plot!
orf_dat <- data.frame(ATG=c(53430, 28059), CTG=c(74227, 38482), 
                      GTG=c(59214, 31571), TTG=c(52403, 27548),
                      row.names = c("predicted", "unique"))
orf_plot <- reshape2::melt(t(orf_dat))

#Fig a
#General number of predicted versus unique ORFs
ggplot(data=orf_plot, aes(x=Var2, y=value, fill=Var1))+
  geom_bar(position="dodge", stat="identity", color="black", alpha=0.8) + 
  labs(fill="Start Codon", y="Number of ORFs", x="ORFs")+
  theme_bw()+scale_fill_brewer(palette = "Set2")+
  theme(axis.text.x = element_text(size=13, colour="black"),
        axis.title=element_text(size=13, color="black", face="bold"), 
        axis.text.y=element_text(size=13, color="black"), 
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="bottom",
        legend.title=element_text(size=13, color="black", face="bold"),
        legend.text = element_text(size=12, colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())

#read ORF data!
setwd("E:/Janki/bifunctional_genes/data/ORFs/uniqueORFs_nmdMatch")
atg_orfs <- read.csv("2024-07-17_unique0RFs_bifunc_ATG_nmd_blast_mRNAs_KSS.csv")
ctg_orfs <- read.csv("2024-07-17_unique0RFs_bifunc_CTG_nmd_blast_mRNAs_KSS.csv")
gtg_orfs <- read.csv("2024-07-17_unique0RFs_bifunc_GTG_nmd_blast_mRNAs_KSS.csv")
ttg_orfs <- read.csv("2024-07-17_unique0RFs_bifunc_TTG_nmd_blast_mRNAs_KSS.csv")


#get orf length (in aa)
num_orfs <- data.frame(categories=c("150-300", "301-1500", "1501-3000", "3001-6000", ">6000"))
num_orfs$ATG <- c(length(which(atg_orfs$len_strip<301)),
                  length(intersect(which(atg_orfs$len_strip>=301), which(atg_orfs$len_strip>=1500))),
                  length(intersect(which(atg_orfs$len_strip>=1501), which(atg_orfs$len_strip>=3000))),
                  length(intersect(which(atg_orfs$len_strip>=3001), which(atg_orfs$len_strip>=6000))),
                  length(which(atg_orfs$len_strip>6000)))
num_orfs$CTG <- c(length(which(ctg_orfs$len_strip<301)),
                  length(intersect(which(ctg_orfs$len_strip>=301), which(ctg_orfs$len_strip>=1500))),
                  length(intersect(which(ctg_orfs$len_strip>=1501), which(ctg_orfs$len_strip>=3000))),
                  length(intersect(which(ctg_orfs$len_strip>=3001), which(ctg_orfs$len_strip>=6000))),
                  length(which(ctg_orfs$len_strip>6000)))
num_orfs$GTG <- c(length(which(gtg_orfs$len_strip<301)),
                  length(intersect(which(gtg_orfs$len_strip>=301), which(gtg_orfs$len_strip>=1500))),
                  length(intersect(which(gtg_orfs$len_strip>=1501), which(gtg_orfs$len_strip>=3000))),
                  length(intersect(which(gtg_orfs$len_strip>=3001), which(gtg_orfs$len_strip>=6000))),
                  length(which(gtg_orfs$len_strip>6000)))
num_orfs$TTG <- c(length(which(ttg_orfs$len_strip<301)),
                  length(intersect(which(ttg_orfs$len_strip>=301), which(ttg_orfs$len_strip>=1500))),
                  length(intersect(which(ttg_orfs$len_strip>=1501), which(ttg_orfs$len_strip>=3000))),
                  length(intersect(which(ttg_orfs$len_strip>=3001), which(ttg_orfs$len_strip>=6000))),
                  length(which(ttg_orfs$len_strip>6000)))

#plot
num_plot <- reshape2::melt(num_orfs)
num_plot$categories <- factor(num_plot$categories, levels = c("150-300", "301-1500", "1501-3000", "3001-6000", ">6000"))

#Fig b
#plot the lengths of ORFs!
ggplot(data=num_plot, aes(y=categories, x=value, fill=variable))+
  geom_bar(position="dodge", stat="identity", color="black", alpha=0.8) + xlim(0,28000)+
  labs(fill="Start Codon", x="Number of unique ORFs", y="Predicted ORF lengths (nucleotides)")+
  theme_bw()+theme(legend.position="bottom")+
  scale_fill_brewer(palette = "Set2")+
  geom_text(aes(label=value), position=position_dodge(width=0.9), hjust=-0.25)+
  theme(axis.text.x = element_text(size=13, colour="black"),
        axis.title=element_text(size=13, color="black", face="bold"), 
        axis.text.y=element_text(size=13, color="black"), 
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="bottom",
        legend.title=element_text(size=13, color="black", face="bold"),
        legend.text = element_text(size=12, colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())
#Fig D
#Kozak Score versus blast status for predicted ORFs!
atg_orfs$start_codon="ATG"
ctg_orfs$start_codon="CTG"
gtg_orfs$start_codon="GTG"
ttg_orfs$start_codon="TTG"
all_orfs <- rbind.data.frame(atg_orfs, ctg_orfs, gtg_orfs, ttg_orfs)
ggplot(data=all_orfs, aes(x=similarity_score, color=blast_status))+ 
  geom_density(linewidth=1.25)+labs(x="Kozak Similarity Score", 
                                    y="Density of predicted unique ORFs", colour="blastp status")+
  facet_grid(vars(start_codon))+ 
  geom_vline(xintercept = 0.64, color = "darkred", linetype = "dashed", size=1.5)+
  theme(axis.text.x = element_text(size=13, colour="black"),
        axis.title=element_text(size=13, color="black", face="bold"), 
        axis.text.y=element_text(size=13, color="black"), 
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="bottom",
        legend.title=element_text(size=13, color="black", face="bold"),
        legend.text = element_text(size=12, colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 13, face="bold"))


#Fig c
##Relative start position in the transcript
all_orfs$start_relative_pos <- all_orfs$start_position/all_orfs$sequence_length
ggplot(data=all_orfs, aes(x=start_relative_pos, color=start_codon))+ 
  geom_density(linewidth=1.25)+theme_bw()+ 
  scale_color_brewer(palette = "Set2")+ theme(legend.position="bottom")+
  labs(x="Start position along the transcript", y="Density", colour="start codon")+
  theme(axis.text.x = element_text(size=13, colour="black"),
        axis.title=element_text(size=13, color="black", face="bold"), 
        axis.text.y=element_text(size=13, color="black"), 
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="bottom",
        legend.title=element_text(size=13, color="black", face="bold"),
        legend.text = element_text(size=12, colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank())


#hydrophobicity
library(Peptides)
all_orfs$miyazawa <- sapply(all_orfs$x, function(x){
  hydrophobicity(x, scale="Miyazawa")
})
all_orfs$miyazawa_30aa <- sapply(all_orfs$x, function(x){
  prot_30aa <- substr(x, nchar(x)-29,nchar(x))
  hydrophobicity(prot_30aa,scale = "Miyazawa")
})
ggplot(data=all_orfs, aes(x=miyazawa_30aa, color=blast_status))+ 
  geom_density(linewidth=1.25)+labs(x="Miyazawa Hydrophobicity Score (of last 30aa)", 
                                    y="Density of predicted unique proteins", colour="blastp status")+
  facet_grid(vars(start_codon))+
  theme(axis.text.x = element_text(size=13, colour="black"),
        axis.title=element_text(size=13, color="black", face="bold"), 
        axis.text.y=element_text(size=13, color="black"), 
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="bottom",
        legend.title=element_text(size=13, color="black", face="bold"),
        legend.text = element_text(size=12, colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 13, face="bold"))

#start position along the transcript 
ggplot(data=all_orfs, aes(x=start_relative_pos, color=blast_status))+ 
  geom_density(linewidth=1.25)+labs(x="Start position along the transcript", 
                                    y="Density of predicted unique ORFs", colour="blastp status")+
  facet_grid(vars(start_codon))+
  theme(axis.text.x = element_text(size=13, colour="black"),
        axis.title=element_text(size=13, color="black", face="bold"), 
        axis.text.y=element_text(size=13, color="black"), 
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="bottom",
        legend.title=element_text(size=13, color="black", face="bold"),
        legend.text = element_text(size=12, colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 13, face="bold"))

#fig g
#denote as high-confidence and all_predicted ORFs!
all_orfs$hi_conf_label  <- "predicted"
all_orfs$hi_conf_label[intersect(which(all_orfs$blast_status=="blast_Cterm"),
                                 which(all_orfs$similarity_score>=0.64))] <- "high_confidence"
freq_orf<- data.frame(table(all_orfs$hi_conf_label))  
colnames(freq_orf) <- c("ORF_type", "Freq")

freq_orf$kozak_hi <- c(length(which(all_orfs$similarity_score>=0.64)),
                       dim(all_orfs)[1] - length(which(all_orfs$similarity_score>=0.64)))
freq_orf$blast_status <- c(length(which(all_orfs$blast_status=="blast_Cterm")),
                           dim(all_orfs)[1] - length(which(all_orfs$blast_status=="blast_Cterm")))

#get number of genes
hi_conf_orfs <- all_orfs[which(all_orfs$hi_conf_label=="high_confidence"),]
genes_hi <- hi_conf_orfs$genes[-grep(",", hi_conf_orfs$genes)]
genes_hi <- unique(unlist(c(genes_hi, sapply(hi_conf_orfs$genes[grep(",", hi_conf_orfs$genes)], function(x){
  unlist(strsplit(x, split=","))
}))))

#make pie plots for the ORFs 
ggplot(data=freq_orf, aes(x=" ", y=kozak_hi, fill=ORF_type)) +
  geom_bar(width =1, stat = "identity", color="white", linewidth=1, alpha=0.8) +
  coord_polar("y", start=0) + labs(fill="", title="Kozak Similarity Score >=0.64")+
  geom_text(aes(label = kozak_hi), position = position_stack(vjust = 0.5), size=6.75)+
  theme(
    plot.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(), 
    plot.title= element_text(size=13, colour = "black", face="bold", hjust=0.5,),
    axis.line.x = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position="bottom",
    strip.text.y = element_text(size=13, colour = "black", face="bold"),
    legend.title = element_text(size=13, colour = "black", face="bold"),
    legend.text = element_text(size=13,colour = "black"))+scale_fill_brewer(palette = "Set2", direction=-1)

ggplot(data=freq_orf, aes(x=" ", y=blast_status, fill=ORF_type)) +
  geom_bar(width =1, stat = "identity", color="white", linewidth=1, alpha=0.8) +
  coord_polar("y", start=0) + labs(fill="ORFs by confidence", title="ORFs by blasting status of resulting peptide")+
  geom_text(aes(label = blast_status), position = position_stack(vjust = 0.5), size=6.75)+
  theme(
    plot.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title= element_text(size=13, colour = "black", face="bold", hjust=0.5),
    axis.line.x = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position="bottom",
    strip.text.y = element_text(size=13, colour = "black", face="bold"),
    legend.title = element_text(size=13, colour = "black", face="bold"),
    legend.text = element_text(size=13,colour = "black"))+scale_fill_brewer(palette = "Set2", direction=-1)

ggplot(data=freq_orf, aes(x=" ", y=Freq, fill=ORF_type)) +
  geom_bar(width =1, stat = "identity", color="white", linewidth=1, alpha=0.8) +
  coord_polar("y", start=0) + labs(title="High-confidence ORF subset")+
  geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), size=6.75)+
  theme(
    plot.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(size=13, colour = "black", face="bold", hjust=0.5),
    axis.line.x = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position="bottom",
    strip.text.y = element_text(size=13, colour = "black", face="bold"),
    legend.title = element_blank(),
    legend.text = element_text(size=13,colour = "black"))+scale_fill_brewer(palette = "Set2", direction=-1)


ggplot(data=all_orfs, aes(x=miyazawa_30aa, color=hi_conf_label))+ 
  geom_density(linewidth=1.25)+labs(x="Miyazawa Hydrophobicity Score (of last 30aa)", 
                                    y="Density of predicted unique proteins", colour="ORF status")+
  facet_grid(vars(start_codon))+
  theme(axis.text.x = element_text(size=13, colour="black"),
        axis.title=element_text(size=13, color="black", face="bold"), 
        axis.text.y=element_text(size=13, color="black"), 
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="bottom",
        legend.title=element_text(size=13, color="black", face="bold"),
        legend.text = element_text(size=12, colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 13, face="bold"))+scale_color_brewer(palette = "Set2", direction=-1)

ggplot(data=all_orfs, aes(x=miyazawa_30aa, color=hi_conf_label))+ 
  geom_density(linewidth=1.25)+labs(x="Miyazawa Hydrophobicity Score (of last 30aa)", 
                                    y="Density of predicted unique proteins", colour="ORF status")+
  theme(axis.text.x = element_text(size=13, colour="black"),
        axis.title=element_text(size=13, color="black", face="bold"), 
        axis.text.y=element_text(size=13, color="black"), 
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position="bottom",
        legend.title=element_text(size=13, color="black", face="bold"),
        legend.text = element_text(size=12, colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 13, face="bold"))+scale_color_brewer(palette = "Set2", direction=-1)


write.csv(all_orfs, file="27052026_all_ORFs_nmb_blst_kss_hydrophobicity_confidence.csv")
