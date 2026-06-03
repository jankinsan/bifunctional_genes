#figure for ORFs!
library("ggplot2")
library("RColorBrewer")

orf_dat <- data.frame(ATG=c(53430, 28059), CTG=c(74227, 38482), 
                      GTG=c(59214, 31571), TTG=c(52403, 27548),
                      row.names = c("predicted", "unique"))
orf_plot <- reshape2::melt(t(orf_dat))

ggplot(data=orf_plot, aes(x=Var2, y=value, fill=Var1))+
  geom_bar(position="dodge", stat="identity", color="black", alpha=0.8) + 
  labs(fill="Start Codon", y="Number of ORFs", x="ORFs")+
  theme_bw()+scale_fill_brewer(palette = "Set2")

#read ORF data!
setwd("E:/MSR BLY/THESIS WORK/NP NR Project/15122023/ORFs/uniqueORFs_nmdMatch")
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
ggplot(data=num_plot, aes(y=categories, x=value, fill=variable))+
  geom_bar(position="dodge", stat="identity", color="black", alpha=0.8) + xlim(0,28000)+
  labs(fill="Start Codon", x="Number of unique ORFs", y="Predicted ORF lengths (nucleotides)")+
  theme_bw()+scale_fill_brewer(palette = "Set2")+geom_text(aes(label=value), position=position_dodge(width=0.9), hjust=-0.25)
ggplot(data=num_plot, aes(y=categories, x=value, fill=variable))+
  geom_bar(position="dodge", stat="identity", color="black", alpha=0.8) + xlim(0,28000)+
  labs(fill="Start Codon", x="Number of unique ORFs", y="Predicted ORF lengths (nucleotides)")+
  theme_bw()+theme(legend.position="bottom")+
  scale_fill_brewer(palette = "Set2")+geom_text(aes(label=value), position=position_dodge(width=0.9), hjust=-0.25)


#KSS versus blast status!
atg_orfs$start_codon="ATG"
ctg_orfs$start_codon="CTG"
gtg_orfs$start_codon="GTG"
ttg_orfs$start_codon="TTG"
all_orfs <- rbind.data.frame(atg_orfs, ctg_orfs, gtg_orfs, ttg_orfs)
ggplot(data=all_orfs, aes(x=similarity_score, color=blast_status))+ 
geom_density(linewidth=1.25)+labs(x="Kozak Similarity Score", y="Density of predicted unique proteins", colour="blastp status")+
  facet_grid(vars(start_codon))+ theme_bw()+theme(legend.position = "bottom")

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
  geom_density(linewidth=1.25)+labs(x="Miyazawa Hydrophobicity Score (of last 30aa)", y="Density of predicted unique proteins", colour="blastp status")+
  facet_grid(vars(start_codon))+ theme_bw()


all_orfs$start_relative_pos <- all_orfs$start_position/all_orfs$sequence_length
ggplot(data=all_orfs, aes(x=start_relative_pos, color=start_codon))+ 
  geom_density(linewidth=1.25)+theme_bw()+ 
  scale_color_brewer(palette = "Set2")+ theme(legend.position="bottom")+
  labs(x="Start position along the transcript", y="Density", colour="start codon")

write.csv(all_orfs, file="24112024_all_ORFs_nmb_blst_kss_hydrophobicity.csv")
