##############################################################################################################
#figure s1 
library(ggplot2)
library(UpSetR)
library(cowplot)

#read transcripts counts and status as per refseq
counts_hg38 <- read.csv("F:/Janki/bifunctional_genes/data/counts/23112024_GRCh38.p14_genewise_status.csv")
counts_hg38$status<- sub("mRNA encoding", "coding", counts_hg38$status)
counts_hg38$status<- sub("ncRNA encoding", "noncoding", counts_hg38$status)

bifunc <- counts_hg38[counts_hg38$status=="bifunctional", ]
bifunc$percent_noncoding_transcripts <- bifunc$NR_XR/bifunc$total

#Percent noncoding transcripts
#Density plot for Fig S1 (G)
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
#for Fig S1 (a)  TO (d)
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

print(sessionInfo())
