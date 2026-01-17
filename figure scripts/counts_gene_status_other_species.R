library("ggplot2")
setwd("F:/Janki/bifunctional_genes")
#counts for other species: analysis!
counts_status <- function(filePath, org){
  counts_dat <- read.csv(filePath)
  counts_dat$status<- "coding"
  counts_dat$status[counts_dat$NR_XR>0]<- "noncoding"
  counts_dat$status[intersect(which(counts_dat$NM_XM>0), which(counts_dat$NR_XR>0))]<- "bifunctional"
  write.csv(counts_dat, file = paste0(Sys.Date(), "_counts_", org, "_genewise_status.csv"))
  write.csv(counts_dat[counts_dat$status=="bifunctional",], file = paste0(Sys.Date(), "_counts_", org, "_bifunc_genes.csv"))
  summary_stat<- c(sum(counts_dat$total), sum(counts_dat$NM)+sum(counts_dat$NR),
                   sum(counts_dat$NM), sum(counts_dat$NR),
                   sum(counts_dat$XM)+sum(counts_dat$XR),
                   sum(counts_dat$XM), sum(counts_dat$XR), dim(counts_dat)[1],
                   length(which(counts_dat$status=="bifunctional")),
                   length(which(counts_dat$status=="coding")),
                   length(which(counts_dat$status=="noncoding")))
  return(summary_stat)
  }


#summarize datasets 
human_summary <- counts_status(filePath = "./data/counts/20062024__GRCh38.p14_genewise_status.csv",
                               org = "human")
mouse_summary <- counts_status(filePath = "./data/counts/other species/05032024_mouse_transcripts_counts.csv",
                              org = "mouse")
chimp_summary <- counts_status(filePath = "./data/counts/other species/05032024_chimpanzee_transcripts_counts.csv",
                               org = "chimpanzee")
zebfish_summary <- counts_status(filePath = "./data/counts/other species/05032024_zebrafish_transcripts_counts.csv",
                               org = "zebrafish")

#18-02-2025
rat_summary <- counts_status(filePath = "./data/counts/other species/18022025_rat_transcripts_counts.csv",
                             org = "Rat")
drosophila_summary <- counts_status(filePath = "./data/counts/other species/20042025_fruitfly_transcripts_counts.csv",
                             org = "Drosophila")
yeast_summary <- counts_status(filePath = "./data/counts/other species/18022025_yeast_transcripts_counts.csv",
                             org = "Yeast")
chicken_summary <- counts_status(filePath = "./data/counts/other species/18022025_chicken_transcripts_counts.csv",
                             org = "Chicken")
xenopus_summary <- counts_status(filePath = "./data/counts/other species/18022025_xenopus_transcripts_counts.csv",
                             org = "Xenopus")
celegans_summary <- counts_status(filePath = "./data/counts/other species/18022025_celegans_transcripts_counts.csv",
                             org = "c_elegans")

summ_species <- rbind.data.frame(c("Organism", "Human", "Chimpanzee", "Mouse", "Rat", "Frog", "Chicken", "Zebrafish", "Drosophila", "Worm", "Yeast"),
  cbind(c("Total transcripts", "Validated Transcripts", "NM_", "NR_", 
             "Predicted Transcripts", "XM_", "XR_", "Transcribed genes", 
             "Bifunctional genes", "Coding genes", "Noncoding genes"),
        human_summary, chimp_summary, mouse_summary, rat_summary, xenopus_summary, chicken_summary, zebfish_summary, drosophila_summary, celegans_summary, yeast_summary))
#add genome size
summ_species[13, 1:11]<- c("Genome size (billion base pairs)", 3.099,	3.050,	2.718,	2.8495,	2.7424,	1.0533,	1.373,	0.1437,	0.1003,	0.0121)
write.csv(summ_species, file=paste0(Sys.Date(), "_summary_other_species.csv"), row.names = FALSE)

#plot
summ_plot_num <- t(summ_species[10:12, -1])
colnames(summ_plot_num)<- summ_species[10:12, 1]
rownames(summ_plot_num)<- summ_species[1,-1]
summ_plot<- reshape2::melt(summ_plot_num)
summ_plot$value<- as.numeric(summ_plot$value)
plot_stacked_annot <- ggplot(data = summ_plot, 
                             aes(y = Var1, x=value, 
                                 fill=Var2, width=0.65)) +
  geom_bar(position="stack", stat="identity", color="black")+
  labs(x="Number of Genes", y="Organism", fill="Annotation") +
  theme(plot.background = element_blank(),
        axis.text= element_text(size = 13, color = 'black'),
        axis.title = element_text(size = 13, color = 'black', face = "bold"),
        panel.border = element_rect(fill=NA), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "lightgrey"),
        panel.grid.minor = element_blank(),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text = element_text(size=12, colour = "black"))+
  scale_fill_manual(values = c("Coding genes"="#e8e288", "Bifunctional genes" = "#82c0cc", "Noncoding genes" = "#ffa5a5"))
  #scale_x_discrete(labels=c("All Genes", "PC Genes", "NPC Genes"))
plot(plot_stacked_annot)


#22-03-2025
cattle_summary <- counts_status(filePath = "./data/counts/other species/22032025_cattle_transcripts_counts.csv",
                                org = "cattle")
arabidposis_summary <- counts_status(filePath = "./data/counts/other species/22032025_arabidopsis_transcripts_counts.csv",
                                     org = "arabidopsis")
rice_summary <- counts_status(filePath = "./data/counts/other species/22032025_rice_transcripts_counts.csv",
                              org = "rice")
aniger_summary <- counts_status(filePath = "./data/counts/other species/22032025_aniger_transcripts_counts.csv",
                                org = "aniger")  
afumi_summary <- counts_status(filePath = "./data/counts/other species/22032025_afumi_transcripts_counts.csv",
                               org = "afumi")
candida_summary <- counts_status(filePath = "./data/counts/other species/22032025_candida_transcripts_counts.csv",
                                 org = "candida")
rhizopus_summary <- counts_status(filePath = "./data/counts/other species/22032025_rhizopus_transcripts_counts.csv",
                                  org = "rhizopus")
spombe_summary <- counts_status(filePath = "./data/counts/other species/22032025_spombe_transcripts_counts.csv",
                                org = "spombe")
wheat_summary <- counts_status(filePath = "./data/counts/other species/22032025_wheat_transcripts_counts.csv",
                               org="wheat")


#summary 
#for plant species
summ_plant <- rbind.data.frame(c("Organism", "Arabidopsis thaliana", "Oryza sativa", "Triticum aestivum"),
cbind(c("Total transcripts", "Validated Transcripts", "NM_", "NR_", 
        "Predicted Transcripts", "XM_", "XR_", "Transcribed genes", 
        "Bifunctional genes", "Coding genes", "Noncoding genes"),
       arabidposis_summary, rice_summary, wheat_summary))
write.csv(summ_plant, file=paste0(Sys.Date(), "_summary_plant_species.csv"), row.names = FALSE)
#plot
summ_plant_num <- t(summ_plant[10:12, -1])
colnames(summ_plant_num)<- summ_plant[10:12, 1]
rownames(summ_plant_num)<- summ_plant[1,-1]
summ_plant_plot<- reshape2::melt(summ_plant_num)
summ_plant_plot$value<- as.numeric(summ_plant_plot$value)
plot_stacked_plant <- ggplot(data = summ_plant_plot, 
                             aes(y = Var1, x=value, 
                                 fill=Var2, width=0.65)) +
  geom_bar(position="stack", stat="identity", color="black")+
  labs(x="Number of Genes", y="Organism", fill="Annotation") +
  theme(plot.background = element_blank(),
        axis.text.y= element_text(size = 12, color = 'black', face = "italic"),
        axis.text.x= element_text(size = 13, color = 'black'),
        axis.title = element_text(size = 13, color = 'black', face = "bold"),
        panel.border = element_rect(fill=NA), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "lightgrey"),
        panel.grid.minor = element_blank(),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text = element_text(size=12, colour = "black"))+
  scale_fill_manual(values = c("Coding genes"="#e8e288", "Bifunctional genes" = "#82c0cc", "Noncoding genes" = "#ffa5a5"))
#scale_x_discrete(labels=c("All Genes", "PC Genes", "NPC Genes"))
svg("./figures/counts/other species/plants_summary.svg")
plot(plot_stacked_plant)
dev.off()

#fungi species plot!
summ_fungal<- rbind.data.frame(c("Organism",  "Aspergillus niger", "Aspergillus fumigatus", "Candida albicans", "Rhizopus microsporus", "Schizosaccharomyces pombe", "Saccharomyces cerevisiae"),
                                 cbind(c("Total transcripts", "Validated Transcripts", "NM_", "NR_", 
                                         "Predicted Transcripts", "XM_", "XR_", "Transcribed genes", 
                                         "Bifunctional genes", "Coding genes", "Noncoding genes"),
                                       aniger_summary, afumi_summary, candida_summary, rhizopus_summary, spombe_summary, yeast_summary))
write.csv(summ_fungal, file=paste0(Sys.Date(), "_summary_fungal_species.csv"), row.names = FALSE)
#plot
summ_fungal_num <- t(summ_fungal[10:12, -1])
colnames(summ_fungal_num)<- summ_fungal[10:12, 1]
rownames(summ_fungal_num)<- summ_fungal[1,-1]
summ_fungal_plot<- reshape2::melt(summ_fungal_num)
summ_fungal_plot$value<- as.numeric(summ_fungal_plot$value)
plot_stacked_fungal <- ggplot(data = summ_fungal_plot, 
                              aes(y = Var1, x=value, 
                                  fill=Var2, width=0.65)) +
  geom_bar(position="stack", stat="identity", color="black")+
  labs(x="Number of Genes", y="Organism", fill="Annotation") +
  theme(plot.background = element_blank(),
        axis.text.y= element_text(size = 12, color = 'black', face = "italic"),
        axis.text.x= element_text(size = 13, color = 'black'),
        axis.title = element_text(size = 13, color = 'black', face = "bold"),
        panel.border = element_rect(fill=NA), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "lightgrey"),
        panel.grid.minor = element_blank(),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text = element_text(size=12, colour = "black"))+
  scale_fill_manual(values = c("Coding genes"="#e8e288", "Bifunctional genes" = "#82c0cc", "Noncoding genes" = "#ffa5a5"))
#scale_x_discrete(labels=c("All Genes", "PC Genes", "NPC Genes"))
svg("./figures/counts/other species/fungals_summary.svg")
plot(plot_stacked_fungal)
dev.off()

#Plotting the number of transcribed, bifunctional, coding and noncoding genes in all organisms
library(ggplot2)
#taking only data needed
summ_t <- t(summ_species[9:12,-1])
colnames(summ_t) <- summ_species[9:12, 1]
#converting to numeric 
summ_t <- apply(summ_t, 2, as.numeric)
rownames(summ_t)<- colnames(summ_species)[-1]
plot_num <- reshape2::melt(summ_t)
#ggplot with facets for each type of gene: transcribed, bifunctional, coding or noncoding?
ggplot(plot_num, aes(x=Var1, y=value, fill=Var1, label=value))+geom_bar(stat="identity")+
  geom_text(vjust=-0.2)+ylim(0,50000)+
  facet_wrap(~Var2)+theme_bw()+theme(axis.text.x = element_text(angle=90, hjust=-0.1))

#linear regression plots
#formatting data to plot! and melting the dataframe!
summ_species[14,] <- c("Transcribed genes (tens of thousands)", as.numeric(summ_species[9,-1])/10000)
summ_species[15,] <- c("Bifunctional genes (thousands)", as.numeric(summ_species[10,-1])/1000)
summ_reg<- t(summ_species[13:15, -1])
colnames(summ_reg)<- summ_species[13:15, 1]
rownames(summ_reg)<- summ_species[1,-1]
summ_reg_plot<- reshape2::melt(summ_reg)
summ_reg_plot$value<- as.numeric(summ_reg_plot$value)

#get organism names separately to use for x-axis labels 
org_names <- levels(factor(summ_reg_plot$Var1)) # Get ordered category names
plot_reg <- ggplot(data = summ_reg_plot, 
                  aes(x = as.numeric(factor(Var1)), color= Var2, 
                      shape=Var2, y=value)) +
  geom_point(size=5, alpha=0.7)+
  geom_smooth(method = "lm", se = FALSE, linewidth=1.5) +
  labs(y="Number", x="Organism", color="", shape="") +
  guides(color = guide_legend(nrow = 3, byrow = TRUE))+
  theme(plot.background = element_blank(),
        axis.text.x= element_text(size = 13, color = 'black', angle=45, vjust = 0.65),
        axis.text.y= element_text(size = 13, color = 'black'),
        axis.title = element_text(size = 13, color = 'black', face = "bold"),
        panel.border = element_rect(fill=NA), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "lightgrey"),
        panel.grid.minor = element_blank(),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text = element_text(size=12, colour = "black")) +
  scale_color_brewer(palette = "Dark2")+
  scale_x_continuous(breaks = seq_along(org_names), labels = org_names)
plot(plot_reg)

ggplot(summ_reg_plot, aes(x = Var1, y = value, color= Var2, 
                          shape=Var2, y=value)) +
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun = mean, geom = "line")
# plot the data


#drosophila redo!#drosophil summary redo!
#NO NEED TO RUN since I was able to modify the python code and write one specifically for
drosophila_gtf <- read.delim("./data/counts/other species/fruitfly/genomic.gtf",
                             skip=4, header=FALSE)
droso_trans <- drosophila_gtf[drosophila_gtf[,3]=="transcript",]
droso_trans$gene <- sapply(droso_trans[,9], function(x){
  split<-strsplit(x, split=";")[[1]]
  return(strsplit(split[grep("gene ", as.character(split))], split="gene ")[[1]][2])})
droso_trans$transcript_id <- sapply(droso_trans[,9], function(x){strsplit(x, split=";")[[1]][2]})
droso_trans$transcript_id <- substr(droso_trans$transcript_id, 15, 42)

