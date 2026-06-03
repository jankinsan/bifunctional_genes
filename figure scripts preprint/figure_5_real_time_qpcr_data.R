# September 17, 2025
library("ggplot2")
library("reshape2")
library("dplyr")
library("tidyr")
library("ggpubr")
library("cowplot")
setwd("F:/Janki/bifunctional_genes")
# Real-time expression data for noncoding transcripts of bifunctional genes!
real_time_dat <- read.csv("./data/expression/cell_lines_real_time_reps.csv")
real_time_dat$Cell_line <- as.factor(real_time_dat$Cell_line)
real_time_dat$Average <- as.numeric(real_time_dat$Average)
real_time_dat$stdev <- as.numeric(real_time_dat$stdev)

real_time_dat_fin <- real_time_dat[-which(real_time_dat$Cell_line=="SHSY5Y"),]


#color with cell line
plot_dat_cell <- function(Gene, id) { # nolint: object_name_linter. # nolint
  plot_dat_cell <- real_time_dat_fin[which(real_time_dat_fin$Gene == Gene), ]
  
  # Reshape replicate columns to long format
  plot_dat_cell_long <- plot_dat_cell %>%
    pivot_longer(cols = starts_with("Rep"), names_to = "Replicate", values_to = "Value")
  
  # plot with p-value labels
  ggplot(plot_dat_cell_long, aes(x = Cell_line, y = Average, fill=Cell_line)) +
    geom_bar(
      stat = "identity", color = "black", position = position_dodge(),
      width = 0.7,) + 
    ylim (0, max(plot_dat_cell$Average + plot_dat_cell$stdev)*1.22)+
    labs(x = "Cell Line", y = "Relative Expression", title=paste0("nc", Gene)) +
    geom_errorbar(aes(ymin = Average, ymax = Average + stdev),
                  width = .2, position = position_dodge(.9)
    ) +
    geom_point(
      data = plot_dat_cell_long, aes(x = Cell_line, y = Value), inherit.aes = TRUE, # nolint
      position = position_jitter(width = 0.15), size = 2, alpha=0.6
    ) +
    stat_compare_means(
      data = plot_dat_cell_long, aes(x = Cell_line, y = Value), method = "t.test",
      method.args = list(var.equal = TRUE, alternative = "two.sided"), 
      ref.group = "A549",
      label = "p.signif", # Show stars instead of p-values
      label.y = max(plot_dat_cell$Average + plot_dat_cell$stdev) * 1.1
    ) +
    theme(
      plot.background = element_blank(),
      axis.text.x = element_text(
        angle = 45, hjust = 1,
        color = "black", face = "bold"),
      axis.text.y = element_text(color = "black", face = "bold"),
      axis.title.y = element_text(color = "black", face = "bold"),
      axis.title.x=element_blank(), 
      plot.title = element_text(
        size = 15, colour = "black",
        face = "bold.italic", hjust = 0.5),
      panel.border = element_rect(linewidth = 1, fill = NA),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(colour = "lightgrey"),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) +
    scale_fill_brewer(palette = "Pastel2")
}
plot_grid(
  plot_dat_cell(Gene = "SNX16", id = 1),
  plot_dat_cell(Gene = "XIAP", id = 2),
  plot_dat_cell(Gene = "TIA1", id = 3),
  plot_dat_cell(Gene = "USP15", id = 4),
  plot_dat_cell(Gene = "FOXP1", id = 5),
  plot_dat_cell(Gene = "TGFBR1", id = 6),
  plot_dat_cell(Gene = "DDX3Y", id = 7),
  plot_dat_cell(Gene = "ACVR1B", id = 8),
  plot_dat_cell(Gene = "FUBP1", id = 9),
  ncol=3)

# thp1
thp_dat <- read.csv("./data/expression/thp_cell_lines_real_time.csv")
thp_dat$Condition <- as.factor(thp_dat$Condition)
thp_dat$Average <- as.numeric(thp_dat$Average)
thp_dat$stdev <- as.numeric(thp_dat$stdev)
plot_dat_thp <- function(Gene) {
  plot_dat <- thp_dat[which(thp_dat$Gene == Gene), ]
  # Reshape replicate columns to long format
  plot_dat_long <- plot_dat %>%
    pivot_longer(cols = starts_with("Rep"), 
                 names_to = "Replicate", values_to = "Value")
  
  # plot with p-value labels
  plot_fin <- ggplot(plot_dat_long, aes(x = Condition, y = Average, fill = Condition)) +
    geom_bar(alpha=0.7,
             stat = "identity", color = "black", position = position_dodge(),
             width = 0.45
    ) + ylim (0, max(plot_dat$Average + plot_dat$stdev)*1.22)+
    #scale_fill_manual(values = c("THP1" = "lightgray", "THP1+PMA" = "#7CAE00", "THP1+PMA+LPS" = "#CD9600")) +
    scale_fill_brewer(palette = "BuPu")+
    labs(x = "Condition", y = "Relative Expression", title = paste0("nc", Gene)) +
    geom_errorbar(aes(ymin = Average, ymax = Average + stdev),
                  width = .2, position = position_dodge(.9)
    ) +
    geom_point(
      data = plot_dat_long, aes(x = Condition, y = Value), # nolint # nolint
      position = position_jitter(width = 0.15), size = 2, alpha=0.6
    ) +
    stat_compare_means(
      data = plot_dat_long, aes(x = Condition, y = Value), method = "t.test",
      paired = TRUE, ref.group = "THP1",
      label = "p.signif", # Show stars instead of p-values
      label.y = max(plot_dat$Average + plot_dat$stdev) * 1.1
    ) +
    theme(
      plot.background = element_blank(),
      title=element_text(face="italic"),
      axis.text.x = element_blank(),
      axis.ticks.x=element_blank(),
      #element_text(size=13, angle = 45, hjust = 1, color = "black", face = "bold"),
      axis.text.y = element_text(size=13, color = "black"),
      axis.title.y = element_text(color = "black", face = "bold"),
      axis.title.x=element_blank(), 
      plot.title = element_text(size = 15, colour = "black", face = "bold.italic", hjust = 0.5),
      panel.background = element_blank(),
      panel.border = element_rect(fill=NA),
      panel.grid.major = element_line(colour = "lightgrey"),
      panel.grid.minor = element_blank(),
      #legend.position = "none"
    )
}

combined <- plot_grid(
  plot_dat_thp(Gene = "SNX16")+theme(legend.position="none"),
  plot_dat_thp(Gene = "XIAP")+theme(legend.position="none"),
  plot_dat_thp(Gene = "TIA1")+theme(legend.position="none"),
  plot_dat_thp(Gene = "USP15")+theme(legend.position="none"),
  plot_dat_thp(Gene = "FOXP1")+theme(legend.position="none"),
  plot_dat_thp(Gene = "TGFBR1")+theme(legend.position="none"),
  plot_dat_thp(Gene = "DDX3Y")+theme(legend.position="none"),
  plot_dat_thp(Gene = "ACVR1B")+theme(legend.position="none"),
  plot_dat_thp(Gene = "FUBP1")+theme(legend.position="none"),
  ncol=3)

legend <- get_legend(plot_dat_thp(Gene = "FUBP1")+
                       theme(legend.position="bottom", 
                             legend.title = element_text(size=13, color = "black", face="bold"), 
                             legend.text = element_text(size=13, color = "black", face="bold")))

plot_grid(combined, legend, ncol = 1, rel_heights = c(1, 0.1))

sessionInfo()

