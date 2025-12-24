##Enrichment/Over-representation Analysis figure
## What processes are bifunctional genes involved in?
## Figure 3
#September 16, 2025 v4
library("GSA")
library("clusterProfiler")  # A popular R package for enrichment analysis
library("org.Hs.eg.db")
library("ggplot2")
library("stringr")

setwd("F:/Janki/bifunctional_genes/")

#READ hg38 genes
hg38_genes <- read.csv("./data/GCF_000001405.40_GRCh38.p14_genes.csv")
hg38_genes$length <- hg38_genes$end - hg38_genes$start

#read counts and status as for bifunctionality!
gene_status <- read.csv("./data/counts/18122023_GRCh38.p14_genewise_counts_status.csv")
#update the symbols so that the gene symbols match!
changed_idx <- which(is.na(match(gene_status$X, hg38_genes$gene)))
gene_status$final_gene <- gene_status$X
gene_status$final_gene[changed_idx] <- c("HSALR1", "GAR1-DT")
hg38_genes$gene_status<-gene_status$status[match(hg38_genes$gene, gene_status$final_gene)]

#get only transcribed genes
#getting the biotypes of transcribed genes
hg38_genes2<-hg38_genes[-which(is.na(hg38_genes$gene_status)),]
hg38_genes2$gene_id <- sapply(hg38_genes2$details, function(x){
  strsplit(strsplit(x, split=";")[[1]][3], split="GeneID:")[[1]][2]
})

#get gene ids into separate lists!
geneids <-  hg38_genes2$gene_id
hybrid_geneids <- hg38_genes2$gene_id[hg38_genes2$gene_status=="hybrid"]
coding_geneids <- hg38_genes2$gene_id[hg38_genes2$gene_status=="only_mRNA"]
noncoding_geneids<- hg38_genes2$gene_id[hg38_genes2$gene_status=="only_ncRNA"]

#Over-representation analysis from clusterprofiler
ora_type <- function(sigdb_type) {
  # assigning gene lists:
  genes_of_interest <- hybrid_geneids # Vector of 4378 gene IDs
  background_genes <- hg38_genes2$gene_id  # Vector of 42987 gene IDs
  
  # Perform GO enrichment analysis: Over-representation Analysis
  print("Performing ORA for all gene categories..")
  if(length(grep("GO", sigdb_type))==1){
    sigdb_end <- strsplit(sigdb_type, "GO_")[[1]][2]
    go_enrichment <- enrichGO(gene = genes_of_interest,
                              universe = background_genes,
                              OrgDb = org.Hs.eg.db,  # For human genes
                              ont = sigdb_end, 
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.01,
                              qvalueCutoff = 0.2)
    
    go_enrichment_mrna <- enrichGO(gene = coding_geneids,
                                   universe = background_genes,
                                   OrgDb = org.Hs.eg.db,  # For human genes
                                   ont = sigdb_end,  
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.01,
                                   qvalueCutoff = 0.2)
    
    go_enrichment_ncrna <- enrichGO(gene = noncoding_geneids,
                                    universe = background_genes,
                                    OrgDb = org.Hs.eg.db,  # For human genes
                                    ont = sigdb_end,  
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.01,
                                    qvalueCutoff = 0.2)
    GO_IDs <- unique(c(go_enrichment@result$ID, 
                       go_enrichment_mrna@result$ID,
                       go_enrichment_ncrna@result$ID))
    
    # combining data from each gene type into one!
    print("Combining all gene categories into one dataset..")
    combined_data <- do.call(rbind, lapply(GO_IDs, function(id) {
      bifunc_id = which(go_enrichment@result$ID == id)
      mrna_id = which(go_enrichment_mrna@result$ID == id)
      ncrna_id = which(go_enrichment_ncrna@result$ID == id)
      
      bifunc_data <- if(length(bifunc_id) > 0) go_enrichment@result[bifunc_id, -c(1, 2)] else as.data.frame(matrix(NA, nrow=1, ncol=ncol(go_enrichment@result) - 2))
      mrna_data <- if(length(mrna_id) > 0) go_enrichment_mrna@result[mrna_id, -c(1, 2)] else as.data.frame(matrix(NA, nrow=1, ncol=ncol(go_enrichment_mrna@result) - 2))
      ncrna_data <- if(length(ncrna_id) > 0) go_enrichment_ncrna@result[ncrna_id, -c(1, 2)] else as.data.frame(matrix(NA, nrow=1, ncol=ncol(go_enrichment_ncrna@result) - 2))
      
      colnames(bifunc_data) <- paste("bifunc", colnames(go_enrichment@result)[-c(1, 2)], sep = ".")
      colnames(mrna_data) <- paste("mrna", colnames(go_enrichment_mrna@result)[-c(1, 2)], sep = ".")
      colnames(ncrna_data) <- paste("ncrna", colnames(go_enrichment_ncrna@result)[-c(1, 2)], sep = ".")
      
      cbind(bifunc_data, mrna_data, ncrna_data)
    }))
    
    # Extract the ID and Description columns once from go_enrichment
    id_description <- go_enrichment@result[match(GO_IDs, go_enrichment@result$ID), c("ID", "Description")]
    
    
    # Combine the ID and Description with the combined data
    final_combined_data <- cbind(id_description, combined_data)
    
  } else if (sigdb_type == "KEGG"){
    kegg_enrichment <- enrichKEGG(gene = genes_of_interest,
                                  universe = background_genes,
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.01,
                                  qvalueCutoff = 0.2)
    
    kegg_enrichment_mrna <- enrichKEGG(gene = coding_geneids,
                                       universe = background_genes,
                                       pAdjustMethod = "BH",
                                       pvalueCutoff = 0.01,
                                       qvalueCutoff = 0.2)
    
    kegg_enrichment_ncrna <- enrichKEGG(gene = noncoding_geneids,
                                        universe = background_genes,
                                        pAdjustMethod = "BH",
                                        pvalueCutoff = 0.01,
                                        qvalueCutoff = 0.2)
    kegg_IDs <- unique(c(kegg_enrichment@result$ID, 
                         kegg_enrichment_mrna@result$ID,
                         kegg_enrichment_ncrna@result$ID))
    
    
    
    # combining data from each gene type into one!
    print("Combining all gene categories into one dataset..")
    combined_data <- do.call(rbind, lapply(kegg_IDs, function(id) {
      bifunc_id = which(kegg_enrichment@result$ID == id)
      mrna_id = which(kegg_enrichment_mrna@result$ID == id)
      ncrna_id = which(kegg_enrichment_ncrna@result$ID == id)
      
      bifunc_data <- if(length(bifunc_id) > 0) kegg_enrichment@result[bifunc_id, -c(1, 2)] else as.data.frame(matrix(NA, nrow=1, ncol=ncol(kegg_enrichment@result) - 2))
      mrna_data <- if(length(mrna_id) > 0) kegg_enrichment_mrna@result[mrna_id, -c(1, 2)] else as.data.frame(matrix(NA, nrow=1, ncol=ncol(kegg_enrichment_mrna@result) - 2))
      ncrna_data <- if(length(ncrna_id) > 0) kegg_enrichment_ncrna@result[ncrna_id, -c(1, 2)] else as.data.frame(matrix(NA, nrow=1, ncol=ncol(kegg_enrichment_ncrna@result) - 2))
      
      colnames(bifunc_data) <- paste("bifunc", colnames(kegg_enrichment@result)[-c(1, 2)], sep = ".")
      colnames(mrna_data) <- paste("mrna", colnames(kegg_enrichment_mrna@result)[-c(1, 2)], sep = ".")
      colnames(ncrna_data) <- paste("ncrna", colnames(kegg_enrichment_ncrna@result)[-c(1, 2)], sep = ".")
      
      cbind(bifunc_data, mrna_data, ncrna_data)
    }))
    
    # Extract the ID and Description columns once from go_enrichment
    id_description <- kegg_enrichment@result[match(kegg_IDs, kegg_enrichment@result$ID), c("ID", "Description")]
    # Combine the ID and Description with the combined data
    final_combined_data <- cbind(id_description, combined_data)
    } else if (sigdb_type == "DAVID"){
    david_enrichment <- enrichDAVID(gene = genes_of_interest,
                                    universe = background_genes,
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.01,
                                    qvalueCutoff = 0.2)
    
    david_enrichment_mrna <- enrichDAVID(gene = coding_geneids,
                                         universe = background_genes,
                                         pAdjustMethod = "BH",
                                         pvalueCutoff = 0.01,
                                         qvalueCutoff = 0.2)
    
    david_enrichment_ncrna <- enrichDAVID(gene = noncoding_geneids,
                                          universe = background_genes,
                                          pAdjustMethod = "BH",
                                          pvalueCutoff = 0.01,
                                          qvalueCutoff = 0.2)
    david_IDs <- unique(c(david_enrichment@result$ID, 
                          david_enrichment_mrna@result$ID,
                          david_enrichment_ncrna@result$ID))
    
    
    
    # combining data from each gene type into one!
    print("Combining all gene categories into one dataset..")
    combined_data <- do.call(rbind, lapply(david_IDs, function(id) {
      bifunc_id = which(david_enrichment@result$ID == id)
      mrna_id = which(david_enrichment_mrna@result$ID == id)
      ncrna_id = which(david_enrichment_ncrna@result$ID == id)
      
      bifunc_data <- if(length(bifunc_id) > 0) david_enrichment@result[bifunc_id, -c(1, 2)] else as.data.frame(matrix(NA, nrow=1, ncol=ncol(david_enrichment@result) - 2))
      mrna_data <- if(length(mrna_id) > 0) david_enrichment_mrna@result[mrna_id, -c(1, 2)] else as.data.frame(matrix(NA, nrow=1, ncol=ncol(david_enrichment_mrna@result) - 2))
      ncrna_data <- if(length(ncrna_id) > 0) david_enrichment_ncrna@result[ncrna_id, -c(1, 2)] else as.data.frame(matrix(NA, nrow=1, ncol=ncol(david_enrichment_ncrna@result) - 2))
      
      colnames(bifunc_data) <- paste("bifunc", colnames(david_enrichment@result)[-c(1, 2)], sep = ".")
      colnames(mrna_data) <- paste("mrna", colnames(david_enrichment_mrna@result)[-c(1, 2)], sep = ".")
      colnames(ncrna_data) <- paste("ncrna", colnames(david_enrichment_ncrna@result)[-c(1, 2)], sep = ".")
      
      cbind(bifunc_data, mrna_data, ncrna_data)
    }))
    
    # Extract the ID and Description columns once from go_enrichment
    id_description <- david_enrichment@result[match(david_IDs, david_enrichment@result$ID), c("ID", "Description")]
    # Combine the ID and Description with the combined data
    final_combined_data <- cbind(id_description, combined_data)
    
  } else if (sigdb_type == "MKEGG"){
    mkegg_enrichment <- enrichMKEGG(gene = genes_of_interest,
                                    universe = background_genes,
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.01,
                                    qvalueCutoff = 0.2)
    
    mkegg_enrichment_mrna <- enrichMKEGG(gene = coding_geneids,
                                         universe = background_genes,
                                         pAdjustMethod = "BH",
                                         pvalueCutoff = 0.01,
                                         qvalueCutoff = 0.2)
    mkegg_IDs <- unique(c(mkegg_enrichment@result$ID, 
                          mkegg_enrichment_mrna@result$ID))

    # combining data from each gene type into one!
    print("Combining all gene categories into one dataset..")
    combined_data <- do.call(rbind, lapply(mkegg_IDs, function(id) {
      bifunc_id = which(mkegg_enrichment@result$ID == id)
      mrna_id = which(mkegg_enrichment_mrna@result$ID == id)
      
      bifunc_data <- if(length(bifunc_id) > 0) mkegg_enrichment@result[bifunc_id, -c(1, 2)] else as.data.frame(matrix(NA, nrow=1, ncol=ncol(mkegg_enrichment@result) - 2))
      mrna_data <- if(length(mrna_id) > 0) mkegg_enrichment_mrna@result[mrna_id, -c(1, 2)] else as.data.frame(matrix(NA, nrow=1, ncol=ncol(mkegg_enrichment_mrna@result) - 2))
    
      colnames(bifunc_data) <- paste("bifunc", colnames(mkegg_enrichment@result)[-c(1, 2)], sep = ".")
      colnames(mrna_data) <- paste("mrna", colnames(mkegg_enrichment_mrna@result)[-c(1, 2)], sep = ".")
      
      cbind(bifunc_data, mrna_data)
    }))
    
    # Extract the ID and Description columns once from go_enrichment
    id_description <- mkegg_enrichment@result[match(mkegg_IDs, mkegg_enrichment@result$ID), c("ID", "Description")]
    final_combined_data <- cbind(id_description, combined_data)
    
  } else if (sigdb_type == "PC"){
    pc_enrichment <- enrichPC(gene = genes_of_interest,
                              universe = background_genes,
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.01,
                              qvalueCutoff = 0.2)
    
    pc_enrichment_mrna <- enrichPC(gene = coding_geneids,
                                   universe = background_genes,
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.01,
                                   qvalueCutoff = 0.2)
    
    pc_enrichment_ncrna <- enrichPC(gene = noncoding_geneids,
                                    universe = background_genes,
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.01,
                                    qvalueCutoff = 0.2)
    pc_IDs <- unique(c(pc_enrichment@result$ID, 
                       pc_enrichment_mrna@result$ID,
                       pc_enrichment_ncrna@result$ID))
    
    
    
    # combining data from each gene type into one!
    print("Combining all gene categories into one dataset..")
    combined_data <- do.call(rbind, lapply(pc_IDs, function(id) {
      bifunc_id = which(pc_enrichment@result$ID == id)
      mrna_id = which(pc_enrichment_mrna@result$ID == id)
      ncrna_id = which(pc_enrichment_ncrna@result$ID == id)
      
      bifunc_data <- if(length(bifunc_id) > 0) pc_enrichment@result[bifunc_id, -c(1, 2)] else as.data.frame(matrix(NA, nrow=1, ncol=ncol(pc_enrichment@result) - 2))
      mrna_data <- if(length(mrna_id) > 0) pc_enrichment_mrna@result[mrna_id, -c(1, 2)] else as.data.frame(matrix(NA, nrow=1, ncol=ncol(pc_enrichment_mrna@result) - 2))
      ncrna_data <- if(length(ncrna_id) > 0) pc_enrichment_ncrna@result[ncrna_id, -c(1, 2)] else as.data.frame(matrix(NA, nrow=1, ncol=ncol(pc_enrichment_ncrna@result) - 2))
      
      colnames(bifunc_data) <- paste("bifunc", colnames(pc_enrichment@result)[-c(1, 2)], sep = ".")
      colnames(mrna_data) <- paste("mrna", colnames(pc_enrichment_mrna@result)[-c(1, 2)], sep = ".")
      colnames(ncrna_data) <- paste("ncrna", colnames(pc_enrichment_ncrna@result)[-c(1, 2)], sep = ".")
      
      cbind(bifunc_data, mrna_data, ncrna_data)
    }))
    
    # Extract the ID and Description columns once from go_enrichment
    id_description <- pc_enrichment@result[match(pc_IDs, pc_enrichment@result$ID), c("ID", "Description")]
    final_combined_data <- cbind(id_description, combined_data)
    
  } else if (sigdb_type == "WP"){
    wp_enrichment <- enrichWP(gene = genes_of_interest,
                              universe = background_genes,
                              organism = "Homo sapiens",
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.01,
                              qvalueCutoff = 0.2)
    
    wp_enrichment_mrna <- enrichWP(gene = coding_geneids,
                                   universe = background_genes,
                                   organism = "Homo sapiens",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.01,
                                   qvalueCutoff = 0.2)
    
    wp_enrichment_ncrna <- enrichWP(gene = noncoding_geneids,
                                    universe = background_genes,
                                    organism = "Homo sapiens",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.01,
                                    qvalueCutoff = 0.2)
    wp_IDs <- unique(c(wp_enrichment@result$ID, 
                       wp_enrichment_mrna@result$ID,
                       wp_enrichment_ncrna@result$ID))
    
    
    
    # combining data from each gene type into one!
    print("Combining all gene categories into one dataset..")
    combined_data <- do.call(rbind, lapply(wp_IDs, function(id) {
      bifunc_id = which(wp_enrichment@result$ID == id)
      mrna_id = which(wp_enrichment_mrna@result$ID == id)
      ncrna_id = which(wp_enrichment_ncrna@result$ID == id)
      
      bifunc_data <- if(length(bifunc_id) > 0) wp_enrichment@result[bifunc_id, -c(1, 2)] else as.data.frame(matrix(NA, nrow=1, ncol=ncol(wp_enrichment@result) - 2))
      mrna_data <- if(length(mrna_id) > 0) wp_enrichment_mrna@result[mrna_id, -c(1, 2)] else as.data.frame(matrix(NA, nrow=1, ncol=ncol(wp_enrichment_mrna@result) - 2))
      ncrna_data <- if(length(ncrna_id) > 0) wp_enrichment_ncrna@result[ncrna_id, -c(1, 2)] else as.data.frame(matrix(NA, nrow=1, ncol=ncol(wp_enrichment_ncrna@result) - 2))
      
      colnames(bifunc_data) <- paste("bifunc", colnames(wp_enrichment@result)[-c(1, 2)], sep = ".")
      colnames(mrna_data) <- paste("mrna", colnames(wp_enrichment_mrna@result)[-c(1, 2)], sep = ".")
      colnames(ncrna_data) <- paste("ncrna", colnames(wp_enrichment_ncrna@result)[-c(1, 2)], sep = ".")
      
      cbind(bifunc_data, mrna_data, ncrna_data)
    }))
    
    # Extract the ID and Description columns once from go_enrichment
    id_description <- wp_enrichment@result[match(wp_IDs, wp_enrichment@result$ID), c("ID", "Description")]
    # Combine the ID and Description with the combined data
    final_combined_data <- cbind(id_description, combined_data)
  }
  
  # Convert the list to a data frame
  final_combined_data <- as.data.frame(final_combined_data)
  
  write.csv(final_combined_data, file=paste0("results_", sigdb_type,".csv"))
  
  # Function to create plot data
  create_plot_data <- function(significant_data, type) {
    plot_dat <- data.frame(Description=character(), type=character(), 
                           FoldEnrichment=numeric(), padj=numeric())
    if(length(grep("GO", sigdb_type))==1){
    for(row in rownames(significant_data)) {
      for (gene_type in c("bifunc", "mrna", "ncrna")) {
        ora_dat <- significant_data[row, grep(gene_type, colnames(significant_data))]
        add_dat <- c(significant_data[row, "Description"], gene_type, ora_dat[c(4, 7)])
        names(add_dat) <- c("Description", "type", "FoldEnrichment", "padj")
        plot_dat <- rbind(plot_dat, add_dat, make.row.names = FALSE)
      }
    }} else if(sigdb_type=="KEGG"){
    for(row in rownames(significant_data)) {
      for (gene_type in c("bifunc", "mrna", "ncrna")) {
        ora_dat <- significant_data[row, grep(gene_type, colnames(significant_data))]
        add_dat <- c(significant_data[row, "Description"], gene_type, ora_dat[c(6, 9)])
        names(add_dat) <- c("Description", "type", "FoldEnrichment", "padj")
        plot_dat <- rbind(plot_dat, add_dat, make.row.names = FALSE)
      }
    }} else if (sigdb_type=="WP"){
      for(row in rownames(significant_data)) {
        for (gene_type in c("bifunc", "mrna", "ncrna")) {
          ora_dat <- significant_data[row, grep(gene_type, colnames(significant_data))]
          add_dat <- c(significant_data[row, "Description"], gene_type, ora_dat[c(4, 7)])
          names(add_dat) <- c("Description", "type", "FoldEnrichment", "padj")
          plot_dat <- rbind(plot_dat, add_dat, make.row.names = FALSE)
        }
      }} else if(sigdb_type=="MKEGG"){
        for(row in rownames(significant_data)) {
          for (gene_type in c("bifunc", "mrna")) {
            ora_dat <- significant_data[row, grep(gene_type, colnames(significant_data))]
            add_dat <- c(significant_data[row, "Description"], gene_type, ora_dat[c(4, 7)])
            names(add_dat) <- c("Description", "type", "FoldEnrichment", "padj")
            plot_dat <- rbind(plot_dat, add_dat, make.row.names = FALSE)
          }
        }}
    
    return(plot_dat)
  } 
  
  # Function to create dot plot
  create_dot_plot <- function(plot_dat, title, type) {
    #wrap description labels since they take up too much space!
    plot_dat$Description_wrapped <- factor(stringr::str_wrap(plot_dat$Description, width = 45), levels = unique(stringr::str_wrap(plot_dat$Description, width = 45)))
    #plot
    dotplot <- ggplot(plot_dat, aes(y = Description_wrapped, colour = as.factor(type), x = FoldEnrichment, size = -log10(padj))) +
      geom_point() + xlim (0, max(plot_dat$FoldEnrichment)*1.2)+ 
      theme(
        plot.background = element_blank(),
        axis.text = element_text(color = 'black', face = "bold", size=12),
        axis.title = element_text(color = 'black', face = "bold", size=13),
        panel.border = element_rect(linewidth = 1, fill = NA),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "lightgrey"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title= element_text(size=13, color="black", face="bold"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.box = "vertical",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.justification = c(1, 0)
      ) +
      scale_y_discrete() +
      scale_colour_manual(labels = c("Bifunctional", "Coding", "Noncoding"),
        values = c( "#82c0cc","#e8e288", "#ffa5a5")) +
      labs(title = title,
        x = "Fold Enrichment",
        y = "Description",
        size = "-log10(p-value)",
        colour = "Gene Type") + 
      guides(color = guide_legend(nrow = 1),
             shape = guide_legend(nrow = 1))
    
    #save the plot
    #svg(paste0(title, "_", type, ".svg"))
    plot(dotplot)
    #dev.off()
    print("plotted another dotplot..")
  }
  
  if(sigdb_type!="MKEGG"){
  # Plot the data!
  significant_data <- final_combined_data[!is.na(final_combined_data$bifunc.p.adjust) & final_combined_data$bifunc.p.adjust < 0.05, ]
  plot_dat <- create_plot_data(significant_data, "bifunc")
  create_dot_plot(plot_dat, title=sigdb_type, type="bifunc")
  
  # mRNA plot
  #keep the sets with p-value < 0.05 for m-RNA enoding genes!
  significant_data2 <- na.omit(final_combined_data[!is.na(final_combined_data$mrna.p.adjust) & final_combined_data$mrna.p.adjust < 0.05, ])
  #keeping only top 20 genesets with the highest enrcihment score for mRNA-coding genes!
  idxs <- order(significant_data2$mrna.FoldEnrichment, decreasing=TRUE)
  significant_data2 <- na.omit(significant_data2[idxs[1:20], ])
  plot_dat2 <- create_plot_data(significant_data2, "mrna")
  create_dot_plot(plot_dat2, sigdb_type, type="mrna")  
  
  # ncRNA plot
  #keep the sets with p-value < 0.05 for ncRNA encoding genes!
  significant_data3 <- na.omit(final_combined_data[!is.na(final_combined_data$mrna.p.adjust) & final_combined_data$ncrna.p.adjust < 0.05, ])
  #keeping only top 20 gene-sets with the highest enrichment score for ncRNA-coding genes!
  idxs <- order(significant_data3$ncrna.FoldEnrichment, decreasing=TRUE)
  significant_data3 <- na.omit(significant_data3[idxs[1:20], ])
  plot_dat3 <- create_plot_data(significant_data3, "ncrna")
  create_dot_plot(plot_dat3, sigdb_type, type="ncrna") 
  
  } else {
    # Plot the data!
    significant_data <- final_combined_data[!is.na(final_combined_data$bifunc.p.adjust) & final_combined_data$bifunc.pvalue < 0.05, ]
    plot_dat <- create_plot_data(significant_data, "bifunc")
    create_dot_plot(plot_dat, title=sigdb_type, type="bifunc")
    
    # mRNA plot
    #keep the sets with p-value < 0.05 for m-RNA enoding genes!
    significant_data2 <- na.omit(final_combined_data[!is.na(final_combined_data$mrna.p.adjust) & final_combined_data$mrna.pvalue < 0.05, ])
    #keeping only top 25 genesets with the highest enrcihment score for mRNA-coding genes!
    idxs <- order(significant_data2$mrna.FoldEnrichment, decreasing=TRUE)
    significant_data2 <- na.omit(significant_data2[idxs[1:20], ])
    plot_dat2 <- create_plot_data(significant_data2, "mrna")
    create_dot_plot(plot_dat2, sigdb_type, type="mrna")
  }
  
}

dir.create("./figures/enrichment")
setwd("./figures/enrichment/")
ora_type(sigdb_type="GO_BP")
gc()
ora_type(sigdb_type="GO_MF")
gc()
ora_type(sigdb_type="GO_CC")
gc()
ora_type(sigdb_type="KEGG") 
gc()
ora_type(sigdb_type="MKEGG") 
gc()
#ora_type(sigdb_type="DAVID") # can't run because one of the required packages have been removed by Bioconductor!
#gc()
ora_type(sigdb_type="WP") 
gc()
#ora_type(sigdb_type="PC") #needs gene symbols?? not doing for now!
#gc()

sessionInfo()
