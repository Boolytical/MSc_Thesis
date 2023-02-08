#Function: iqi_source
#Takes: data
#Returns: plot
#Calls: NA
#Job: create a historgram of the IQI values per source, make a boxplot of the IQI values
#       per sample site coloured by source, and make a barplot of the mean IQI values
#       per sample site coloured by source.
iqi_source = function(data){
  IQI_per_site = aggregate(data$IQI, list(data$Site, data$Source), mean)
  colnames(IQI_per_site) = c("Site", "Source", "IQI")
  
  hist_IQI = ggplot(data, aes(x = IQI, group = Source, fill = Source, colour = Source)) +
    geom_density(alpha = 0.5) +
    ggtitle("(a)") +
    labs(x = 'IQI', y = 'Density') +
    scale_color_brewer(palette = 'Paired', aesthetics = c("colour", "fill")) +
    geom_vline(xintercept = 0.75, color = "green",  linetype = 'dashed') +
    geom_vline(xintercept = 0.64, color = "green",  linetype = 'dashed') +
    geom_vline(xintercept = 0.44, color = "green",  linetype = 'dashed') +
    geom_vline(xintercept = 0.24, color = "green",  linetype = 'dashed') +
    theme(plot.title = element_text(size = 12)) 
  
  box_IQI = ggplot(data = data, aes(x = reorder(Site, as.numeric(factor(Source))), y = IQI, fill = Source)) +
    geom_boxplot() +
    ggtitle("(b)") +
    labs(x = 'Site', y = 'IQI') +
    guides(fill = guide_legend(nrow = 8, byrow = T)) +
    scale_color_brewer(palette = 'Paired', aesthetics = c("colour", "fill")) +
    geom_hline(yintercept = 0.75, color = "green",  linetype = 'dashed') +
    geom_hline(yintercept = 0.64, color = "green",  linetype = 'dashed') +
    geom_hline(yintercept = 0.44, color = "green",  linetype = 'dashed') +
    geom_hline(yintercept = 0.24, color = "green",  linetype = 'dashed') +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.title = element_text(size = 12)) 
  
  bar_IQI = ggplot(data = IQI_per_site, aes(x = reorder(Site, as.numeric(factor(Source))), y = IQI, fill = Source)) +
    geom_bar(stat = "identity", position =  position_dodge(width = 0.9, preserve = "total")) +
    ggtitle("(c)") +
    labs(x = 'Site', y = 'IQI') +
    scale_color_brewer(palette = 'Paired', aesthetics = c("colour", "fill")) +
    geom_hline(yintercept = 0.75, color = "green",  linetype = 'dashed') +
    geom_hline(yintercept = 0.64, color = "green",  linetype = 'dashed') +
    geom_hline(yintercept = 0.44, color = "green",  linetype = 'dashed') +
    geom_hline(yintercept = 0.24, color = "green",  linetype = 'dashed') +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.title = element_text(size = 12))
  
  ggarrange(hist_IQI, box_IQI, bar_IQI, nrow = 3, ncol = 1, common.legend = T, legend = "bottom") %>%
    ggexport(filename = "figures/IQI_per_source.pdf", width = 7, height = 12)
}


#Function: rarefaction
#Takes: NA
#Returns: plot, summary on taxa reads
#Calls: NA
#Job: make a histogram of the read count per site coloured by source and return the diagnostics
rarefaction = function(MDS4){
  # Determine the number of reads per sample site
  data_taxa = MDS4[ , grepl("Bacteria", names(MDS4))]
  taxa_reads = rowSums(data_taxa)
  taxa_reads = data.frame("Reads" = taxa_reads, "Source" = MDS4$Source, "IQI" = MDS4$IQI)
  taxa_reads$Source = factor(taxa_reads$Source, levels=rev(sort(unique(taxa_reads$Source))))
  
  # Make a summary data frame
  reads_summary = taxa_reads %>%
    group_by(Source) %>%
    dplyr::summarise(min = min(Reads),
                     max = max(Reads))
  
  # Make a violin plot to show the frequency of reads per sample site over the sources
  violin_reads = ggplot(taxa_reads, aes(Reads, Source)) +
    geom_violin(aes(group = Source, fill = Source, colour = Source, alpha = 0.5)) +
    geom_text(data = reads_summary, aes(x = max, y = Source, label = max),
              size = 3, hjust = -0.5) +
    geom_text(data = reads_summary, aes(x = min, y = Source, label = min),
              size = 3, hjust = 1.5) +
    scale_color_brewer(palette = 'Paired', aesthetics = c("colour", "fill"), direction = -1) + 
    theme(legend.position = "none") +
    scale_x_continuous(limits = c(- 10000, max(taxa_reads$Reads) + 25000))

  violin_reads %>%
    ggexport(filename = paste("figures/rarefaction_read_counts_", toString(unique(MDS4$Source)), ".pdf", sep = ''))
  
  # Make rarefaction curves of the ASV reads per sample coloured by Source
  data_taxa$Source = as.numeric(as.factor(MDS4$Source))
  
  pdf(paste("figures/rarefaction_curves_", toString(unique(MDS4$Source)), ".pdf", sep = ''), height=12)
  par(mfrow = c(2, 1), mar = c(4, 4, 4, 1))
  rarecurve(subset(data_taxa, select = -Source), step = 50, col = data_taxa$Source, label = FALSE, xlim = c(0, 100000))
  title("(a)", adj = 0)
  rarecurve(subset(data_taxa, select = -Source), step = 50, col = data_taxa$Source, label = FALSE, xlim = c(0, 15000))
  title("(b)", adj = 0)
  dev.off()
  
  write.csv(taxa_reads %>% arrange(Reads) %>% head(n = 15),
            file = paste("figures/rarefaction_read_counts_head_Source_", toString(unique(MDS4$Source)), ".csv", sep = ''), row.names = TRUE)
  
  return(taxa_reads %>% arrange(Reads) %>% head(n = 15))
}

#Function: alpha_diversity
#Takes: data
#Returns: plot
#Calls: NA
#Job: create boxplots of the Chao1 and Shannon diversity indexes per site coloured by source
alpha_diversity = function(data){
  data_taxa = data[ , grepl("Bacteria", names(data))]
  
  # Calculate Chao1 and Shannon index per sample site
  chao1 = c()
  for (i in 1:nrow(data)){
    estimate = estimateR(as.integer(data_taxa[i, ]))
    chao1[i] = estimate[2]
  }
  shannon = diversity(data_taxa, index = "shannon")
  indices = data.frame("Chao1" = chao1,
                       "Shannon" = shannon,
                       "Source" = data$Source,
                       "Site" = data$Site)
  
  # Make a boxplot per index coloured by site
  box_chao1 = ggplot(data = indices, aes(x = reorder(Site, as.numeric(factor(Source))), y = Chao1, fill = Source)) +
    geom_boxplot() +
    labs(x = "Site", y = "Chao1") +
    guides(fill = guide_legend(nrow = 1, byrow = T)) +
    scale_color_brewer(palette = "Paired", aesthetics = c("colour", "fill")) +
    theme(legend.position = "bottom",
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 12))
  
  box_shannon = ggplot(data = indices, aes(x = reorder(Site, as.numeric(factor(Source))), y = Shannon, fill = Source)) +
    geom_boxplot() +
    labs(x = "Site", y = "Shannon") +
    guides(fill = guide_legend(nrow = 1, byrow = T)) +
    scale_color_brewer(palette = "Paired", aesthetics = c("colour", "fill")) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.title = element_text(size = 12))
  
  ggarrange(box_chao1, box_shannon, nrow = 2, ncol = 1, common.legend = T, legend = "bottom",
            widths = c(7, 7), heights = c(3.5, 6)) %>%
    ggexport(filename = "figures/alpha_diversity.pdf")
}


#Function: beta_diversity
#Takes: data
#Returns: plot
#Calls: NA
#Job: create a Bray-Curtis pairwise comparison graph to show the beta diversity per site.
#     make a PCoA plot based on this matrix to display the site differences including IQI indicator.
beta_diversity = function(data){
  # Select the bacterial taxa and add the IQI value
  taxa_columns = data[ , grepl("Bacteria", names(data))]
  taxa_columns$IQI = as.numeric(data$IQI)
  
  # Generate data set for Bray-Curtis
  data_bc = aggregate(taxa_columns, list(data$Site), mean) # Take bacterial and IQI mean per Site
  rownames(data_bc) = data_bc$Group.1; data_bc$Group.1 = NULL
  
  ## Bray-Curtis pairwise distancing test on sites ##
  bray_curtis = bcdist(data_bc[, -ncol(data_bc)]) # Exclude IQI for Bray-Curtis analysis
  
  pdf("figures/beta_bray_curtis.pdf")
  par(mar = c(5, 4, 1, 1))
  distogram(bray_curtis, colFn = colorRampPalette(brewer.pal(n = 9, 'Blues')))
  dev.off()
  
  ## Perform PCoA analysis based on Bray-Curtis distances; per Site ##
  pcoa = cmdscale(bray_curtis)
  colnames(pcoa) = c("PCoA1", "PCoA2")
  pcoa_data = as.data.frame(pcoa)
  
  # Determine effect of IQI
  fit_iqi = envfit(pcoa, data_bc[, ncol(data_bc)])
  fit_data = as.data.frame(scores(fit_iqi, display = 'vectors'))
  fit_data = cbind(fit_data, Taxa = colnames(data_bc)[ncol(data_bc)])
  
  # Make the PCoA plot with a vector for the IQI effect
  fit_plot_site = ggplot(pcoa_data, aes(x = PCoA1, y = PCoA2, label = rownames(pcoa_data))) +
    coord_fixed(xlim = c(-0.50, 0.50), ylim = c(-0.50, 0.50)) +
    geom_text_repel(max.overlaps = 1000, show.legend = FALSE) +
    geom_segment(aes(x = 0, xend = fit_data$PCoA1, y = 0, yend = fit_data$PCoA2),
                 arrow = arrow(length = unit(0.25, "cm")), colour = "green") +
    geom_text(aes(x = fit_data$PCoA1 - 0.02, y = fit_data$PCoA2 + 0.02,
                  label = fit_data$Taxa), size = 5, col = "green")
  
  fit_plot_site %>%
    ggexport(filename = "figures/beta_pcoa.pdf")
}

