# for simple heatmaps
library(gplots)
library(viridis)
library(vegan)

library(dendextend) # for cutting tree

library("ComplexHeatmap")

library(ggrepel) # for labeling with geom_text_repel
library(ggforce) # for geom_mark_ellipse

library(gridExtra)
library(RColorBrewer)

# build Amundsen ASV by running AMUN_BIOINFO_ASV

meta.amun.small_nonseq_below100 <- meta.amun.small %>%
  filter (seq.filt == 0, `Depth [m]` > 100)

proptab.small.taxan.saman_nonseq_below100 <- proptab.small.taxan.saman %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  mutate(plate.geno = sapply(strsplit(Sample, " ; "), tail, 1)) %>%     # Extract last part as plate.geno
  filter(plate.geno %in% meta.amun.small_nonseq_below100$plate.geno) %>%  # Keep only matching samples
  select(-plate.geno) %>%                                               # Remove plate.geno column
  column_to_rownames("Sample") %>%                                      # Restore row names after filtering
  .[colSums(. > 0) > 2]                                                 # select species present in more than 2 samples

#
### What does your ASV table look like in terms of data distribution? What normalization should you apply?
hist(as.matrix(proptab.small.taxan.saman_nonseq_below100))         # looks pretty bad
hist(as.matrix(proptab.small.taxan.saman_nonseq_below100)^0.5)     # looks better
hist(as.matrix(proptab.small.taxan.saman_nonseq_below100)^0.25)    # looks optimal: a lot of 0s but normal distribution to the left


#####################
##### Simple Heatmap

#
### build a simple heatmap to explore all the data

## open pdf
pdf(file="OUTPUT/amun_simple_hm_nonseq_below100.pdf", width=35, height=20)

## simple heatmap
hm = heatmap.2(as.matrix(proptab.small.taxan.saman_nonseq_below100)^.25, trace="none", col=viridis,
               distfun=function(x) vegdist(x, method="bray"),
               hclustfun=function(x) hclust(x, method="ward.D2"))

## see how many clusters you may have for samples (rows)
clusnum=7
oldclus = cutree(hm$rowDendrogram, k=clusnum)
oldorder = unname(rle(oldclus[as.hclust(hm$rowDendrogram)$order])$values)
neworder = (1:clusnum)
names(neworder) = oldorder
clus = unname(neworder[as.character(oldclus)])
names(clus) = names(oldclus)

## redraw the map with clusters colored
hm = heatmap.2(as.matrix(proptab.small.taxan.saman_nonseq_below100)^.25, trace="none", key=F, margins=c(20,30), col=viridis,
               distfun=function(x) vegdist(x, method="bray"),
               hclustfun=function(x) hclust(x, method="ward.D2"),
               RowSideColors = rainbow(max(clus))[clus])

## close the pdf
dev.off()


#######################
##### ComplexHeatmap

# colnames(meta.amun.small) # see which meta parameters you can investigate

#
# Compute Bray-Curtis distance and hierarchical clustering
tdist <- vegdist(as.matrix(proptab.small.taxan.saman_nonseq_below100)^.25, method = "bray")
ttree <- hclust(tdist, method = "ward.D2")

# Assign cluster numbers
tclus <- cutree(ttree, k = 4)

# Order clusters so they appear in order from top to bottom
ordered_samples <- names(sort(tclus))  # Get sample names in correct order

# Reorder the data matrix by the new order
ordered_matrix <- (as.matrix(proptab.small.taxan.saman_nonseq_below100)^.25)[ordered_samples, ]
ordered_clusters <- tclus[ordered_samples]  # Match clusters to reordered samples

#
## open pdf
pdf(file="OUTPUT/amun_complex_hm_nonseq_below100.pdf", width=50, height=20)

ht_opt$ROW_ANNO_PADDING = unit(2, "mm")                                          # Padding between the heatmap and annotations

#
## create a complex heatmap
amun_complex_hm_nonseq_below100 <- 
  
  Heatmap(ordered_matrix,
          
          # clustering_distance_rows = bray_dist,
          clustering_distance_columns = bray_dist,
          # clustering_method_rows = "ward.D2",                                  #passed to hclust
          clustering_method_columns = "ward.D2",                                 #passed to hclust
          
          name = "amun_complex_hm_nonseq_below100",
          col = blues9,
          
          right_annotation = rowAnnotation(bar1 = anno_points(log2(meta.amun.small_nonseq_below100$"Depth [m]"),
                                                              width = unit(3, "cm"),
                                                              axis_param = list(
                                                                gp = gpar(fontsize = 8))),
                                           bar2 = anno_points(meta.amun.small_nonseq_below100$"Temperature [deg C]",
                                                              width = unit (3, "cm"),
                                                              axis_param = list(
                                                                gp = gpar(fontsize = 8))),
                                           bar3 = anno_points(meta.amun.small_nonseq_below100$"Salinity",
                                                              width = unit (3, "cm"),
                                                              axis_param = list(
                                                                gp = gpar(fontsize = 8))),
                                           bar4 = anno_points(meta.amun.small_nonseq_below100$"Dissolved Oxygen [uM]",
                                                              width = unit (3, "cm"),
                                                              axis_param = list(
                                                                gp = gpar(fontsize = 8))),
                                           bar5 = anno_points(meta.amun.small_nonseq_below100$"Fluorescence [mg/m^3]",
                                                              width = unit (3, "cm"),
                                                              axis_param = list(
                                                                gp = gpar(fontsize = 8))),
                                           bar6 = anno_points(meta.amun.small_nonseq_below100$"Nitrate [uM]",
                                                              width = unit (3, "cm"),
                                                              axis_param = list(
                                                                gp = gpar(fontsize = 8))),
                                           gap = unit(2, "mm"),                  # Gap between annotations
                                           show_annotation_name = FALSE),        # Assign label later through decoration
          
          cluster_rows = FALSE,                                                  # To do clustering yourself (see ttree)
          row_split = ordered_clusters,                                          # Split in the right order
          # row_dend_reorder = TRUE,                                             # Rotate tree (if you don't do clustering yourself)
          
          column_dend_reorder = TRUE,
          
          heatmap_height = unit(40, "cm"),
          heatmap_width = unit(65, "cm"),
          row_dend_width = unit(4, "cm"),
          column_dend_height = unit(2, "cm"),
          
          row_title = "Samples",
          row_title_gp = gpar(fontsize = 25, fontface = "bold"),
          # row_labels = rep("", nrow(proptab.small.taxan.saman_nonseq_below100)),    # If you want to remove row names (aka samples)
          
          column_title = "Taxa",
          column_title_gp = gpar(fontsize = 25, fontface = "bold"),
          
          row_names_gp = gpar(fontsize = 25, fontface = "bold"),
          column_names_gp = gpar(fontsize = 1.2, fontface = "bold"),
          column_names_rot = 45,
          
          heatmap_legend_param = list(
            title = "Fourth-root transformed relative abundances of microbial taxa during the 2023 Amundsen expedition
(non-sequential filtration, below 100 meters)\n",
            title_gp = gpar(fontsize = 30, fontface = "bold"),
            legend_direction = "horizontal",
            legend_width = unit(2, "cm"),
            color_bar = "continuos"),
          
          # show_heatmap_legend = FALSE                                # Turn off heatmap legend b/c numbers are meaningless  
  )

#
## Draw the complex heatmap
draw(amun_complex_hm_nonseq_below100,                                       # First draw, then decorate, place legend in a good spot
     heatmap_legend_side="top")                                           # Padding between row annotations

#
## Add decorations

# Define parameter labels for the annotations (see bar1, bar2, etc.)
annotation_labels <- c(
  "Depth, m",
  "Temperature, deg C",
  "Salinity",
  "Oxygen, uM",
  "Fluorescence, ug/l",
  "Nitrate, uM"
)

# Loop over annotation bars
for (i in 1:6) {
  decorate_annotation(paste0("bar", i), {
    grid.text(annotation_labels[i],
              y = unit(1, "npc") + unit(2, "mm"), just = "bottom",               # Position label above the annotation
              gp = gpar(fontsize = 8, fontface = "bold"))
  })
}

#
##

# Generate k number of mild colors for your clusters
cluster_palette <- brewer.pal(min(4, 4), "Set1")  # set 4 distinct colors to be generated

# Loop over cluster slices
for (i in 1:4) {
  decorate_row_names("amun_complex_hm_nonseq_below100", {
    grid.rect(gp = gpar(fill = paste0(cluster_palette[i], "80")),  # Apply mild transparency
              width = unit(13.3, "cm"),                              # Width of the colored box
              x = unit(1, "npc") + unit(0.6, "cm"),                  # Adjust alignment
              just = "center")
  }, slice = i)
}

## close the pdf
dev.off()


######## Check out stats for the cluster order
proptab.small.taxan.saman_nonseq_below100 %>%
  mutate(Cluster = tclus[rownames(.)]) %>%  # Add cluster column
  cbind(meta.amun.small_nonseq_below100) %>%
  group_by(Cluster) %>%
  summarize(min_depth = min(`Depth [m]`, na.rm = TRUE),
            max_depth = max(`Depth [m]`, na.rm = TRUE),
            mean_depth = mean(`Depth [m]`, na.rm = TRUE))


########################################
##### db-RDA (based on Bray Curtis Distance matrix)
# Load microbial ASV data (response variables)
microbial_data <- as.matrix(proptab.small.taxan.saman_nonseq_below100)^0.25

# Check how your ASV data looks like
# If the distribution is not normal - go back and modify.
# Here we have already done prop.table (by sample) and ^0.25 transformations on the original ASV matrix
hist(microbial_data)    # looks optimal: a lot of 0s but normal distribution to the left

# Load metadata (explanatory variables)
metadata <- meta.amun.small_nonseq_below100

# Select numeric environmental variables correctly
env_vars <- metadata %>%
  select(where(is.numeric)) %>%
  
  # interpolate NO3 w/ CTD readings if you are planning to use all the samples at the expense of removing other nutrients
  # not used here as the variance explained and adjusted R2 drops substantially
  mutate(`Nitrate [uM]` = if_else(is.na(`Nitrate [uM]`), `Nitrate CTD [uM]`, `Nitrate [uM]`)) %>%
  
  select (-"Cast",                                        # meaningless
          -"Latitude [deg]", -"Longitude [deg]",          # meaningless + colored by location
          - "Year", -"Month", -"Day",                     # work with DOY instead if needed
          -"Bot. Depth [m]",                              # meaningless + colored by location
          -"Sea Pressure [db]",                           # correlates with depth
          -"Bottle",                                      # meaningless
          -"prefilter.size", -"seq.filt", -"pore.size",   # will display seq.filt and pore.size on the plot
          -"volume.l", -"dna.conc",                       # do not matter for relative abundance in a sample
          -"plate",                                       # meaningless
          -"Sigma-T [kg/m^3]", -"Sigma-Theta [kg/m^3]",   # correlates with salinity
          -"Sound Velocity [m/s]",                        # meaningless
          -"PAR [uE/s/m^2]",                              # include light attenuation instead
          
          -"DOY",                                         # it has a big effect b/c two cruises were temporally separate
          
          -"Oxygen solubility [ml/l]",                    # correlates with temperature based on ggpair()
          -"Nitrate CTD [uM]",                            # correlates with bottle nitrate based on ggpair()
          
          -"Dissolved Oxygen [uM]",                       # has high VIF and correlates with nitrate (ggpair())
          -"Silicate [uM]",                               # has high VIF and correlates with nitrate and phosphate (ggpair())
          
          # -"Light Attenuation Coefficient [1/m]",         # has weak significance (*) based on anova.cca
          -"Light Transmission [%]",                      # has weak significance based on anova.cca
          -"Fluorescence [mg/m^3]",                       # has weak significance based on anova.cca (can use chl a later?)
          # -"CDOM [mg/m^3]",                               # has weak significance (**) based on anova.cca
          
          -"Phosphate [uM]",                              # high correlation with nitrate
          -"Nitrite [uM]"                                 # insignificant, too many NAs and can't interpolate
  ) %>%
  rename_with(~ str_replace_all(., "\\s*\\[.*?\\]", "")) %>%                    # Remove brackets
  rename_with(~ gsub(" ", "_", .))                                              # Replace spaces with underscores

# Identify the variables that strongly correlate with each other according to Pearson cor. coef. (R >= ~0.9) - remove one of them
# If you leave two strongly correlating environmental variables, it does not add much to the total variation explained b/c 
# the corr. var-le appearing in the model formula first already explains the variance the other corr. var-le also tries to explain
# It leads to inflated variance inflation factors (VIF) for both corr. var-les: it is hard to know which var-le drives the variation
ggpairs(env_vars, 
        lower = list(continuous = wrap("smooth", method = "lm", se = TRUE, color = "blue", size = 0.8)), 
        diag = list(continuous = wrap("barDiag", bins = 20, color = "cyan")),  # Histogram on the diagonal
        upper = list(continuous = wrap("cor", size = 3, color = "black")))     # Correlation in the upper panel

# Leave only the data where all parameters are non-NA (complete cases), unless you can find a way to harmlessly interpolate
# Check if the data was evenly removed
complete_cases <- complete.cases(env_vars)                      # here, rivers were removed
microbial_data <- microbial_data[complete_cases, , drop = FALSE]
env_vars <- env_vars[complete_cases, , drop = FALSE]

# The environmental/explanatory variables are all in different units and need to be standardized prior to computing distance measures to perform most ordination analyses
env_stand <- decostand(env_vars, method="standardize")  # subtract the mean and divide by the standard deviation
# this centers each variable around 0 (mean = 0) and scales it so that its standard deviation is 1
apply(env_stand, 2, mean) # check if the data are now centered (means~0)
apply(env_stand, 2, sd)   # check if the data are now scaled (standard deviations=1)

# Before applying db-RDA, you can check which dissimilarity index works best for you
rankindex(env_stand, microbial_data, indices = c("euc", "man", "gow", "bra", "kul"),
          stepacross= FALSE, method = "spearman")

# Perform db-RDA (distance-based redundancy analysis) with all variables
rda_model <- dbrda(microbial_data ~ ., data = env_stand, distance = "bray", dfun = vegdist)

# Check if any variables significantly correlate (VIF>10~20) in your constrained ordination (e.g., RDA)
# If yes, remove the ones that makes sense to remove above and rerun
vif.cca(rda_model) %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%  # Convert row names to a proper column
  rename(VIF = ".") %>%               # Rename the unnamed column to "VIF"
  mutate(VIF = round(VIF, 0)) %>%     # Round to 0 decimal places
  arrange(desc(VIF)) %>%              # Sort in descending order for easy review)
  print()

# Perform per-variable ANOVA to extract significant constraints:
# anova.cca performs an ANOVA like permutation test for cca, rda, or dbRDA to assess the significance of constraints
# Leave only significant environmental variables (high F value and p < 0.01** or p < 0.001***), remove the rest above and rerun
# High F value: the values within each group do not differ significantly from one another, while those across different groups do
anova_per_var <- anova.cca(rda_model, permutations = 1000, by = "term") %>%
  rownames_to_column("Variable") %>%                        # Convert rownames to a column
  select(Variable, F, `Pr(>F)`) %>%                         # Keep relevant columns
  rename(`P-value` = `Pr(>F)`) %>%                          # Rename Pr(>F) to p-value
  filter(Variable != "Residual") %>%                        # Remove "Residual" row
  mutate(
    F = round(F, 2),  
    `P-value` = sprintf("%.6f", as.numeric(`P-value`)),     # Ensure 6 decimal places
    `Significance code` = case_when(
      as.numeric(`P-value`) <= 0.001 ~ "***",
      as.numeric(`P-value`) <= 0.01  ~ "**",
      as.numeric(`P-value`) <= 0.05  ~ "*",
      as.numeric(`P-value`) <= 0.1   ~ ".",
      TRUE                           ~ "")) %>%             # Assign stars based on p-value significance
  arrange(desc(F)) %>%                                      # Sort by descending F-value
  as.data.frame() %>%                                       # Convert to data frame
  print()                                                   # Explicitly print the table

# Print variance explained by axes
perc <- round(100 * (summary(rda_model)$cont$importance[2, 1:2]), 2) %>%
  print()

# See if your percentages are OK, maybe you excluded too many variables after anova.cca

# Print proportion of variance explained by environmental variables
variance_partitioning <- (summary(rda_model)$constr.chi / summary(rda_model)$tot.chi) %>%
  {paste("variance explained:", round(. * 100, 2), "%")} %>%
  print()

# Print final adjusted R²
rda_adj_r2 <- RsquareAdj(rda_model)$adj.r.squared %>%
  {paste("adjusted R²:", round(. * 100, 1), "%")} %>%
  print()

# Make sure the model is significant (p-value < 0.001), and check if variable placement matters

# Print significance of the model.
anova.cca(rda_model, permutations = 1000)

# Print marginal effects of terms
anova.cca(rda_model, permutations = 1000, by = "margin")

# Extract RDA scores for plotting
rda_scores <- scores(rda_model, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Sample")

# Extract environmental variable loadings (arrows)
env_scores <- scores(rda_model, display = "bp") %>%
  as.data.frame() %>%
  rownames_to_column("Variable")

# Ensure metadata order matches RDA scores
rda_scores <- rda_scores %>%
  mutate(plate.geno = sapply(strsplit(Sample, " ; "), tail, 1))

# Ensure order consistency between dendrogram and db-RDA plot
ordered_samples_clean <- sapply(strsplit(ordered_samples, " ; "), tail, 1)  # Extract plate.geno
ordered_samples_clean <- trimws(ordered_samples_clean)  # Remove any extra spaces
ordered_samples_clean <- toupper(ordered_samples_clean)  # Standardize case

# Map clusters correctly using `match` to ensure proper indexing
ordered_clusters <- tclus[match(ordered_samples, names(tclus))]  # Ensure correct mapping

# Ensure metadata matches clustering
metadata_filtered <- metadata %>%
  filter(plate.geno %in% ordered_samples_clean) %>%
  arrange(match(plate.geno, ordered_samples_clean)) %>%
  mutate(Cluster = as.factor(ordered_clusters))  # Assign corrected clusters

# Define cluster colors
cluster_levels <- levels(as.factor(metadata_filtered$Cluster))  # Get unique cluster numbers
cluster_colors <- setNames(brewer.pal(min(length(cluster_levels), 4), "Set1"), cluster_levels)

# Compute dynamic arrow scaling factor
arrow_scaling_factor <- max(abs(c(env_scores$dbRDA1, env_scores$dbRDA2)), na.rm = TRUE) * 3

# Extract RDA scores for plotting
# Extract plate.geno from Sample names in RDA scores
sc_si <- scores(rda_model, display = "sites", choices = c(1,2), scaling = 1) %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  mutate(plate.geno = sapply(strsplit(Sample, " ; "), tail, 1)) %>%  # Extract plate.geno
  left_join(select(metadata_filtered, plate.geno, Cluster), by = "plate.geno") %>%  # Join on plate.geno
  mutate(Cluster = factor(Cluster, levels = cluster_levels)) %>%  # Ensure correct order
  select(Sample, dbRDA1, dbRDA2, Cluster)

# Extract scores for environmental variables (arrows)
sc_bp <- scores(rda_model, display = "bp", choices = c(1,2), scaling = 1) %>%
  as.data.frame() %>%
  rownames_to_column("Variable")

# Define proper axis limits
xlim_vals <- range(c(sc_si$dbRDA1, sc_bp$dbRDA1 * arrow_scaling_factor), na.rm = TRUE) * 1.3
ylim_vals <- range(c(sc_si$dbRDA2, sc_bp$dbRDA2 * arrow_scaling_factor), na.rm = TRUE) * 1.3

## Open PDF
pdf(file="OUTPUT/amun_db-RDA_nonseq_below100.pdf", width = 35, height = 20)

# Prepare the db-RDA plot (colored by cluster)
rda_plot <- ggplot() +
  geom_point(data = sc_si, aes(x = dbRDA1, y = dbRDA2, fill = Cluster, color = Cluster), 
             shape = 21, size = 8, stroke = 2) +  # No sequential filtration or pore size
  geom_segment(data = sc_bp, 
               aes(x = 0, y = 0, xend = dbRDA1 * arrow_scaling_factor, yend = dbRDA2 * arrow_scaling_factor), 
               arrow = arrow(length = unit(0.5, "cm")), size = 1.5, color = "black") + 
  geom_text_repel(data = sc_bp, 
                  aes(x = dbRDA1 * arrow_scaling_factor, 
                      y = dbRDA2 * arrow_scaling_factor, 
                      label = Variable),
                  color = "black", size = 10, fontface = 'bold',                 # Increased label size
                  max.overlaps = Inf, box.padding = 1.5, point.padding = 0.7,
                  segment.color = "grey50") +                                    
  scale_fill_manual(name = "Cluster", values = cluster_colors) +  
  scale_color_manual(name = "Cluster", values = cluster_colors) +  
  labs(
    title = paste0(
      "db-RDA plot for Amundsen 2023 (non-sequential filtration, below 100 m)\n",
      "(", variance_partitioning, "; ",
      rda_adj_r2, ")"),
    x = paste0("dbRDA1 (", perc[1], "%)"),
    y = paste0("dbRDA2 (", perc[2], "%)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 50, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 35, face = "bold"),
    axis.title.y = element_text(size = 35, face = "bold"),
    plot.margin = margin(t = 40, r = 60, b = 40, l = 60),
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 35, face = "bold"),
    legend.spacing.y = unit(0.5, 'cm'),
    legend.key.size = unit(2, "cm")
  )

# Convert plots to grobs
rda_grob <- ggplotGrob(rda_plot)
grid.draw(rda_grob)  # Print plot on first page

# Print a table with ANOVA results for the significance of environmental variables
grid.newpage()
grid.draw(tableGrob(anova_per_var, rows = NULL,
                    theme = ttheme_minimal(
                      base_size = 60,
                      padding = unit(c(10, 10), "mm"),
                      core = list(bg_params = list(fill = c("white", "grey95"), col = "black")), 
                      colhead = list(bg_params = list(fill = "grey80", col = "black")))))

# Close PDF
dev.off()
