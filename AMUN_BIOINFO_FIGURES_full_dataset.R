# for simple heatmaps
library(gplots)
library(viridis)
library(vegan)

library(dendextend) # for cutting tree

library("ComplexHeatmap")

library(ggrepel) # for labeling with geom_text_repel
library(ggforce) # for geom_mark_ellipse

library(gridExtra) # for grid.draw
library(GGally)    # for ggpairs

# build Amundsen ASV by running AMUN_BIOINFO_ASV

#
### What does your ASV table look like in terms of data distribution? What normalization should you apply?
hist(as.matrix(proptab.small.taxan.saman))         # looks pretty bad
hist(as.matrix(proptab.small.taxan.saman)^0.5)     # looks better
hist(as.matrix(proptab.small.taxan.saman)^0.25)    # looks optimal: a lot of 0s but normal distribution to the left


#####################
##### Simple Heatmap

#
### build a simple heatmap to explore all the data

## open pdf
pdf(file="OUTPUT/amun_simple_hm.pdf", width=35, height=20)

## simple heatmap
hm = heatmap.2(as.matrix(proptab.small.taxan.saman)^.25, trace="none", col=viridis,
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
hm = heatmap.2(as.matrix(proptab.small.taxan.saman)^.25, trace="none", key=F, margins=c(20,30), col=viridis,
               distfun=function(x) vegdist(x, method="bray"),
               hclustfun=function(x) hclust(x, method="ward.D2"),
               RowSideColors = rainbow(max(clus))[clus])

## close the pdf
dev.off()


#######################
#######################
##### ComplexHeatmap

# colnames(meta.amun.small) # see which meta parameters you can investigate

#
# Compute Bray-Curtis distance and hierarchical clustering
tdist <- vegdist(as.matrix(proptab.small.taxan.saman)^.25, method = "bray")
ttree <- hclust(tdist, method = "ward.D2")

# Assign cluster numbers
tclus <- cutree(ttree, k = 7)

# Order clusters so they appear in order from top to bottom
ordered_samples <- names(sort(tclus))  # Get sample names in correct order

# Reorder the data matrix by the new order
ordered_matrix <- (as.matrix(proptab.small.taxan.saman)^.25)[ordered_samples, ]
ordered_clusters <- tclus[ordered_samples]  # Match clusters to reordered samples

#
## open pdf
pdf(file="OUTPUT/amun_complex_hm.pdf", width=33, height=20)

ht_opt$ROW_ANNO_PADDING = unit(2, "mm")                                          # Padding between the heatmap and annotations

#
## create a complex heatmap
amun_complex_hm <- 
  
  Heatmap(ordered_matrix,
          
          # clustering_distance_rows = bray_dist,
          clustering_distance_columns = bray_dist,
          # clustering_method_rows = "ward.D2",                                  #passed to hclust
          clustering_method_columns = "ward.D2",                                 #passed to hclust
          
          name = "amun_complex_hm",
          col = blues9,
          
          right_annotation = rowAnnotation(bar1 = anno_points(log2(meta.amun.small$"Depth [m]"),
                                                              width = unit(3, "cm"),
                                                              axis_param = list(
                                                                gp = gpar(fontsize = 8))),
                                           bar2 = anno_points(meta.amun.small$"Temperature [deg C]",
                                                              width = unit (3, "cm"),
                                                              axis_param = list(
                                                                gp = gpar(fontsize = 8))),
                                           bar3 = anno_points(meta.amun.small$"Salinity",
                                                              width = unit (3, "cm"),
                                                              axis_param = list(
                                                                gp = gpar(fontsize = 8))),
                                           bar4 = anno_points(meta.amun.small$"Dissolved Oxygen [uM]",
                                                              width = unit (3, "cm"),
                                                              axis_param = list(
                                                                gp = gpar(fontsize = 8))),
                                           bar5 = anno_points(meta.amun.small$"Fluorescence [mg/m^3]",
                                                              width = unit (3, "cm"),
                                                              axis_param = list(
                                                                gp = gpar(fontsize = 8))),
                                           bar6 = anno_points(meta.amun.small$"Nitrate [uM]",
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
          # row_labels = rep("", nrow(proptab.small.taxan.saman)),    # If you want to remove row names (aka samples)
          
          column_title = "Taxa",
          column_title_gp = gpar(fontsize = 25, fontface = "bold"),
          
          row_names_gp = gpar(fontsize = 5, fontface = "bold"),
          column_names_gp = gpar(fontsize = 1, fontface = "bold"),
          column_names_rot = 45,
          
          heatmap_legend_param = list(
            title = "Fourth-root transformed relative abundances of microbial taxa during the 2023 Amundsen expedition\n",
            title_gp = gpar(fontsize = 30, fontface = "bold"),
            legend_direction = "horizontal",
            legend_width = unit(2, "cm"),
            color_bar = "continuos"),
          
          # show_heatmap_legend = FALSE                                # Turn off heatmap legend b/c numbers are meaningless  
  )

#
## Draw the complex heatmap
draw(amun_complex_hm,                                                     # First draw, then decorate, place legend in a good spot
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
colors <- substr(rainbow(7, alpha = 0.25), 1, 7)                                 # Adding transparency with alpha

# Loop over 7 slices
for (i in 1:7) {
  decorate_row_names("amun_complex_hm", {
    grid.rect(gp = gpar(fill = paste0(colors[i], "40")),                         # Apply mild transparency
              width = unit(13.3, "cm"),                                          # Width of the colored box
              x = unit(1, "npc") + unit(0.6, "cm"),                              # Play with this to shift the box onto the labels
              just = "center")
  }, slice = i)
}

## close the pdf
dev.off()


#### Additional
proptab.small.taxan.saman %>%
  mutate(Cluster = tclus[rownames(.)]) %>%  # Add cluster column
  cbind(meta.amun.small) %>%
  group_by(Cluster) %>%
  summarize(min_depth = min(`Depth [m]`, na.rm = TRUE),
            max_depth = max(`Depth [m]`, na.rm = TRUE),
            mean_depth = mean(`Depth [m]`, na.rm = TRUE))


########################################
##### db-RDA (distance-based Redundancy Analysis; in this case, we will base it on the Bray Curtis Distance matrix)
# to identify how environmental variables influence microbial relative abundance
# the vector of the arrow points in the direction of its greatest influence
# PCA is an unconstrained form of ordination - searches for any variables best explaining the species composition
# RDA is constrained - searches for the best explanatory (environmental) variables

# Info that helped me to put together a db-RDA pipeline below:
# https://youtu.be/Wqa7_18jf5A
# https://youtu.be/_bpQ16uUqs4
# https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html
# https://wiki.qcbs.ca/r_workshop10
# https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/rda-and-dbrda/
# https://sites.ualberta.ca/~ahamann/teaching/renr690/seminars/dbRDA.pdf
# https://vegandevs.github.io/vegan/reference/goodness.cca.html


# Load microbial ASV data (response variables)
microbial_data <- as.matrix(proptab.small.taxan.saman)^0.25

# Check how your ASV data looks like
# If the distribution is not normal - go back and modify.
# Here we have already done prop.table (by sample) and ^0.25 transformations on the original ASV matrix
hist(microbial_data)    # looks optimal: a lot of 0s but normal distribution to the left

# Load metadata (explanatory variables)
metadata <- meta.amun.small

# Select numeric environmental variables correctly
env_vars <- metadata %>%
  select(where(is.numeric)) %>%
  
  # # interpolate NO3 w/ CTD readings if you are planning to use all the samples at the expense of removing other nutrients
  # # not used here as the variance explained and adjusted R2 drops substantially
  # mutate(`Nitrate [uM]` = if_else(is.na(`Nitrate [uM]`), `Nitrate CTD [uM]`, `Nitrate [uM]`)) %>%
  
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
          
          -"Dissolved Oxygen [uM]",                       # has high VIF (14) and correlates with nitrate (ggpair())
          -"Silicate [uM]",                               # has high VIF (14) and correlates with nitrate and phosphate (ggpair())

          -"Light Attenuation Coefficient [1/m]",         # has weak significance (*) based on anova.cca
          -"Light Transmission [%]",                      # has weak significance (*) based on anova.cca
          -"Fluorescence [mg/m^3]",                       # has weak significance (*) based on anova.cca (can use chl a later?)
          # -"CDOM [mg/m^3]",                               # has weak significance (**) based on anova.cca
          
          # -"Phosphate [uM]",                              # too many NAs and can't interpolate + weak significance (**, anova.cca)
          # -"Nitrite [uM]"                                 # too many NAs and can't interpolate
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

metadata_filtered <- metadata %>%
  filter(plate.geno %in% rda_scores$plate.geno) %>%
  arrange(match(plate.geno, rda_scores$plate.geno))

# Add metadata variables to rda_scores
rda_scores <- rda_scores %>%
  mutate(
    Sequential_Filtration = ifelse(metadata_filtered$seq.filt == 1, "Yes", "No"),
    Pore_Size = as.character(metadata_filtered$pore.size),
    geo_loc_name = metadata_filtered$geo_loc_name) %>%
  mutate(Pore_Size = factor(Pore_Size, levels = c("0.22", "3", "20")),
         Sequential_Filtration = factor(Sequential_Filtration, levels = c("Yes", "No")))

# Custom plot settings
region_colors <- c("The North Western Passages" = "red1",
                   "Gulf of Boothia" = "magenta1",
                   "Fury and Hecla Strait" = "green",
                   "Foxe Basin" = "darkgoldenrod1",
                   "Hudson Strait" = "aquamarine4",
                   "Baffin Bay" = "chocolate1",
                   "Davis Strait" = "mediumpurple1",
                   "Labrador Sea" = "blue",
                   "The Saint Lawrence River" = "cyan",
                   "Nunavut" = "#4F4F4F")

# Compute dynamic arrow scaling factor
arrow_scaling_factor <- max(abs(c(env_scores$dbRDA1, env_scores$dbRDA2)), na.rm = TRUE) * 2

# Extract site scores
# Scaling = 1 if you are more interested in the relationship between samples and species (distance between points)
# Scaling = 2 if you are more interested in the correlation between env. variables
sc_si <- scores(rda_model, display = "sites", choices = c(1,2), scaling = 1) %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  left_join(select(rda_scores, Sample, geo_loc_name, Sequential_Filtration, Pore_Size), by = "Sample") %>%
  select(Sample, dbRDA1, dbRDA2, geo_loc_name, Sequential_Filtration, Pore_Size)

# Extract scores for environmental variables (arrows)
sc_bp <- scores(rda_model, display = "bp", choices = c(1,2), scaling = 1) %>%
  as.data.frame() %>%
  rownames_to_column("Variable")

# Define proper axis limits
xlim_vals <- range(c(sc_si$dbRDA1, sc_bp$dbRDA1 * arrow_scaling_factor), na.rm = TRUE) * 1.3
ylim_vals <- range(c(sc_si$dbRDA2, sc_bp$dbRDA2 * arrow_scaling_factor), na.rm = TRUE) * 1.3

### Plot using ggplot2

## Open PDF
pdf(file="OUTPUT/amun_db-RDA.pdf", width = 35, height = 20)

## Prepare the distance based RDA plot (db-RDA)
rda_plot <- ggplot() +
  geom_point(data = sc_si, aes(x = dbRDA1, y = dbRDA2, 
                               shape = Pore_Size, 
                               fill = geo_loc_name, 
                               color = geo_loc_name, 
                               alpha = Sequential_Filtration), 
             size = 8, stroke = 5) +                                             # Increased point size and stroke
  
  geom_segment(data = sc_bp, 
               aes(x = 0, y = 0, xend = dbRDA1 * arrow_scaling_factor, yend = dbRDA2 * arrow_scaling_factor), 
               arrow = arrow(length = unit(0.5, "cm")),                          # Increased arrow size
               size = 1.5, color = "black") +                                    # Changed arrows to black
  
  geom_text_repel(data = sc_bp, 
                  aes(x = dbRDA1 * arrow_scaling_factor, 
                      y = dbRDA2 * arrow_scaling_factor, 
                      label = Variable),
                  color = "black", size = 10, fontface = 'bold',                 # Increased label size
                  max.overlaps = Inf, box.padding = 1.5, point.padding = 0.7,
                  segment.color = "grey50") +                                    
  
  # Region colors now properly shown in the legend
  scale_fill_manual(name = "Region", values = region_colors) +  
  scale_color_manual(name = "Region", values = region_colors, guide = guide_legend(override.aes = list(shape = 21, size = 8))) +  
  
  scale_shape_manual(name = "Pore Size", values = c("0.22" = 21, "3" = 22, "20" = 23)) +
  scale_alpha_manual(name = "Sequential Filtration", values = c("Yes" = 0.4, "No" = 1)) +  # Semitransparent vs filled points
  
  labs(
    title = paste0(
      "db-RDA plot for Amundsen 2023\n",
      "(", variance_partitioning, "; ",
      rda_adj_r2, ")"),
    x = paste0("dbRDA1 (", perc[1], "%)"),
    y = paste0("dbRDA2 (", perc[2], "%)")
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(size = 50, face = "bold", hjust = 0.5, margin = margin(b = 20)),  # Added margin for spacing
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 35, face = "bold"),
    axis.title.y = element_text(size = 35, face = "bold"),
    plot.margin = margin(t = 40, r = 60, b = 40, l = 60),            # Add left and right margins
    legend.text = element_text(size = 30),                           # Make legend text large
    legend.title = element_text(size = 35, face = "bold"),           # Larger legend title
    legend.spacing.y = unit(0.5, 'cm'),                              # Reduce space between legend items
    legend.key.size = unit(2, "cm")                                  # Increase legend key size
  )

rda_grob <- ggplotGrob(rda_plot)                                     # Convert ggplot object to grob for compatibility
grid.draw(rda_grob)                                                  # Print plot on first page


## Print a table with ANOVA results for the significance of environmental variables
grid.newpage()                                                                              # Start a new page for the table
grid.draw(tableGrob(anova_per_var, rows = NULL,                                             # Print the table
                    theme = ttheme_minimal(
                      base_size = 60,   # Adjust text size for readability
                      padding = unit(c(10, 10), "mm"),  # Adjust padding for spacing
                      core = list(bg_params = list(fill = c("white", "grey95"), col = "black")),    # Alternating row colors
                      colhead = list(bg_params = list(fill = "grey80", col = "black")))))           # Header background
           

## Close PDF
dev.off()


# ###### ODV for Amundsen DNA
# write.csv(
#   proptab.small.taxan.saman %>%
#     mutate(Cluster = tclus[rownames(.)]) %>%  # Add cluster column
#     cbind(meta.amun.small),                   # Combine with the meta table
#   "OUTPUT/AMUN_DNA_ODV.csv",
#   row.names = FALSE)

