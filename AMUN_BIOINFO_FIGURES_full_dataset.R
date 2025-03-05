# for simple heatmaps
library(gplots)

library(viridis)
library(viridisLite)  # for turbo coloring

library(vegan)        # for RDA, standardization, etc.

library(dendextend)   # for cutting tree

library("ComplexHeatmap")

library(ggrepel)      # for labeling with geom_text_repel
library(ggforce)      # for geom_mark_ellipse

library(gridExtra)    # for grid.draw
library(GGally)       # for ggpairs

library(ggtree)       # for dendogram drawing
library(ape)          # for as.phylo

# build Amundsen ASV by running AMUN_BIOINFO_ASV
source("AMUN_BIOINFO_ASV.R")

#
### What does your ASV table look like in terms of data distribution? What normalization should you apply?
# Apply one of these transformation for optimal clustering, choose based on visual inspection of the elbow/silhouette and heatmaps
hist(as.matrix(proptab.small.taxan.saman))         # looks pretty bad
hist(as.matrix(proptab.small.taxan.saman)^0.5)     # looks better
hist(as.matrix(proptab.small.taxan.saman)^0.25)    # looks good: a lot of 0s but normal distribution to the left
hist(as.matrix(proptab.small.taxan.saman)^0.2)     # looks good
hist(as.matrix(proptab.small.taxan.saman)^0.125)   # looks good


#############################################################
### Statistically discover the optimal number of clusters ###

kr <- 3:20

# --------------------------------------------
# Step 1: Compute Bray-Curtis distance matrix
# --------------------------------------------
bray_dist <- vegdist(as.matrix(proptab.small.taxan.saman)^0.2, method = "bray")

# Convert to matrix for calculations
bray_matrix <- as.matrix(bray_dist)

# --------------------------------------------
# Step 2: Hierarchical clustering with ward.D2
# --------------------------------------------
hclust_result <- hclust(bray_dist, method = "ward.D2")

# --------------------------------------------
# Step 3: Elbow Plot (WCSS)
# --------------------------------------------
# Calculate Within-Cluster Sum of Squares (WCSS)
wcss <- sapply(kr, function(k) {
  clusters <- cutree(hclust_result, k)
  sum(sapply(1:k, function(i) {
    cluster_indices <- which(clusters == i)
    if (length(cluster_indices) > 0) {
      sum(bray_matrix[cluster_indices, cluster_indices]^2) / (2 * length(cluster_indices))
    } else 0
  }))
})

# Plot elbow curve
elbow_plot <- ggplot(data.frame(k = kr, wcss = wcss), aes(k, wcss)) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_point(color = "steelblue", size = 3) +
  labs(title = "Elbow Plot for Optimal Clusters",
       x = "Number of Clusters (k)",
       y = "Within-Cluster Sum of Squares (WCSS)") +
  theme_minimal()

# --------------------------------------------
# Step 4: Silhouette Analysis
# --------------------------------------------
# Calculate average silhouette width
library(cluster)
sil_width <- sapply(kr, function(k) {
  clusters <- cutree(hclust_result, k)
  sil <- silhouette(clusters, bray_dist)
  mean(sil[, 3])
})

# Plot silhouette widths
silhouette_plot <- ggplot(data.frame(k = kr, sil = sil_width), aes(k, sil)) +
  geom_line(color = "firebrick", size = 1.2) +
  geom_point(color = "firebrick", size = 3) +
  labs(title = "Silhouette Analysis",
       x = "Number of Clusters (k)",
       y = "Average Silhouette Width") +
  theme_minimal()

# --------------------------------------------
# Combine Plots and Print out optimal number of clusters
# --------------------------------------------
library(patchwork)
combined_plot <- (elbow_plot | silhouette_plot)
print(combined_plot)


#####################
##### Simple Heatmap

#
### build a simple heatmap to explore all the data and see how many clusters you want

## open pdf
pdf(file="OUTPUT/AMUN/amun_simple_hm.pdf", width=35, height=20)

## simple heatmap
hm = heatmap.2(as.matrix(proptab.small.taxan.saman)^.2, trace="none", col=viridis,
               distfun=function(x) vegdist(x, method="bray"),
               hclustfun=function(x) hclust(x, method="ward.D2"))

## see how many clusters you may have for samples (rows)
clusnum=9
oldclus = cutree(hm$rowDendrogram, k=clusnum)
oldorder = unname(rle(oldclus[as.hclust(hm$rowDendrogram)$order])$values)
neworder = (1:clusnum)
names(neworder) = oldorder
clus = unname(neworder[as.character(oldclus)])
names(clus) = names(oldclus)

## redraw the map with clusters colored
hm = heatmap.2(as.matrix(proptab.small.taxan.saman)^.2, trace="none", key=F, margins=c(20,30), col=viridis,
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
bray_dist <- function(x) vegdist(x,method = "bray")

k = 9 # number of clusters you decided on

tdist <- vegdist(as.matrix(proptab.small.taxan.saman)^.2, method = "bray")
ttree <- hclust(tdist, method = "ward.D2")

# Assign cluster numbers
tclus <- cutree(ttree, k = k)

# Order clusters so they appear in order from top to bottom
ordered_samples <- names(sort(tclus))  # Get sample names in correct order

# Reorder the data matrix by the new order
ordered_matrix <- (as.matrix(proptab.small.taxan.saman)^.2)[ordered_samples, ]
ordered_clusters <- tclus[ordered_samples]  # Match clusters to reordered samples

# Reorder the environmental variables by the new order
ordered_meta <- meta.amun.small %>%
  arrange(match(plate.geno, sapply(strsplit(rownames(ordered_matrix), " ; "), tail, 1)))  # Reorder to match heatmap

#
## open pdf
pdf(file="OUTPUT/AMUN/amun_complex_hm.pdf", width=33, height=20)

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
          col = viridis(256),
          
          right_annotation = rowAnnotation(bar1 = anno_points(log2(ordered_meta$"Depth [m]"),
                                                              width = unit(3, "cm"),
                                                              axis_param = list(
                                                                at = log2(c(2, 10, 60, 200, 600, 2000)), # Log2-transformed depth
                                                                labels = c(2, 10, 60, 200, 600, 2000),   # Original depth labels
                                                                gp = gpar(fontsize = 8))),
                                           bar2 = anno_points(ordered_meta$"Temperature [deg C]",
                                                              width = unit (3, "cm"),
                                                              axis_param = list(
                                                                gp = gpar(fontsize = 8))),
                                           bar3 = anno_points(ordered_meta$"Salinity",
                                                              width = unit (3, "cm"),
                                                              axis_param = list(
                                                                gp = gpar(fontsize = 8))),
                                           bar4 = anno_points(ordered_meta$"Dissolved Oxygen [uM]",
                                                              width = unit (3, "cm"),
                                                              axis_param = list(
                                                                gp = gpar(fontsize = 8))),
                                           bar5 = anno_points(ordered_meta$"Fluorescence [mg/m^3]",
                                                              width = unit (3, "cm"),
                                                              axis_param = list(
                                                                gp = gpar(fontsize = 8))),
                                           bar6 = anno_points(ordered_meta$"Nitrate [uM]",
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
            title = "Fifth-root transformed relative abundances of microbial taxa during the 2023 Amundsen expedition\n",
            title_gp = gpar(fontsize = 30, fontface = "bold"),
            legend_direction = "horizontal",
            legend_width = unit(2, "cm"),
            color_bar = "continuos"),
          
          # show_heatmap_legend = FALSE                                   # Turn off heatmap legend b/c numbers are meaningless  
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

# Generate Turbo colors dynamically for clusters
cluster_palette <- turbo(k, alpha = 0.5)  # k colors with transparency

# Loop over cluster slices
for (i in 1:k) {
  decorate_row_names("amun_complex_hm", {
    grid.rect(gp = gpar(fill = cluster_palette[i]),  # Apply transparency
              width = unit(13.3, "cm"),
              x = unit(1, "npc") + unit(0.6, "cm"),
              just = "center")
  }, slice = i)
}

## close the pdf
dev.off()

################################################################################

#### Some cluster stats
proptab.small.taxan.saman %>%
  as.data.frame() %>%
  mutate(Cluster = tclus[rownames(.)]) %>%  # Add cluster column
  cbind(meta.amun.small) %>%
  group_by(Cluster) %>%
  summarize(
    most_freq_geo_loc = names(which.max(table(geo_loc_name))),
    most_freq_cruise = names(which.max(table(Cruise))),
    depth_range = paste(round(min(`Depth [m]`, na.rm = TRUE)), "-", round(max(`Depth [m]`, na.rm = TRUE))),  # Min-max depth
    second_depth_range = paste(round(sort(na.omit(`Depth [m]`)))[2], "-", round(sort(na.omit(`Depth [m]`), decreasing = TRUE))[2]),
    median_pore = median(`pore.size`, na.rm = TRUE),
    mean_seq = mean(seq.filt))

# Cluster most_freq_geo_loc          most_freq_cruise depth_range second_depth_range median_pore mean_seq
# <int> <chr>                      <chr>            <chr>       <chr>                    <dbl>    <dbl>
# 1       1 Hudson Strait              FOXSIPP          2 - 198     2 - 178                   3       0.702
# 2       2 Hudson Strait              FOXSIPP          118 - 518   119 - 430                 3       0.778
# 3       3 Foxe Basin                 FOXSIPP          2 - 94      2 - 74                    3       0.85 
# 4       4 The North Western Passages FOXSIPP          Inf - -Inf  NA - NA                   0.22    1    
# 5       5 Foxe Basin                 KEBABB           3 - 157     5 - 156                   0.22    0.375
# 6       6 Baffin Bay                 KEBABB           3 - 1845    136 - 693                 0.22    0.05 
# 7       7 Baffin Bay                 KEBABB           2 - 539     2 - 51                    0.22    0    
# 8       8 The North Western Passages KEBABB           2 - 116     2 - 54                    0.22    0    
# 9       9 Davis Strait               KEBABB           2 - 33      2 - 31                    0.22    0 

################################################################################

### Dendrogram
# Create cluster data frame for annotation
cluster_data <- data.frame(label = names(tclus),
                           cluster = tclus,
                           stringsAsFactors = FALSE)
phylo <- as.phylo(ttree)

# Open PDF with optimized dimensions
pdf(file = "OUTPUT/AMUN/amun_dend.pdf", width = 8, height = 8)

cluster_colors <- setNames(viridisLite::turbo(k), as.character(1:k))  # adjust cluster colors, the same one is defined in db-RDA

# Create base circular tree with compact layout
p <- ggtree(phylo, layout = 'circular', ladderize = TRUE, size = 0.8) %<+% cluster_data +
  geom_tree(aes(color = factor(cluster))) +                           # color branches by cluster
  geom_tiplab(aes(color = factor(cluster)),                           # label spacing
              size = 1,                        # label size
              offset = 0.5,                    # pushing labels outward
              align = TRUE,                    # prevent overlap
              show.legend = TRUE) +
  scale_color_manual(name = "Cluster", values = cluster_colors,
                     guide = guide_legend(override.aes = list(size = 5))) +     # your predefined colors for clusters
  expand_limits(x = 6)                         # expands 

print(p)
dev.off()

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

#
## Select and transform microbial and environmental 

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
          
          # only for KEBABB
          -"Chla [mg/m^3]",                # average total chlorophyll a (mg/m^3), similar to the unit of µg/L
          -"Chla >5 um [mg/m^3]",          # average >5 um chlorophyll a in water (mg/m^3), similar to the unit of µg/L
          -"LNA bacteria [1/ml]",          # low nucleic acid content
          -"HNA bacteria [1/ml]",          # high nucleic acid content
          -"Pico PE-cyans [1/ml]",         # Pico- (<2 µm), phycoerythrin-containing cyanobacteria
          -"Nano PE-cyans [1/ml]",         # Nano- (2-20 µm), phycoerythrin-containing cyanobacteria
          -"Pico PC-cyans [1/ml]",         # Pico- (<2 µm), phycocyanin-containing cyanobacteria
          -"Nano PC-cyans [1/ml]",         # Nano- (2-20 µm), phycocyanin-containing cyanobacteria
          -"Picoeukaryotes [1/ml]",
          -"Nanoeukaryotes [1/ml]",
          
          -"DOY",                                         # it has a big effect b/c two cruises were temporally separate
          
          -"Oxygen solubility [ml/l]",                    # correlates with temperature based on ggpair() + useless
          -"Nitrate CTD [uM]",                            # correlates with bottle nitrate based on ggpair() + we have bottle nitrate

          -"Dissolved Oxygen [uM]",                       # has strong correlation with nitrate (ggpair())
          -"Silicate [uM]",                               # has strong correlation with depth and other nutrients (ggpair())

          -"Light Attenuation Coefficient [1/m]",         # has weak significance (*) based on anova.cca + using fluo
          -"Light Transmission [%]"                       # has weak significance (**) based on anova.cca + using fluo
          # -"Fluorescence [mg/m^3]",                       
          # -"CDOM [mg/m^3]",                               

          # -"Phosphate [uM]",                              
          # -"Nitrite [uM]"
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
complete_cases <- complete.cases(env_vars)          # 177/186: 80 FOXSIPP and 97 KEBABB
microbial_data <- microbial_data[complete_cases, , drop = FALSE]
env_vars <- env_vars[complete_cases, , drop = FALSE]

# The environmental/explanatory variables are all in different units and need to be standardized prior to computing distance measures to perform most ordination analyses
env_stand <- decostand(env_vars, method="standardize")  # subtract the mean and divide by the standard deviation
# this centers each variable around 0 (mean = 0) and scales it so that its standard deviation is 1
apply(env_stand, 2, mean) # check if the data are now centered (means~0)
apply(env_stand, 2, sd)   # check if the data are now scaled (standard deviations=1)

#
## Build db-RDA model and remove/add environmental variables above + rerun, as needed based on the results of the checks below

# Before applying db-RDA, you can check which dissimilarity index works best for you
# Usually, we choose Bray anyway because it works best with microbial data, calculates distances better
rankindex(env_stand, microbial_data, indices = c("euc", "man", "gow", "bra", "kul"),
          stepacross= FALSE, method = "pearson")

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

# Perform ANOVA checks for the model, environmental variables an axes.
# anova.cca performs an ANOVA like permutation test for cca, rda, or dbRDA to assess the significance of constraints
# Also, extract and look at environmental variable scores, % variance explained, R²

# Print significance of the model (has to be p < 0.001)
anova.cca(rda_model, permutations = 1000)

# Extract scores for environmental variables (arrows)
# Scaling = 1 if interested in clustering of samples (distance between points)
# Scaling = 2 if interested in correlations between env. variables (arrows) and how they influence community structure
# Look at the variables that impact each of your dbRDA axis
# interpret how environmental variables influence community structure
sc_bp <- scores(rda_model, display = "bp", choices = c(1:ncol(env_stand)), scaling = 2) %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  print()

# Print significance of your environmental variables
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

# Print marginal effects of terms: how well each variable adds to explaining the total variance
# Yet another check as to whether or not you should remove some variables due to their strong correlation
# While "terms" assesses the sequential effects of var-les (order matters), "margin" evaluates each var-le independently
anova.cca(rda_model, permutations = 1000, by = "margin") %>%
  arrange(desc(F))                                          # Sort by descending F-value

# Print significance of the first two dbRDA axis (has to be p < 0.001)
anova.cca(rda_model, permutations = 1000, by = "axis")

# Print variance explained by axes
perc <- round(100 * (summary(rda_model)$cont$importance[2, 1:2]), 2) %>%
  print()

# Print proportion of variance explained by environmental variables
variance_partitioning <- (summary(rda_model)$constr.chi / summary(rda_model)$tot.chi) %>%
  {paste("variance explained:", round(. * 100, 2), "%")} %>%
  print()

# Print final adjusted R²
rda_adj_r2 <- RsquareAdj(rda_model)$adj.r.squared %>%
  {paste("adjusted R²:", round(. * 100, 1), "%")} %>%
  print()

#
## Prepare data for plotting

# Ensure order consistency between dendrogram and db-RDA plot
ordered_samples_clean <- sapply(strsplit(ordered_samples, " ; "), tail, 1)  # Extract plate.geno

# Map clusters correctly using `match` to ensure proper indexing
ordered_clusters <- tclus[match(ordered_samples, names(tclus))]  # Ensure correct mapping

# Ensure metadata matches clustering
metadata_filtered <- meta.amun.small %>%
  arrange(match(plate.geno, sapply(strsplit(ordered_samples, " ; "), tail, 1))) %>%
  mutate(Cluster = as.factor(tclus[match(ordered_samples, names(tclus))]))    # assign corrected clusters

# Extract RDA scores for plotting and add metadata to it
sc_si <- scores(rda_model, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  mutate(plate.geno = sapply(strsplit(Sample, " ; "), tail, 1)) %>%
  left_join(metadata_filtered, by = "plate.geno") %>%
  mutate(Sequential_Filtration = ifelse(seq.filt == 1, "Yes", "No"),
         Pore_Size = factor(as.character(pore.size), levels = c("0.22", "3", "20")),
         Sequential_Filtration = factor(Sequential_Filtration, levels = c("Yes", "No")))

# Define cluster colors
cluster_levels <- levels(as.factor(metadata_filtered$Cluster))  # Get unique cluster numbers
cluster_colors <- setNames(viridisLite::turbo(length(cluster_levels)), cluster_levels) # use turbo: better and can also use in ODV
# turbo: an improved Jet rainbow color map for reducing false detail, banding and color blindness ambiguity

# Compute dynamic arrow scaling factor
arrow_scaling_factor <- max(abs(c(sc_bp$dbRDA1, sc_bp$dbRDA2)), na.rm = TRUE) * 2

# Define proper axis limits
xlim_vals <- range(c(sc_si$dbRDA1, sc_bp$dbRDA1 * arrow_scaling_factor), na.rm = TRUE) * 1.3
ylim_vals <- range(c(sc_si$dbRDA2, sc_bp$dbRDA2 * arrow_scaling_factor), na.rm = TRUE) * 1.3

## Open PDF
pdf(file="OUTPUT/AMUN/amun_db-RDA.pdf", width = 35, height = 20)

# Loop through all numeric variables in `sc_si` to generate additional db-RDA plots
numeric_vars <- c("DOY", "Depth [m]", "volume.l", "dna.conc",
                  "Temperature [deg C]", "Salinity",
                  "Dissolved Oxygen [uM]", "Oxygen solubility [ml/l]",
                  "PAR [uE/s/m^2]", "Light Attenuation Coefficient [1/m]", "Light Transmission [%]",
                  "CDOM [mg/m^3]", "Fluorescence [mg/m^3]",              
                  "Nitrate CTD [uM]", "Nitrate [uM]", "Nitrite [uM]", "Phosphate [uM]", "Silicate [uM]")

for (param in numeric_vars) {
  
  param_min <- min(sc_si[[param]], na.rm = TRUE)
  param_max <- max(sc_si[[param]], na.rm = TRUE)
  
  param_plot <- ggplot() +
    
    # Base layer: Ensures ALL shapes appear in the legend
    geom_point(
      data = sc_si %>%
        mutate(Pore_Size = factor(Pore_Size, levels = c("0.22", "3", "20"))),  # Explicit levels to force display
      aes(x = dbRDA1, y = dbRDA2, shape = Pore_Size),
      alpha = 0, size = 0  # Invisible layer to enforce legend
    ) +
    
    # Sequential Filtration ("Yes") - Filled Points
    geom_point(
      data = sc_si %>% filter(Sequential_Filtration == "Yes"),
      aes(x = dbRDA1, y = dbRDA2,
          fill = Cluster,
          shape = Pore_Size,
          size = !!sym(param)),
      alpha = 0.4, color = "black", stroke = 3  # Strong outline for better visibility
    ) +
    
    # Non-Sequential Filtration ("No") - Hollow Points with Colored Outline
    geom_point(
      data = sc_si %>% filter(Sequential_Filtration == "No"),
      aes(x = dbRDA1, y = dbRDA2,
          color = Cluster,
          shape = Pore_Size,
          size = !!sym(param)),
      fill = "white", alpha = 1, stroke = 3    # Proper stroke for visibility
    ) +
    
    geom_segment(data = sc_bp, 
                 aes(x = 0, y = 0, xend = dbRDA1 * arrow_scaling_factor, yend = dbRDA2 * arrow_scaling_factor), 
                 arrow = arrow(length = unit(0.5, "cm")), size = 1.5, color = "black") +
    
    geom_text_repel(data = sc_bp, 
                    aes(x = dbRDA1 * arrow_scaling_factor, 
                        y = dbRDA2 * arrow_scaling_factor, 
                        label = Variable),
                    color = "black", size = 10, fontface = 'bold',
                    max.overlaps = Inf, box.padding = 1.5, point.padding = 0.7,
                    segment.color = "grey50") + 
    
    # Ensure Cluster and Pore Size Legends are Correct
    scale_fill_manual(name = "Cluster", values = cluster_colors) +
    scale_color_manual(name = "Cluster", values = cluster_colors) +
    scale_shape_manual(name = "Pore Size [μm]", 
                       values = c("0.22" = 21, "3" = 22, "20" = 23), 
                       drop = FALSE) +         # Prevent ggplot from removing unused shapes
    
    scale_size(range = scales::rescale(range(sc_si[[param]], na.rm = TRUE), to = c(2, 20))) +
  
    # Create custom legend for filtration type
    guides(
      fill = guide_legend(
        title = "Cluster (Filled = Sequential)",
        override.aes = list(alpha = 0.4, shape = 21, color = "black", size = 15)           # Make sure filled points appear
      ),
      color = guide_legend(
        title = "Cluster (Hollow = Non-Sequential)",
        override.aes = list(alpha = 1, shape = 21, fill = "white", size = 15, stroke = 3)  # Ensure proper stroke
      ),
      shape = guide_legend(title = "Pore Size [μm]", override.aes = list(size = 15)),  
      size = guide_legend(order = 2)  
    ) +
    
    labs(title = paste0(
      "db-RDA plot for Amundsen 2023: ", param,
      "\n(", variance_partitioning, "; ", rda_adj_r2, ")"),
      x = paste0("dbRDA1 (", perc[1], "%)"),
      y = paste0("dbRDA2 (", perc[2], "%)")) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(size = 50, face = "bold", hjust = 0.5, margin = margin(b = 20)),
      axis.text = element_text(size = 16),
      axis.title.x = element_text(size = 35, face = "bold"),
      axis.title.y = element_text(size = 35, face = "bold"),
      plot.margin = margin(t = 40, r = 60, b = 40, l = 60),
      legend.text = element_text(size = 30),
      legend.title = element_text(size = 35, face = "bold"),
      legend.key.size = unit(2, "cm"),    # Reduce legend key size
      legend.spacing.y = unit(0.1, "cm")  # Reduce vertical spacing between legend items
    )
  
  print(param_plot)
}

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

################################################################################

# Additional material

# ## Open PDF
# pdf(file="OUTPUT/amun_db-RDA.pdf", width = 35, height = 20)
# 
# # Prepare the db-RDA plot (colored by cluster)
# rda_plot <- ggplot() +
#   
#   geom_point(data = sc_si, aes(x = dbRDA1, y = dbRDA2,
#                                fill = Cluster,
#                                color = Cluster,
#                                shape = Pore_Size, 
#                                alpha = Sequential_Filtration,
#                                size = `Fluorescence [mg/m^3]`), 
#              stroke = 2) +
#   
#   geom_segment(data = sc_bp, 
#                aes(x = 0, y = 0, xend = dbRDA1 * arrow_scaling_factor, yend = dbRDA2 * arrow_scaling_factor), 
#                arrow = arrow(length = unit(0.5, "cm")), size = 1.5, color = "black") +
#   
#   geom_text_repel(data = sc_bp, 
#                   aes(x = dbRDA1 * arrow_scaling_factor, 
#                       y = dbRDA2 * arrow_scaling_factor, 
#                       label = Variable),
#                   color = "black", size = 10, fontface = 'bold',                 # Increased label size
#                   max.overlaps = Inf, box.padding = 1.5, point.padding = 0.7,
#                   segment.color = "grey50") + 
#   
#   scale_fill_manual(name = "Cluster", values = cluster_colors, 
#                     guide = guide_legend(override.aes = list(size = 15))) +      # Larger legend points
#   
#   scale_color_manual(name = "Cluster", values = cluster_colors, 
#                      guide = guide_legend(override.aes = list(size = 15))) +     # Larger legend points
#   
#   scale_shape_manual(name = "Pore Size", 
#                      values = c("0.22" = 21, "3" = 22, "20" = 23), 
#                      guide = guide_legend(override.aes = list(size = 15))) +     # Larger shape size
#   
#   scale_alpha_manual(name = "Sequential Filtration", values = c("Yes" = 0.4, "No" = 1), 
#                      guide = guide_legend(override.aes = list(size = 15))) +     # Larger alpha legend
#   
#   scale_size_continuous(range = c(0, 15)) +
#   
#   labs(
#     title = paste0(
#       "db-RDA plot for Amundsen 2023",
#       " (", variance_partitioning, "; ",
#       rda_adj_r2, ")"),
#     x = paste0("dbRDA1 (", perc[1], "%)"),
#     y = paste0("dbRDA2 (", perc[2], "%)")
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 50, face = "bold", hjust = 0.5, margin = margin(b = 20)),
#     axis.text = element_text(size = 16),
#     axis.title.x = element_text(size = 35, face = "bold"),
#     axis.title.y = element_text(size = 35, face = "bold"),
#     plot.margin = margin(t = 40, r = 60, b = 40, l = 60),
#     legend.text = element_text(size = 30),
#     legend.title = element_text(size = 35, face = "bold"),
#     legend.spacing.y = unit(0.5, 'cm'),
#     legend.key.size = unit(2.5, "cm")  # Increase legend key size
#   )
# 
# # Convert plots to grobs
# rda_grob <- ggplotGrob(rda_plot)
# grid.draw(rda_grob)  # Print plot on first page
# 
# # Print a table with ANOVA results for the significance of environmental variables
# grid.newpage()
# grid.draw(tableGrob(anova_per_var, rows = NULL,
#                     theme = ttheme_minimal(
#                       base_size = 60,
#                       padding = unit(c(10, 10), "mm"),
#                       core = list(bg_params = list(fill = c("white", "grey95"), col = "black")), 
#                       colhead = list(bg_params = list(fill = "grey80", col = "black")))))
# 
# # Close PDF
# dev.off()


# ###### ODV for Amundsen DNA
# write.csv(
#   proptab.small.taxan.saman %>%
#     mutate(Cluster = tclus[rownames(.)]) %>%  # Add cluster column
#     cbind(meta.amun.small),                   # Combine with the meta table
#   "OUTPUT/AMUN_DNA_ODV.csv",
#   row.names = FALSE)

