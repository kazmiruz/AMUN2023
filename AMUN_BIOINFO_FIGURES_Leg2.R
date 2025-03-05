# for simple heatmaps
library(gplots)

library(viridis)
library(viridisLite)  # for turbo coloring

library(vegan)        # for RDA, standardization, etc.

library(dendextend)   # for cutting tree

library("ComplexHeatmap") # for complexheatmap

library(ggrepel)      # for labeling with geom_text_repel
library(ggforce)      # for geom_mark_ellipse

library(gridExtra)    # for grid.draw and grid.table
library(GGally)       # for ggpairs

library(ggtree)       # for dendrogram drawing
library(ape)          # for as.phylo

# build Amundsen ASV by running AMUN_BIOINFO_ASV
source("AMUN_BIOINFO_ASV.R")

#
### filter data
<<<<<<< HEAD

# Leg 2&4: there is a difference in filt. protocol (Leg 2: only 0.22 non-seq) and nut. data avail-ty (Leg 2: has nuts for all samples). 
# While it is true that some 0.22 non-sequential filtration samples are from Leg 4, they are few and interfere w/ Leg 2 clustering.
# While deep Amundsen samples cluster separately, they're included in the analysis due to lack of clear depth sep-n and for comparison.

# Hence, here, Leg 2 is analysed completely and separately from Leg 4.
=======
# Since there is a big difference in filtration protocol and env. data avail-ty btwn Leg 2 and Leg 4, Leg 2 is looked at separately.
# While it is true that some 0.22 non-seqiential filtration samples are from Leg 4, they are very few of those.
# It was also found that deeper Amundsen samples (~100 m +) cluster separately from other the surface ones.
>>>>>>> 1266ca155daf646c4adce3d7719e777e811d11ee

meta_leg2 <- meta.amun.small %>%
  filter(Cruise == "KEBABB") %>%              # 97 of Leg 2 samples from 186 in total
  # filter (`Depth [m]` <= 80) %>%            # top 65 meters leaves 73/97 Leg 2 samples, next deepest sample is at 115 m
  mutate("Bacteria [1/ml]" = `LNA bacteria [1/ml]` + `HNA bacteria [1/ml]`,
         "Cyanobacteria PE [1/ml]" = `Pico PE-cyans [1/ml]` + `Nano PE-cyans [1/ml]`,
         "Cyanobacteria PC [1/ml]" = `Pico PC-cyans [1/ml]` + `Nano PC-cyans [1/ml]`,
         "Eukaryotes [1/ml]" = `Picoeukaryotes [1/ml]` + `Nanoeukaryotes [1/ml]`) %>%
  select(-"Chla >5 um [mg/m^3]",          # average >5 um chlorophyll a in water (mg/m^3), similar to the unit of µg/L
         -"LNA bacteria [1/ml]",          # low nucleic acid content
         -"HNA bacteria [1/ml]",          # high nucleic acid content
         -"Pico PE-cyans [1/ml]",         # Pico- (<2 µm), phycoerythrin-containing cyanobacteria
         -"Nano PE-cyans [1/ml]",         # Nano- (2-20 µm), phycoerythrin-containing cyanobacteria
         -"Pico PC-cyans [1/ml]",         # Pico- (<2 µm), phycocyanin-containing cyanobacteria
         -"Nano PC-cyans [1/ml]",         # Nano- (2-20 µm), phycocyanin-containing cyanobacteria
         -"Picoeukaryotes [1/ml]",
         -"Nanoeukaryotes [1/ml]")

proptab_leg2 <- proptab.small.taxan.saman %>%                         # the column names: "Phylum/Division; Order; Genus; ASV #"
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  mutate(plate.geno = sapply(strsplit(Sample, " ; "), tail, 1)) %>%   # Extract last part as plate.geno
  filter(plate.geno %in% meta_leg2$plate.geno) %>%                    # Keep only matching samples
  arrange(match(plate.geno, meta_leg2$plate.geno)) %>%                # Ensure order matches meta_leg2
  select(-plate.geno) %>%                                             # Remove plate.geno column
  column_to_rownames("Sample") %>%                                    # Restore row names after filtering
  .[colSums(. > 0) > 2] %>%                                           # select species present in more than 2 samples
  `rownames<-`(gsub("^KEBABB; ", "", gsub("; 0\\.22 to 350 um ; ", "; ", rownames(.))))      # Remove unnecessary info from rownames

# ggplot(meta_leg2, aes(x = `Fluorescence [mg/m^3]`, y = `Chla [mg/m^3]`)) + geom_point() + geom_smooth(method=lm)

################################################################################

#
### What does your ASV table look like in terms of data distribution? What normalization should you apply?
# Apply one of these transformation for optimal clustering, choose based on visual inspection of the elbow/silhouette and heatmaps
hist(as.matrix(proptab_leg2))        # looks pretty bad
hist(as.matrix(proptab_leg2)^.5)     # looks better
hist(as.matrix(proptab_leg2)^.25)    # looks good: a lot of 0s but normal distribution to the left
hist(as.matrix(proptab_leg2)^.2)     # looks good
hist(as.matrix(proptab_leg2)^.125)   # looks good


#############################################################
### Statistically discover the optimal number of clusters ###

kr <- 3:20

# --------------------------------------------
# Step 1: Compute Bray-Curtis distance matrix
# --------------------------------------------
bray_dist <- vegdist(as.matrix(proptab_leg2)^0.2, method = "bray")

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


######################
### Simple Heatmap ### 
# build a simple heatmap to explore all the data and visually assess how many clusters you want

## open pdf
pdf(file="OUTPUT/Leg2/leg2_simple_hm.pdf", width=35, height=20)

## simple heatmap
hm = heatmap.2(as.matrix(proptab_leg2)^.2, trace="none", col=viridis,
               distfun=function(x) vegdist(x, method="bray"),
               hclustfun=function(x) hclust(x, method="ward.D2"))

## see how many clusters you may have for samples (rows)
clusnum=6
oldclus = cutree(hm$rowDendrogram, k=clusnum)
oldorder = unname(rle(oldclus[as.hclust(hm$rowDendrogram)$order])$values)
neworder = (1:clusnum)
names(neworder) = oldorder
clus = unname(neworder[as.character(oldclus)])
names(clus) = names(oldclus)

## redraw the map with clusters colored
hm = heatmap.2(as.matrix(proptab_leg2)^.2, trace="none", key=F, margins=c(20,30), col=viridis,
               distfun=function(x) vegdist(x, method="bray"),
               hclustfun=function(x) hclust(x, method="ward.D2"),
               RowSideColors = rainbow(max(clus))[clus])

## close the pdf
dev.off()


######################
### ComplexHeatmap ###

# colnames(meta_leg2) # see which meta parameters you want to investigate

# Compute Bray-Curtis distance and hierarchical clustering
bray_dist <- function(x) vegdist(x, method = "bray")

k = 6 # number of clusters you decided on

tdist <- vegdist(as.matrix(proptab_leg2)^.2, method = "bray")
ttree <- hclust(tdist, method = "ward.D2")

# Make a dendrogram variable and color branches with Turbo
row_dend <- as.dendrogram(ttree) %>%
  dendextend::color_branches(., k = k, col = turbo(k))

#
## open pdf
pdf(file="OUTPUT/Leg2/leg2_complex_hm.pdf", width=33, height=20)

leg2_complex_hm <- function() {

ht_opt$ROW_ANNO_PADDING = unit(2, "mm")                                          # Padding between the heatmap and annotations

# Create heatmap with dendrogram
amun_complex_hm <- Heatmap(as.matrix(proptab_leg2)^.2,
                           name = "amun_complex_hm",
                           col = viridis(256),
                           
                           # Row clustering with pre-computed dendrogram
                           cluster_rows = row_dend %>%
                             set("branches_lwd", 3),                           # set branch line width for better visibility
                           row_split = k,
                           
                           # Column clustering
                           clustering_distance_columns = bray_dist,
                           clustering_method_columns = "ward.D2",
                           column_dend_reorder = TRUE,
                           
                           # Annotations
                           right_annotation = rowAnnotation(
                             bar1 = anno_points(log2(meta_leg2$"Depth [m]"),
                                                size = unit(3, "mm"),
                                                width = unit(3, "cm"),
                                                axis_param = list(
                                                  at = log2(c(2, 10, 50, 250, 2000)),
                                                  labels = c(2, 10, 50, 250, 2000),
                                                  gp = gpar(fontsize = 10))),
                             bar2 = anno_points(meta_leg2$"Temperature [deg C]",
                                                size = unit(3, "mm"),
                                                width = unit(3, "cm"),
                                                axis_param = list(gp = gpar(fontsize = 10))),
                             bar3 = anno_points(meta_leg2$"Salinity",
                                                size = unit(3, "mm"),
                                                width = unit(3, "cm"),
                                                axis_param = list(gp = gpar(fontsize = 10))),
                             bar4 = anno_points(meta_leg2$"Dissolved Oxygen [uM]",
                                                size = unit(3, "mm"),
                                                width = unit(3, "cm"),
                                                axis_param = list(gp = gpar(fontsize = 10))),
                             bar5 = anno_points(meta_leg2$"Fluorescence [mg/m^3]",
                                                size = unit(3, "mm"),
                                                width = unit(3, "cm"),
                                                axis_param = list(gp = gpar(fontsize = 10))),
                             bar6 = anno_points(meta_leg2$"Nitrate [uM]",
                                                size = unit(3, "mm"),
                                                width = unit(3, "cm"), 
                                                axis_param = list(gp = gpar(fontsize = 10))),
                             gap = unit(2, "mm"),
                             show_annotation_name = FALSE
                           ),
                           
                           # Rest of your parameters
                           heatmap_height = unit(40, "cm"),
                           heatmap_width = unit(65, "cm"),
                           
                           row_dend_width = unit(6, "cm"),
                           column_dend_height = unit(2, "cm"),
                           
                           row_title = "Samples",
                           row_title_gp = gpar(fontsize = 25, fontface = "bold"),
                           
                           column_title = "Taxa",
                           column_title_gp = gpar(fontsize = 25, fontface = "bold"),
                           
                           row_names_gp = gpar(fontsize = 11, fontface = "bold"),
                           
                           column_names_gp = gpar(fontsize = 1.2, fontface = "bold"), 
                           column_names_rot = 45,
                           
                           heatmap_legend_param = list(
                             title = "Fifth-root transformed relative abundances of microbial taxa during Leg 2 of the 2023 Amundsen expedition\n",
                             title_gp = gpar(fontsize = 30, fontface = "bold"),
                             legend_direction = "horizontal",
                             legend_width = unit(2, "cm"),
                             color_bar = "continuous")
)

#
## Draw the complex heatmap
draw(amun_complex_hm, heatmap_legend_side = "top")

#
## Add decorations

# Fix annotation bars
annotation_labels <- c( "Depth, m",                                              # Define parameter labels in order (bar1, bar2, etc.)
                        "Temp, oC", "Sal",
                        "O2, uM", "Fluo, ug/l", "NO3, uM")
for (i in 1:6) {
  decorate_annotation(paste0("bar", i), {
    grid.text(annotation_labels[i],
              y = unit(1, "npc") + unit(2, "mm"), just = "bottom",               # Position label above the annotation
              gp = gpar(fontsize = 17, fontface = "bold"))
  })
}

# Color the cluster row labels (samples) from top/first slice (later named as "Cluster 1") to the bottom/last (later named"Cluster k")
cluster_palette <- turbo(k, alpha = 0.5)            # generate Turbo colors for clusters to color from top to bottom
for (i in 1:k) {                                    # loop over cluster slices
  decorate_row_names("amun_complex_hm", {
    grid.rect(gp = gpar(fill = cluster_palette[i]),
              width = unit(13.3, "cm"),
              x = unit(1, "npc") + unit(0.6, "cm"),
              just = "center")
  }, slice = i)
}
}

leg2_complex_hm()

## close the pdf
dev.off()


########################################
### Cluster order + add to metadata ####
# reorder clusters in accordance with the heatmap order
tclus_reordered <- cutree(ttree, k = k) %>%
  setNames(rownames(proptab_leg2)) %>%  # Assign sample names
  (\(tclus) tclus %>%
     factor(levels = unique(tclus[order.dendrogram(row_dend)]), labels = 1:k) %>%
     as.numeric() %>%
     setNames(names(tclus)))()

meta_leg2$Cluster <- tclus_reordered[match(meta_leg2$plate.geno, sub(".*; ", "", names(tclus_reordered)))]


###########################
### Circular Dendrogram ###

# Create cluster data frame for annotation
cluster_data <- data.frame(label = names(tclus_reordered),
                           cluster = tclus_reordered,
                           stringsAsFactors = FALSE)
phylo <- as.phylo(ttree)

# Open PDF with optimized dimensions
pdf(file = "OUTPUT/Leg2/leg2_circ_dend.pdf", width = 8, height = 8)

cluster_colors <- setNames(turbo(k), as.character(1:k))  # adjust cluster colors, the same one is defined in the heatmap and db-RDA

# Create base circular tree with compact layout
p <- ggtree(phylo, layout = 'circular', ladderize = TRUE, size = 0.8) %<+% cluster_data +
  geom_tree(aes(color = factor(cluster))) +                           # color branches by cluster
  geom_tiplab(aes(color = factor(cluster)),                           # label spacing
              size = 2,                        # label size
              offset = 0.5,                      # pushing labels outward
              align = TRUE,                      # prevent overlap
              show.legend = TRUE) +
  scale_color_manual(name = "Cluster", values = cluster_colors,
                     guide = guide_legend(override.aes = list(size = 5))) +     # your predefined colors for clusters
  expand_limits(x = 5)                           # expands - smaller number means greater exppansion

print(p)
dev.off()


#####################################################################################################################
### db-RDA (distance-based Redundancy Analysis; in this case, we will base it on the Bray Curtis Distance matrix) ###
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
microbial_data <- as.matrix(proptab_leg2)^0.2

# Check how your ASV data looks like
# If the distribution is not normal - go back and modify.
# Here we have already done prop.table (by sample) and ^0.25 transformations on the original ASV matrix
hist(microbial_data)    # looks optimal: a lot of 0s but normal distribution to the left

# Load metadata (explanatory variables)
metadata <- meta_leg2

# Select numeric environmental variables correctly
env_vars <- metadata %>%
  select(where(is.numeric)) %>%
  select (-"Cluster",
          
          -"Cast",                                        # meaningless
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
          -"Chla [mg/m^3]",
          -"Bacteria [1/ml]",
          -"Cyanobacteria PE [1/ml]",
          -"Cyanobacteria PC [1/ml]",
          -"Eukaryotes [1/ml]",
          
          -"DOY",                                         # it has a big effect, but it is likely due to locations differing in taxa

          -"Oxygen solubility [ml/l]",                    # correlates with temperature based on ggpair() + useless
          -"Nitrate CTD [uM]",                            # correlates with bottle nitrate based on ggpair() + we have bottle nitrate

          -"Dissolved Oxygen [uM]",                       # has strong correlation with nitrate (ggpair())

          -"Light Attenuation Coefficient [1/m]",         # has weak significance based on anova.cca, specifically margin effect
          -"Light Transmission [%]",                      # has weak significance based on anova.cca
          -"Fluorescence [mg/m^3]",                       # has weak significance based on anova.cca
          -"CDOM [mg/m^3]",                               # has weak significance based on anova.cca
          
          -"Silicate [uM]"                                # has strong correlation with depth and other nutrients (ggpair())
          # -"Phosphate [uM]",                              # has strong correlation with nitrate (ggpair()) but acceptable (anova.cca)
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
complete_cases <- complete.cases(env_vars)                      # 97/97: all have nuts and we don't employ KEBABB data here
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

# Add ALL metadata parameters to plot by size
# metadata_filtered <- meta_leg2 %>%
#   mutate(Cluster = as.factor(tclus_reordered[match(ordered_samples, names(tclus_reordered))]))    # assign corrected clusters

# Extract RDA scores for plotting and add metadata to it
sc_si <- scores(rda_model, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  mutate(plate.geno = sapply(strsplit(Sample, "; "), tail, 1)) %>%
  left_join(metadata, by = "plate.geno") %>%
  mutate(Cluster = as.factor(Cluster))

# Define cluster colors
cluster_levels <- levels(as.factor(metadata$Cluster))  # Get unique cluster numbers
cluster_colors <- setNames(turbo(length(cluster_levels)), cluster_levels) # use turbo: better and can also use in ODV
# turbo: an improved Jet rainbow color map for reducing false detail, banding and color blindness ambiguity

# Compute dynamic arrow scaling factor
arrow_scaling_factor <- max(abs(c(sc_bp$dbRDA1, sc_bp$dbRDA2)), na.rm = TRUE) * 2

# Define proper axis limits
xlim_vals <- range(c(sc_si$dbRDA1, sc_bp$dbRDA1 * arrow_scaling_factor), na.rm = TRUE) * 1.3
ylim_vals <- range(c(sc_si$dbRDA2, sc_bp$dbRDA2 * arrow_scaling_factor), na.rm = TRUE) * 1.3

## Open PDF
pdf(file="OUTPUT/Leg2/leg2_db-RDA.pdf", width = 35, height = 20)

# Loop through all numeric variables in `sc_si` to generate additional db-RDA plots
numeric_vars <- c("DOY", "Bot. Depth [m]",
                  "Depth [m]", "volume.l", "dna.conc",
                  "Temperature [deg C]", "Salinity",
                  "Dissolved Oxygen [uM]", "Oxygen solubility [ml/l]",
                  "PAR [uE/s/m^2]", "Light Attenuation Coefficient [1/m]", "Light Transmission [%]",
                  "CDOM [mg/m^3]", "Fluorescence [mg/m^3]",              
                  "Nitrate CTD [uM]", "Nitrate [uM]", "Nitrite [uM]", "Phosphate [uM]", "Silicate [uM]")

for (param in numeric_vars) {
  
  param_min <- min(sc_si[[param]], na.rm = TRUE)
  param_max <- max(sc_si[[param]], na.rm = TRUE)
  
  param_plot <- ggplot() +
    
    geom_point(
      data = sc_si,
      aes(x = dbRDA1, y = dbRDA2,
          color = Cluster,
          size = !!sym(param)),
      alpha = 0.5, stroke = 3  # Strong outline for better visibility
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
    
    scale_color_manual(name = "Cluster", values = cluster_colors) +
    scale_size(range = scales::rescale(range(sc_si[[param]], na.rm = TRUE), to = c(2, 20))) +
    
    # Create custom legend for filtration type
    guides(
      color = guide_legend(
        # title = "Cluster",
        override.aes = list(size = 15)),
      size = guide_legend(order = 2)) +
    
    labs(title = paste0(
      "db-RDA plot for Leg 2 of the 2023 Amundsen expedition: ", param,
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


###################################
### Environmental cluster stats ###
# make a summary table for clusters
cluster_summary <- proptab_leg2 %>%
  as.data.frame() %>%
  cbind(meta_leg2) %>%
  group_by(Cluster) %>%
  summarize(
    Location_freq = names(which.max(table(geo_loc_name))),
    "Bottom depth, m" = round(mean(`Bot. Depth [m]`)),
    "Target depth" = names(which.max(table(target.depth))),
    "Depth, m" = paste(round(min(`Depth [m]`, na.rm = TRUE)), "-", round(max(`Depth [m]`, na.rm = TRUE))),  # Min-max depth
    "Depth_cor, m" = paste(round(sort(na.omit(`Depth [m]`)))[2], "-", round(sort(na.omit(`Depth [m]`), decreasing = TRUE))[2]),
    "Temp, oC" = round(mean(`Temperature [deg C]`), 2),
    "Sal" = round(mean(`Salinity`), 1),
    "Fluo, ug/l"  = round(mean(`Fluorescence [mg/m^3]`), 2),
    "Bacteria, 1/ml" = round(mean(`Bacteria [1/ml]`, na.rm=TRUE)),
    "Eukaryota, 1/ml" = round(mean(`Eukaryotes [1/ml]`, na.rm=TRUE)),
    "NO3, uM" = round(mean(`Nitrate [uM]`), 2),
    "NO2, uM" = round(mean(`Nitrite [uM]`), 2)) %>%
  print()


########################
### Cluster analysis ###

# Compute Spearman correlation matrix between ASVs
# Spearman is useful because it uses ranks rather than absolute values
adjacency_matrix <- cor(proptab_leg2, method = "spearman")
adjacency_matrix[abs(adjacency_matrix) < 0.6] <- 0            # Replace weak correlations with zero (set threshold)
diag(adjacency_matrix) <- 0                                   # Remove self-connections

# Create a network graph object
library(igraph)
network <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", weighted = TRUE)
E(network)$width <- E(network)$weight * 0.5                   # Set edge widths

library(Rtsne)
network_layout <- as.matrix(Rtsne(t(proptab_leg2)^0.2, perplexity=30, dims=2, check_duplicates=FALSE, 
                                  pca=TRUE, max_iter=2500, verbose=FALSE)$Y)

# Aggregate ASV abundances by cluster
asv_clus <- aggregate(proptab_leg2, by=list(meta_leg2$Cluster), FUN = mean) # mean ASV abundance by cluster
rownames(asv_clus) = asv_clus$Group.1                                       # row names set to cluster labels
asv_clus$Group.1 = NULL

# For each ASV, determine in which cluster it is the most dominant (aka the highest mean abundance)
# This helps determine which cluster each ASV is most associated with - the nodes (ASVs) are colored accordingly later
asv_dominant_clus <- rownames(asv_clus)[apply(asv_clus, 2, which.max)]

# Define colors for cluster
vertex_colors <- setNames(turbo(k), as.character(1:k))[asv_dominant_clus]
vertex_colors[is.na(vertex_colors)] <- "gray"                            # default color for unassigned ASVs

# Out of the ASVs most associated with each cluster, extract top 10 by relative abundance
# Notice: this is different from the most abundant ASVs per cluster!
# Some ASVs may be very abundant in the cluster but NOT associated with it as they are more abundant in other clusters
# Such ASVs may have bigger nods than Cluster-spec. ASVs, but the color is dif. - it is of the cluster they are more associated w/
top_asv_clus <- asv_clus %>%                                             # mean ASV abundance by cluster
  as.data.frame() %>%
  rownames_to_column("Cluster") %>%                                      # cluster as column
  pivot_longer(-Cluster, names_to = "ASV", values_to = "Abundance") %>%  # 2 columns: ASV and its mean Abundance per cluster
  left_join(
    tibble(ASV = colnames(asv_clus),                                     # a tibble mapping ASV names to their dominant cluster
           Cluster = as.character(asv_dominant_clus),                    # assigns dominant cluster for each ASV
           ASV_Number = gsub(".*; ASV_", "ASV_", colnames(asv_clus))),   # extract ASV numbers from names
    by = "ASV"                                                           
  ) %>%
  filter(Cluster.x == Cluster.y,                                         # keep only ASVs that are assigned to the correct cluster
         !grepl("Unknown|Domain", ASV),                                  # remove ASVs with "Unknown" or "Domain" in the name
         Abundance >= 0.01) %>%                                          # retain only ASVs with at least 1% relative abundance
  select(Cluster = Cluster.x, ASV, Abundance) %>%                        # remove intermediate columns
  mutate(Abundance = signif(Abundance, 2)) %>%                           # round to 2 significant figures
  arrange(Cluster, desc(Abundance)) %>%                                  # sort clusters & abundances in order
  print()

#
## Open PDF
pdf(file="OUTPUT/Leg2/leg2_cluster_analysis.pdf", width = 35, height = 20)

#
## Regenerate the complex heatmap
leg2_complex_hm()

#
## Generate the network map

# Loop through clusters 1 to k
# Circle/nod size - relative abundance of the ASV
# All ASVs are included in the layout, regardless of their primary cluster
# If sizable nods of other colors appear in the ASV network graph for a specific cluster, it is due to cross-cluster relationships
# If the cluster nodes are closely connected, it is highly self-contained, but if dif. colors are mixed - cross-cluster association

for (i in 1:k) {
  node_size <- sqrt(asv_clus[as.character(i),])

  # Extract only part of the taxonomic name
  # asv_labels <- gsub(".*; (ASV_\\d+)", "\\1", colnames(asv_clus))          # Extracts ASV_#"
  asv_labels <- gsub(".*; ([^;]+); ASV_\\d+$", "\\1", colnames(asv_clus))  # Extract Genus

  plot(network,
       layout = network_layout,
       vertex.size = node_size * 30,
       vertex.color = vertex_colors,
       # vertex.frame.color = NA,
       edge.width = E(network)$width,
       edge.color = "gray50",
       vertex.label = asv_labels,      # show only part of the name
       vertex.label.cex = 0.2)         # make the name smaller

  title(main = paste("Square-root transformed relative abundances of microbial taxa in Cluster", i, "(ASVs are colored by Cluster they are most associated with)"), cex.main = 3)
}

#
## Add the table with environmental parameter values specific to each cluster
grid.arrange(textGrob("Cluster-specific values of defining environmental parameters", 
                      gp = gpar(fontsize = 35, fontface = "bold"),
                      vjust = 18),                                                            # move down the title
             grid.table(cluster_summary, 
                        rows = NULL,                                                          # remove row numbers
                        theme = ttheme_default(
                          core = list(fg_params = list(cex = 2.2)),                           # cell text size
                          colhead = list(fg_params = list(cex = 2.2, fontface = "bold")))),   # column label size
             nrow = 2, heights = c(0.2, 1))

#
# Add the table with the most abundant ASVs specific to each cluster
grid.arrange(textGrob("Most abundant Cluster-specific ASVs\n(>1% mean relative abundance)", 
                      gp = gpar(fontsize = 28, fontface = "bold")),
             grid.table(top_asv_clus, 
                        rows = NULL,                                                          # remove row numbers
                        theme = ttheme_default(
                          core = list(fg_params = list(cex = 1.5)),                           # cell text size
                          colhead = list(fg_params = list(cex = 1.8, fontface = "bold")))),   # column label size
             nrow = 2, heights = c(0.2, 1))


## close PDF
dev.off()


##############################
### ODV for Leg 2 Amundsen ###

# should have ctd.nut.amun or ctd.nut.kebabb off of AMUN_CTD_NUT_KEBABB.R
ctd.nut.kebabb <- read.csv("OUTPUT/ctd.nut.kebabb.csv", stringsAsFactors = FALSE, check.names = FALSE)

write.csv(ctd.nut.kebabb %>% filter(Cruise == "KEBABB") %>%
            left_join(meta_leg2 %>%
                        # Remove overlapping columns except join keys
                        select(-any_of(names(ctd.nut.kebabb)[!names(ctd.nut.kebabb) %in% c("Station", "Cast", "Depth [m]", "Bottle")])),
      join_by(Station, Cast, `Depth [m]`, Bottle)),
      "OUTPUT/Leg2/leg2_ODV.csv", row.names = FALSE)


################################################################################
### ADDITIONAL ###

# ####################################
# ### ComplexHeatmap (alternative) ###
# # clusters are ordered properly here but the dendrogram is absent
# 
# colnames(meta_leg2) # see which meta parameters you can investigate
# 
# # #
# # # Compute Bray-Curtis distance and hierarchical clustering
# # #
# # Compute Bray-Curtis distance and hierarchical clustering
# bray_dist <- function(x) vegdist(x,method = "bray")
# 
# k = 6 # number of clusters you decided on
# 
# tdist <- vegdist(as.matrix(proptab_leg2)^.2, method = "bray")
# ttree <- hclust(tdist, method = "ward.D2")
# 
# # Assign cluster numbers
# tclus <- cutree(ttree, k = k)
# 
# # Order clusters so they appear in order from top to bottom
# ordered_samples <- names(sort(tclus))  # Get sample names in correct order
# 
# # Reorder the data matrix by the new order
# ordered_matrix <- (as.matrix(proptab_leg2)^.2)[ordered_samples, ]
# ordered_clusters <- tclus[ordered_samples]  # Match clusters to reordered samples
# #
# # Reorder the environmental variables by the new order
# ordered_meta <- meta_leg2 %>%
#   arrange(match(plate.geno, sapply(strsplit(rownames(ordered_matrix), "; "), tail, 1)))  # Reorder to match heatmap
# 
# #
# ## open pdf
# pdf(file="OUTPUT/Leg2/additional/leg2_complex_hm_no_dend.pdf", width=33, height=20)
# 
# ht_opt$ROW_ANNO_PADDING = unit(2, "mm")                                          # Padding between the heatmap and annotations
# 
# #
# ## create a complex heatmap
# amun_complex_hm <-
# 
#   Heatmap(ordered_matrix,
# 
#           # clustering_distance_rows = bray_dist,
#           clustering_distance_columns = bray_dist,
#           # clustering_method_rows = "ward.D2",                                  #passed to hclust
#           clustering_method_columns = "ward.D2",                                 #passed to hclust
# 
#           name = "amun_complex_hm",
#           col = viridis(256),
# 
#           right_annotation = rowAnnotation(bar1 = anno_points(log2(ordered_meta$"Depth [m]"),
#                                                               width = unit(3, "cm"),
#                                                               axis_param = list(
#                                                                 at = log2(c(2, 10, 60, 200, 600, 2000)), # Log2-transformed depth
#                                                                 labels = c(2, 10, 60, 200, 600, 2000),   # Original depth labels
#                                                                 gp = gpar(fontsize = 8))),
#                                            bar2 = anno_points(ordered_meta$"Temperature [deg C]",
#                                                               width = unit (3, "cm"),
#                                                               axis_param = list(
#                                                                 gp = gpar(fontsize = 8))),
#                                            bar3 = anno_points(ordered_meta$"Salinity",
#                                                               width = unit (3, "cm"),
#                                                               axis_param = list(
#                                                                 gp = gpar(fontsize = 8))),
#                                            bar4 = anno_points(ordered_meta$"Dissolved Oxygen [uM]",
#                                                               width = unit (3, "cm"),
#                                                               axis_param = list(
#                                                                 gp = gpar(fontsize = 8))),
#                                            bar5 = anno_points(ordered_meta$"Fluorescence [mg/m^3]",
#                                                               width = unit (3, "cm"),
#                                                               axis_param = list(
#                                                                 gp = gpar(fontsize = 8))),
#                                            bar6 = anno_points(ordered_meta$"Nitrate [uM]",
#                                                               width = unit (3, "cm"),
#                                                               axis_param = list(
#                                                                 gp = gpar(fontsize = 8))),
#                                            gap = unit(2, "mm"),                  # Gap between annotations
#                                            show_annotation_name = FALSE),        # Assign label later through decoration
# 
#           cluster_rows = FALSE,                                                  # To do clustering yourself (see ttree)
#           row_split = ordered_clusters,                                          # Split in the right order
#           # row_dend_reorder = TRUE,                                             # Rotate tree (if you don't do clustering yourself)
# 
#           column_dend_reorder = TRUE,
# 
#           heatmap_height = unit(40, "cm"),
#           heatmap_width = unit(65, "cm"),
#           row_dend_width = unit(4, "cm"),
#           column_dend_height = unit(2, "cm"),
# 
#           row_title = "Samples",
#           row_title_gp = gpar(fontsize = 25, fontface = "bold"),
#           # row_labels = rep("", nrow(proptab_kebabb)),    # If you want to remove row names (aka samples)
# 
#           column_title = "Taxa",
#           column_title_gp = gpar(fontsize = 25, fontface = "bold"),
# 
#           row_names_gp = gpar(fontsize = 10, fontface = "bold"),
#           column_names_gp = gpar(fontsize = 1.2, fontface = "bold"),
#           column_names_rot = 45,
# 
#           heatmap_legend_param = list(
#             title = "Fourth-root transformed relative abundances of microbial taxa during the 2023 Amundsen expedition\n",
#             title_gp = gpar(fontsize = 30, fontface = "bold"),
#             legend_direction = "horizontal",
#             legend_width = unit(2, "cm"),
#             color_bar = "continuos"),
# 
#           # show_heatmap_legend = FALSE                                   # Turn off heatmap legend b/c numbers are meaningless
#   )
# 
# #
# ## Draw the complex heatmap
# draw(amun_complex_hm,                                                     # First draw, then decorate, place legend in a good spot
#      heatmap_legend_side="top")                                           # Padding between row annotations
# 
# #
# ## Add decorations
# 
# # Define parameter labels for the annotations (see bar1, bar2, etc.)
# annotation_labels <- c(
#   "Depth, m",
#   "Temp-re, deg C",
#   "Salinity",
#   "Oxygen, uM",
#   "Fluorescence, ug/l",
#   "Nitrate, uM"
# )
# 
# # Loop over annotation bars
# for (i in 1:6) {
#   decorate_annotation(paste0("bar", i), {
#     grid.text(annotation_labels[i],
#               y = unit(1, "npc") + unit(2, "mm"), just = "bottom",               # Position label above the annotation
#               gp = gpar(fontsize = 12, fontface = "bold"))
#   })
# }
# 
# #
# ##
# 
# # Generate Turbo colors dynamically for clusters
# cluster_palette <- turbo(k, alpha = 0.5)  # k colors with transparency
# 
# # Loop over cluster slices
# for (i in 1:k) {
#   decorate_row_names("amun_complex_hm", {
#     grid.rect(gp = gpar(fill = cluster_palette[i]),  # Apply transparency
#               width = unit(13.3, "cm"),
#               x = unit(1, "npc") + unit(0.6, "cm"),
#               just = "center")
#   }, slice = i)
# }
# 
# ## close the pdf
# dev.off()

# ###########################################################
# ### Network analysis w/ the rest of the ASVs grayed out ###
# pdf(file="OUTPUT/Leg2/additional/leg2_network_highlight.pdf", width = 35, height = 20)
# 
# for (i in 1:k) {
#   node_size <- sqrt(asv_clus[as.character(i),])
# 
#   # Highlight only the current cluster, gray out others
#   cluster_mask <- asv_dominant_clus == as.character(i)
#   node_colors <- ifelse(cluster_mask, vertex_colors, "gray70")  # Gray out other clusters
# 
#   plot(network,
#        layout = network_layout,
#        vertex.size = node_size * 30,
#        vertex.color = node_colors,
#        vertex.frame.color = NA,
#        edge.width = E(network)$width,
#        edge.color = "gray50",  # Keep all edges visible
#        vertex.label = NA)
#   title(main = paste("Network Graph Highlighting Cluster", i), cex.main = 3)
# }
# 
# dev.off()


