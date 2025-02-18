# for simple heatmaps
library(gplots)
library(viridis)
library(vegan)

library(dendextend) # for cutting tree

library("ComplexHeatmap")

library(ggrepel) # for labeling with geom_text_repel
library(ggforce) # for geom_mark_ellipse

library(gridExtra)

# build Amundsen ASV by running AMUN_BIOINFO_ASV

meta.amun.small_exp <- meta.amun.small %>%
  filter (seq.filt == 0)

proptab.small.taxan.saman_exp <- proptab.small.taxan.saman %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  mutate(plate.geno = sapply(strsplit(Sample, " ; "), tail, 1)) %>%   # Extract last part as plate.geno
  filter(plate.geno %in% meta.amun.small_exp$plate.geno) %>%                 # Keep only matching samples
  select(-plate.geno) %>%                                             # Remove plate.geno column
  column_to_rownames("Sample") %>%                                    # Restore row names after filtering
  .[colSums(. > 0) > 2]                                               # select species present in more than 2 samples


#
### What does your ASV table look like in terms of data distribution? What normalization should you apply?
hist(as.matrix(proptab.small.taxan.saman_exp))         # looks pretty bad
hist(as.matrix(proptab.small.taxan.saman_exp)^0.5)     # looks better
hist(as.matrix(proptab.small.taxan.saman_exp)^0.25)    # looks optimal: a lot of 0s but normal distribution to the left


#####################
##### Simple Heatmap

#
### build a simple heatmap to explore all the data

## open pdf
pdf(file="OUTPUT/amun_simple_hm_exp.pdf", width=35, height=20)

## simple heatmap
hm = heatmap.2(as.matrix(proptab.small.taxan.saman_exp)^.25, trace="none", col=viridis,
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
hm = heatmap.2(as.matrix(proptab.small.taxan.saman_exp)^.25, trace="none", key=F, margins=c(20,30), col=viridis,
               distfun=function(x) vegdist(x, method="bray"),
               hclustfun=function(x) hclust(x, method="ward.D2"),
               RowSideColors = rainbow(max(clus))[clus])

## close the pdf
dev.off()


#######################
##### ComplexHeatmap

# colnames(meta.amun.small) # see which meta parameters you can investigate

#
## set the clustering
bray_dist <- function(x) vegdist(x,method = "bray")
tdist = vegdist(as.matrix(proptab.small.taxan.saman_exp)^.25, method="bray")
ttree = hclust(tdist, method = "ward.D2")
tclus = cutree(ttree, k = 7)

#
## open pdf
pdf(file="OUTPUT/amun_complex_hm_exp.pdf", width=33, height=20)

ht_opt$ROW_ANNO_PADDING = unit(2, "mm")                                          # Padding between the heatmap and annotations

#
## create a complex heatmap
amun_complex_hm_exp <- 
  
  Heatmap(as.matrix(proptab.small.taxan.saman_exp)^.25,
          
          # clustering_distance_rows = bray_dist,
          clustering_distance_columns = bray_dist,
          # clustering_method_rows = "ward.D2",                                  #passed to hclust
          clustering_method_columns = "ward.D2",                                 #passed to hclust
          
          name = "amun_complex_hm_exp",
          col = blues9,
          
          right_annotation = rowAnnotation(bar1 = anno_points(log2(meta.amun.small_exp$"Depth [m]"),
                                                              width = unit(3, "cm"),
                                                              axis_param = list(
                                                                at = log2(c(2, 10, 60, 200)), # Log2-transformed depth
                                                                labels = c(2, 10, 60, 200),   # Original depth labels
                                                                gp = gpar(fontsize = 8))),
                                           bar2 = anno_points(meta.amun.small_exp$"Temperature [deg C]",
                                                              width = unit (3, "cm"),
                                                              axis_param = list(
                                                                gp = gpar(fontsize = 8))),
                                           bar3 = anno_points(meta.amun.small_exp$"Salinity",
                                                              width = unit (3, "cm"),
                                                              axis_param = list(
                                                                gp = gpar(fontsize = 8))),
                                           bar4 = anno_points(meta.amun.small_exp$"Dissolved Oxygen [uM]",
                                                              width = unit (3, "cm"),
                                                              axis_param = list(
                                                                gp = gpar(fontsize = 8))),
                                           bar5 = anno_points(meta.amun.small_exp$"Fluorescence [mg/m^3]",
                                                              width = unit (3, "cm"),
                                                              axis_param = list(
                                                                gp = gpar(fontsize = 8))),
                                           bar6 = anno_points(meta.amun.small_exp$"Nitrate [uM]",
                                                              width = unit (3, "cm"),
                                                              axis_param = list(
                                                                gp = gpar(fontsize = 8))),
                                           gap = unit(2, "mm"),                  # Gap between annotations
                                           show_annotation_name = FALSE),        # Assign label later through decoration
          
          
          cluster_rows = as.dendrogram(ttree),                                   # To do clustering yourself (see ttree)
          row_split = 7,                                                         # Number of clusters (if you do clustering yourself)
          # row_dend_reorder = TRUE,                                             # Rotate tree (if you don't do clustering yourself)
          
          column_dend_reorder = TRUE,
          
          heatmap_height = unit(40, "cm"),
          heatmap_width = unit(65, "cm"),
          row_dend_width = unit(4, "cm"),
          column_dend_height = unit(2, "cm"),
          
          row_title = "Samples",
          row_title_gp = gpar(fontsize = 25, fontface = "bold"),
          # row_labels = rep("", nrow(proptab.small.taxan.saman_exp)),               # If you want to remove row names (aka samples)
          
          column_title = "Taxa",
          column_title_gp = gpar(fontsize = 25, fontface = "bold"),
          
          row_names_gp = gpar(fontsize = 8.5, fontface = "bold"),
          column_names_gp = gpar(fontsize = 1.2, fontface = "bold"),
          column_names_rot = 45,
          
          heatmap_legend_param = list(
            title = "Fourth-root transformed relative abundances of microbial taxa during the 2023 Amundsen expedition
(non-sequential filtration)\n",
            title_gp = gpar(fontsize = 30, fontface = "bold"),
            legend_direction = "horizontal",
            legend_width = unit(2, "cm"),
            color_bar = "continuos"),
          
          # show_heatmap_legend = FALSE                                            # Turn off heatmap legend b/c numbers are meaningless  
  )

#
## Draw the complex heatmap
draw(amun_complex_hm_exp,                                                        # First draw, then decorate, place legend in a good spot
     heatmap_legend_side="top")                                                  # Padding between row annotations

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
  decorate_row_names("amun_complex_hm_exp", {
    grid.rect(gp = gpar(fill = paste0(colors[i], "40")),                         # Apply mild transparency
              width = unit(13.3, "cm"),                                           # Width of the colored box
              x = unit(1, "npc") + unit(0.6, "cm"),                              # Play with this to shift the box onto the labels
              just = "center")
  }, slice = i)
}

## close the pdf
dev.off()


########################################
##### db-RDA (based on Bray Curtis Distance matrix)
# to identify how environmental variables influence microbial relative abundance

# Load metadata (explanatory variables)
metadata <- meta.amun.small_exp

# Load microbial ASV data (response variables)
microbial_data <- as.matrix(proptab.small.taxan.saman_exp)^0.25

# Select numeric environmental variables correctly
env_vars <- metadata %>%
  select(where(is.numeric)) %>%
  select(all_of(c("Depth [m]",
                  # "Month",
                  # "seq.filt", "pore.size", "volume.l", "dna.conc",
                  "Temperature [deg C]", "Salinity", "Dissolved Oxygen [uM]", "PAR [uE/s/m^2]",
                  "CDOM [mg/m^3]", "Fluorescence [mg/m^3]", "Nitrate [uM]", "Nitrite [uM]",
                  "Phosphate [uM]", "Silicate [uM]"))) %>%
  rename_with(~ str_replace_all(., "\\s*\\[.*?\\]", "")) %>%                              # Remove brackets
  rename(Oxygen = `Dissolved Oxygen`)                                                     # Rename Oxygen

# Identify complete cases (rows with no NA values)
complete_cases <- complete.cases(env_vars)

# Filter both datasets to include only complete cases
microbial_data <- microbial_data[complete_cases, , drop = FALSE]
env_vars <- env_vars[complete_cases, , drop = FALSE]

# Convert column names to valid R variable names
colnames(env_vars) <- make.names(colnames(env_vars))

# Compute Bray-Curtis dissimilarity
bray_curtis_dist <- vegdist(microbial_data, method = "bray")

# Perform db-RDA (distance-based redundancy analysis) with all variables
rda_model <- dbrda(bray_curtis_dist ~ ., data = env_vars)
summary(rda_model)

# Print proportion of variance explained
variance_partitioning <- summary(rda_model)$constr.chi / summary(rda_model)$tot.chi
print(paste("Variance explained by environmental variables:", round(variance_partitioning * 100, 2), "%"))

# Compute adjusted R²
rda_adj_r2 <- RsquareAdj(rda_model)$adj.r.squared
print(paste("Adjusted R²:", round(rda_adj_r2, 2)))

# Perform per-variable ANOVA to extract significant constraints
anova_per_var <- anova.cca(rda_model, by = "terms", permutations = 10000)

# Extract significant environmental variables (P < 0.001) and fix names
significant_env <- rownames(anova_per_var)[anova_per_var$`Pr(>F)` < 0.001]
significant_env <- significant_env[!is.na(significant_env)]                     # Remove NAs
valid_vars <- make.names(significant_env)  # Convert to valid R names
valid_vars <- valid_vars[valid_vars %in% colnames(env_vars)]                    # Ensure they exist in env_vars

# Run final db-RDA with only significant variables
final_formula <- as.formula(paste("bray_curtis_dist ~", paste(valid_vars, collapse = " + ")))
final_rda_model <- dbrda(final_formula, data = env_vars %>% select(all_of(valid_vars)))

# Print final db-RDA summary
summary(final_rda_model)

# Print proportion of variance explained
final_variance_partitioning <- summary(final_rda_model)$constr.chi / summary(rda_model)$tot.chi
print(paste("Variance explained by environmental variables:", round(final_variance_partitioning * 100, 2), "%"))

# Print final adjusted R²
final_rda_adj_r2 <- RsquareAdj(final_rda_model)$adj.r.squared
print(paste("Adjusted R²:", round(final_rda_adj_r2, 3)))

# Print significance of each variable
final_anova_results <- anova.cca(final_rda_model, permutations = 10000, by = "term") %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  select(Variable, F, `Pr(>F)`) %>%
  mutate(
    F = round(F, 2), 
    `Pr(>F)` = sprintf("%.4f", pmin(`Pr(>F)`, 0.9999))) %>%                   # Ensures 4 decimals & prevents 1.0000
  print()


# # Print significance of the final model
# anova.cca(final_rda_model, permutations = 10000) %>%
#   as.data.frame() %>%
#   slice(1) %>%
#   mutate(
#     F = round(F, 2), 
#     `Pr(>F)` = sprintf("%.4f", pmin(`Pr(>F)`, 0.9999))) %>%                   # Ensures 4 decimals & prevents 1.0000
#   transmute(Result = paste("F-value:", F, "| p-value:", `Pr(>F)`)) %>%
#   pull(Result) %>%
#   print()

# # Print significance of the first two dbRDA axis
# anova.cca(final_rda_model, permutations = 10000, by = "axis") %>%
#   as.data.frame() %>%
#   mutate(
#     Percentage = round((SumOfSqs / sum(SumOfSqs)) * 100, 2),                   # Calculate variance explained
#     F = round(F, 2),
#     `Pr(>F)` = sprintf("%.4f", pmin(`Pr(>F)`, 0.9999))) %>%                    # Ensures 4 decimals & prevents 1.0000
#   filter(grepl("dbRDA", rownames(.))) %>%                                      # Keep only dbRDA axes
#   transmute(
#     Result = paste(
#       "F-value:", F,
#       "| p-value:", `Pr(>F)`,
#       "| Variance explained:", Percentage, "%")) %>%
#   pull(Result) %>%
#   paste(collapse = "\n") %>%
#   cat("\n")

# Extract RDA scores for plotting
rda_scores <- scores(final_rda_model, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Sample")

# Extract environmental variable loadings (arrows)
env_scores <- scores(final_rda_model, display = "bp") %>%
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
    geo_loc_name = str_remove(metadata_filtered$geo_loc_name, "^Canada: ")
  ) %>%
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

# Compute variance explained by axes
perc <- round(100 * (summary(final_rda_model)$cont$importance[2, 1:2]), 2)

# Extract site scores
sc_si <- scores(final_rda_model, display = "sites", choices = c(1,2), scaling = 1) %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  left_join(select(rda_scores, Sample, geo_loc_name, Sequential_Filtration, Pore_Size), by = "Sample") %>%
  select(Sample, dbRDA1, dbRDA2, geo_loc_name, Sequential_Filtration, Pore_Size)

# Extract scores for environmental variables (arrows)
sc_bp <- scores(final_rda_model, display = "bp", choices = c(1,2), scaling = 1) %>%
  as.data.frame() %>%
  rownames_to_column("Variable")

# Define proper axis limits
xlim_vals <- range(c(sc_si$dbRDA1, sc_bp$dbRDA1 * arrow_scaling_factor), na.rm = TRUE) * 1.3
ylim_vals <- range(c(sc_si$dbRDA2, sc_bp$dbRDA2 * arrow_scaling_factor), na.rm = TRUE) * 1.3

### Plot using ggplot2

## Open PDF
pdf(file="OUTPUT/amun_db-RDA_exp.pdf", width = 35, height = 20)

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
      "db-RDA plot for Amundsen 2023: non-sequential filtration\n",
      "(variance explained by the environmental variables: ", round(final_variance_partitioning * 100, 2), "%; ",
      "adjusted R²: ", round(final_rda_adj_r2, 2), ")"
    ),
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

rda_grob <- ggplotGrob(rda_plot) # Convert ggplot object to grob for compatibility
grid.draw(rda_grob)              # Print plot on first page


## Prepare a table with ANOVA results
# Here, we show how well environmental variables explain the Bray-Curtis distance matrix of the taxonomic relative abundance
anova_table <- final_anova_results %>%
  select(Variable, F, `Pr(>F)`) %>%
  filter(Variable != "Residual") %>%  # Only keep environmental variables
  mutate(
    F = round(F, 2),
    `Pr(>F)` = sprintf("%.4f", pmin(as.numeric(`Pr(>F)`), 0.9999))  # Ensure 4 decimal places
  )


# Define custom theme with grid lines for a professional look
table_theme <- ttheme_minimal(
  base_size = 70,   # Adjust text size for readability
  core = list(fg_params = list(fontface = "plain")),   # Regular font
  colhead = list(fg_params = list(fontface = "bold")), # Bold headers
  padding = unit(c(10, 10), "mm"),  # Adjust padding for spacing
  core = list(bg_params = list(fill = c("white", "grey95"), col = "black")),  # Alternating row colors
  colhead = list(bg_params = list(fill = "grey80", col = "black"))  # Header background
)

grid.newpage()                                                                           # Start a new page for the table
grid.draw(tableGrob(anova_table, rows = NULL, theme = table_theme))                      # Print the table

## Close PDF
dev.off()














