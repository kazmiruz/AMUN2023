# 3A_CLUS.R --------------------------------------------------------------------------------------
#
# Purpose:
# Cluster paired-marker samples using fourth-root RA, Bray-Curtis dissimilarity and Ward.D2.
# Retain seven broad clusters and one supported subcluster; identify prevalent cluster and cluster-pair bioindicators.
# Save reusable clustering and indicator results and draw the supplementary heatmap.
#
# Inputs (paths are relative to the project working directory containing EXT, OUT and DNA):
# EXT/meta.amun.rds; EXT/proptab.amun.rds (from 2A; loaded here)
# EXT/indicator_res.rds (cached 9,999-permutation multipatt result; optional regeneration block below)
#
# Outputs:
# EXT/amun_clustering.rds -> cluster assignments, tree and colours for 4A-6A
# EXT/amun_bioindicators.rds -> tested patterns and prevalent ASVs for 4A-6A
# OUT/3A_CLUS/simple_hm.pdf; amun_complex_hm.pdf / .png (Figure S1); amun_dendro_pure.pdf / .png
#
# Run sections in order within this script. Later scripts load saved inputs; they do not need earlier scripts in memory.
# Heavy calculations remain optional or cached. Existing EXT and OUT folders are assumed to exist.
#---------------------------------------------------------------------------------------------------

library(tidyverse)
library(gplots)                          # diagnostic heatmap.2
# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)                  # annotated heatmap and pure dendrogram
library(viridis)                         # heatmap colours
library(cluster)                         # silhouette widths
library(dendextend)                      # rotate and colour dendrograms
library(vegan)                           # Bray-Curtis dissimilarity
library(indicspecies)                    # indicator taxa analysis
options("lifecycle_verbosity"="warning") # show lifecycle warnings every time they occur
# -------------------------------------------------------------------------

### ### Interpretation notes

#   Since microbes are the main drivers of biogeochemical cycles, community composition data can directly explain
#   and predict environmental shifts better than biochemical parameters alone. This mirrors medical diagnostics:
#   doctors would rather test the patient to identify the parasite than assign treatment based on symptoms alone.
#   Moreover, microbes act as sensitive environmental sensors - much like clinical biomarkers for early cancer detection.
#   Finally, microbes can be used as tracers of water masses, complementing traditional parameters like T and S.

# - EXPLORE the data with simple figures and tables before deciding which parameters contribute to the story.
# - Prioritize NOVEL or UNEXPECTED patterns over well-established differences.
# - Focus on findings with IMPLICATIONS for researchers, modelers, and northern communities.
# - Use each figure to COMMUNICATE a clear message and support the overall manuscript narrative.
# -------------------------------------------------------------------------


##############################
### LOAD CLUSTERING INPUTS ###
# Both markers contribute their non-rare ASV abundances to sample clustering.
# Uses: Saved outputs of 2A; no objects from 1A or 2A are needed in memory.
# Produces: meta.amun and proptab.amun, used below.


meta.amun <- readRDS("EXT/meta.amun.rds")
proptab.amun <- readRDS("EXT/proptab.amun.rds")


#######################
### SAMPLE CLUSTERS ###

### ### Select an abundance transformation
# Transformations reduce the influence of dominant ASVs; Bray–Curtis does not require normally distributed abundances
# hist(as.matrix(proptab.amun))
# hist(as.matrix(proptab.amun)^0.5)
# hist(as.matrix(proptab.amun)^0.25)
# hist(as.matrix(proptab.amun)^0.125)

# Fourth root selected after comparing the distributions
proptab.trans <- as.matrix(proptab.amun)^0.25

# Construct the sample dendrogram
dist_row <- vegan::vegdist(proptab.trans, method = "bray")
tree_row <- hclust(dist_row, method = "ward.D2")
dend_row <- as.dendrogram(tree_row)
sample_order <- attr(dist_row, "Labels")

# Rotate branches toward decreasing latitude; this changes only their display orientation
ord_row <- rownames(meta.amun)[order(meta.amun$`Latitude [deg]`, decreasing = TRUE)]
dend_row_order <- dendextend::rotate(dend_row, ord_row)
row_ordered <- labels(dend_row_order)


### ### Select the broad cluster solution
# Compare mean silhouette widths of solutions with 5–20 clusters
silhouette_k <- map_dfr(5:20, \(k) {
  groups <- cutree(tree_row, k); sil <- cluster::silhouette(groups, dist_row)[, 3]
  tibble(k, Mean_sil = mean(sil), Min_group_sil = min(tapply(sil, groups, mean)),
         Neg_frac = mean(sil < 0), Min_n = min(table(groups)))
})

# Plot and print the results
p_silhouette <- ggplot(silhouette_k, aes(k, Mean_sil)) + geom_line(linewidth = 0.8, colour = "#3157E2") +
  geom_point(size = 2.8, colour = "#3157E2") + scale_x_continuous(breaks = 5:20) +
  labs(x = "Number of clusters", y = "Mean silhouette width") + theme_bw(base_size = 14)
print(silhouette_k, n = Inf)
# print(p_silhouette)  # inspect interactively when choosing the broad cut

# Seven main clusters selected after inspecting the silhouette profile
base_k <- 7L

# Create the temporary seven-cluster solution and identify clusters containing at least 20 samples
base_cluster <- cutree(tree_row, k = base_k)
base_n <- table(base_cluster); print(base_n)
big_parents <- as.integer(names(base_n[base_n >= 20]))


### ### Test the next split of each large cluster
score_split <- function(parent) {
  members <- names(base_cluster)[base_cluster == parent]
  
  # Find the first finer dendrogram cut that divides this parent
  split_k <- purrr::detect(seq.int(base_k + 1, nrow(proptab.trans)),
                           \(k) dplyr::n_distinct(cutree(tree_row, k)[members]) > 1)
  child <- cutree(tree_row, split_k)[members]
  
  # Treat the smaller child as the candidate split-off branch
  branch <- members[as.character(child) == names(sort(table(child)))[1]]
  remainder <- setdiff(members, branch)
  
  # Evaluate the resulting base_k + 1 partition
  proposed <- setNames(paste0("Parent_", base_cluster), names(base_cluster))
  proposed[members] <- paste0("Child_", child)
  sil <- setNames(cluster::silhouette(as.integer(factor(proposed[sample_order])), dist_row)[, 3], sample_order)
  
  tibble(Parent = parent, Split_k = split_k, Parent_n = length(members), Branch_n = length(branch), Remain_n = length(remainder),
         Branch_sil = mean(sil[branch]), Remain_sil = mean(sil[remainder]), Members = list(branch))
}

split_scores <- map_dfr(big_parents, score_split) %>% arrange(desc(Branch_sil))
print(split_scores %>% select(-Members), n = Inf)
# Parent - original cluster being examined for subdivision; Split_k- first dendrogram cut at which that parent divides into two branches;
# Parent_n - number of samples in the original parent cluster; Branch_n -	samples in the smaller, candidate split-off branch;
# Remain_n - samples remaining in the parent cluster after removing that branch; Branch_sil	- mean silhouette width of the candidate split-off branch;
# Remain_sil - mean silhouette width of the remaining samples.

# Begin with the initial seven-cluster cut; retain the single best-supported candidate branch when its mean silhouette is >=0.50.
selected_splits <- split_scores %>% filter(Branch_sil >= 0.50) %>% slice_head(n = 1)
if (nrow(selected_splits) != 1L) stop("The specified eight-cluster refinement was not supported.")

final_id <- setNames(paste0("Parent_", base_cluster), names(base_cluster))
final_id[selected_splits$Members[[1]]] <- "Selected_1"
custom_num <- match(final_id, unique(final_id[row_ordered]))
k <- n_distinct(custom_num)
clus <- setNames(factor(custom_num, levels = seq_len(k), labels = as.character(seq_len(k))), names(final_id))
cl_labels <- paste0("CL", seq_len(k))
meta.amun$Cluster <- clus[rownames(meta.amun)]

# Report final cluster sizes and mean silhouette widths
final_sil <- setNames(cluster::silhouette(as.integer(clus[sample_order]), dist_row)[, 3], sample_order)
cluster_summary <- tibble(Cluster = as.character(clus[sample_order]), Sil = unname(final_sil)) %>%
  group_by(Cluster) %>% summarise(n = n(), Mean_sil = mean(Sil), .groups = "drop")
print(cluster_summary, n = Inf)

# What just happened:
# Samples were hierarchically clustered using Bray–Curtis dissimilarity and Ward.D2 linkage. Among solutions with 5-20 clusters,
# the seven-cluster dendrogram cut had the highest mean silhouette width (0.19). Clusters containing at least 20 samples were tested
# for one additional dendrogram split, resulting in the identification of an eighth cluster with a high mean silhouette width (0.50).

# Colour the final sample dendrogram
cluster_palette <- setNames(c("#0066FF", "#00DFFF", "#E60000", "#ED8A00", "#FFD500", "#D17DFF", "#00A823", "#00FF66"), levels(clus))
dend_row_order <- dend_row_order %>% dendextend::color_branches(clusters = as.integer(clus[row_ordered]), col = unname(cluster_palette)) %>%
  dendextend::set("branches_lwd", 3)

# Construct and rotate the ASV dendrogram by abundance-weighted latitude
dist_col <- vegan::vegdist(t(proptab.trans), method = "bray"); tree_col <- hclust(dist_col, method = "ward.D2")
taxon_latitude <- apply(proptab.trans, 2, \(x) weighted.mean(meta.amun$`Latitude [deg]`, w = x))
dend_col_order <- dendextend::rotate(as.dendrogram(tree_col), names(sort(taxon_latitude, decreasing = TRUE))) %>%
  dendextend::set("branches_lwd", 3)
col_ordered <- labels(dend_col_order)

# Store separator positions for ComplexHeatmap and heatmap.2
row_cl <- as.character(clus[row_ordered]); row_sep <- which(head(row_cl, -1) != tail(row_cl, -1))
row_sep_y <- 1 - row_sep / length(row_ordered)
hm_rows <- rev(row_ordered); hm_cl <- as.character(clus[hm_rows])
hm_rowsep <- which(head(hm_cl, -1) != tail(hm_cl, -1))


#########################################
### SAVE CLUSTERING FOR LATER SCRIPTS ###
# Keep the fitted tree, final cluster assignments, order and colours together.
# Uses: The clustering calculated above.
# Produces: EXT/amun_clustering.rds -> 4A, 5A and 6A.


saveRDS(list(clus = clus, cl_labels = cl_labels, cluster_palette = cluster_palette, k = k,
             dend_row_order = dend_row_order, row_ordered = row_ordered,
             tree_row = tree_row, cluster_summary = cluster_summary, silhouette_k = silhouette_k,
             split_scores = split_scores %>% select(-Members)), "EXT/amun_clustering.rds")

# Diagnostic heatmap for inspection
pdf("OUT/3A_CLUS/simple_hm.pdf", width = 35, height = 20)
heatmap.2(proptab.trans, Rowv = dend_row_order, Colv = dend_col_order, trace = "none", key = FALSE, margins = c(20, 30),
          col = viridis(256), rowsep = hm_rowsep, sepcolor = "white", sepwidth = c(0.003, 0.003),
          RowSideColors = unname(cluster_palette[as.character(clus[rownames(proptab.trans)])]))
dev.off()


##########################################################
### BIOINDICATORS OF SINGLE CLUSTERS AND CLUSTER PAIRS ###
# Use untransformed non-rare RA; save the tested patterns for Figures S1 and 5 and Results 3.3.

# Use meta.amun$Cluster because its sample names align with proptab.amun
# Run on untransformed relative abundances, allow single clusters and cluster pairs, and save the result for later sessions.
# set.seed(161193)
# indicator_res <- multipatt(proptab.amun, meta.amun$Cluster, func = "IndVal.g",
#                            duleg = FALSE, # If TRUE, Cluster combos aren't considered, only the original Clusters (1 to k, Dufrêne & Legendre, 1997)
#                            max.order = 2, control = how(nperm = 9999))
# attr(indicator_res, "input_hash") <- rlang::hash(list(proptab.amun, as.character(meta.amun$Cluster)))
# write_rds(indicator_res, "EXT/indicator_res.rds")       # eight clusters and their pairs
indicator_res <- read_rds("EXT/indicator_res.rds")
# Check abundance values as well as ASV/sample names and cluster assignments.
if (!identical(attr(indicator_res, "input_hash"), rlang::hash(list(proptab.amun, as.character(meta.amun$Cluster)))))
  stop("Indicator inputs changed or the cache has no input fingerprint; regenerate EXT/indicator_res.rds using the optional block above.")

# # Explore the indicspecies output to understand the concept
# indicator_res$cluster   # sample clusters
# indicator_res$comb      # candidate single-cluster and cluster-pair patterns
#
# # A (specificity) = sum of the ASV's mean relative abundances in the target cluster(s) /
# #                   sum of its mean relative abundances across all individual clusters.
# # B (fidelity) = samples containing the ASV in the target cluster(s) /
# #                all samples in the target cluster(s).
# # stat = sqrt(A * B)
#
# plot(density(indicator_res$A, from = 0, na.rm = TRUE), col = "blue", main = "A vs B",
#      xlim = range(c(indicator_res$A, indicator_res$B), na.rm = TRUE),
#      ylim = c(0, max(density(indicator_res$A, na.rm = TRUE)$y, density(indicator_res$B, na.rm = TRUE)$y)))
# lines(density(indicator_res$B, from = 0, na.rm = TRUE), col = "red")
#
# # $sign gives some important values: s.1, .. - group preferred by each ASV; index - cluster comb-n; stat - sq.root(IndVal); p.value - sign-ce of stat
# indicator_res$sign
# plot(colMeans(proptab.amun), indicator_res$sign$stat, log = 'xy')     # for each ASV, plot Mean RelAb (cols = ASVs) vs stat (same order of ASVs)
# abline(h = 0.6, col = "red", lwd = 2)                                 # most highly abundant bioindicators have stat >0.6
# plot(indicator_res$sign$p.value, indicator_res$sign$stat, log = 'xy') # p-value and stat correlate
# abline(h = 0.6, col = "red", lwd = 2)
#
# # Because the clusters were derived from these ASVs, p-values describe associations within this clustering.
# hist(indicator_res$sign$p.value, breaks = 1000)

### ### Select good indicators
rn <- rownames(indicator_res$A); cn <- colnames(indicator_res$A); idx <- indicator_res$sign$index

# # Dufrêne & Legendre (1997) classic threshold is IndVal > 0.25, which in our case means stat > 0.5 (remember: stat = sqrt(IndVal) = sqrt(A*B)).
# # One can also be more specific (and stringent) and set up A >= 0.5 and B >= 0.5, guaranteeing IndVal >= 0.25 and stopping A or B from dropping too low.
# # For example:
# valid <- indicator_res$sign$p.value < 0.01 & indicator_res$A[cbind(seq_along(rn), idx)] >= 0.5 & indicator_res$B[cbind(seq_along(rn), idx)] >= 0.5

# For this dataset, a more stringent stat threshold can be used: retain bioindicators with stat > 0.6 and p < 0.01.
# Beware that p-values generally decrease as stat increases.
# I do not impose separate A and B cutoffs b/c some bioindicators may not show up in many samples but still be relatively abundant and important in the cluster.
# Jacquemot et al. (2024, https://doi.org/10.1002/edn3.70028) also used an stat >=0.6 and p <=0.01. They called stat "IndVal score". 
valid <- !is.na(indicator_res$sign$p.value) & indicator_res$sign$p.value < 0.01 & indicator_res$sign$stat > 0.6

specific_asvs <- rn[valid]
chosen_group <- setNames(rep(NA_character_, length(rn)), rn)
chosen_group[valid] <- cn[idx[valid]]

# For pair-associated indicators, assign the ASV for display to the cluster w/ the highest mean relative abundance; retain the tested pair for interpretation.
single_cols <- cn[!str_detect(cn, "[-+]")]
sample_cluster <- as.character(meta.amun$Cluster[match(rownames(proptab.amun), rownames(meta.amun))])
mean_rel_ab <- sapply(single_cols, \(g) colMeans(proptab.amun[sample_cluster == g, , drop = FALSE], na.rm = TRUE))

chosen_cluster <- map2_chr(chosen_group, rn, \(group, asv) {
  if (is.na(group) || !str_detect(group, "[-+]")) return(group)
  pair <- intersect(str_split(group, "[-+]")[[1]], single_cols)
  if (!length(pair)) return(NA_character_)
  pair[which.max(mean_rel_ab[asv, pair])]
}) %>% set_names(rn)


# Identify prevalent ASVs (>=2% RelAb in at least 2 samples), including prevalent bioindicators
# This is a simple and unbiased logic, compared to gatekeeping ASVs according to their average abundance in the cluster or group of clusters they peak in
thr_display <- 0.02; min_samples_display <- 2L
display_n <- colSums(as.matrix(proptab.amun) >= thr_display, na.rm = TRUE)
asv_in_dotplot <- intersect(names(display_n)[display_n >= min_samples_display], specific_asvs)
proptab.display <- proptab.amun[, asv_in_dotplot, drop = FALSE]

selected_rows <- match(asv_in_dotplot, rn)
selected_patterns <- idx[selected_rows]
selected_asv_tbl <- tibble(ASV = asv_in_dotplot, Pattern = unname(chosen_group[asv_in_dotplot]), cluster = unname(chosen_cluster[asv_in_dotplot]),
                           A = indicator_res$A[cbind(selected_rows, selected_patterns)], B = indicator_res$B[cbind(selected_rows, selected_patterns)],
                           stat = indicator_res$sign$stat[selected_rows], p_value = indicator_res$sign$p.value[selected_rows],
                           number_of_samples_passing_2pct = as.integer(display_n[asv_in_dotplot])) %>% arrange(cluster, desc(stat), p_value)
options(scipen = 999)
print(selected_asv_tbl, n = Inf)
# Save the tested pattern separately from the cluster used to order colours and labels in the figures.
saveRDS(list(specific_asvs = specific_asvs, chosen_group = chosen_group, chosen_cluster = chosen_cluster,
             asv_in_dotplot = asv_in_dotplot, selected_asv_tbl = selected_asv_tbl), "EXT/amun_bioindicators.rds")


#################################
### ASV HEATMAP FOR FIGURE S1 ###
# Uses the fitted sample and ASV trees and indicator selection above; creates the supplementary heatmap.

draw_amun_complex_hm <- function() {
  ht_opt$DIMNAME_PADDING <- unit(0.2, "mm")
  ht_opt$ROW_ANNO_PADDING <- unit(2.5, "mm")
  Heatmap(proptab.trans,
          name = "amun_complex_hm",
          col = viridis(256),
          cluster_rows = dend_row_order,
          cluster_columns = dend_col_order,
          right_annotation = rowAnnotation(
            # bar1 = anno_points(log1p(meta.amun$`Distance to coast [km]`), size = unit(3, "mm"), width = unit(2.6, "cm"),
            #                    axis_param = list(at = log1p(c(0, 25, 250)), labels = c(0, 25, 250), labels_rot = 0, gp = gpar(fontsize = 10))),
            # bar2 = anno_points(log10(meta.amun$`Depth [m]`), size = unit(3, "mm"), width = unit(2.6, "cm"),
            #                    axis_param = list(at = log10(c(1, 10, 75)), labels = c(1, 10, 75), labels_rot = 0, gp = gpar(fontsize = 10))),
            bar1 = anno_points(meta.amun$`Conservative temperature [deg C]`, size = unit(3, "mm"), width = unit(2.6, "cm"),
                               axis_param = list(labels_rot = 0, gp = gpar(fontsize = 10))),
            bar2 = anno_points(meta.amun$`Salinity`, size = unit(3, "mm"), width = unit(2.6, "cm"),
                               axis_param = list(labels_rot = 0, gp = gpar(fontsize = 10))),
            bar3 = anno_points(meta.amun$`fChl-a [mg/m^3]`, size = unit(3, "mm"), width = unit(2.6, "cm"),
                               axis_param = list(labels_rot = 0, gp = gpar(fontsize = 10))),
            bar4 = anno_points(meta.amun$`fCDOM [mg/m^3]`, size = unit(3, "mm"), width = unit(2.6, "cm"),
                               axis_param = list(labels_rot = 0, gp = gpar(fontsize = 10))),
            bar5 = anno_points(meta.amun$`Nitrate [uM]`, size = unit(3, "mm"), width = unit(2.6, "cm"),
                               axis_param = list(labels_rot = 0, gp = gpar(fontsize = 10))),
            bar6 = anno_points(meta.amun$`Nitrite [uM]`, size = unit(3, "mm"), width = unit(2.6, "cm"),
                               axis_param = list(labels_rot = 0, gp = gpar(fontsize = 10))),
            gap = unit(3, "mm"), show_annotation_name = FALSE),
          
          row_labels = paste0("  ", rownames(proptab.trans)),  # move rownames a bit to the right
          
          top_annotation = {
            asv_ids <- colnames(proptab.trans); all_indic <- intersect(asv_ids, specific_asvs); focal_indic <- intersect(asv_ids, asv_in_dotplot)
            
            top_vec <- setNames(rep(NA_character_, length(asv_ids)), asv_ids)
            top_vec[all_indic] <- as.character(chosen_cluster[all_indic])
            
            bottom_vec <- setNames(rep(NA_character_, length(asv_ids)), asv_ids)
            bottom_vec[focal_indic] <- as.character(chosen_cluster[focal_indic])
            
            top <- factor(top_vec, levels = as.character(seq_len(k)))
            bottom <- factor(bottom_vec, levels = as.character(seq_len(k)))
            
            columnAnnotation("All bioindicators" = top, "Prevalent bioindicators" = bottom,
                             col = list("All bioindicators" = cluster_palette, "Prevalent bioindicators" = cluster_palette),
                             show_annotation_name = TRUE, annotation_name_side = "left",
                             annotation_name_gp = gpar(fontsize = 16, fontface = "bold"), show_legend = FALSE)
          },
          
          heatmap_height = unit(65, "cm"), heatmap_width = unit(75, "cm"),
          row_dend_width = unit(4, "cm"), column_dend_height = unit(2, "cm"),
          row_title = "Samples", row_title_gp = gpar(fontsize = 34, fontface = "bold"),
          column_title = "ASVs", column_title_gp = gpar(fontsize = 34, fontface = "bold"),
          row_names_gp = gpar(fontsize = 13, fontface = "bold"),
          row_names_max_width = unit(7.6, "cm"),
          column_names_gp = gpar(fontsize = 1, fontface = "bold"), column_names_rot = 80,
          
          show_heatmap_legend = FALSE,
          heatmap_legend_param = list(
            title = NULL,
            # title = "Fourth-root transformed relative abundances\n",
            title_gp = gpar(fontsize = 30, fontface = "bold"),
            legend_direction = "horizontal", legend_width = unit(8, "cm"),
            color_bar = "continuous")) %>%
    ComplexHeatmap::draw(., heatmap_legend_side = "top", padding = unit(c(10, 10, 10, 10), "mm"))  # top, right, bottom, left
  
  decorate_heatmap_body("amun_complex_hm", {
    grid.segments(x0 = unit(0, "npc"), x1 = unit(1, "npc"), y0 = unit(row_sep_y, "npc"), y1 = unit(row_sep_y, "npc"),
                  gp = gpar(col = "white", lwd = 4))
  })
  
  # Decorate slices (colors first, then labels on top)
  # for (i in 1:k) {
  #   decorate_row_names("amun_complex_hm", {
  #     grid.rect(gp = gpar(col = cluster_palette[i], lwd = 3.2, fill = NA),
  #               width = unit(7.6, "cm"), x = unit(1, "npc") + unit(-0.4, "cm"), just = "center")
  #   }, slice = i)
  # }
  # for (i in 1:k) {
  #   decorate_row_names("amun_complex_hm", {
  #     grid.text(cl_labels[i], x = unit(2, "npc") + unit(-4, "cm"), just = "left",
  #               gp = gpar(fontsize = 27, fontface = "bold"))
  #   }, slice = i)
  # }
  
  annotation_labels <- list(
    list(name = "Θ", unit = expression(bold("(" * degree * "C)"))),
    list(name = "S", unit = NULL),
    list(name = "fChl-a", unit = expression(bold("(mg m"^{-3} * ")"))),
    list(name = "fCDOM", unit = expression(bold("(mg m"^{-3} * ")"))),
    list(name = "Nitrate", unit = expression(bold("(µM)"))),
    list(name = "Nitrite", unit = expression(bold("(µM)"))))
  
  for (i in seq_along(annotation_labels)) {
    decorate_annotation(paste0("bar", i), {
      if (is.null(annotation_labels[[i]]$unit)) {
        grid.text(annotation_labels[[i]]$name,
                  y = unit(1, "npc") + unit(5.5, "mm"), just = "bottom",
                  gp = gpar(fontsize = 18, fontface = "bold"))
      } else {
        grid.text(annotation_labels[[i]]$name,
                  y = unit(1, "npc") + unit(11.0, "mm"), just = "bottom",
                  gp = gpar(fontsize = 18, fontface = "bold"))
        
        grid.text(annotation_labels[[i]]$unit,
                  y = unit(1, "npc") + unit(1.5, "mm"), just = "bottom",
                  gp = gpar(fontsize = 18, fontface = "bold"))
      }
    })
  }
}

cairo_pdf(file = "OUT/3A_CLUS/amun_complex_hm.pdf", width = 32, height = 27)
draw_amun_complex_hm()
dev.off()

png("OUT/3A_CLUS/amun_complex_hm.png", width = 32, height = 27, units = "in", res = 1000, type = "cairo")
draw_amun_complex_hm()                       # calls the function which includes all decorations
dev.off()


################################
### SAMPLE DENDROGRAM EXPORT ###
# Uses the fitted sample tree; creates a separate dendrogram for inspection and figure assembly.

# Create a completely empty heatmap structure
empty_mat <- matrix(0, nrow = nrow(proptab.trans), ncol = 0)
rownames(empty_mat) <- rownames(proptab.trans)

# Create the heatmap w/ a name
ht <- Heatmap(empty_mat,
              name = "dendro_plot",
              cluster_rows = dend_row_order,
              row_labels = rownames(proptab.trans),
              row_names_side = "right",
              row_dend_side = "left",
              show_row_names = TRUE,
              show_column_names = FALSE,
              show_heatmap_legend = FALSE,
              row_title = NULL,
              column_title = NULL,
              heatmap_height = unit(72, "cm"),
              row_dend_width = unit(10, "cm"),
              row_names_gp = gpar(fontsize = 20, fontface = "bold"))

cairo_pdf(file = "OUT/3A_CLUS/amun_dendro_pure.pdf", width = 12, height = 30)
ComplexHeatmap::draw(ht)
dev.off()

png("OUT/3A_CLUS/amun_dendro_pure.png", width = 12, height = 30, units = "in", res = 300)
ComplexHeatmap::draw(ht)
dev.off()
