# 4A_PLANKT.R ------------------------------------------------------------------------------------
#
# Purpose:
# Calculate major group composition from the original eukaryotic and prokaryotic ASV tables and assemble Figure 4.
# Use cached rarefied Shannon indices to summarize and compare diversity among clusters.
# Save taxonomy, group abundances and statistics for subsequent figures and manuscript exploration.
#
# Inputs (paths are relative to the project working directory containing EXT, OUT and DNA):
# EXT/meta.amun.rds; EXT/proptab.amun.rds; EXT/proptab.amun.original.rds (2A)
# EXT/amun_taxonomy.rds (2A); EXT/amun_clustering.rds; EXT/amun_bioindicators.rds (3A)
# EXT/amun_shannon_diversity.rds (cached indices); EXT/amun_original_counts.rds (input check and recalculation)
# Optional EXT/amun_shannon_depth_assessment.rds records the chosen depths and stability tolerance
#
# Outputs:
# OUT/4A_PLANKT/plankton_communities.pdf / .png (Figure 4)
# EXT/amun_composition.rds -> original-ASV taxonomy and abundance summaries for 5A and 6A
# EXT/amun_shannon_statistics.rds -> diversity summaries and tests
# Optional original-count recalculation updates EXT/amun_shannon_diversity.rds and amun_shannon_depth_assessment.rds
#
# Run sections in order within this script. Later scripts load saved inputs; they do not need earlier scripts in memory.
# Heavy calculations remain optional or cached. Existing EXT and OUT folders are assumed to exist.
#---------------------------------------------------------------------------------------------------

library(tidyverse)
library(vegan)
library(patchwork)


#############################################
### LOAD COMPOSITION AND DIVERSITY INPUTS ###
# Read the saved inputs so 4A can run in a fresh R session.
# Uses: 2A taxonomy and original RA; 3A clustering and bioindicators; cached Shannon values.
# Produces: The objects used in this script.

meta.amun <- readRDS("EXT/meta.amun.rds")
proptab.amun.original <- readRDS("EXT/proptab.amun.original.rds")
proptab.amun <- readRDS("EXT/proptab.amun.rds")
clustering <- readRDS("EXT/amun_clustering.rds")
clus <- clustering$clus; cl_labels <- clustering$cl_labels; k <- clustering$k
cluster_palette <- clustering$cluster_palette; dend_row_order <- clustering$dend_row_order; row_ordered <- clustering$row_ordered
meta.amun$Cluster <- clus[rownames(meta.amun)]
taxonomy <- readRDS("EXT/amun_taxonomy.rds")
asv16.original <- taxonomy$asv16.original; asv18.original <- taxonomy$asv18.original
bioindicators <- readRDS("EXT/amun_bioindicators.rds")
asv_in_dotplot <- bioindicators$asv_in_dotplot

# Keep labels, axis numbers and tick marks black; individual panels still control font sizes and hide unused axes.
black_theme <- theme(text = element_text(colour = "black"), axis.text = element_text(colour = "black"),
                     axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
                     axis.title = element_text(colour = "black"), axis.ticks = element_line(colour = "black"),
                     legend.text = element_text(colour = "black"), legend.title = element_text(colour = "black"),
                     strip.text = element_text(colour = "black"), plot.title = element_text(colour = "black"),
                     plot.subtitle = element_text(colour = "black"), plot.caption = element_text(colour = "black"))

community_theme <- theme_bw() + black_theme +
  theme(axis.text.y = element_text(face = "bold", color = "black", size = 16),
        axis.text.x = element_text(color = "black", size = 15),
        axis.title = element_text(color = "black", size = 16, face = "bold"),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(color = "black", size = 16),
        legend.title = element_blank())

####################################################################################################
### SHANNON INDEX by cluster ###
# Run steps 1–4 to recompute; otherwise load the cached values using step 5.

# Original integer counts are also used to check that cached Shannon values match the current data.
shannon_counts <- readRDS("EXT/amun_original_counts.rds")
shannon_counts <- lapply(shannon_counts, \(tab) tab[rownames(meta.amun), , drop = FALSE])

# # 1. Define the read totals and depth-selection tolerance.
# shannon_reads <- list(Eukaryota = meta.amun$`Eukaryotic Reads`, Prokaryota = meta.amun$`Prokaryotic Reads`)
# shannon_tolerance <- 0.05                                              # maximum acceptable change in mean Shannon
# 
# # 2. Calculate expected Shannon under rarefaction without replacement
# # This gives the mean over possible subsamples directly, avoiding simulations at every candidate depth.
# expected_shannon <- function(x, depth) {
#   total <- sum(x); frequencies <- table(x[x > 0])
#   sum(as.numeric(frequencies) * vapply(as.integer(names(frequencies)), \(n) {
#     k <- seq.int(max(1L, depth - total + n), min(n, depth)); p <- k / depth
#     -sum(p * log(p) * stats::dhyper(k, n, total - n, depth))
#   }, numeric(1)))
# }
# 
# # Evaluate observed read counts using the same better-sequenced half of samples at every candidate depth
# assess_shannon <- function(tab, reads) {
#   evaluation <- reads >= median(reads); reference_depth <- min(reads[evaluation])
#   tab <- tab[evaluation, , drop = FALSE]; depths <- sort(unique(reads[reads <= reference_depth]))
#   reference <- apply(tab, 1, expected_shannon, depth = reference_depth)
#   change <- vapply(depths, \(d) max(abs(apply(tab, 1, expected_shannon, depth = d) - reference)), numeric(1))
#   chosen <- which(change <= shannon_tolerance)[1]                       # smallest depth meeting the tolerance in every evaluation sample
#   tibble(Minimum_reads = depths[chosen], Max_Shannon_change = change[chosen])
# }
# 
# # Select cutoffs independently, then retain samples meeting BOTH cutoffs
# depth_choice <- map2(shannon_counts, shannon_reads, assess_shannon) %>% bind_rows(.id = "Domain")
# cutoffs <- setNames(depth_choice$Minimum_reads, depth_choice$Domain)
# keep <- shannon_reads$Eukaryota >= cutoffs["Eukaryota"] & shannon_reads$Prokaryota >= cutoffs["Prokaryota"]
# final_depths <- vapply(shannon_reads, \(x) min(x[keep]), numeric(1))      # marker-specific minima among retained samples
# depth_choice <- depth_choice %>%
#   mutate(Final_rarefaction_depth = unname(final_depths[Domain]), Paired_samples_retained = sum(keep))
# print(depth_choice, width = Inf)
# 
# # Compare reference-depth and selected-depth Shannon with Shannon from all available reads
# full_depth_check <- imap_dfr(shannon_counts, \(tab, domain) {
#   reads <- shannon_reads[[domain]]; evaluation <- reads >= median(reads)
#   reference_depth <- min(reads[evaluation]); tab <- tab[evaluation, , drop = FALSE]
#   full <- vegan::diversity(tab, index = "shannon")
#   tibble(Domain = domain,
#          Reference_to_full = max(abs(full - apply(tab, 1, expected_shannon, depth = reference_depth))),
#          Cutoff_to_full = max(abs(full - apply(tab, 1, expected_shannon, depth = cutoffs[domain]))))
# })
# print(full_depth_check, width = Inf)
# 
# # 3. Average Shannon across 1,000 random rarefactions at the final depths
# mean_shannon <- function(tab, depth, reps = 1000L) {
#   rowMeans(replicate(reps, vegan::diversity(vegan::rrarefy(tab, sample = depth), index = "shannon")))
# }
# set.seed(161193)
# shannon_tbl <- tibble(Sample = rownames(meta.amun)[keep],
#                       Shannon_Euk = mean_shannon(shannon_counts$Eukaryota[keep, , drop = FALSE], final_depths["Eukaryota"]),
#                       Shannon_Prok = mean_shannon(shannon_counts$Prokaryota[keep, , drop = FALSE], final_depths["Prokaryota"]))
# 
# # 4. Cache Shannon values and the sequencing-depth summary
# attr(shannon_tbl, "input_hash") <- rlang::hash(shannon_counts)
# saveRDS(shannon_tbl, "EXT/amun_shannon_diversity.rds")
# saveRDS(list(Selected_depths = depth_choice, Tolerance = shannon_tolerance), "EXT/amun_shannon_depth_assessment.rds")

# 5. Attach cached Shannon values to metadata; the cache contains 107 paired samples at 1,420 / 7,160 reads.
shannon_tbl <- readRDS("EXT/amun_shannon_diversity.rds")
if (!identical(attr(shannon_tbl, "input_hash"), rlang::hash(shannon_counts)))
  stop("Shannon counts changed or the cache has no input fingerprint; regenerate it using steps 1-4 above.")
meta.amun <- meta.amun %>% select(-any_of(c("Shannon_Euk", "Shannon_Prok"))) %>% rownames_to_column("Sample") %>%
  left_join(shannon_tbl, by = "Sample") %>% column_to_rownames("Sample")


### ### Plot Shannon diversity and test differences among clusters
domain_order <- c("Eukaryota", "Prokaryota")
shannon_long <- meta.amun %>% rownames_to_column("Sample") %>%
  pivot_longer(c(Shannon_Euk, Shannon_Prok), names_to = "Domain", values_to = "Shannon", values_drop_na = TRUE) %>%
  mutate(Domain = recode_values(Domain, "Shannon_Euk" ~ "Eukaryota", "Shannon_Prok" ~ "Prokaryota"),
         Domain = factor(Domain, levels = domain_order), Cluster = factor(cl_labels[as.integer(as.character(Cluster))], levels = rev(cl_labels), ordered = TRUE))
by_dom <- split(shannon_long, shannon_long$Domain)

# Plot Shannon diversity for eukaryotes and prokaryotes
shannon_plot <- ggplot(shannon_long, aes(Shannon, Cluster, color = Domain)) +
  geom_boxplot(fill = "white", width = 0.6, coef = Inf, outlier.shape = NA, position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8, seed = 161193), size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c(Eukaryota = "#D55E00", Prokaryota = "#0072B2")) +
  scale_x_continuous(expand = expansion(mult = 0.02)) + scale_y_discrete(drop = FALSE) +
  labs(x = "Shannon Diversity Index", y = NULL) +
  community_theme + theme(plot.margin = margin(l = 10, r = 20))

# Inspect residual normality and group variances before interpreting ANOVA and Tukey results
fits <- lapply(by_dom, \(d) aov(Shannon ~ Cluster, data = d))
imap_dfr(fits, \(m, dom) tibble(Domain = dom, `Shapiro p-value` = shapiro.test(residuals(m))$p.value)) %>%
  mutate(Domain = factor(Domain, levels = domain_order)) %>% arrange(Domain) %>% print()
imap_dfr(by_dom, \(d, dom) tibble(Domain = dom, `Levene p-value` = car::leveneTest(Shannon ~ Cluster, data = d, center = median)$`Pr(>F)`[1])) %>%
  mutate(Domain = factor(Domain, levels = domain_order)) %>% arrange(Domain) %>% print()

# Medians describe the plotted distributions; means and SD accompany the ANOVA and Tukey comparisons.
shannon_summary <- shannon_long %>% group_by(Domain, Cluster) %>%
  summarise(n = n(), Median = median(Shannon), Mean = mean(Shannon), SD = sd(Shannon), .groups = "drop") %>% arrange(Domain, desc(Cluster))

# Overall differences in mean Shannon diversity among clusters
shannon_anova <- imap_dfr(fits, \(m, dom) {
  a <- summary(m)[[1]]
  tibble(Domain = dom, Cluster_df = a["Cluster", "Df"], Residual_df = a["Residuals", "Df"],
         F = a["Cluster", "F value"], `p-value` = a["Cluster", "Pr(>F)"])
}) %>% mutate(Domain = factor(Domain, levels = domain_order)) %>% arrange(Domain)

# Significant Tukey HSD pairwise comparisons for domains with a significant overall ANOVA
# No BH adjustment b/c unlike Dunn and Spearman, Tukey already supplies its own multiple-comparison adjustment
shannon_tukey <- imap_dfr(fits, \(m, dom) {
  tk <- TukeyHSD(m, "Cluster")$Cluster %>% as.data.frame() %>% rownames_to_column("Comparison")
  groups <- str_split_fixed(tk$Comparison, "-", 2)
  tibble(Domain = dom,
         `Direction of difference` = if_else(tk$diff > 0, paste(groups[, 1], ">", groups[, 2]), paste(groups[, 2], ">", groups[, 1])),
         `Mean difference` = abs(tk$diff), `Tukey-adjusted p-value` = tk$`p adj`)
}) %>%
  filter(Domain %in% as.character(shannon_anova$Domain[shannon_anova$`p-value` < 0.01]),
         `Tukey-adjusted p-value` < 0.01) %>%
  mutate(Domain = factor(Domain, levels = domain_order)) %>% arrange(Domain, `Tukey-adjusted p-value`)


##########################################
### COMMUNITY COMPOSITION AND FIGURE 4 ###

### ### Reduce the original sample dendrogram to one tip per final cluster
cluster_reps <- row_ordered[!duplicated(as.character(clus[row_ordered]))]
dendro_clusters <- dendextend::prune(dend_row_order, setdiff(row_ordered, cluster_reps))
dendro_seg <- dendextend::as.ggdend(dendro_clusters)$segments

dendro_plot <- ggplot(dendro_seg) +
  geom_segment(aes(x = y, y = (k + 1) - x, xend = yend, yend = (k + 1) - xend),
               colour = "grey20", linewidth = 1, lineend = "round") +
  scale_y_continuous(expand = expansion(add = 0.6)) +
  scale_x_reverse(expand = expansion(mult = c(0.05, 0.02))) +
  theme_void() +
  theme(plot.margin = margin(l = 10, r = 5, t = 5, b = 5))


### ### Prepare major groups for Figure 4 and detailed taxonomy for further exploration

# Read the taxonomic names attached to each ASV so related organisms can be combined into major groups below.
parse_tax <- function(asv) {
  tx <- asv %>% str_remove(";\\s*ASV_\\d+$") %>% stringr::str_remove_all("\\s*\\[\\d+\\]"); parts <- strsplit(tx, ";\\s*")[[1]]
  tibble(Taxonomy = tx, Domain = if (length(parts) >= 1) parts[1] else NA_character_, Phylum = if (length(parts) >= 2) parts[2] else NA_character_,
         Order = if (length(parts) >= 3) parts[3] else NA_character_, Genus = if (length(parts) >= 4) parts[4] else NA_character_)
}

# Collect class and family names from the original 16S and 18S taxonomy tables.
# The class names help determine which ASVs belong to groups such as diatoms and non-diatom Ochrophyta.
class_family <- bind_rows(asv16.original %>% transmute(ASV_id = ASV, Domain, Class, Family), asv18.original %>% transmute(ASV_id = ASV, Domain, Class, Family))

# Keep the detailed taxonomy of each ASV for identifying genera and other taxa in the optional plots and summaries at the end of 6A_ENV.R.
genus_taxonomy <- tibble(ASV = colnames(proptab.amun.original), ASV_id = str_extract(ASV, "ASV_\\d+$"), Domain = str_extract(ASV, "^[^;]+")) %>%
  left_join(bind_rows(asv16.original %>% transmute(ASV_id = ASV, Domain_group = "Prokaryote", Domain, Supergroup = NA_character_,
                                                   Phylum_or_Division = Phylum, Subdivision = NA_character_, Class, Order, Family, Genus),
                      asv18.original %>% transmute(ASV_id = ASV, Domain_group = "Eukaryote", Domain, Supergroup,
                                                   Phylum_or_Division = Division, Subdivision, Class, Order, Family, Genus)),
            by = c("ASV_id", "Domain")) %>% select(-ASV_id)

# Assign each ASV to a major taxonomic group using the rules below.
# These assignments determine which organisms are combined in Figure 4 and in major group summaries at the end of 6A_ENV.R.
# Ecological notes below refer to the earlier non-rare summary; the optional group audit at the end of 6A recalculates them for original ASVs.
all_asvs_tax <- tibble(ASV_full = colnames(proptab.amun.original)) %>% rownames_to_column() %>% rowwise() %>%
  mutate(bind_rows(parse_tax(ASV_full))) %>%
  ungroup() %>% mutate(ASV_id = str_extract(ASV_full, "ASV_\\d+$")) %>%
  left_join(class_family, by = c("ASV_id", "Domain")) %>%                               # join Class and Family
  transmute(ASV = ASV_full, Domain2 = if_else(Domain %in% c("Bacteria", "Archaea"), "Prokaryote", "Eukaryote"),
            Group = if_else(Domain == "Archaea", "Archaea",
                            case_when(
                              # in these samples, Coscinodiscophyceae were completely absent from diatoms
                              # 98 % of diatom abundance is Chaetoceros but Minidiscus and Thalassiosira also attain high relab of 10-25% in some samples
                              Class %in% c("Bacillariophyceae", "Mediophyceae", "Coscinodiscophyceae") ~ "Diatoms",
                              
                              # non-diatom Ochrophyta = selected non-diatom Gyrista classes; these include heterotrophs:
                              # in these samples, only Pelagophyceae, Chrysophyceae, Dictyochophyceae, Phaeophyceae, Bolidophyceae, and MOCH-2 were present,
                              # but other potential members are also listed for consistency if other studies would reuse the code
                              # Pedinellales, the only order detected in Dictyochophyceae, is 71 % of all "non-diatom Ochrophyta"
                              # However, Chrysophyceae (golden algae) and Phaeophyceae (brown algae) also appear on the dotplot
                              Class %in% c("Pelagophyceae", "Chrysophyceae", "Dictyochophyceae", "Phaeophyceae","Bolidophyceae", "MOCH-2", 
                                           "Eustigmatophyceae", "Raphidophyceae", "Xanthophyceae") ~ "non-diatom Ochrophyta",
                              
                              Class == "Dinophyceae" ~ "Dinophyceae",                   # includes photosynthetic, heterotrophic and parasitic genera
                              Class == "Syndiniales" ~ "Syndiniales",                   # parasitic dinoflagellates
                              Order == "Eugregarinorida" ~ "Eugregarinorida",           # 99.7% of Gregarinomorphea
                              
                              Class == "Alphaproteobacteria" ~ "Alphaproteobacteria",   # split Proteobacteria into one class + 2 orders
                              Order == "Pseudomonadales" ~ "Pseudomonadales", Order == "Enterobacterales" ~ "Enterobacterales",
                              Order == "Flavobacteriales" ~ "Flavobacteriales",         # 99.3 % of Bacteroidota and often mentioned in lit-re as such
                              TRUE ~ str_replace(Phylum, "^uncl.\\s+", ""))))

# Combine each sample's ASV abundances with its cluster membership and taxonomic assignments.
# This table supplies the major group bars in Figure 4, the optional taxon plots and abundance summaries at the end of 6A_ENV.R.
# Keep zero abundances so averages include samples in which a taxon was absent.
sample_asvs <- meta.amun %>% rownames_to_column("Sample") %>% select(Sample, Cluster) %>%
  inner_join(proptab.amun.original %>% as_tibble(rownames = "Sample") %>% pivot_longer(-Sample, names_to = "ASV", values_to = "RelAb"), by = "Sample") %>%
  left_join(all_asvs_tax, by = "ASV") %>%
  left_join(genus_taxonomy %>% select(ASV, Domain_group, Phylum_or_Division, Class, Order, Family, Genus), by = "ASV")

# Add together ASVs belonging to the same major group within each sample, then average those totals within each cluster.
# These mean percentages are used for the bars in Figure 4 and for inspecting group contributions below.
cluster_group <- sample_asvs %>% group_by(Sample, Cluster, Domain2, Group) %>% summarise(GroupAbundSample = sum(RelAb), .groups = "drop") %>%
  group_by(Cluster, Domain2, Group) %>% summarise(AvgAbund = 100 * mean(GroupAbundSample), .groups = "drop") %>%
  mutate(Cluster = factor(cl_labels[as.integer(as.character(Cluster))], levels = rev(cl_labels), ordered = TRUE))

# Inspect how much rare ASVs add to each group's mean RA within each cluster.
# Rare_ASVs_pp is the contribution in percentage points; zeros remain included in every mean.
# View(sample_asvs %>% group_by(Sample, Cluster, Domain2, Group) %>%
#        summarise(Original_pct = 100 * sum(RelAb), Nonrare_pct = 100 * sum(RelAb[ASV %in% colnames(proptab.amun)]), .groups = "drop") %>%
#        group_by(Cluster, Domain2, Group) %>%
#        summarise(across(c(Original_pct, Nonrare_pct), mean), .groups = "drop") %>%
#        mutate(Rare_ASVs_pp = Original_pct - Nonrare_pct, across(c(Original_pct, Nonrare_pct, Rare_ASVs_pp), ~ round(.x, 3)),
#               Cluster = factor(cl_labels[as.integer(as.character(Cluster))], levels = cl_labels)) %>% arrange(desc(Rare_ASVs_pp)))

# Inspect which classes, orders and genera were combined into each major group.
# Use this to review the grouping rules before making the plots.
# View(sample_asvs %>% distinct(Domain2, Group, Phylum_or_Division, Class, Order, Genus) %>% arrange(Domain2, Group, Phylum_or_Division, Class, Order, Genus))

# Each plotted group reaches >2% mean relative abundance in at least one cluster.
# Pseudomonadales and Enterobacterales are shown separately; residual Pseudomonadota is omitted from the bars.

# Choose which major groups appear in Figure 4 and assign their colours.
# Groups absent from these lists are omitted from the bars.
euk_colors <- c("Diatoms" = "#478214", "non-diatom Ochrophyta" = "#D6C200", "Chlorophyta" = "#00BC00",                 # Includes phototrophs and heterotrophs
                "Dinophyceae" = "#1E5AFF", "Ciliophora" = "#B07C96", "Cercozoa" = "#00CDDA", "Radiolaria" = "#A000FF", # Mixotrophs and protist consumers
                "Eugregarinorida" = "#FF1A2E", "Syndiniales" = "#FE8A1F",                                              # Parasites
                "Metazoa" = "#F978BE")                                                                                 # Other

prok_colors <- c("Cyanobacteriota" = "#00BC00",                                                                     # Primary producers
                 "Alphaproteobacteria" = "#FF1A2E", "Pseudomonadales" = "#FE8A1F", "Enterobacterales" = "#D6C200",  # Pseudomonadota/Proteobacteria (majority)
                 "Flavobacteriales" = "#F978BE",                                                                    # Complex OM degraders
                 "Planctomycetota" = "#1E5AFF", "Verrucomicrobiota" = "#00CDDA",                                    # PVC sisters
                 "Actinomycetota" = "#A000FF", "Archaea" = "#7B7F55")                                               # Other

# Inspect each group's mean percentage in every cluster and whether it will appear in Figure 4.
# Use these tables to check the >2% criterion and review which groups are omitted.
# View(cluster_group %>%
#        mutate(Plotted = (Domain2 == "Eukaryote" & Group %in% names(euk_colors)) | (Domain2 == "Prokaryote" & Group %in% names(prok_colors))) %>%
#        transmute(Domain2, Group, Cluster, Mean_pct = round(AvgAbund, 2), Plotted) %>% arrange(Domain2, Group, Cluster))
# View(cluster_group %>%
#        mutate(Plotted = (Domain2 == "Eukaryote" & Group %in% names(euk_colors)) | (Domain2 == "Prokaryote" & Group %in% names(prok_colors))) %>%
#        transmute(Domain2, Group, Cluster, Mean_pct = round(AvgAbund, 2), Plotted) %>%
#        slice_max(Mean_pct, n = 1, by = Group) %>%
#        arrange(desc(Mean_pct)))

### ### Make a bar plots for major groups (once you have decided what to highlight)

# Formatting shared by all panels
# Filter data to include ONLY the chosen groups in the stacked bar plot
major_group_plot <- function(data, domain, colors, xlab, left, right, legend_space) {
  data %>% filter(Domain2 == domain, Group %in% names(colors)) %>%
    mutate(Group = factor(Group, levels = rev(names(colors)), ordered = TRUE)) %>%
    ggplot(aes(AvgAbund, Cluster, fill = Group)) +
    geom_col(orientation = "y", position = position_stack(reverse = FALSE),
             width = 0.7, color = "black", linewidth = 0.1) +
    scale_fill_manual(values = colors, breaks = names(colors), drop = FALSE) + scale_y_discrete(drop = FALSE) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
    labs(x = xlab, y = NULL) +
    community_theme +
    theme(plot.margin = margin(l = left, r = right), legend.key.spacing.x = grid::unit(legend_space, "cm")) +
    guides(fill = guide_legend(ncol = 2, byrow = FALSE))
}

# Plot major groups for eukaryotes and prokaryotes
euk_plot <- major_group_plot(cluster_group, "Eukaryote", euk_colors, "Eukaryotic Relative Abundance (%)", 15, 15, 2)
prok_plot <- major_group_plot(cluster_group, "Prokaryote", prok_colors, "Prokaryotic Relative Abundance (%)", 20, 10, 1.5)


### ### Combine plots and output them
final_plot <- (dendro_plot + shannon_plot + euk_plot + prok_plot) + plot_layout(nrow = 1, widths = c(0.2, 0.8, 1, 1))
cairo_pdf("OUT/4A_PLANKT/plankton_communities.pdf", width = 21, height = 9.5, onefile = FALSE); print(final_plot); dev.off()
png("OUT/4A_PLANKT/plankton_communities.png", width = 21, height = 9.5, units = "in", res = 300); print(final_plot); dev.off()


############################################
### SAVE COMPOSITION AND SHANNON RESULTS ###
# Provide the exact group assignments, colours and abundances used in Figure 4 to 5A and 6A.
# Uses: Original ASVs, full taxonomy, the colour lists and Shannon results above.
# Produces: EXT/amun_composition.rds; EXT/amun_shannon_statistics.rds. No taxonomic summaries use the rare-ASV filter.

# These marker counts support the Results 3.3 introduction and Figure 5.
bioindicator_counts <- sample_asvs %>% filter(ASV %in% asv_in_dotplot) %>% distinct(Domain_group, ASV) %>% count(Domain_group)
print(bioindicator_counts)
# Eukaryote: 36 prevalent bioindicators; Prokaryote: 59 prevalent bioindicators.

saveRDS(list(sample_asvs = sample_asvs, all_asvs_tax = all_asvs_tax, genus_taxonomy = genus_taxonomy,
             cluster_group = cluster_group, euk_colors = euk_colors, prok_colors = prok_colors), "EXT/amun_composition.rds")
saveRDS(list(summary = shannon_summary, anova = shannon_anova, tukey = shannon_tukey), "EXT/amun_shannon_statistics.rds")

# Across all clusters, the printed ranges show total RA represented by the prevalent bioindicators (Results 3.3).
# These are sums of ASV percentages within samples, followed by cluster means; zeros remain included.
bioindicator_coverage <- sample_asvs %>% group_by(Domain_group, Cluster) %>%
  summarise(Pct = 100 * sum(RelAb[ASV %in% asv_in_dotplot]) / n_distinct(Sample), .groups = "drop_last") %>%
  summarise(Min = min(Pct), Max = max(Pct), .groups = "drop")
print(bioindicator_coverage)
# Domain_group   Min   Max
# Eukaryote     63.0  81.6
# Prokaryote    60.0  90.2
