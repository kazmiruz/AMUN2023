# 7A_dbRDA.R -------------------------------------------------------------------------------------
#
# Purpose:
# Assess environmental associations with eukaryotic and prokaryotic composition in CL3-CL6.
# Assess additional variation in prokaryotic composition associated with prevalent eukaryotic ASVs after environmental adjustment.
# Draw Figure 6: environmental dbRDAs (A-B) and partitioning between selected environmental and eukaryotic ASV predictors (C).
#
# Inputs (paths are relative to the project working directory containing EXT, OUT and DNA):
# EXT/proptab.amun.rds (2A; non-rare combined ASV table)
# EXT/amun_environment_setup.rds (6A; metadata, transformations, colours and sampling-phase symbols)
#
# Outputs:
# OUT/7A_dbRDA/amun_CL3-6_dbRDA.pdf / .png / .svg (Figure 6, panels A-C)
# OUT/7A_dbRDA/amun_CL3-6_dbRDA_statistics.xlsx (model tests, explained fractions, predictor contributions and candidate ASVs)
#
# Variation partitioning: Borcard et al. (1992); adjusted R2: Peres-Neto et al. (2006).
# Panels B and C use exactly the same selected environmental predictors; ASV selection conditions on them from the start.
# ASV eligibility depends only on eukaryotic abundance; bioindicator status is not used in this analysis.
#
# Run sections in order; models and tests are recalculated each time. Earlier scripts need not be in memory.
# Existing EXT and OUT/7A_dbRDA folders are assumed to exist.
#---------------------------------------------------------------------------------------------------

library(tidyverse)
library(vegan)
library(permute)
library(patchwork)
library(ggrepel)
library(openxlsx)

### ### Save sortable Excel tables
# Stored values remain unrounded. Display counts as integers, other values to one decimal, and small values scientifically.
save_tables <- function(tables, file) {
  wb <- createWorkbook()
  for (sheet in names(tables)) {
    d <- as.data.frame(tables[[sheet]])
    addWorksheet(wb, sheet)
    writeDataTable(wb, sheet, d, tableStyle = "TableStyleMedium2", keepNA = TRUE, na.string = "NA")
    freezePane(wb, sheet, firstRow = TRUE, firstCol = TRUE)
    setColWidths(wb, sheet, cols = seq_along(d), widths = "auto")
    addStyle(wb, sheet, createStyle(wrapText = TRUE), rows = 1, cols = seq_along(d), gridExpand = TRUE)
    setRowHeights(wb, sheet, rows = 1, heights = 32)
    for (j in which(vapply(d, is.numeric, logical(1)))) {
      fmt <- if (str_detect(names(d)[j], "^P($| )")) "0.0E+00" else if (is.integer(d[[j]])) "0" else "0.0"
      if (nrow(d)) addStyle(wb, sheet, createStyle(numFmt = fmt), rows = seq_len(nrow(d)) + 1L, cols = j, stack = TRUE)
      small <- which(is.finite(d[[j]]) & d[[j]] != 0 & abs(d[[j]]) < .1)
      if (length(small)) addStyle(wb, sheet, createStyle(numFmt = "0.0E+00"), rows = small + 1L, cols = j, stack = TRUE)
    }
  }
  saveWorkbook(wb, file, overwrite = TRUE)
}

### ### Shared plot appearance
black_theme <- theme(text = element_text(colour = "black", face = "bold", family = "sans"), axis.text = element_text(colour = "black"),
                     axis.title = element_text(colour = "black"), axis.ticks = element_line(colour = "black"),
                     legend.text = element_text(colour = "black"), legend.title = element_text(colour = "black"),
                     strip.text = element_text(colour = "black"), plot.title = element_text(colour = "black"))
cluster_key <- \(data, params, size) {
  if (is.na(data$shape)) grid::rectGrob(width = .9, height = .55, gp = grid::gpar(fill = data$fill, col = NA))
  else ggplot2::draw_key_point(data, params, size)
}

################################################
### LOAD INPUTS AND SELECT THE FOCAL SAMPLES ###
# Use the same CL3-CL6 samples, environmental transformations and symbols as 6A.
setup <- readRDS("EXT/amun_environment_setup.rds")
proptab.amun <- readRDS("EXT/proptab.amun.rds")
metadata <- setup$meta_grouped %>% filter(as.character(Cluster) %in% as.character(3:6))
cl_labels <- setup$cl_labels; cluster_palette <- setup$cluster_palette
sampling_periods <- setup$sampling_periods; sampling_shapes <- setup$sampling_shapes; sampling_sizes <- setup$sampling_sizes
plot_sizes <- .35 * sampling_sizes        # Preserve the shape-size compensation used in 6A, at this figure's smaller scale.
seed <- setup$seed; nperm <- 9999; alpha <- .01; ra_power <- .25
out <- "OUT/7A_dbRDA/amun_CL3-6_dbRDA"
clean_env <- \(x) str_remove(x, "\\s*\\[[^]]+\\]$") %>% str_replace_all("[^[:alnum:]]+", "_") %>% str_remove("_$")
env_labels_db <- setNames(str_remove(unname(setup$env_labels), "\\s*\\([^)]*\\)$"), clean_env(names(setup$env_labels)))

### ### Community responses and environmental candidates
# Environmental dbRDAs retain the established Bray-Curtis dissimilarities of fourth-root-transformed RA.
comm <- as.matrix(proptab.amun[metadata$Sample, , drop = FALSE])
euk <- comm[, str_detect(colnames(comm), "^Eukaryota;"), drop = FALSE]
euk <- euk[, colSums(euk) > 0, drop = FALSE]
prok <- comm[, str_detect(colnames(comm), "^(Bacteria|Archaea);"), drop = FALSE]
d_euk <- vegdist(euk^ra_power, "bray"); d_prok <- vegdist(prok^ra_power, "bray")
env <- metadata %>% select(all_of(setup$env_vars)) %>%
  mutate(across(all_of(intersect(setup$log1p_env, setup$env_vars)), log1p)) %>%
  rename_with(clean_env) %>% decostand("standardize") %>% as.data.frame()
rownames(env) <- metadata$Sample

### ### Prevalent eukaryotic ASVs as candidate predictors
# Require at least 2% RA in at least two CL3-CL6 samples. This uses eukaryotic abundance, not cluster associations.
# Fourth-root-transform and standardize the eligible ASVs before conditional forward selection.
candidate_ra <- euk[, colSums(euk >= .02) >= 2, drop = FALSE]
candidate_key <- tibble(ASV_full = colnames(candidate_ra), ASV = str_extract(ASV_full, "ASV_[0-9]+$"),
                        Taxon = str_trim(str_extract(ASV_full, "[^;]+(?=;\\s*ASV_[0-9]+$)")))
colnames(candidate_ra) <- candidate_key$ASV
asv_candidates <- as.data.frame(decostand(candidate_ra^ra_power, "standardize"))

###############################################
### FIT MODELS AND TEST EXPLAINED VARIATION ###
# A-B select environmental predictors. C uses exactly the environmental predictors retained in B.
# Select ASVs for additional explanation after environmental adjustment, then partition using these same two predictor sets.

### ### Sample and whole-station permutations
set.seed(seed); perm <- shuffleSet(nrow(metadata), nset = nperm)
station_rows <- lapply(split(seq_len(nrow(metadata)), interaction(metadata$Cruise, metadata$Station, drop = TRUE)),
                       \(i) i[order(metadata$`Depth [m]`[i])])
set.seed(seed)
station_order <- shuffleSet(length(station_rows), nset = nperm, control = how(blocks = factor(lengths(station_rows))))
perm_station <- t(apply(station_order, 1, \(i) {
  p <- seq_len(nrow(metadata))
  p[unlist(station_rows, use.names = FALSE)] <- unlist(station_rows[i], use.names = FALSE); p
}))

### ### Forward selection for the two environmental models and the partial ASV model
# First test all candidates together; then select at p <0.01, bounded by the full model's adjusted R2.
select_model <- function(d, x, z = NULL) {
  full <- if (is.null(z)) dbrda(d ~ ., data = x) else dbrda(d ~ . + Condition(as.matrix(z)), data = x)
  null <- if (is.null(z)) dbrda(d ~ 1, data = x) else dbrda(d ~ 1 + Condition(as.matrix(z)), data = x)
  full_p <- anova.cca(full, permutations = perm)$`Pr(>F)`[1]
  if (!is.finite(full_p) || full_p >= alpha || RsquareAdj(full)$adj.r.squared <= 0)
    stop("Forward selection stopped: the full candidate model must have p <0.01 and positive adjusted R2.")
  set.seed(seed)
  selected <- ordiR2step(null, scope = formula(full), Pin = alpha - .Machine$double.eps, R2scope = TRUE, permutations = perm, trace = FALSE)
  if (is.null(selected$CCA) || selected$CCA$rank == 0) stop("No explanatory predictors passed forward selection.")
  attr(selected, "full_model_p") <- full_p
  selected
}
m_euk <- select_model(d_euk, env); m_prok <- select_model(d_prok, env)
env_selected <- env[, attr(terms(m_prok), "term.labels"), drop = FALSE]
m_partial <- select_model(d_prok, asv_candidates, env_selected)
selected_asvs <- intersect(attr(terms(m_partial), "term.labels"), colnames(asv_candidates))
asv_selected <- asv_candidates[, selected_asvs, drop = FALSE]
models <- list("Eukaryotic dbRDA" = m_euk, "Prokaryotic dbRDA" = m_prok, "Partial prokaryotic dbRDA" = m_partial)

### ### Model significance and selected predictor contributions
# Final model and marginal p-values are post-selection tests; they do not repeat selection within each permutation.
test_p <- \(m, permutations) anova.cca(m, permutations = permutations, model = "reduced")$`Pr(>F)`[1]
model_summary <- imap(models, \(m, nm) tibble(Model = nm, Samples = nrow(metadata), Stations = length(station_rows),
                                              `Explained variation (%)` = 100 * RsquareAdj(m)$adj.r.squared,
                                              `P (samples)` = test_p(m, perm), `P (whole stations)` = test_p(m, perm_station),
                                              `P (full candidates, samples)` = attr(m, "full_model_p"))) %>% list_rbind()
# Report marginal p-values under sample and whole-station permutations, and VIFs, for the selected predictors in each fitted model.
predictor_summary <- imap(models, \(m, nm) {
  p_station <- anova.cca(m, permutations = perm_station, by = "margin")
  anova.cca(m, permutations = perm, by = "margin") %>% as.data.frame() %>% rownames_to_column("Predictor") %>%
    filter(Predictor != "Residual") %>%
    transmute(Model = nm, Parameter = coalesce(unname(env_labels_db[Predictor]), Predictor),
              `Marginal contribution (%)` = 100 * SumOfSqs / m$tot.chi,
              `P (samples)` = `Pr(>F)`, `P (whole stations)` = p_station[Predictor, "Pr(>F)"], VIF = unname(vif.cca(m)[Predictor]))
}) %>% list_rbind()
selection_steps <- imap(models, \(m, nm) as.data.frame(m$anova) %>% rownames_to_column("Step") %>%
                          transmute(Model = nm, Step, `Explained variation (%)` = 100 * R2.adj, `Entry p-value` = `Pr(>F)`)) %>% list_rbind()

### ### Unique and shared fractions of variation in prokaryotic composition
# Use the environmental predictors from B and the ASVs retained by the partial model; do not reselect either set.
# vegan two-set indfract order: [a] X1|X2, [b] X2|X1, [c] shared, [d] unexplained.
vp <- varpart(d_prok, env_selected, asv_selected)
m_env_unique <- dbrda(d_prok ~ . + Condition(as.matrix(asv_selected)), data = env_selected)
fractions <- tibble(Fraction = c("Environment only", "Eukaryotes only", "Shared", "Unexplained"),
                    `Explained variation (%)` = 100 * vp$part$indfract$Adj.R.square,
                    `P (samples)` = c(test_p(m_env_unique, perm), model_summary$`P (samples)`[3], NA, NA),
                    `P (whole stations)` = c(test_p(m_env_unique, perm_station), model_summary$`P (whole stations)`[3], NA, NA))

###########################################################
### FIGURE 6: ENVIRONMENTAL dbRDAs AND UNIQUE FRACTIONS ###
# A-B display community scores and environmental arrows at scaling 2. Signs only set the familiar figure orientation.
# C displays the two selected predictor sets; it does not assign ASV effects to individual clusters.
present_clusters <- names(cluster_palette)[names(cluster_palette) %in% unique(as.character(metadata$Cluster))]
present_phases <- sampling_periods[sampling_periods %in% unique(as.character(metadata$Sampling_period))]
site_scores <- \(m) scores(m, display = "sites", choices = 1:2, scaling = 2)
# Axis limits are +/-2, extended in steps of 0.5 if any sample score falls outside.
plot_limit <- max(2, ceiling((max(abs(c(site_scores(m_euk), site_scores(m_prok)))) + .05) * 2) / 2)
reference <- site_scores(m_euk)
f1 <- ifelse(c(mean(reference[as.character(metadata$Cluster) == "6", 1]),
               mean(reference[as.character(metadata$Cluster) == "5", 2])) >= 0, 1, -1)
reference <- sweep(reference, 2, f1, "*")
f2 <- ifelse(vapply(1:2, \(i) cor(reference[, i], site_scores(m_prok)[rownames(reference), i]), numeric(1)) >= 0, 1, -1)
make_db_plot <- function(m, title, flips) {
  si <- as.data.frame(sweep(site_scores(m), 2, flips, "*")) %>% rownames_to_column("Sample") %>%
    left_join(metadata %>% select(Sample, Cluster, Sampling_period), by = "Sample") %>%
    mutate(Cluster = factor(as.character(Cluster), levels = present_clusters),
           Sampling_period = factor(as.character(Sampling_period), levels = present_phases))
  bp <- scores(m, display = "bp", choices = 1:2, scaling = 2) %>% sweep(2, flips, "*") %>%
    as.data.frame() %>% rownames_to_column("Variable") %>%
    mutate(Variable = str_remove_all(Variable, "`"), Label = str_replace(unname(env_labels_db[Variable]),
                                                                         "Conservative temperature", "Conservative\ntemperature"))
  pct <- 100 * as.numeric(eigenvals(m, model = "constrained")[1:2]) / m$tot.chi
  arrow_scale <- 1.55 / max(abs(as.matrix(bp[, c("dbRDA1", "dbRDA2")])))
  label_nudge <- if ("Silicate" %in% bp$Variable) if_else(bp$Variable == "Silicate", .35,
                                                          if_else(bp$Variable == "SWHmax", -.2, 0)) else rep(0, nrow(bp))
  ggplot(si, aes(dbRDA1, dbRDA2, fill = Cluster, shape = Sampling_period, size = Sampling_period)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = .25) + geom_vline(xintercept = 0, colour = "grey85", linewidth = .25) +
    geom_point(colour = "black", stroke = .25, key_glyph = cluster_key) +
    geom_segment(data = bp, aes(x = 0, y = 0, xend = dbRDA1 * arrow_scale, yend = dbRDA2 * arrow_scale),
                 inherit.aes = FALSE, linewidth = .35, arrow = grid::arrow(length = grid::unit(.12, "cm"))) +
    geom_text_repel(data = bp, aes(dbRDA1 * arrow_scale, dbRDA2 * arrow_scale, label = Label),
                    inherit.aes = FALSE, size = 2, fontface = "bold", seed = 1, max.overlaps = Inf, lineheight = .9,
                    box.padding = .2, point.padding = .1, nudge_y = label_nudge, force = 1.5, max.iter = 30000, max.time = Inf,
                    segment.colour = "grey40", segment.size = .2, xlim = c(-plot_limit, plot_limit), ylim = c(-plot_limit, plot_limit)) +
    scale_size_manual(values = plot_sizes[present_phases], guide = "none") +
    scale_fill_manual("Cluster", values = cluster_palette[present_clusters], breaks = present_clusters,
                      labels = cl_labels[match(present_clusters, names(cluster_palette))],
                      guide = guide_legend(order = 1, nrow = 1, title.position = "top", title.hjust = .5,
                                           override.aes = list(shape = NA))) +
    scale_shape_manual("Phase", values = sampling_shapes[present_phases], breaks = present_phases,
                       guide = guide_legend(order = 2, nrow = 1, title.position = "top", title.hjust = .5,
                                            override.aes = list(size = 1.8 * unname(plot_sizes[present_phases]), fill = "grey75", stroke = .35))) +
    labs(tag = substr(title, 1, 2), title = sprintf("%s\n(%.1f%% of variation explained)", substring(title, 4), 100 * RsquareAdj(m)$adj.r.squared),
         x = sprintf("dbRDA1 (%.1f%% of total)", pct[1]), y = sprintf("dbRDA2 (%.1f%% of total)", pct[2])) +
    coord_equal(xlim = c(-plot_limit, plot_limit), ylim = c(-plot_limit, plot_limit), expand = FALSE) + theme_bw(base_size = 8) + black_theme +
    theme(panel.grid = element_blank(), plot.title = element_text(size = 8, face = "bold", lineheight = .9, margin = margin(b = 2)),
          axis.title = element_text(size = 6.4, face = "bold"), axis.text = element_text(size = 6, face = "plain"),
          axis.ticks = element_line(linewidth = .25),
          plot.margin = margin(4, 5, 4, 3))
}
p1 <- make_db_plot(m_euk, "A. dbRDA of eukaryotic composition", f1)
p2 <- make_db_plot(m_prok, "B. dbRDA of prokaryotic composition", f2)

### ### Partitioning panel
# Circle areas are schematic. Labels show adjusted percentages of variation in prokaryotic composition; the overlap is shared.
fraction_pct <- setNames(fractions$`Explained variation (%)`, fractions$Fraction)
angle <- seq(0, 2 * pi, length.out = 361)
venn <- tibble(Block = rep(c("Environmental predictors", "Eukaryotic ASVs"), each = length(angle)),
               x = rep(c(-.8, .8), each = length(angle)) + 1.3 * rep(cos(angle), 2), y = 1.3 * rep(sin(angle), 2))
predictor_colours <- c("Environmental predictors" = "#0072B2", "Eukaryotic ASVs" = "#009E73")
p3 <- ggplot(venn, aes(x, y, group = Block)) +
  geom_polygon(aes(colour = Block), fill = NA, linewidth = .5) +
  scale_colour_manual(values = predictor_colours, guide = "none") +
  annotate("text", x = c(-1.3, 0, 1.3), y = 0, size = 3.2, fontface = "bold",
           label = sprintf("%.1f%%", fraction_pct[c("Environment only", "Shared", "Eukaryotes only")])) +
  annotate("text", x = c(-1.05, 1.05), y = 1.7, size = 2.7, fontface = "bold", lineheight = .95,
           colour = unname(predictor_colours), label = c("Environmental\npredictors", "Eukaryotic ASV\npredictors")) +
  annotate("text", x = 0, y = -1.65, size = 3.2, fontface = "bold",
           label = sprintf("Unexplained: %.1f%%", fraction_pct["Unexplained"])) +
  labs(tag = "C.", title = "Variation partitioning\nof prokaryotic composition") +
  coord_equal(xlim = c(-2.2, 2.2), ylim = c(-2, 2.35), expand = FALSE) + theme_void(base_size = 8) +
  theme(plot.title = element_text(size = 8, face = "bold", lineheight = .9, margin = margin(b = 2)),
        plot.margin = margin(4, 5, 4, 3))
p <- (p1 | p2 | p3) + plot_layout(guides = "collect", widths = c(1, 1, 1)) &
  theme(legend.position = "bottom", legend.box = "horizontal", legend.spacing.x = grid::unit(.3, "cm"),
        legend.title = element_text(size = 7.2, face = "bold"), legend.text = element_text(size = 6.5, face = "plain"),
        legend.key.size = grid::unit(.4, "cm"), legend.margin = margin(0, 0, 0, 0),
        plot.title.position = "plot", plot.tag.position = c(0, 1),
        plot.title = element_text(size = 8, face = "bold", hjust = 0, lineheight = .95, margin = margin(l = 11, b = 6)),
        plot.tag = element_text(size = 8, face = "bold", hjust = 0, vjust = 1, margin = margin()))
ggsave(paste0(out, ".pdf"), p, width = 170, height = 80, units = "mm", device = cairo_pdf)
ggsave(paste0(out, ".png"), p, width = 170, height = 80, units = "mm", dpi = 300, bg = "white")
ggsave(paste0(out, ".svg"), p, width = 170, height = 80, units = "mm")

#####################################
### SAVE NUMBERS USED IN THE TEXT ###
# Model, fraction and selection percentages use adjusted R2; predictor contributions and displayed axes are unadjusted.
# ASV marginal contributions condition on the environment and the other selected ASVs; they are not causal effects.
axis_table <- imap(models[1:2], \(m, nm)
                   tibble(Model = nm, Axis = names(eigenvals(m, model = "constrained"))[1:2],
                          `Total variation (%)` = 100 * as.numeric(eigenvals(m, model = "constrained")[1:2]) / m$tot.chi)) %>% list_rbind()
candidate_table <- candidate_key %>% select(-ASV_full) %>%
  mutate(`Samples with RA >=2%` = as.integer(colSums(candidate_ra >= .02)[ASV]),
         Selected = ASV %in% selected_asvs)
# Include the combined explanation quoted in Results as the final row of the Fractions sheet.
save_tables(list("Models" = model_summary,
                 "Fractions" = bind_rows(fractions, tibble(Fraction = "Total explained",
                                                           `Explained variation (%)` = sum(fractions$`Explained variation (%)`[1:3]))),
                 "Predictor contributions" = predictor_summary,
                 "Selection steps" = selection_steps, "Candidate ASVs" = candidate_table, "Displayed axes" = axis_table),
            paste0(out, "_statistics.xlsx"))
# View(fractions); View(predictor_summary); View(candidate_table %>% filter(Selected))
