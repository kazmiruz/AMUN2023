# 5A_INDIC.R -------------------------------------------------------------------------------------
#
# Purpose:
# Draw Figure 5: prevalent eukaryotic and prokaryotic bioindicators across all samples.
# Retain the tested cluster or cluster pair in each label; use the higher-mean member only for display order.
#
# Inputs (paths are relative to the project working directory containing EXT, OUT and DNA):
# EXT/meta.amun.rds; EXT/proptab.amun.rds (2A)
# EXT/amun_clustering.rds; EXT/amun_bioindicators.rds (3A)
# EXT/amun_composition.rds (4A; major group assignments and colours)
#
# Outputs:
# OUT/5A_INDIC/amun_bioindicators.pdf / .png (Figure 5)
# No analysis inputs or cluster assignments are modified by this plotting script.
#
# Run sections in order within this script. Later scripts load saved inputs; they do not need earlier scripts in memory.
# Heavy calculations remain optional or cached. Existing EXT and OUT folders are assumed to exist.
#---------------------------------------------------------------------------------------------------

library(tidyverse)
library(ggdendro)
library(patchwork)


#######################################
### LOAD BIOINDICATOR FIGURE INPUTS ###
# Display the prevalent ASVs selected in 3A with the major group colours used in Figure 4.
# Uses: Saved metadata, non-rare ASVs, clustering, indicator results and composition assignments.
# Produces: Data and display settings for Figure 5.

meta.amun <- readRDS("EXT/meta.amun.rds")
proptab.amun <- readRDS("EXT/proptab.amun.rds")
clustering <- readRDS("EXT/amun_clustering.rds")
cl_labels <- clustering$cl_labels; dend_row_order <- clustering$dend_row_order
meta.amun$Cluster <- clustering$clus[rownames(meta.amun)]
bioindicators <- readRDS("EXT/amun_bioindicators.rds")
chosen_cluster <- bioindicators$chosen_cluster; chosen_group <- bioindicators$chosen_group
proptab.display <- proptab.amun[, bioindicators$asv_in_dotplot, drop = FALSE]
composition <- readRDS("EXT/amun_composition.rds")
all_asvs_tax <- composition$all_asvs_tax; euk_colors <- composition$euk_colors; prok_colors <- composition$prok_colors


###############################
### FIGURE SIZE AND SPACING ###

# Scale the complete bioindicator figure for the paper
S <- (24/52)

# File resolution
width_inches <- 64 * S
height_inches <- 45 * S # export height (inches)
res_dpi <- 300          # export resolution (DPI)

# Horizontal widths (plots only)
TREE_SCALE <- 0.6   # tree width
W_TREE <- 0.14      # the border around the tree

W_EUK <- 0.75
W_PROK <- 1.25

W_READS <- 0.065
W_SAMPLE <- 0.075

# Spacer widths (relative). GAP_SMALL can be negative to tighten.
GAP_MED <- -0.002   # gap between Euk and Prok dotplots
GAP_SMALL <- -0.007 # small gaps between oceanography columns on the right side

# ASV label placement (labels drawn *inside* dot-panels at y = 0.5)
LABEL_SIZE <- 30 * S
LABEL_ANGLE <- 60    # the smaller - the more to the left the tail
LABEL_Y_AT <- 0.50   # anchor at panel bottom edge
LABEL_VJUST <- 1.4   # >1 places label just below the anchor
LABEL_XNUDGE <- -0.4 # small left shift (counters rightward drift at angle)
LABEL_DNUDGE <- 0.00 # keep 0 so labels never dip below y = 0.5 (prevents clipping)

BOTTOM_MARGIN_PTS <- 86 * S # bottom-only padding so rotated tails aren't clipped
LEGEND_GAP_HT <- 2.5        # gap between plots row and legends row

euk_colors_dot <- euk_colors
prok_colors_dot <- prok_colors


### ### Functions for labels and abundance symbols
normalize_sample <- function(x) x %>% str_replace_all("[\u2212\u2013\u2014]", "-") %>% str_squish()
ensure_other_color <- function(cols) if (!"Other" %in% names(cols)) c(cols, Other = "black") else cols
to_hex <- function(cols) {m <- col2rgb(cols); apply(m, 2, \(v) rgb(v[1], v[2], v[3], maxColorValue = 255))}
label_total_reads_ticks <- function(breaks) { sapply(breaks, \(b) if (is.na(b)) "" else if (b == 0) "0" else if (b == 50000) "5×10⁴"
                                                     else if (b == 100000) "10⁵" else "") }

get_visual_mapping <- function(RelAb) {
  list(alpha = rep(1, length(RelAb)),                                # all symbols fully opaque
       size = S * if_else(RelAb <= 0.25, 2 + 8 * sqrt(RelAb / 0.25),
                          case_when(RelAb <= 0.50 ~ 8.6, RelAb <= 0.75 ~ 12.4, TRUE ~ 11.4)),
       shape = case_when(RelAb <= 0.25 ~ 16L, RelAb <= 0.50 ~ 17L,
                         RelAb <= 0.75 ~ 18L, TRUE ~ 15L))
}

add_bottom_margin <- function(p, b_pts = BOTTOM_MARGIN_PTS) p + theme(plot.margin = margin(0, 0, b_pts, 0))

### ### Order ASVs by displayed cluster and taxonomic group
prepare_domain_data <- function(data_plot, data_all, group_order, group_colors, asv_mean_info) {
  asv_ordering <- data_plot %>% distinct(ASV) %>% left_join(data_all %>% distinct(ASV, Group), by = "ASV") %>% left_join(asv_mean_info, by = "ASV") %>%
    mutate(IndicatorCluster = unname(chosen_cluster[ASV]), Cluster_order = as.integer(IndicatorCluster),
           Group_order = match(as.character(Group), group_order)) %>%
    arrange(Cluster_order, Group_order, as.character(chosen_group[ASV]), desc(MaxRelAb), ASV) %>% mutate(asv_x = row_number())
  
  extended_colors <- ensure_other_color(group_colors)
  d <- data_plot %>% select(-any_of("Group")) %>% left_join(asv_ordering %>% select(ASV, asv_x, Group), by = "ASV") %>%
    mutate(Group = as.character(Group), Group = if_else(!is.na(Group) & Group %in% names(extended_colors), Group, "Other"),
           Group = factor(Group, levels = names(extended_colors)))
  
  gl <- asv_ordering %>% left_join(data_plot %>% distinct(ASV, label), by = "ASV") %>%
    distinct(asv_x, ASV, label, Group, IndicatorCluster) %>%
    mutate(Group = as.character(Group), Group = if_else(!is.na(Group) & Group %in% names(extended_colors), Group, "Other"),
           Group = factor(Group, levels = names(extended_colors))) %>% arrange(asv_x)
  
  list(d = d, gl = gl, extended_colors = extended_colors)
}

add_indicator_info <- function(obj) {                                           # cluster labels
  pattern <- map2_chr(unname(chosen_group[obj$gl$ASV]), unname(chosen_cluster[obj$gl$ASV]),
                      \(g, first) paste0("CL", c(first, setdiff(str_split(g, "[-+]")[[1]], first)), collapse = "+"))
  obj$gl <- obj$gl %>% mutate(label = paste0(label, " [", pattern, "]"))
  obj
}

build_tick_df <- function(gl_tbl, ph_cols) {                                    # community labels
  ph_cols_ext <- ensure_other_color(ph_cols)
  gl_tbl %>% distinct(asv_x, ASV, label, Group) %>%
    mutate(ph = if_else(as.character(Group) %in% names(ph_cols_ext), as.character(Group), "Other"), col_hex = to_hex(unname(ph_cols_ext[ph]))) %>%
    arrange(asv_x)
}


### ### Draw abundance panels and ASV labels
create_dot_panel <- function(obj, title, ph_cols, cluster_ranges, n_samples, label_angle = LABEL_ANGLE, label_size = LABEL_SIZE, label_y_at = LABEL_Y_AT, label_vjust = LABEL_VJUST) {
  d <- obj$d
  vm <- get_visual_mapping(d$RelAb)
  d <- d %>% mutate(mapped_alpha = vm$alpha, mapped_size = vm$size, shape_var = vm$shape)
  xmx <- max(d$asv_x, na.rm = TRUE) + 0.5
  ticks <- build_tick_df(obj$gl, ph_cols)
  xgrid_df <- tibble(x = sort(unique(ticks$asv_x)))
  ygrid_df <- tibble(y = seq_len(n_samples))
  y_top <- max(ygrid_df$y) + 0.5
  p <- ggplot(d) +
    geom_segment(data = xgrid_df, aes(x = x, xend = x, y = 0.5, yend = max(ygrid_df$y) + 0.5), inherit.aes = FALSE, color = "grey85", linewidth = 0.25) +
    geom_segment(data = ygrid_df, aes(y = y, yend = y, x = 0.5, xend = xmx), inherit.aes = FALSE, color = "grey85", linewidth = 0.25) +
    geom_rect(data = cluster_ranges, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
              fill = NA, color = "black", linewidth = 0.25, linetype = "solid") +
    geom_point(aes(x = asv_x, y = ix, color = Group, alpha = mapped_alpha, size = mapped_size, shape = shape_var), na.rm = TRUE) +
    geom_text(data = ticks, aes(x = asv_x, y = label_y_at, label = label, colour = I(col_hex)),
              angle = label_angle, hjust = 1, vjust = label_vjust, size = label_size / ggplot2::.pt,
              inherit.aes = FALSE, position = position_nudge(x = LABEL_XNUDGE, y = LABEL_DNUDGE)) +
    # geom_text(data = ticks, aes(x = asv_x, y = y_top, label = Community),
    #           fontface = "bold.italic", size = label_size / ggplot2::.pt, inherit.aes = FALSE, vjust = -0.3) +       # community labels if you want
    scale_shape_identity() +
    scale_color_manual(values = ph_cols, limits = names(ph_cols), drop = FALSE, guide = "none") +
    scale_alpha_identity() + scale_size_identity() +
    scale_x_continuous(limits = c(0.5, xmx), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, y_top), expand = c(0, 0)) +
    labs(title = title) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 10) +
    theme(text = element_text(color = "black"),
          plot.title.position = "panel",
          plot.title = element_text(hjust = 0.5, size = LABEL_SIZE + 8, face = "bold", margin = margin(t = 5, r = 0, b = 7, l = 0)),
          axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
          legend.position = "none",
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  add_bottom_margin(p)
}

### ### Draw taxonomic group and abundance legends
make_col_phylum_legend <- function(phylum_names, color_map, title_text, num_cols) {
  n <- length(phylum_names)
  items_per_col <- ceiling(n / num_cols)
  legend_data <- tibble(Phylum = phylum_names,
                        x_dot   = rep(seq(1, by = 1, length.out = num_cols), each = items_per_col)[1:n] - 0.25,
                        x_label = rep(seq(1, by = 1, length.out = num_cols), each = items_per_col)[1:n] - 0.20,
                        y = rep(items_per_col:1, length.out = n))
  ggplot(legend_data, aes(x = x_dot, y = y)) +
    geom_point(aes(color = Phylum), size = 10 * S, show.legend = FALSE) +
    geom_text(aes(x = x_label, y = y, label = Phylum), hjust = 0, vjust = 0.5, size = 13 * S) +
    scale_color_manual(values = color_map, na.translate = FALSE) +
    labs(title = title_text) +
    theme_void() +
    theme(plot.title = element_text(size = LABEL_SIZE + 6, face = "bold", hjust = 0.5), plot.margin = margin(0, 60 * S, 0, 60 * S)) + #incr. margin to part
    scale_x_continuous(limits = c(0.5, max(legend_data$x_label) + 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, max(legend_data$y) + 0.5), expand = c(0, 0)) +
    coord_cartesian(clip = "off")
}

create_clean_legend_with_dots <- function() {
  # Only show phyla that actually occur in the dotplot for this figure
  euk_names  <- intersect(c(names(euk_colors_dot), "Other"),  euk_phyla_used)
  prok_names <- intersect(c(names(prok_colors_dot), "Other"), prok_phyla_used)
  euk_legend_plot  <- make_col_phylum_legend(euk_names, euk_colors_dot_extended[euk_names], "Major Eukaryotic Groups", 2)
  prok_legend_plot <- make_col_phylum_legend(prok_names, prok_colors_dot_extended[prok_names], "Major Prokaryotic Groups", 2)
  vm <- get_visual_mapping(c(0.01, 0.05, 0.10, 0.25))
  ab_df <- tibble(RelAb = c("≤1%", "5%", "10%", "25%", ">25-50%", ">50-75%", ">75%"),
                  alpha_val = 1, size_val = c(vm$size, S * c(8.6, 12.4, 11.4)), shape_val = c(vm$shape, 17L, 18L, 15L),
                  x = c(1.01, 1.26, 1.47, 1.68, 1.99, 2.30, 2.61), y = 2)
  
  abundance_legend_plot <- ggplot(ab_df, aes(x = x, y = y)) +
    geom_point(aes(alpha = alpha_val, size = size_val, shape = shape_val), color = "black") +
    geom_text(aes(label = RelAb), nudge_y = -2.5, size = 12 * S, vjust = 0.5) +
    scale_alpha_identity() + scale_size_identity() + scale_shape_identity() +
    labs(title = "Relative Abundance") +
    theme_void() +
    theme(plot.title.position = "panel",
          plot.title = element_text(size = LABEL_SIZE + 6, face = "bold", hjust = 0.5, vjust = 1, margin = margin(t = 0, b = 0)),
          plot.margin = margin(0, 0, 0, 0)) +
    scale_x_continuous(limits = c(0.8, 2.8), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1.0, 2.2), expand = c(0, 0)) +
    coord_cartesian(clip = "off")
  centers <- c(0.19, 0.63); y_mid <- 0.825; pad <- 0.03; box_sz <- 9 * S
  tr_df <- tibble(kind = c("Eukaryota", "Prokaryota"), fill = c("black", "lightgrey"), cx = centers, cy = y_mid, xtext = centers + pad)
  
  total_reads_plot <- ggplot(tr_df) +
    geom_point(aes(x = cx, y = cy, fill = fill), shape = 22, size = box_sz, stroke = 1, color = "black") +
    geom_text(aes(x = xtext, y = cy, label = kind), hjust = 0, vjust = 0.4, size = 13 * S) +
    labs(title = "Total Reads") +
    scale_fill_identity() +
    theme_void() +
    theme(plot.title = element_text(size = LABEL_SIZE + 6, face = "bold", hjust = 0.5, vjust = 0.5, margin = margin(t = 15)),
          plot.margin = margin(0, 40, 0, 40)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.0, 1.5), expand = c(0, 0)) +
    coord_cartesian(clip = "off")
  right_column <- abundance_legend_plot / total_reads_plot + plot_layout(heights = c(1, 1))
  euk_legend_plot + plot_spacer() + prok_legend_plot + plot_spacer() + right_column +
    plot_layout(nrow = 1, widths = c(1.35, 0, 1.4, 0.25, 1.7)) # legend, spacer, legend, spacer, legend
}


####################################
### ASSEMBLE BIOINDICATOR PANELS ###

build_main_row <- function(proptab.amun) {
  
  # Full data for plotting (NOT restricted to affinity clusters)
  dat_all <- proptab.amun %>% as_tibble(rownames = "Sample") %>%
    pivot_longer(-Sample, names_to = "ASV", values_to = "RelAb") %>%
    left_join(all_asvs_tax %>% select(ASV, Domain2, Group), by = "ASV") %>%
    left_join(meta.amun %>% rownames_to_column("Sample") %>%
                transmute(Sample, Cluster = as.character(Cluster)), by = "Sample") %>%
    mutate(Genus = str_extract(ASV, "(?<=; )[^;]*(?=; ASV_\\d+$)"))
  
  dat_plot <- dat_all %>% filter(RelAb > 0)
  
  # Compute per-ASV cluster with maximum mean abundance
  asv_max_clus <- dat_all %>%
    group_by(ASV, Cluster) %>%
    summarise(mean_rel = mean(RelAb, na.rm = TRUE), .groups = "drop") %>%
    group_by(ASV) %>%
    summarise(
      MaxCluster = Cluster[which.max(mean_rel)],
      MaxRelAb   = max(mean_rel, na.rm = TRUE),
      .groups = "drop")
  
  dd <- dendro_data(dend_row_order, type = "rectangle")
  n <- nrow(dd$labels)
  
  # Build idx_lut: get positions from dendrogram, then join with cluster from metadata
  idx_lut <- tibble(Sample = dd$labels$label, ix = max(dd$labels$x) + 1 - dd$labels$x) %>%
    left_join(meta.amun %>% rownames_to_column("Sample") %>% select(Sample, Cluster), by = "Sample")
  
  dat_all <- dat_all %>% left_join(idx_lut %>% select(-Cluster), by = "Sample") %>% mutate(label = paste0(Genus, " (", str_extract(ASV, "ASV_\\d+$"), ")"))
  dat_plot <- dat_plot %>% left_join(idx_lut %>% select(-Cluster), by = "Sample") %>% mutate(label = paste0(Genus, " (", str_extract(ASV, "ASV_\\d+$"), ")"))
  cluster_ranges <- idx_lut %>% group_by(Cluster) %>%
    summarise(ymin = min(ix) - 0.5, ymax = max(ix) + 0.5, y_center = (min(ix) + max(ix)) / 2, .groups = "drop")
  
  # Panels
  E <- prepare_domain_data(dat_plot %>% filter(Domain2 == "Eukaryote"), dat_all %>% filter(Domain2 == "Eukaryote"),
                           names(euk_colors_dot), euk_colors_dot, asv_max_clus) %>% add_indicator_info()
  P <- prepare_domain_data(dat_plot %>% filter(Domain2 == "Prokaryote"), dat_all %>% filter(Domain2 == "Prokaryote"),
                           names(prok_colors_dot), prok_colors_dot, asv_max_clus) %>% add_indicator_info()
  
  # Which phyla actually appear in this figure (incl. 'Other' only if used)
  used_euk_phyla  <- sort(unique(as.character(E$d$Group)))
  used_prok_phyla <- sort(unique(as.character(P$d$Group)))
  
  euk_cols_ext <- ensure_other_color(E$extended_colors)
  prok_cols_ext <- ensure_other_color(P$extended_colors)
  p_euk <- create_dot_panel(E, "Eukaryota", ph_cols = euk_cols_ext, cluster_ranges = cluster_ranges, n_samples = n)
  p_prok <- create_dot_panel(P, "Prokaryota", ph_cols = prok_cols_ext, cluster_ranges = cluster_ranges, n_samples = n)
  
  # Total Reads (stacked bars)
  bar_wide <- idx_lut %>% mutate(Sample_key = normalize_sample(Sample)) %>%
    left_join(meta.amun %>% rownames_to_column("Sample") %>% mutate(Sample_key = normalize_sample(Sample)) %>%
                select(Sample_key, `Eukaryotic Reads`, `Prokaryotic Reads`), by = "Sample_key") %>%
    mutate(`Eukaryotic Reads` = as.numeric(coalesce(`Eukaryotic Reads`, 0)),
           `Prokaryotic Reads` = as.numeric(coalesce(`Prokaryotic Reads`, 0)),
           y0 = ix - 0.4, y1 = ix + 0.4)
  bar_segments <- bind_rows(transmute(bar_wide, ix, xmin = 0, xmax = `Eukaryotic Reads`, ymin = y0, ymax = y1, read_type = "Eukaryotes"),
                            transmute(bar_wide, ix, xmin = `Eukaryotic Reads`, xmax = `Eukaryotic Reads` + `Prokaryotic Reads`,
                                      ymin = y0, ymax = y1, read_type = "Prokaryotes"))
  p_stacked_bar <- ggplot(bar_segments) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = read_type), color = "black") +
    scale_fill_manual(values = c("Eukaryotes" = "black", "Prokaryotes" = "lightgrey")) +
    geom_rect(data = cluster_ranges, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = NA, color = "black", linewidth = 0.25, linetype = "solid") +
    scale_x_continuous(breaks = c(0, 50000, 100000), labels = label_total_reads_ticks,
                       expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    labs(x = "Total Reads", y = NULL) +
    guides(fill = "none") +
    theme_minimal() +
    theme(text = element_text(color = "black"),
          axis.line.x = element_line(linewidth = 0.25, color = "black"),
          axis.ticks.x = element_line(linewidth = 0.25, color = "black"),
          panel.grid.major.x = element_line(linewidth = 0.25, color = "grey85"),
          panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(size = 27 * S, color = "black"),
          axis.title.x = element_text(size = 29 * S, face = "bold", margin = margin(t = -60, r = 0, b = 0, l = 0), vjust = -0.25))
  p_stacked_bar <- add_bottom_margin(p_stacked_bar)
  
  # Dendrogram
  dd_segs <- dd$segments
  xmax <- max(dd_segs$y, dd_segs$yend, na.rm = TRUE)
  seg_tree <- dd_segs %>% mutate(leaf_x = x, leaf_xend = xend) %>%
    transmute(x = TREE_SCALE * (xmax - y), xend = TREE_SCALE * (xmax - yend), y = n + 1 - leaf_x, yend = n + 1 - leaf_xend)
  p_tree <- ggplot(seg_tree) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.5, lineend = "round", na.rm = TRUE, color = "black") +
    geom_rect(data = cluster_ranges, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
              fill = NA, color = "black", linewidth = 0.25, linetype = "solid", inherit.aes = FALSE) +
    geom_text(data = cluster_ranges, aes(x = -0.7, y = y_center, label = cl_labels[as.integer(as.character(Cluster))]),
              color = "black", size = 14 * S, fontface = "bold", inherit.aes = FALSE, hjust = 0.5, vjust = 0.5) +
    # annotate("text",  x = xmax - 0.2, y = n + 0.4, label = "Community", hjust = 1, vjust = -0.3, fontface = "bold.italic", size = 10) + # left/top/right-adj
    scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    scale_x_continuous(limits = c(-1.5, TREE_SCALE * xmax), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void() + theme(plot.margin = margin(t = 20 * S, r = 0, b = BOTTOM_MARGIN_PTS, l = 0))
  
  # Side-panels
  panel_df <- idx_lut
  
  create_component_plot <- function(component, title_txt, box_right = 0.8, pad_rel = 0.05) {
    pad_data <- box_right * pad_rel
    ggplot(panel_df) +
      geom_rect(data = cluster_ranges, aes(xmin = -Inf, xmax = box_right, ymin = ymin, ymax = ymax),
                fill = NA, color = "black", linewidth = 0.25, linetype = "solid", inherit.aes = FALSE) +
      geom_text(aes(y = ix, label = !!sym(component)), x = box_right - pad_data, hjust = 1, size = 7.8 * S, color = "black") +
      scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
      scale_x_continuous(limits = c(0, box_right), expand = c(0, 0)) +
      coord_cartesian(clip = "on") +
      labs(title = title_txt) +
      theme_void() +
      theme(plot.title = element_text(size = LABEL_SIZE-1, face = "bold", hjust = 0.5, vjust = 0, margin = margin(b = 5)))
  }
  
  side_specs <- "Sample"
  width_map <- c("Sample" = W_SAMPLE)
  
  side_plots <- side_specs %>% set_names() %>% map(~ add_bottom_margin(create_component_plot(.x, .x)))
  plots_row <- list(p_tree, p_euk, plot_spacer(), p_prok, plot_spacer(), p_stacked_bar)
  widths_row <- c(W_TREE, W_EUK, GAP_MED, W_PROK, GAP_SMALL, W_READS)
  for (nm in names(side_plots)) {
    plots_row <- append(plots_row, list(plot_spacer(), side_plots[[nm]]))
    widths_row <- c(widths_row, GAP_SMALL, unname(width_map[nm]))
  }
  list(
    plot_row        = wrap_plots(plots_row, nrow = 1) + plot_layout(widths = widths_row),
    euk_cols_ext    = euk_cols_ext,
    prok_cols_ext   = prok_cols_ext,
    euk_phyla_used  = used_euk_phyla,
    prok_phyla_used = used_prok_phyla,
    E               = E,
    P               = P
  )
}


#####################
### SAVE FIGURE 5 ###

# --- Saving the plots ---
save_plot <- function(p, path, w, h, dpi){
  if (tolower(tools::file_ext(path)) == "pdf")
    cairo_pdf(path, width = w, height = h, onefile = FALSE)
  else
    png(path, width = w, height = h, units = "in", res = dpi, bg = "white", type = "cairo")
  print(p); dev.off()
}

### ### Combine panels and legends
row_display <- build_main_row(proptab.display)
E <- row_display$E; P <- row_display$P
euk_colors_dot_extended <- row_display$euk_cols_ext; prok_colors_dot_extended <- row_display$prok_cols_ext
euk_phyla_used <- row_display$euk_phyla_used; prok_phyla_used <- row_display$prok_phyla_used

legend_plot <- create_clean_legend_with_dots()
final_display <- row_display$plot_row / plot_spacer() / legend_plot + plot_layout(heights = c(10, LEGEND_GAP_HT, 1)) +
  plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))

save_plot(final_display, "OUT/5A_INDIC/amun_bioindicators.pdf", width_inches, height_inches, res_dpi)
save_plot(final_display, "OUT/5A_INDIC/amun_bioindicators.png", width_inches, height_inches, res_dpi)
