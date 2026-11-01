# 6A_ENV.R ---------------------------------------------------------------------------------------
#
# Purpose:
# Map sampling stations and clusters (Figures 1-2) and compare environmental conditions (Figure 3).
# Assemble Tables S1-S7 in one supplementary workbook and save abundance summaries for Results and Discussion.
# Save the metadata and shared settings loaded by 7A_dbRDA.R.
#
# Inputs (paths are relative to the project working directory containing EXT, OUT and DNA):
# EXT/meta.amun.rds; EXT/proptab.amun.original.rds; EXT/proptab.amun.rds (2A)
# EXT/amun_clustering.rds; EXT/amun_bioindicators.rds (3A)
# EXT/amun_shannon_diversity.rds; EXT/amun_shannon_statistics.rds; EXT/amun_composition.rds (4A)
# EXT/ne_10m/ shapefiles; EXT/marmap/marmap_coord_-96;59;-52;78_res_4.csv
#
# Outputs:
# OUT/6A_ENV/amun_station_map.png; amun_clusters_combined.pdf / .png; amun_environment_by_cluster.png (Figures 1-3)
# OUT/6A_ENV/Supplementary materials Tables S1-S7.xlsx (one tab per supplementary table)
# OUT/6A_ENV/amun_group_summaries.xlsx; amun_genera_summaries.xlsx; amun_bioindicator_summaries.xlsx
# EXT/amun_environment_setup.rds -> metadata, scales, colours and symbols used by 7A_dbRDA.R
#
# Run sections in order within this script. Saved inputs allow a fresh R session without running earlier scripts first.
# Existing EXT and OUT folders are assumed to exist.
#---------------------------------------------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(sf)
library(ggspatial)
library(marmap)
library(dunn.test)
library(openxlsx)

### ### Save Excel tables with consistent formatting
# One tab per table; missing entries are NA. Measurements use one decimal, and counts/identifiers use whole numbers.
# P-values and nonzero measurements below 0.1 use scientific notation with one decimal in the mantissa.
# RA fractions in Tables S2-S3 use six decimals (scientific notation below 0.000001); correlations use two.
# Formatting changes only the display; stored values retain full precision.
save_tables <- function(tables, file) {
  wb <- createWorkbook()
  decimal <- createStyle(numFmt = "0.0"); integer <- createStyle(numFmt = "0"); scientific <- createStyle(numFmt = "0.0E+00")
  for (sheet in names(tables)) {
    d <- as.data.frame(tables[[sheet]]) %>% mutate(across(where(is.character), ~ na_if(str_trim(.x), "")))
    addWorksheet(wb, sheet)
    writeDataTable(wb, sheet, d, tableStyle = "TableStyleMedium2", withFilter = TRUE, keepNA = TRUE, na.string = "NA")
    freezePane(wb, sheet, firstRow = TRUE, firstCol = TRUE)
    addStyle(wb, sheet, createStyle(wrapText = TRUE), rows = 1, cols = seq_along(d), gridExpand = TRUE, stack = TRUE)
    setRowHeights(wb, sheet, rows = 1, heights = 32)
    if (nrow(d)) {
      rows <- seq_len(nrow(d)) + 1L; numeric_cols <- which(vapply(d, is.numeric, logical(1)))
      counts <- which(names(d) %in% c("Samples", "Stations", "n", "Cast", "Bottle", "Eukaryotic Reads", "Prokaryotic Reads"))
      pcols <- which(str_detect(names(d), "p-value$"))
      if (length(numeric_cols)) addStyle(wb, sheet, decimal, rows = rows, cols = numeric_cols, gridExpand = TRUE, stack = TRUE)
      if (length(counts)) addStyle(wb, sheet, integer, rows = rows, cols = counts, gridExpand = TRUE, stack = TRUE)
      small <- map_dfr(setdiff(numeric_cols, counts), function(j) {
        r <- which(is.finite(d[[j]]) & d[[j]] != 0 & abs(d[[j]]) < 0.1)
        tibble(Row = r + 1L, Column = rep(j, length(r)))
      })
      if (nrow(small)) addStyle(wb, sheet, scientific, rows = small$Row, cols = small$Column, stack = TRUE)
      if (length(pcols)) addStyle(wb, sheet, scientific, rows = rows, cols = pcols, gridExpand = TRUE, stack = TRUE)
      if (sheet %in% c("Table S2", "Table S3")) {
        addStyle(wb, sheet, createStyle(numFmt = "0.000000"), rows = rows, cols = numeric_cols, gridExpand = TRUE, stack = TRUE)
        for (j in numeric_cols) {
          tiny <- which(is.finite(d[[j]]) & d[[j]] != 0 & abs(d[[j]]) < 1e-6)
          if (length(tiny)) addStyle(wb, sheet, scientific, rows = tiny + 1L, cols = j, stack = TRUE)
        }
      }
      if ("ρ" %in% names(d)) addStyle(wb, sheet, createStyle(numFmt = "0.00"), rows = rows, cols = match("ρ", names(d)), stack = TRUE)
      if (any(vapply(d, function(x) is.character(x) && any(str_detect(x, "\n"), na.rm = TRUE), logical(1)))) {
        addStyle(wb, sheet, createStyle(wrapText = TRUE), rows = rows, cols = seq_along(d), gridExpand = TRUE, stack = TRUE)
        setRowHeights(wb, sheet, rows = rows, heights = 34)
      }
    }
    setColWidths(wb, sheet, cols = seq_along(d), widths = if (sheet %in% c("Table S2", "Table S3")) 22 else "auto")
    if (sheet %in% c("Table S2", "Table S3")) setRowHeights(wb, sheet, rows = 1, heights = 84)
  }
  saveWorkbook(wb, file, overwrite = TRUE)
}

### ### Black theme for plots
# Keep labels, axis numbers and tick marks black; individual panels still control font sizes and hide unused axes.
black_theme <- theme(text = element_text(colour = "black"), axis.text = element_text(colour = "black"),
                     axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
                     axis.title = element_text(colour = "black"), axis.ticks = element_line(colour = "black"),
                     legend.text = element_text(colour = "black"), legend.title = element_text(colour = "black"),
                     strip.text = element_text(colour = "black"), plot.title = element_text(colour = "black"),
                     plot.subtitle = element_text(colour = "black"), plot.caption = element_text(colour = "black"))


#############################################
### LOAD METADATA AND CLUSTER ASSIGNMENTS ###
# Load saved metadata from 2A and cluster assignments from 3A for maps and environmental comparisons.

meta.amun <- readRDS("EXT/meta.amun.rds")
clustering <- readRDS("EXT/amun_clustering.rds")
cl_labels <- clustering$cl_labels; cluster_palette <- clustering$cluster_palette
dend_row_order <- clustering$dend_row_order
meta.amun$Cluster <- clustering$clus[rownames(meta.amun)]


############################################################
### SHARED DEPTH LAYERS, PHASES AND ENVIRONMENTAL SCALES ###
# Define the four depth layers, sampling phases and display scales used by the figures and saved for 7A.

seed <- 161193
depth_levels <- c("Surface (1–6 m)", "Subsurface (9–55 m)", "Intermediate (56–200 m)", "Deep (>225 m)")
# Cutoffs separate the four sampled bands; their display names reflect the observed ranges (there are no 201-225 m samples).
depth_cutoffs <- c(6, 55, 200)
depth_axis_breaks <- c(1, 6, 25, 55, 100, 200, 500, 1000, 2000)
sampling_periods <- c("Leg 2.1", "Leg 2.2", "Leg 2.3", "Leg 4")
sampling_shapes <- c("Leg 2.1" = 21, "Leg 2.2" = 24, "Leg 2.3" = 23, "Leg 4" = 22)
# Base sample-symbol sizes used throughout; dbRDA legend symbols are enlarged for readability.
sampling_sizes <- c("Leg 2.1" = 4.2, "Leg 2.2" = 3.6, "Leg 2.3" = 3.8, "Leg 4" = 4.6)
# Fill legends use borderless rectangles; phase legends retain their sample shapes.
cluster_key <- \(data, params, size) {
  if (is.na(data$shape)) grid::rectGrob(width = .9, height = .55, gp = grid::gpar(fill = data$fill, col = NA))
  else ggplot2::draw_key_point(data, params, size)
}
sampling_period_labels <- c("Leg 2.1" = "Leg 2.1 (13–18 Aug)", "Leg 2.2" = "Leg 2.2 (20 Aug–1 Sep)",
                            "Leg 2.3" = "Leg 2.3 (2–5 Sep)", "Leg 4" = "Leg 4 (9–21 Oct)")

meta_grouped <- meta.amun %>% rownames_to_column("Sample") %>%
  mutate(`Longitude [deg]` = as.numeric(`Longitude [deg]`), `Latitude [deg]` = as.numeric(`Latitude [deg]`),
         Depth_group = cut(`Depth [m]`, c(0, depth_cutoffs, Inf), depth_levels, include.lowest = TRUE),
         Sampling_period = factor(case_when(Cruise == "AM2023_04" ~ "Leg 4",
                                            Cruise == "AM2023_02" & DOY <= 231 ~ "Leg 2.1",
                                            Cruise == "AM2023_02" & DOY < 245 ~ "Leg 2.2",
                                            Cruise == "AM2023_02" & DOY >= 245 ~ "Leg 2.3"), levels = sampling_periods))

### ### Environmental variables and readable labels
# Nine sample properties for cluster comparisons and Shannon correlations; 7A also uses two station properties.
sample_env_vars <- c("Depth [m]", "Conservative temperature [deg C]", "Salinity", "fChl-a [mg/m^3]", "fCDOM [mg/m^3]", "Nitrate [uM]",
                     "Nitrite [uM]", "Phosphate [uM]", "Silicate [uM]")

env_vars <- c(sample_env_vars, "Pycnocline depth [m]", "SWHmax [m]")

env_headers <- c("Depth\n(m)" = "Depth [m]", "Conservative temperature\n(°C)" = "Conservative temperature [deg C]", "Salinity" = "Salinity",
                 "fChl-a\n(mg m⁻³)" = "fChl-a [mg/m^3]", "fCDOM\n(mg m⁻³)" = "fCDOM [mg/m^3]", "Nitrate\n(µM)" = "Nitrate [uM]",
                 "Nitrite\n(µM)" = "Nitrite [uM]", "Phosphate\n(µM)" = "Phosphate [uM]", "Silicate\n(µM)" = "Silicate [uM]",
                 "Pycnocline depth\n(m)" = "Pycnocline depth [m]", "SWHmax\n(m)" = "SWHmax [m]")

env_labels <- setNames(str_replace_all(names(env_headers), "\n", " "), unname(env_headers))
plot_headers <- env_headers[env_headers %in% sample_env_vars]

skewness <- \(x) {
  x <- x[is.finite(x)]; n <- length(x); s <- sd(x)
  n / ((n - 1) * (n - 2)) * sum(((x - mean(x)) / s)^3)
}

# One scale rule for the nine sample properties, assessed across all 115 samples.
# Pycnocline depth and SWHmax intentionally remain untransformed, including in 7A.
# A variable is shown on log1p scale if it is non-negative, raw skew is >=1, and log1p lowers absolute skew by >=0.5.
# Skewness values are rounded to two decimal places before applying these cutoffs.
transform_check <- map_dfr(sample_env_vars, \(v) {
  x <- meta_grouped[[v]][is.finite(meta_grouped[[v]])]
  tibble(
    Variable = v, N = length(x), Unique = n_distinct(x), `Zero (%)` = round(100 * mean(x == 0), 1), Min = min(x), Median = median(x), Max = max(x),
    Skew = round(skewness(x), 2), `Skew log1p` = if (all(x >= 0)) round(skewness(log1p(x)), 2) else NA_real_
  )
}) %>% mutate(Chosen_scale = if_else(!is.na(`Skew log1p`) & Skew >= 1.0 & (abs(Skew) - abs(`Skew log1p`)) >= 0.5, "log1p", "raw"))

log1p_env <- transform_check$Variable[transform_check$Chosen_scale == "log1p"]
display_env <- \(x, variable) if (variable %in% log1p_env) log1p(x) else x

### ### Save the shared settings used by 7A_dbRDA.R
# 7A loads this file to use the same samples, environmental scales, cluster colours and phase symbols.
saveRDS(list(meta_grouped = meta_grouped %>% select(Sample, Cruise, Station, Cluster, Sampling_period, all_of(env_vars)),
             env_vars = env_vars, env_labels = env_labels, log1p_env = log1p_env,
             cl_labels = cl_labels, cluster_palette = cluster_palette, sampling_periods = sampling_periods,
             sampling_shapes = sampling_shapes, sampling_sizes = sampling_sizes, seed = seed), "EXT/amun_environment_setup.rds")


####################################
### SAMPLE METADATA FOR TABLE S1 ###
# Keep only manuscript fields in this publication copy. The analysis metadata remain unchanged.
# Sample links Table S1 to Tables S2-S3; Original station retains the name recorded during sampling.
meta_public <- meta_grouped %>%
  mutate(`Sampling date` = as.Date(sprintf("%04d-%02d-%02d", Year, Month, Day)),
         Cluster = factor(paste0("CL", as.character(Cluster)), levels = cl_labels)) %>%
  select(Sample, Station, `Original station`, `Station group`, Region = geo_loc_name,
         Cruise, `Sampling phase` = Sampling_period, `Sampling date`, `Time [UTC]` = Time,
         Cast, Bottle, `eDNA sample label` = `Sample label`, Target, Layer = Depth_group, Cluster,
         `Latitude [deg]`, `Longitude [deg]`, `Depth [m]`, `Bottom depth [m]`, `Volume filtered [L]`,
         `Conservative temperature [deg C]`, Salinity, `fChl-a [mg/m^3]`, `fCDOM [mg/m^3]`,
         `Nitrate [uM]`, `Nitrite [uM]`, `Phosphate [uM]`, `Silicate [uM]`,
         `Pycnocline depth [m]`, `Chl-a maximum depth [m]`, `SWHmax [m]`, `Eukaryotic Reads`, `Prokaryotic Reads`)


#######################################
### MAP BACKGROUNDS AND STATION MAP ###
# Locate the sampling stations and provide the common geographic background for Figures 1 and 2.
# Uses: meta_grouped, Natural Earth shapefiles and cached bathymetry.
# Produces: map_base and p_inset for Figure 2; OUT/6A_ENV/amun_station_map.png for Figure 1.

sf::sf_use_s2(FALSE)
bb <- st_bbox(c(xmin = -96, xmax = -52, ymin = 59, ymax = 78), crs = 4326)
xlim <- range(meta_grouped$`Longitude [deg]`, na.rm = TRUE) + c(-2, 2)
ylim <- range(meta_grouped$`Latitude [deg]`, na.rm = TRUE) + c(-0.5, 0.5)
rd <- function(path, bounds = bb) {
  suppressWarnings(st_read(path, quiet = TRUE) %>% st_transform(4326) %>% st_make_valid() %>% st_crop(bounds)) %>%
    dplyr::select(geometry, any_of(c("name", "name_en")))
}
ne <- "EXT/ne_10m"
land <- rd(file.path(ne, "ne_10m_land.shp"))
coast <- rd(file.path(ne, "ne_10m_coastline.shp"))
isles <- rd(file.path(ne, "ne_10m_minor_islands.shp"))
islesc <- rd(file.path(ne, "ne_10m_minor_islands_coastline.shp"))

# Densified bounding-box polygon for the Arctic projection
outline <- \(b, pad = c(0, 0, 0, 0), n = 200) {
  x <- c(b["xmin"] - pad[1], b["xmax"] + pad[2])
  y <- c(b["ymin"] - pad[3], b["ymax"] + pad[4])
  xy <- rbind(cbind(seq(x[1], x[2], length.out = n), y[1]), cbind(x[2], seq(y[1], y[2], length.out = n)),
              cbind(seq(x[2], x[1], length.out = n), y[2]), cbind(x[1], seq(y[2], y[1], length.out = n)))
  st_sfc(st_polygon(list(rbind(xy, xy[1, ]))), crs = 4326)
}

# WK-style Arctic inset; the red outline shows the displayed Amundsen map extent
crs_arctic <- 3573
bb_arctic <- st_bbox(c(xmin = -180, xmax = 180, ymin = 35, ymax = 90), crs = 4326)
bb_view <- st_bbox(c(xmin = -179.9, xmax = -10, ymin = 47.5, ymax = 60), crs = 4326)

A <- map(c(land = "ne_10m_land.shp", isles = "ne_10m_minor_islands.shp", coast = "ne_10m_coastline.shp"),
         \(f) rd(file.path(ne, f), bounds = bb_arctic) %>% st_transform(crs_arctic))
view <- st_bbox(st_transform(outline(bb_view, n = 300), crs_arctic))
bb_locator <- st_bbox(c(xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2]), crs = 4326)
amun_box <- st_transform(st_sf(geometry = outline(bb_locator)), crs_arctic)

p_inset <- ggplot() +
  geom_sf(data = A$land, fill = "grey90", colour = NA) +
  geom_sf(data = A$isles, fill = "grey90", colour = NA) +
  geom_sf(data = A$coast, colour = "black", linewidth = 0.25) +
  geom_sf(data = amun_box, fill = NA, colour = "red2", linewidth = 1) +
  coord_sf(crs = crs_arctic, xlim = c(view["xmin"], view["xmax"] - 1000000), ylim = view[c("ymin", "ymax")], expand = FALSE, datum = NA) +
  theme_void() +
  theme(panel.background = element_rect(fill = "#dbeaf4", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 1.2), plot.margin = margin(0, 0, 0, 0))

### ### Load cached bathymetry
# Use 0.5 arc minutes to provide high resolution maps for the manuscript only; use 4 arc minutes to generate maps quickly.
bathy_res <- 4
# Run only if the bathymetry file at the selected resolution has not already been downloaded.
# marmap::getNOAA.bathy(lon1 = as.numeric(bb["xmin"]), lon2 = as.numeric(bb["xmax"]),
#                       lat1 = as.numeric(bb["ymin"]), lat2 = as.numeric(bb["ymax"]), resolution = bathy_res, keep = TRUE, path = "EXT/marmap")
bathy_labs <- c("0–50", "50–200", "200–500", "500–1000", "≥1000")
bathy_cols <- setNames(c("#F1F7FB", "#E6F0F7", "#DBEAF4", "#D0E3EF", "#C5DCE9"), bathy_labs)
bathy_file <- file.path("EXT/marmap", paste0("marmap_coord_", bb["xmin"], ";", bb["ymin"], ";", bb["xmax"], ";", bb["ymax"], "_res_", bathy_res, ".csv"))
bathy_df <- marmap::read.bathy(bathy_file, header = TRUE) %>% fortify.bathy() %>% filter(z <= 0) %>%
  mutate(depth = -z, depth_class = cut(depth, c(0, 50, 200, 500, 1000, Inf), bathy_labs, include.lowest = TRUE, right = FALSE))

### ### Common map background
map_base <- ggplot() +
  geom_tile(data = bathy_df, aes(x, y, fill = depth_class)) +
  scale_fill_manual(values = bathy_cols, breaks = levels(bathy_df$depth_class), drop = TRUE, na.translate = FALSE,
                    name = "Depth (m)", guide = guide_legend(order = 3, ncol = 1, title.position = "top",
                                                             title.hjust = 0, keyheight = unit(1, "cm"), keywidth = unit(1, "cm"))) +
  geom_sf(data = land, fill = "grey90", colour = NA) + geom_sf(data = isles, fill = "grey90", colour = NA) +
  geom_sf(data = coast, colour = "black", linewidth = 0.4) + geom_sf(data = islesc, colour = "black", linewidth = 0.4) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) + theme_minimal() + black_theme +
  theme(panel.grid.major = element_line(colour = "grey85", linewidth = 0.2), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#DBEAF4", colour = NA), legend.box.just = "left",
        axis.title = element_blank(), axis.text = element_text(size = 12, colour = "black"),
        axis.ticks = element_line(colour = "black"), strip.text = element_text(colour = "black", size = 25, face = "bold"),
        plot.title = element_text(colour = "black", size = 27, face = "bold", hjust = 0.5), plot.margin = margin(2, 2, 2, 2))

### ### Draw and save Figure 1
# Labels remain close to stations and retain white backgrounds; connector segments are suppressed.
# Rasterise the bathymetry, land and coastlines; station points and labels remain vectors.
station_labels <- meta_grouped %>% distinct(Cruise, Station, .keep_all = TRUE)
p_locator <- ggrastr::rasterise(map_base, layers = c("Tile", "Sf"), dpi = 300) +
  geom_point(data = station_labels, aes(`Longitude [deg]`, `Latitude [deg]`, shape = Sampling_period, size = Sampling_period),
             colour = "black", fill = "red4", stroke = 0.7) +
  scale_size_manual(values = sampling_sizes, guide = "none") +
  ggrepel::geom_label_repel(data = station_labels, aes(`Longitude [deg]`, `Latitude [deg]`, label = Station),
                            direction = "y", nudge_x = 0.12, hjust = 0, size = 6.5, fontface = "bold", colour = "red4",
                            fill = scales::alpha("white", 0.92), box.padding = 0.08, point.padding = 0.05,
                            label.padding = unit(0.10, "lines"), force = 0.2, force_pull = 20, max.iter = 30000,
                            max.overlaps = Inf, segment.colour = NA, seed = seed) +
  scale_shape_manual(values = sampling_shapes, breaks = sampling_periods, drop = FALSE, name = "Phase",
                     guide = guide_legend(order = 1, title.position = "top", title.hjust = 0, keyheight = unit(0.7, "cm"), keywidth = unit(0.7, "cm"),
                                          override.aes = list(size = 1.2*unname(sampling_sizes), fill = "red4", stroke = 0.7))) +
  annotation_scale(location = "bl", text_cex = 1.4) + annotation_north_arrow(location = "tr", which_north = "true",
                                                                             height = unit(2, "cm"), width = unit(2, "cm")) +
  theme(legend.position = "inside", legend.position.inside = c(0.99, 0.005), legend.justification.inside = c(1, 0),
        legend.direction = "vertical", legend.box = "vertical",
        legend.title = element_text(colour = "black", size = 24, face = "bold"), legend.text = element_text(colour = "black", size = 22),
        legend.background = element_blank(), legend.box.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
        axis.text.x = element_text(size = 20, colour = "black"), axis.text.y = element_text(size = 20, colour = "black")) +
  inset_element(p_inset, left = 0.005, bottom = 0.015, right = 0.335, top = 0.210, align_to = "panel")

ggsave("OUT/6A_ENV/amun_station_map.png", p_locator, width = 14, height = 14.5, units = "in", dpi = 300, bg = "white")
# ggsave("OUT/6A_ENV/amun_station_map.svg", p_locator, device = svglite::svglite, width = 14, height = 14.5, units = "in", bg = "white") # adjust in Inkscape


##############################################
### SAMPLE CLUSTERS ACROSS SPACE AND DEPTH ###
# Relate the fitted communities to geography, sampling phases, profile features and recent waves (Results 3.1).
# Uses: Saved sample dendrogram; meta_grouped; map_base and p_inset from the preceding map section.
# Produces: OUT/6A_ENV/amun_clusters_combined.png / .pdf (Figure 2).

figure2_sizes <- 1.4 * sampling_sizes

## Figure 2A: sample dendrogram with sampling-phase symbols
# Use the existing tree, colours and sample order saved by 3A_CLUS.R.
dd <- dendextend::prepare.ggdend(dendextend::as.ggdend(dend_row_order))
n_samples <- nrow(dd$labels); tree_height <- max(dd$segments$y)
dendro_samples <- dd$labels %>% transmute(Sample = label, Tip = n_samples + 1 - x) %>%
  left_join(meta_grouped %>% select(Sample, Cluster, Sampling_period), by = "Sample")

p_dendro <- ggplot() +
  geom_segment(data = dd$segments, aes(-0.8 * y, n_samples + 1 - x, xend = -0.8 * yend, yend = n_samples + 1 - xend, colour = I(col)), linewidth = 0.7) +
  geom_point(data = dendro_samples, aes(tree_height * 0.035, Tip, fill = factor(Cluster),
                                        shape = Sampling_period, size = Sampling_period), colour = "black", stroke = 0.7) +
  geom_text(data = dendro_samples, aes(tree_height * 0.09, Tip, label = Sample), hjust = 0, size = 7, fontface = "bold") +
  scale_fill_manual(values = cluster_palette, guide = "none") +
  scale_shape_manual(values = sampling_shapes, guide = "none") + scale_size_manual(values = figure2_sizes, guide = "none") +
  scale_x_continuous(limits = c(-0.8*tree_height, 0.9*tree_height), expand = expansion(mult = c(0.002, 0.02))) +
  scale_y_continuous(limits = c(0.5, n_samples + 0.5), expand = c(0, 0)) +
  labs(tag = "A") + theme_void() +
  theme(plot.margin = margin(1, 0, 0, 1), plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 34, face = "bold", colour = "black", hjust = 0, vjust = 1))

## Figure 2B: community clusters in four depth ranges
# Every microbial sample uses the same moderate black outline; chl-max samples are not highlighted separately.
p_depth_maps <- map_base + ggnewscale::new_scale_fill() +
  geom_point(data = meta_grouped, aes(`Longitude [deg]`, `Latitude [deg]`, fill = factor(Cluster),
                                      shape = Sampling_period, size = Sampling_period), stroke = 0.7, colour = "black", key_glyph = cluster_key) +
  scale_size_manual(values = figure2_sizes, guide = "none") +
  scale_fill_manual(values = cluster_palette, guide = "none") +
  scale_shape_manual(values = sampling_shapes, guide = "none") +
  facet_wrap(~ Depth_group, ncol = 2) + annotation_scale(location = "bl", text_cex = 1.6) +
  annotation_north_arrow(location = "tr", which_north = "true", height = unit(1.7, "cm"), width = unit(1.7, "cm")) +
  labs(tag = "B") + theme(panel.spacing = unit(0.6, "cm"), legend.position = "none", axis.text = element_text(size = 18, colour = "black"),
                          strip.text = element_text(size = 30, face = "bold", colour = "black"),
                          plot.margin = margin(1, 0, 18, 10), plot.tag.position = c(0, 1))

## Figure 2C: community clusters across station order and depth
# SWH was calculated in 1A_META.R and is already present in meta.amun; it is not recalculated here.
station_depth <- meta_grouped %>% arrange(Sampling_period, DOY) %>%
  mutate(Station_plot = factor(Station, levels = unique(Station)), Cluster = factor(as.character(Cluster), levels = names(cluster_palette)))

station_structure <- station_depth %>% group_by(Sampling_period, Station_plot) %>%
  summarise(Pycnocline = median(`Pycnocline depth [m]`, na.rm = TRUE), `Chl-a maximum` = median(`Chl-a maximum depth [m]`, na.rm = TRUE),
            Bottom = median(`Bottom depth [m]`, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(c(Pycnocline, `Chl-a maximum`, Bottom), names_to = "Structure", values_to = "Structure_depth")

p_station_depth <- ggplot(station_depth, aes(Station_plot, -log1p(`Depth [m]`), fill = Cluster, shape = Sampling_period, size = Sampling_period)) +
  geom_hline(yintercept = -log1p(c(1, 6, 55, 200, 2000)), colour = "black", linewidth = 0.3) +
  geom_line(data = station_structure, aes(Station_plot, -log1p(Structure_depth), group = Structure, colour = Structure,
                                          linetype = Structure, linewidth = Structure), inherit.aes = FALSE) +
  geom_point(stroke = 0.7, colour = "black", key_glyph = cluster_key) +
  scale_size_manual(values = figure2_sizes, guide = "none") +
  scale_y_continuous(breaks = -log1p(depth_axis_breaks), labels = depth_axis_breaks,
                     expand = expansion(mult = c(0.01, 0.01))) + coord_cartesian(ylim = -log1p(c(2500, 0.5))) +
  scale_fill_manual(values = cluster_palette, breaks = as.character(seq_along(cl_labels)), labels = cl_labels, drop = FALSE, name = "Cluster",
                    guide = guide_legend(order = 1, ncol = 1, title.position = "top", title.hjust = 0,
                                         keyheight = unit(1.075, "cm"), keywidth = unit(1.075, "cm"), override.aes = list(shape = NA))) +
  scale_colour_manual(values = c("Pycnocline" = "blue4", "Chl-a maximum" = "darkgreen", "Bottom" = "black"),
                      limits = c("Pycnocline", "Chl-a maximum", "Bottom"), name = "Depth features",
                      guide = guide_legend(order = 3, ncol = 1, title.position = "top", title.hjust = 0, keywidth = unit(2.5, "cm"),
                                           override.aes = list(linetype = c("dashed", "dotted", "solid"), linewidth = c(1.3, 2, 1.3)))) +
  scale_linewidth_manual(values = c("Pycnocline" = 1.3, "Chl-a maximum" = 2, "Bottom" = 1.3), guide = "none") +
  scale_linetype_manual(values = c("Pycnocline" = "dashed", "Chl-a maximum" = "dotted", "Bottom" = "solid"),
                        limits = c("Pycnocline", "Chl-a maximum", "Bottom"), guide = "none") +
  scale_shape_manual(values = sampling_shapes, breaks = sampling_periods, drop = FALSE, name = "Phase",
                     guide = guide_legend(order = 2, ncol = 1, title.position = "top", title.hjust = 0,
                                          keyheight = unit(1.075, "cm"), keywidth = unit(1.075, "cm"),
                                          override.aes = list(size = 1.1625 * unname(figure2_sizes), fill = "grey75", stroke = 0.7))) +
  facet_grid(cols = vars(Sampling_period), scales = "free_x", space = "free_x",
             labeller = labeller(Sampling_period = as_labeller(sampling_period_labels))) +
  labs(x = NULL, y = "Depth (m)") + theme_bw() + black_theme +
  theme(panel.grid.major.x = element_line(colour = "grey80", linewidth = 0.4), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), axis.text.x = element_text(colour = "black", angle = 60, hjust = 1, vjust = 1, size = 22, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.y = element_text(colour = "black", size = 22, face = "bold", margin = margin(r = 6)),
        strip.text = element_blank(), strip.background = element_blank(), panel.spacing.x = unit(0.15, "lines"), legend.position = "right",
        plot.margin = margin(0, 0, 0, 16))

## Maximum SWH during the 48 h preceding each sample
wave48 <- station_depth %>% group_by(Sampling_period, Station_plot) %>%
  summarise(Max_SWH = median(`SWHmax [m]`, na.rm = TRUE), .groups = "drop") %>%
  mutate(SWH_plot = display_env(Max_SWH, "SWHmax [m]"))
wave_breaks <- 0:ceiling(max(wave48$Max_SWH, na.rm = TRUE))
p_wave48 <- ggplot(wave48, aes(Station_plot, SWH_plot)) +
  geom_col(width = 0.72, fill = "grey80", colour = "black", linewidth = 0.25) +
  facet_grid(cols = vars(Sampling_period), scales = "free_x", space = "free_x",
             labeller = labeller(Sampling_period = as_labeller(sampling_period_labels))) +
  scale_y_continuous(limits = c(0, NA), breaks = display_env(wave_breaks, "SWHmax [m]"),
                     labels = wave_breaks, expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "SWHmax (m)") + theme_bw() + black_theme +
  theme(panel.grid.major.x = element_line(colour = "grey80", linewidth = 0.4), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.y = element_text(colour = "black", size = 22, face = "bold", margin = margin(r = 6)),
        strip.text = element_text(colour = "black", size = 30, face = "bold"),
        strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.6),
        panel.spacing.x = unit(0.15, "lines"), plot.margin = margin(0, 0, 0, 16))

## Assemble Figure 2: A above the legends, B directly above C
# Build the legends separately so Cluster is left of Phase and Depth features.
legend_style <- theme(legend.position = "right", legend.direction = "vertical", legend.box = "vertical",
                      legend.title = element_text(size = 32, face = "bold", colour = "black"),
                      legend.text = element_text(size = 28.75, colour = "black"), legend.key = element_rect(fill = "white", colour = NA))
legend_grob <- \(p) {
  g <- ggplotGrob(p + legend_style)
  g <- g$grobs[[which(vapply(g$grobs, \(x) identical(x$name, "guide-box"), logical(1)))]]
  grid::grobTree(g, vp = grid::viewport(x = grid::unit(0.12, "in"), y = grid::unit(1, "npc") - grid::unit(0.95, "in"),
                                        just = c("left", "top"), width = sum(g$widths), height = sum(g$heights)))
}
p_legends <- wrap_elements(full = legend_grob(p_station_depth + guides(shape = "none", colour = "none"))) +
  wrap_elements(full = legend_grob(p_station_depth + guides(fill = "none", colour = "none"))) +
  wrap_elements(full = legend_grob(p_station_depth + guides(fill = "none", shape = "none"))) +
  plot_layout(design = "AB\nCC", widths = c(0.95, 1.35), heights = c(4.25, 3.75))

p_station <- (p_wave48 + labs(tag = "C") + theme(plot.tag.position = c(-0.005, 1.10))) / (p_station_depth + theme(legend.position = "none")) +
  plot_layout(heights = c(1, 3))
p_left <- wrap_elements(full = p_dendro) / p_legends + plot_layout(heights = c(28, 9))
maps_grob <- grid::grobTree(ggplotGrob(p_depth_maps + labs(tag = NULL)), vp = grid::viewport(y = grid::unit(0.5, "npc") + grid::unit(0.35, "in")))
p_right <- (wrap_elements(full = maps_grob, clip = FALSE) + labs(tag = "B") +
              theme(plot.tag.position = c(-0.005, 1.002), plot.tag = element_text(margin = margin(0, 0, 0, 0)))) / p_station +
  plot_layout(heights = grid::unit(c(-1, 7.0), c("null", "in")))

p_figure2 <- (p_left | p_right) + plot_layout(widths = c(4.5, 22.5)) &
  theme(plot.tag = element_text(size = 34, face = "bold", colour = "black", hjust = 0, vjust = 1))

ggsave("OUT/6A_ENV/amun_clusters_combined.png", p_figure2, width = 30.5, height = 37, units = "in", dpi = 300, bg = "white")


###########################################
### ENVIRONMENTAL PROPERTIES BY CLUSTER ###
# Compare physical, optical and nutrient properties of the eight microbial communities (Results 3.2).
# Uses: All 115 samples in meta_grouped and the nine sample_env_vars; shared display scales.
# Produces: Figure 3, Table S4 environmental summaries and Table S5 Dunn comparisons; tables are exported below.

# Axis positions follow the shared 115-sample scale rule; tick labels remain in measured units.
# Boxes show quartiles and medians; whiskers span the full measured range in each cluster.
plot_data <- meta_grouped %>% mutate(Cluster = factor(as.character(Cluster), levels = names(cluster_palette)))
make_env_box <- \(variable, label) {
  raw <- plot_data[[variable]]
  scale_type <- if (variable == "Depth [m]") "depth" else if (variable %in% log1p_env) "log1p" else "raw"
  y <- display_env(raw, variable); if (scale_type == "depth") y <- -y
  box_data <- plot_data %>% group_by(Cluster) %>% group_modify(\(d, key) {
    q <- display_env(quantile(d[[variable]], c(0, .25, .5, .75, 1), na.rm = TRUE, names = FALSE), variable)
    if (scale_type == "depth") q <- -rev(q)
    as_tibble_row(set_names(q, c("ymin", "lower", "middle", "upper", "ymax")))
  })
  p <- ggplot(plot_data %>% mutate(.y = y), aes(Cluster, .y, fill = Cluster)) +
    geom_boxplot(data = box_data, aes(x = Cluster, ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax, fill = Cluster),
                 stat = "identity", inherit.aes = FALSE, width = .65, alpha = .75, colour = "black", linewidth = .5, na.rm = TRUE) +
    geom_point(aes(shape = Sampling_period, size = Sampling_period), stroke = .4, colour = "black", alpha = .65,
               position = position_jitter(width = .12, height = 0, seed = seed), na.rm = TRUE) +
    scale_size_manual(values = sampling_sizes, guide = "none") + scale_shape_manual(values = sampling_shapes, guide = "none") +
    scale_x_discrete(labels = \(x) paste0("CL", x), drop = FALSE) + scale_fill_manual(values = cluster_palette, drop = FALSE, guide = "none") +
    labs(x = NULL, y = str_replace_all(label, "\n", " ")) + theme_bw() + black_theme +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(colour = "black", size = 13, face = "bold"), axis.text.y = element_text(colour = "black", size = 12),
          axis.title.y = element_text(colour = "black", size = 14, face = "bold", margin = margin(r = 8)),
          axis.ticks = element_line(colour = "black"), plot.margin = margin(6, 7, 6, 7), legend.position = "none")
  if (scale_type == "log1p") {
    ticks <- pretty(range(raw, na.rm = TRUE), n = 5)
    ticks <- ticks[ticks >= 0 & ticks <= max(raw, na.rm = TRUE)]
    p <- p + scale_y_continuous(breaks = log1p(ticks), labels = ticks, expand = expansion(mult = c(.03, .05)))
  }
  if (scale_type == "depth")
    p <- p + scale_y_continuous(breaks = -display_env(depth_axis_breaks, variable), labels = depth_axis_breaks,
                                expand = expansion(mult = c(.02, .02))) +
      coord_cartesian(ylim = -display_env(c(2000, 1), variable))
  p
}

p_environment <- wrap_plots(imap(plot_headers, \(v, label) make_env_box(v, label)), ncol = 3) + plot_layout(guides = "collect") &
  guides(shape = guide_legend(title = "Phase", nrow = 1, title.position = "top", title.hjust = 0.5,
                              override.aes = list(fill = "grey75", size = 1.5 * unname(sampling_sizes), alpha = 1, stroke = 0.7))) &
  theme(legend.position = "bottom", legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 20), legend.key.size = unit(1, "cm"))
ggsave("OUT/6A_ENV/amun_environment_by_cluster.png", p_environment, width = 15, height = 13.5, units = "in", dpi = 300, bg = "white")


######################################################
### TABLE S4 - ENVIRONMENTAL CONDITIONS BY CLUSTER ###
# Summarize physical, optical and nutrient properties and test overall differences among clusters.
# Table S4 contains cluster medians, intervals, main station groups/layers and overall Kruskal-Wallis p-values.

cluster_names <- setNames(as.character(cl_labels), as.character(seq_along(cl_labels)))

### ### Check residual normality
# Check ANOVA residuals on the chosen environmental scales.
# Most were non-normal, so environmental comparisons use Kruskal-Wallis and Dunn consistently.
environment_normality <- map_dfr(sample_env_vars, \(v) {
  d <- na.omit(meta_grouped[, c("Cluster", v)])
  m <- aov(display_env(d[[v]], v) ~ factor(d$Cluster))
  tibble(Variable = v, `Shapiro p-value` = shapiro.test(residuals(m))$p.value)
})

### ### Test overall differences among clusters
# Kruskal-Wallis p-values appear in Table S4 and determine which variables proceed to Dunn's tests for Table S5.
kw_p <- map_dbl(set_names(sample_env_vars), \(var) {
  d <- meta_grouped %>% transmute(Cluster, Value = .data[[var]]) %>% drop_na()
  kruskal.test(Value ~ Cluster, data = d)$p.value
})

### ### Prepare Table S4
# Show median [minimum to maximum] with one decimal, * for missing observations and NA if fewer than half are valid.
# Main station group and Main layer identify the categories with most samples in each cluster; list all ties.
# Keep small nonzero values visible in formatted medians and intervals as well as numeric Excel cells.
number_txt <- function(x) if_else(x != 0 & abs(x) < 0.1, formatC(x, format = "E", digits = 1),
                                  formatC(round(x, 1) + 0, format = "f", digits = 1))
med_range_txt <- function(x) {
  n_total <- length(x); x <- x[is.finite(x)]
  if (!length(x) || length(x) < n_total / 2) return(NA_character_)
  value <- number_txt(c(median(x), min(x), max(x)))
  paste0(value[1], if (length(x) < n_total) "*" else "",
         if (min(x) != max(x)) paste0("\n[", value[2], " to ", value[3], "]") else "")
}
main_category <- function(x) { n <- table(x); paste(names(n)[n == max(n)], collapse = ", ") }
environment_summary <- meta_grouped %>% mutate(Cluster = as.character(Cluster)) %>% group_by(Cluster) %>%
  summarise(Samples = n(), `Main station group` = main_category(`Station group`), `Main layer` = main_category(Depth_group),
            across(all_of(sample_env_vars), med_range_txt), .groups = "drop") %>%
  arrange(as.numeric(Cluster)) %>% mutate(Cluster = unname(cluster_names[Cluster])) %>%
  bind_rows(tibble(Cluster = "Kruskal-Wallis p") %>%
              bind_cols(as_tibble_row(set_names(map_chr(sample_env_vars, function(v) formatC(kw_p[[v]], format = "E", digits = 1)), sample_env_vars)))) %>%
  dplyr::rename(all_of(plot_headers)) %>% mutate(across(where(is.character), function(x) replace_na(x, "NA")))


#############################################################
### TABLE S5 - ENVIRONMENTAL DIFFERENCES BETWEEN CLUSTERS ###
# Apply Dunn's tests after significant Kruskal-Wallis tests; adjust comparisons within each property using BH.
# Save only the parameter, direction based on cluster medians and BH-adjusted p-value for significant comparisons.
environment_dunn <- map_dfr(names(kw_p)[kw_p < 0.01], \(var) {
  d <- meta_grouped %>% transmute(Cluster = as.integer(as.character(Cluster)), Value = .data[[var]]) %>% drop_na()
  invisible(capture.output(dt <- dunn.test(d$Value, d$Cluster, method = "bh", kw = FALSE, altp = TRUE,
                                           table = FALSE, list = FALSE, interpret = FALSE)))
  med <- tapply(d$Value, d$Cluster, median); pair <- str_split_fixed(dt$comparisons, "\\s+-\\s+", 2)
  g1 <- unname(cluster_names[pair[, 1]]); g2 <- unname(cluster_names[pair[, 2]])
  m1 <- as.numeric(med[pair[, 1]]); m2 <- as.numeric(med[pair[, 2]])
  tibble(Parameter = env_labels[[var]],
         `Direction of difference` = case_when(m1 > m2 ~ paste(g1, ">", g2), m2 > m1 ~ paste(g2, ">", g1), TRUE ~ paste(g1, "=", g2)),
         `BH-adjusted p-value` = p.adjust(dt$altP, method = "BH")) # Calculate final BH-adjusted p-values before filtering significant comparisons
}) %>% filter(`BH-adjusted p-value` < 0.01) %>% arrange(Parameter, `BH-adjusted p-value`)


####################################################
### TABLE S7 - SHANNON DIVERSITY AND ENVIRONMENT ###
# Load Shannon estimates from 4A and test associations with the same nine sample properties using available pairs.
# Adjust all nine correlations separately for eukaryotes and prokaryotes before retaining BH-adjusted p <0.01.
shannon_environment <- meta_grouped %>% select(Sample, all_of(sample_env_vars)) %>%
  left_join(readRDS("EXT/amun_shannon_diversity.rds"), by = "Sample") %>%
  pivot_longer(c(Shannon_Euk, Shannon_Prok), names_to = "Domain", values_to = "Shannon") %>%
  pivot_longer(all_of(sample_env_vars), names_to = "Variable", values_to = "Value") %>%
  filter(is.finite(Shannon), is.finite(Value)) %>% group_by(Domain, Variable) %>%
  summarise(n = n(), rho = cor(Shannon, Value, method = "spearman"),
            P = cor.test(Shannon, Value, method = "spearman", exact = FALSE)$p.value, .groups = "drop") %>%
  group_by(Domain) %>% mutate(Q = p.adjust(P, method = "BH")) %>% ungroup() %>% filter(Q < 0.01) %>%                        # BH-adjusted p-values
  transmute(Domain = if_else(Domain == "Shannon_Euk", "Eukaryota", "Prokaryota"), Parameter = unname(env_labels[Variable]),
            n, `ρ` = rho, `BH-adjusted p-value` = Q) %>% arrange(Domain, `BH-adjusted p-value`)


###################################################
### SAVE SUPPLEMENTARY MATERIALS TABLES S1-S7 ###
# This is the only supplementary-table export; 2A and 4A retain their RDS outputs for reuse here.
# Tables S2-S3 retain the original RA fractions (0-1), not percentages, with no renormalization.
supplementary_tables <- list(
  "Table S1" = meta_public,
  "Table S2" = readRDS("EXT/proptab.amun.original.rds") %>% as_tibble(rownames = "Sample"),
  "Table S3" = readRDS("EXT/proptab.amun.rds") %>% as_tibble(rownames = "Sample"),
  "Table S4" = environment_summary,
  "Table S5" = environment_dunn,
  "Table S6" = readRDS("EXT/amun_shannon_statistics.rds")$tukey %>%
    select(Domain, `Direction of difference`, `Tukey-adjusted p-value`),
  "Table S7" = shannon_environment)
save_tables(supplementary_tables, "OUT/6A_ENV/Supplementary materials Tables S1-S7.xlsx")


###################################################
### ABUNDANCE TABLES FOR RESULTS AND DISCUSSION ###
# Use original abundances for Figure 4 groups and genus-level labels, including unclassified taxa; ASV tables contain prevalent bioindicators.
# RA is percent of original eukaryotic or prokaryotic reads. Means include zero abundances.
# Cluster, Station group, Layer and Station group-Layer sheets give mean RA; Samples gives individual RA and conditions.
# Cutoff sheets give sample/station counts and depth, temperature and salinity intervals.
# Cutoffs use unrounded RA and are strict (>); intervals use one decimal and "to" between endpoints.

composition <- readRDS("EXT/amun_composition.rds")
bioindicators <- readRDS("EXT/amun_bioindicators.rds")
indicator_key <- bioindicators$selected_asv_tbl %>%
  transmute(ASV = str_extract(ASV, "ASV_[0-9]+$"), Pattern = str_replace_all(paste0("CL", Pattern), fixed("+"), "+CL"))
sample_asvs <- composition$sample_asvs %>%
  transmute(Sample, Domain = Domain2, Group, Genus, ASV = str_extract(ASV, "ASV_[0-9]+$"), RelAb)
groupings <- list(Cluster = "Cluster", `Station group` = "Station group", Layer = "Layer", `Station group-Layer` = c("Station group", "Layer"))
range_txt <- function(x) if (!any(is.finite(x))) NA_character_ else paste(number_txt(range(x[is.finite(x)])), collapse = " to ")

### ### Calculate sample abundances, spatial means and conditions above each cutoff
taxon_tables <- function(rank) {
  d <- sample_asvs %>% filter(!is.na(.data[[rank]]), .data[[rank]] != "")
  if (rank == "Group") d <- d %>% filter((Domain == "Eukaryote" & Group %in% names(composition$euk_colors)) |
                                           (Domain == "Prokaryote" & Group %in% names(composition$prok_colors)))
  if (rank == "ASV") d <- d %>% inner_join(indicator_key, by = "ASV")
  keys <- c("Domain", "Group", if (rank != "Group") "Genus", if (rank == "ASV") c("ASV", "Pattern"))
  d <- d %>% group_by(across(all_of(c("Sample", keys)))) %>% summarise(RA = 100 * sum(RelAb), .groups = "drop") %>%
    left_join(meta_public %>% select(Sample, Station, Cluster, `Station group`, Layer, `Depth [m]`,
                                     `Conservative temperature [deg C]`, Salinity, `Eukaryotic Reads`, `Prokaryotic Reads`), by = "Sample")
  spatial <- map(groupings, function(by) d %>% group_by(across(all_of(c(by, keys)))) %>%
                   summarise(Samples = n(), `Mean RA` = mean(RA), .groups = "drop") %>%
                   arrange(across(all_of(c(by, "Domain"))), desc(`Mean RA`)))
  cutoffs <- map(set_names(c(0, 25, 50, 75), paste0("Above ", c(0, 25, 50, 75), "% RA")), function(cutoff)
    d %>% filter(RA > cutoff) %>% group_by(across(all_of(keys))) %>%
      summarise(Samples = n(), Stations = n_distinct(Station),
                across(c(`Depth [m]`, `Conservative temperature [deg C]`, Salinity), range_txt), .groups = "drop") %>%
      arrange(Domain, desc(Samples)))
  c(spatial, list(Samples = d %>% select(-Station) %>% arrange(Sample, Domain, desc(RA))), cutoffs)
}

### ### Save sortable tables and retain them in R
group_summaries <- taxon_tables("Group"); save_tables(group_summaries, "OUT/6A_ENV/amun_group_summaries.xlsx")
genera_summaries <- taxon_tables("Genus"); save_tables(genera_summaries, "OUT/6A_ENV/amun_genera_summaries.xlsx")
bioindicator_summaries <- taxon_tables("ASV"); save_tables(bioindicator_summaries, "OUT/6A_ENV/amun_bioindicator_summaries.xlsx")

# View(group_summaries$Layer)
# View(genera_summaries[["Station group-Layer"]])
# View(bioindicator_summaries[["Above 75% RA"]])
