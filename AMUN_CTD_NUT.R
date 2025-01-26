library(tidyverse)

source("natural_earth/get_region_name.R") # for regions
# get_region_name(48.88017, -68.16883) # test

library(factoextra) # for PCA
library(ggrepel)

#####################
### CTD compiling ###

ctd <- "AMUN_CTD" # define the folder containing the metadata CSV files
csv_files <- list.files(ctd, pattern = "\\.csv$", full.names = TRUE) # list all CSV files in the folder

meta_list <- list() # initialize an empty list to store data frames

# Loop through each CSV file (aka cast)
for (file in csv_files) {
  # Read the first 9 rows to extract basic metadata (cruise, cast, date, etc.)
  basic_meta_raw <- read.csv(file, nrows = 9, header = FALSE, stringsAsFactors = FALSE)
  
  # Extract metadata keys (column names) and values
  basic_meta_keys <- gsub("^% |:.*", "", basic_meta_raw$V1)  # Remove "% " prefix and content after ":"
  basic_meta_values <- gsub("^.*: ", "", basic_meta_raw$V1)  # Remove everything before ": "
  
  # Create a data frame of the basic metadata (it would have only 1 row)
  basic_meta <- setNames(as.data.frame(t(basic_meta_values), stringsAsFactors = FALSE), basic_meta_keys)
  
  # Read the main metadata from row 27 onward and select columns A to O
  # Calculate how many rows you got, and repeat appending the basic metadata that many rows 
  meta <- read.csv(file, skip = 26, check.names = FALSE) %>%
    select(1:15) %>%  # Columns A to O
    bind_cols(basic_meta[rep(1, nrow(.)), ])
  
  # Append the resulting metadata for the cast to the list
  meta_list[[file]] <- meta
}

# Combine all casts into a single data frame and wrangle

ctd.amun <- bind_rows(meta_list) %>%
  rename("Bottle" = "BOPO",
         "Sigma-Theta [kg/m^3]" = "SIGP",
         "Sea Pressure [db]" = "PRES",
         "Depth [m]" = "DEPH",
         "Temperature [deg C]" = "TE90",
         "Oxygen solubility [ml/l]" = "OSOL",
         "Salinity" = "PSAL",
         "Sigma-T [kg/m^3]" = "SIGT",
         "Fluorescence [mg/m^3]" = "FLOR",
         "CDOM [mg/m^3]" = "CDOM",
         "PAR [uE/s/m^2]" = "PSAR",
         "Light Transmission [%]" = "TRAN",
         "Dissolved Oxygen [uM]" = "OXYM",
         "Sound Velocity [m/s]" = "SVEL",
         "Nitrate CTD [uM]" = "NTRA",
         "Cruise" = "Cruise_Number",
         "Cast" = "Cast_Number",
         "Latitude [deg]" = "Initial_Latitude [deg]",
         "Longitude [deg]"  = "Initial_Longitude [deg]",
         "Bot. Depth [m]" = "Sounding [m]") %>%
  filter(Station != "test") %>%
  mutate(Cruise = case_when(Cruise == "2023_02" ~ "KEBABB",
                            Cruise == "2023_04" ~ "FOXSIPP"),
         Cast = as.numeric(Cast),   # removing zeroes in front of numbers
         `Latitude [deg]` = as.numeric(`Latitude [deg]`),
         `Longitude [deg]` = as.numeric (`Longitude [deg]`),
         Year = as.numeric(substring(`Start_Date_Time [UTC]`, 1, 4)),  # Extract Year
         Month = as.numeric(substring(`Start_Date_Time [UTC]`, 6, 7)), # Extract Month
         Day = as.numeric(substring(`Start_Date_Time [UTC]`, 9, 10)),  # Extract Day
         Time = substring(`Start_Date_Time [UTC]`, 12),                # Extract Time
         `Bot. Depth [m]` = as.numeric(`Bot. Depth [m]`),
         `Bot. Depth [m]` = if_else(`Bot. Depth [m]` < `Depth [m]`, NA_real_, `Bot. Depth [m]`))  %>%           
  select (- "Cruise_Name", - "Original_Filename", - "Start_Date_Time [UTC]") %>%
  mutate(Station = Station %>%
           str_replace("D1''", "D1-2") %>%
           str_replace("A6/195", "A6") %>%
           str_replace("A7/196", "A7") %>%
           str_replace("A8/197", "A8") %>%
           str_replace("A9/198", "A9") %>%
           str_replace("Sal-01", "SAL-01")) %>%
  mutate(geo_loc_name = map2_chr(`Latitude [deg]`, `Longitude [deg]`, get_region_name)) # region based on the natural_earth data

#####################
### NUT wrangling ###
nut.amun <- read_csv("AMUN_NUT/nutrients_2023_leg_2&4.csv") %>%
  rename ("Bot. Depth [m]" = "Bot.Depth [m]") %>%
  mutate(Cruise = case_when(Cruise == "AN2023_02" ~ "KEBABB",
                            Cruise == "AN2023_04" ~ "FOXSIPP"),
         Year = as.numeric(str_sub(`yyyy-mm-ddThh:mm:ss.sss`, 1, 4)),       # Extract Year (characters 1–4)
         Month = as.numeric(str_sub(`yyyy-mm-ddThh:mm:ss.sss`, 6, 7)),      # Extract Month (characters 6–7)
         Day = as.numeric(str_sub(`yyyy-mm-ddThh:mm:ss.sss`, 9, 10)),       # Extract Day (characters 9–10)
         Time = str_sub(`yyyy-mm-ddThh:mm:ss.sss`, 12, 19),                 # Extract Time (characters 12–19)
         `Bot. Depth [m]` = as.numeric(`Bot. Depth [m]`)) %>%             
  select(-contains("QF"), - "Longitude [degrees_east]", - "Latitude [degrees_north]", - "Type", -"Time", -contains("yyyy")) %>%
  mutate(Station = Station %>%
           str_replace("D1''", "D1-2") %>%
           str_replace("HSS22-003", "HS22-003") %>%
           str_replace("HSS22-005", "HS22-005") %>%
           str_replace("195", "A6") %>%
           str_replace("196", "A7") %>%
           str_replace("197", "A8") %>%
           str_replace("198", "A9")) %>%
  rename ("Nitrate [uM]" = "Nitrate",
          "Nitrite [uM]" = "Nitrite",
          "Phosphate [uM]" = "Phosphate",
          "Silicate [uM]" = "Silicate")


#########################################
### Putting CTD and NUT data together ###

setdiff(ctd.amun$Station, nut.amun$Station) # check which stations don't match
setdiff(nut.amun$Station, ctd.amun$Station) # should be 0 as all nut stations had ctd done on them

intersect(names(ctd.amun), names(nut.amun)) # check common columns

ctd.nut.amun <- ctd.amun %>%
  left_join(nut.amun, by = c(intersect(names(ctd.amun), names(nut.amun)))) %>%
  mutate(`yyyy-mm-ddThh:mm:ss.sss` = paste0(sprintf("%04d", Year), "-",                   # Format Year as 4 digits
                                            sprintf("%02d", Month), "-",                  # Format Month as 2 digits
                                            sprintf("%02d", Day), "T",                  # Format Day as 2 digits and add "T"
                                            ifelse(is.na(Time), "00:00:00", Time))) %>%   # Use "00:00:00" if Time is NA
  select("Cruise", "Station", "Cast", "Latitude [deg]", "Longitude [deg]", "geo_loc_name",
         "yyyy-mm-ddThh:mm:ss.sss", "Year", "Month", "Day", "Time", "Bot. Depth [m]", 
         "Depth [m]", "Sea Pressure [db]", "Bottle", "Temperature [deg C]", "Salinity", "Sigma-T [kg/m^3]", "Sigma-Theta [kg/m^3]",
         "Sound Velocity [m/s]", "Dissolved Oxygen [uM]", "Oxygen solubility [ml/l]",
         "PAR [uE/s/m^2]", "Light Transmission [%]", "CDOM [mg/m^3]", 
         "Fluorescence [mg/m^3]",  "Nitrate CTD [uM]", "Nitrate [uM]", "Nitrite [uM]", "Phosphate [uM]", "Silicate [uM]")


##########################
### EXPLORE and OUTPUT ###

###
### Display the distribution of numerical parameters (histogram and min/max/mean+-sd)
# Purpose: see if the values make sense

colnames(ctd.nut.amun)[sapply(ctd.nut.amun, function(col) !is.numeric(col))] # columns that are not numeric
# Looks like "Cruise", "Station", "geo_loc_name", and "Time"
# "geo_loc_name" out if these is interesting and we might display our numerical parameters at different regions

# Open a PDF to save the histograms
pdf("OUTPUT/amun_parameters_hist.pdf", width = 8.5, height = 12)

# Set up 4 plots per page
par(mfrow = c(3, 2), mar = c(5, 6, 4, 3))  # 2 rows and 2 columns per page; wider left margin

# Loop through numeric columns
for (col_name in names(ctd.nut.amun %>% select(-Cast, -Bottle))) {
  # Check if the column is numeric
  if (is.numeric(ctd.nut.amun[[col_name]])) {
    # Extract the column, removing NA values
    column_data <- na.omit(ctd.nut.amun[[col_name]])
    
    # Skip if the column is empty after removing NAs
    if (length(column_data) == 0) next
    
    # Calculate basic statistics
    col_min <- min(column_data)
    col_max <- max(column_data)
    col_mean <- mean(column_data)
    col_sd <- sd(column_data)
    
    # Generate x-axis range based on the data
    x_range <- range(column_data, na.rm = TRUE)
    
    # Create a simple histogram
    hist(column_data,
         main = NULL,
         xlab = bquote(bold(.(col_name))),
         xlim = x_range,  # Use calculated range for x-axis
         col = "lightblue",
         border = "black")
    
    # Add summary statistics above the graph
    mtext(paste("Min:", round(col_min, 2)), side = 3, line = 2, adj = 0, cex = 0.8)
    mtext(paste("Max:", round(col_max, 2)), side = 3, line = 1, adj = 0, cex = 0.8)
    mtext(paste("Mean ± SD:", round(col_mean, 2), "±", round(col_sd, 2)), side = 3, line = 0, adj = 0, cex = 0.8)
  }
}

# Close the PDF
dev.off()


###
### Display the regions sampled (map)
### Display distributions of numeric parameters by region (boxplot)
### Display how main phys/biochem parameters relate and how similar the regions sampled are based on the behaviour of the parameters (PCA)

# Purpose: assess how parameters behave between regions and relative to each other
         # which parameters correlate? which regions are similar in their phys/biochem state?

# Open pdf
pdf("OUTPUT/amun_parameters_by_region.pdf", width = 12, height = 8.5) # A4 proportions

#
## Create a regional map

# Assign colors per region
# unique(str_remove(ctd.nut.amun$geo_loc_name, "Canada: "))
region_names <- c("The North Western Passages", "Gulf of Boothia", "Fury and Hecla Strait",
                  "Foxe Basin", "Hudson Strait", "Baffin Bay", "Davis Strait",
                  "Labrador Sea", "The Saint Lawrence River", "Nunavut")
custom_colors <- c("The North Western Passages" = "red1",
                   "Gulf of Boothia" = "magenta1",
                   "Fury and Hecla Strait" = "green",
                   "Foxe Basin" = "darkgoldenrod1",
                   "Hudson Strait" = "aquamarine4",
                   "Baffin Bay" = "chocolate1",
                   "Davis Strait" = "mediumpurple1",
                   "Labrador Sea" = "blue",
                   "The Saint Lawrence River" = "cyan",
                   "Nunavut" = "#4F4F4F")

# Load and filter datasets
marine_polys <- st_read("natural_earth/ne_10m_geography_marine_polys/ne_10m_geography_marine_polys.shp") %>% st_transform(4326) %>%
  mutate(name = ifelse(str_detect(name, "River"), paste("The", name), name))
provinces <- ne_states(returnclass = "sf") %>% st_transform(4326)
rivers <- st_read("natural_earth/ne_10m_rivers_lake_centerlines_scale_rank/ne_10m_rivers_lake_centerlines_scale_rank.shp") %>%
  st_transform(4326)
land <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(admin %in% c("Canada", "United States", "Greenland"))
bounding_box <- st_bbox(c(xmin = -125, xmax = -35, ymin = 45, ymax = 85), crs = st_crs(4326))

# Filter regions
marine_regions <- marine_polys %>% filter(name %in% region_names)
province_regions <- provinces %>% filter(name %in% region_names)
river_regions <- rivers %>% filter(name %in% region_names)
all_regions <- bind_rows(st_make_valid(marine_regions), st_make_valid(province_regions), st_make_valid(river_regions))

land_cropped <- st_crop(land, bounding_box)
marine_background <- st_difference(st_union(marine_polys), st_union(all_regions)) %>% st_crop(bounding_box)
all_regions_cropped <- st_crop(all_regions, bounding_box)

# Draw the regional map
map <- ggplot() +
  geom_sf(data = land_cropped, fill = "grey90", color = NA) +
  geom_sf(data = marine_background, fill = "grey80", color = NA) +
  geom_sf(data = all_regions_cropped %>% mutate(name = factor(name, levels = region_names)),
          aes(fill = name), color = "black", lwd = 0.3) +
  scale_fill_manual(values = custom_colors, name = "Regions") +
  labs(title = "Regions sampled during Leg 2 and 4 of Amundsen 2023", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 30)),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, margin = margin(t = 15), face = "bold"),
        axis.title.y = element_text(size = 14, margin = margin(r = 15), face = "bold"),
        plot.margin = margin(t = 40, r = 60, b = 40, l = 60),            # Add left and right margins
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold"))

print(map) # push the map out to the pdf


#
## Make a boxplot for the Bottom Depth
bot_dep_data <- ctd.nut.amun %>%
  group_by(Station) %>%
  reframe(geo_loc_name = unique(str_remove(geo_loc_name, "Canada: ")),
          `Bot. Depth [m]` = unique(`Bot. Depth [m]`)) %>%
  filter(!is.na(`Bot. Depth [m]`)) %>%
  mutate(geo_loc_name = factor(geo_loc_name, levels = region_names))

bot_dep <- ggplot(
  bot_dep_data,
  aes(x = geo_loc_name, y = `Bot. Depth [m]`)) +
  geom_boxplot(aes(color = geo_loc_name), outlier.shape = NA, lwd = 0.3) +
  geom_jitter(width = 0.3, size = 0.8, aes(color = geo_loc_name)) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Bottom depth [m] in the regions sampled", x = NULL, y = "Bottom Depth [m]") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black", face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 15)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 30)),
        plot.margin = margin(t = 40, r = 60, b = 40, l = 80),
        legend.position = "none")

print(bot_dep) # push the boxplot for the station bottom depths to the pdf


#
## Loop over numerical columns and create boxplots for them
for (col_name in colnames(ctd.nut.amun %>%
                          select(-all_of(c("Cast", "Year", "Month", "Day",
                                           "Depth [m]", "Sea Pressure [db]", "Bottle", "Bot. Depth [m]"))) %>%
                          select(where(is.numeric)))){
  
  # Group by geo_loc_name and retain all data without summarizing
  data_to_plot <- ctd.nut.amun %>%
    mutate(geo_loc_name = str_remove(geo_loc_name, "Canada: ")) %>%
    filter(!is.na(.data[[col_name]])) %>% # Exclude rows where the current column is NA
    mutate(geo_loc_name = factor(geo_loc_name, levels = region_names)) # Ensure consistent region order
  
  # Create the boxplot with jitter
  num_var_plot <- ggplot(data_to_plot, aes(x = geo_loc_name, y = .data[[col_name]])) +
    geom_boxplot(aes(color = geo_loc_name), outlier.shape = NA, lwd = 0.3) +
    geom_jitter(width = 0.3, size = 0.7, aes(color = geo_loc_name)) +          # Add data points
    scale_color_manual(values = custom_colors) +
    scale_x_discrete(limits = region_names) +                                  # Consistent x-axis order
    labs(title = paste(col_name, "in the regions sampled"), 
         x = NULL, y = col_name) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black", face = "bold"),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 15)),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 30)),
          plot.margin = margin(t = 40, r = 60, b = 40, l = 80),  # Add left and right margins
          legend.position = "none")
  
  # Print the boxplot to the PDF
  print(num_var_plot)
}


#
## Make a boxplot for the chlorophyll maximum depths

# Prepare the data
max_fluorescence_data <- ctd.nut.amun %>%
  group_by(Station) %>%
  # Calculate the ratio of max fluorescence to the mean fluorescence
  mutate(max_to_mean_ratio = max(`Fluorescence [mg/m^3]`, na.rm = TRUE) / mean(`Fluorescence [mg/m^3]`, na.rm = TRUE)) %>%
  # Filter out stations where the ratio is less than 2
  filter(first(max_to_mean_ratio) >= 2) %>%
  # Select depth with maximum fluorescence
  filter(`Fluorescence [mg/m^3]` == max(`Fluorescence [mg/m^3]`, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(geo_loc_name = str_remove(geo_loc_name, "Canada: ")) %>%
  filter(!is.na(`Depth [m]`)) %>%  # Ensure no NA depths
  mutate(geo_loc_name = factor(geo_loc_name, levels = region_names))

# Create the boxplot
max_fluorescence_plot <- ggplot(max_fluorescence_data,
                                aes(x = geo_loc_name, y = `Depth [m]`)) +
  geom_boxplot(aes(color = geo_loc_name), outlier.shape = NA, lwd = 0.3) +
  geom_jitter(width = 0.3, size = 0.8, aes(color = geo_loc_name)) +
  scale_color_manual(values = custom_colors) +
  labs(title = expression(atop(bold("Chl maximum depths in the regions sampled as inferred from fluorescence"),
                               "(condition for fluorescence: max" >= 2 * " × mean)")),
       x = NULL,
       y = "Depth [m]") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black", face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 15)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 30)),
        plot.margin = margin(t = 40, r = 60, b = 40, l = 80),  # Add left and right margins
        legend.position = "none")

print (max_fluorescence_plot)


#
## Create an ANOVA (y~x) plot for a couple of parameters of interest

# Filter data to exclude rows with NA in the relevant columns
data <- ctd.nut.amun %>%
  filter(!is.na(`Nitrate CTD [uM]`) & !is.na(`Nitrate [uM]`)) %>%
  mutate(geo_loc_name = str_remove(geo_loc_name, "Canada: ")) %>%
  mutate(geo_loc_name = factor(geo_loc_name, levels = region_names))

# Fit a linear regression model
model <- lm(`Nitrate [uM]` ~ `Nitrate CTD [uM]`, data = data)

# Perform ANOVA
anova_result <- anova(model)

# Extract key statistics
r_squared <- summary(model)$r.squared
f_stat <- anova_result$`F value`[1]
p_value <- anova_result$`Pr(>F)`[1]

# Make the graph
nitrate_plot <- ggplot(data, aes(x = `Nitrate CTD [uM]`, y = `Nitrate [uM]`, color = geo_loc_name)) +
  geom_point(size = 1.8) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  scale_color_manual(values = custom_colors, name = "Regions") +
  labs(title = "Relationship between Nitrate CTD [µM] and Nitrate [µM] by region",
       x = "Nitrate CTD [µM]",
       y = "Nitrate [µM]") +
  annotate("text",
           x = max(data$`Nitrate CTD [uM]`) * 0.7,  # Position annotation dynamically
           y = max(data$`Nitrate [uM]`) * 0.9,
           label = paste0(
             "R² = ", signif(r_squared, 5), "\n",
             "F = ", signif(f_stat, 5), "\n",
             "P = ", signif(p_value, 5)),
           hjust = 0, size = 4, color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 15)),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 15)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 30)),
        plot.margin = margin(t = 40, r = 60, b = 40, l = 60),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"))

print(nitrate_plot) # add the ANOVA plot to the pdf


#
## Create a PCA plot with regional centroids

# Get appropriate parameters for the PCA
pca_data <- ctd.nut.amun %>%
  select(where(is.numeric)) %>%              # Keep only numeric columns
  select(-all_of(c("Cast", "Longitude [deg]", "Latitude [deg]",  "Year", "Month", "Day",
                   "Sea Pressure [db]", "Bottle", "Bot. Depth [m]",
                   "Sigma-T [kg/m^3]", "Sigma-Theta [kg/m^3]", "Sound Velocity [m/s]",
                   "Oxygen solubility [ml/l]", "Light Transmission [%]",
                   "Nitrate CTD [uM]"))) %>%
  na.omit() %>%                              # Remove rows with NA values
  select(where(~ var(.) > 0))                # remove constant or zero-variance columns

# Get the associated geo_loc_name
group_data <- ctd.nut.amun %>%
  select(geo_loc_name) %>%  # Replace `geo_loc_name` with your actual grouping column
  mutate (geo_loc_name = str_remove(geo_loc_name, "Canada: ")) %>%
  .[rownames(pca_data), , drop = FALSE] # Ensure `group_data` matches the rows in `numerical_data`

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Compute group centroids for PCA scores
pca_scores <- as.data.frame(pca_result$x) %>%
  mutate(Group = group_data$geo_loc_name) %>%
  group_by(Group) %>%
  summarize(across(starts_with("PC"), mean), .groups = "drop")

# Calculate variance explained by each principal component
explained_variance <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

# Make a PCA plot
pca_plot <- ggplot() +
  geom_point(data = pca_scores %>% mutate(Group = factor(Group, levels = region_names)),
             aes(x = PC1, y = PC2, color = Group), size = 4) +
  geom_segment(data = as.data.frame(pca_result$rotation),
               aes(x = 0, y = 0, xend = PC1*3, yend = PC2*3),  # lengthen the parameter arrows b/c the region points take too much space 
               arrow = arrow(length = unit(0.2, "cm")), color = "grey40") +
  geom_text_repel(data = as.data.frame(pca_result$rotation) %>%
                    mutate(parameter = str_remove(rownames(pca_result$rotation), "\\[.*\\]")),  # Clean parameter names
                  aes(x = PC1*3, y = PC2*3, label = parameter),
                  color = "black", size = 3, fontface = 'bold',
                  max.overlaps = Inf,                # all labels are displayed
                  box.padding = 1.6,                 # padding around labels
                  point.padding = 0,                 # padding around arrowheads
                  force = 15,                        # repelling force for label spacing
                  force_pull = 1,                    # keep labels associated with arrows
                  segment.color = "grey70",          # connectors
                  segment.size = 0.3,                # connectors
                  seed = 7) +
  scale_color_manual(values = custom_colors, name = "Regions") +
  coord_cartesian(xlim = c(-3, 2)) +  # Restrict x-axis range
  labs(title = "PCA plot of CTD and nutrient parameters for the regions sampled",
       x = paste0("PC1 (", explained_variance[1], "%)"),
       y = paste0("PC2 (", explained_variance[2], "%)")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.margin = margin(t = 40, r = 60, b = 40, l = 60),            # Add left and right margins
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold"))

print(pca_plot) # Add the PCA plot to the pdf


#
## Create a PCA plot with individual points colored by region

# Step 1: Filter numeric data and remove unwanted columns
pca_data <- ctd.nut.amun %>%
  select(where(is.numeric)) %>%
  select(-all_of(c("Cast", "Longitude [deg]", "Latitude [deg]", "Year", "Month", "Day",
                   "Sea Pressure [db]", "Bottle", "Bot. Depth [m]",
                   "Sigma-T [kg/m^3]", "Sigma-Theta [kg/m^3]", "Sound Velocity [m/s]",
                   "Oxygen solubility [ml/l]", "Light Transmission [%]",
                   "Nitrate CTD [uM]"))) %>%
  na.omit() %>%
  select(where(~ var(.) > 0))  # Remove constant/zero-variance columns

# Step 2: Perform PCA on the filtered data
pca_result <- prcomp(pca_data, scale. = TRUE)

# Step 3: Perform ANOVA to find significant variables
anova_results <- data.frame(
  parameter = colnames(pca_data),
  p_value = sapply(colnames(pca_data), function(param) {
    anova(lm(pca_result$x[, 1] ~ pca_data[[param]]))$`Pr(>F)`[1]
  })
)

# Step 4: Filter significant parameters based on p-value < 0.05
significant_params <- anova_results %>%
  filter(p_value < 0.05) %>%
  pull(parameter)

# Step 5: Re-filter PCA data to include only significant parameters
pca_data_filtered <- pca_data %>%
  select(all_of(significant_params))

# Recalculate PCA with filtered data
pca_result_filtered <- prcomp(pca_data_filtered, scale. = TRUE)

# Step 6: Add region grouping
group_data <- ctd.nut.amun %>%
  select(geo_loc_name) %>%
  mutate(geo_loc_name = str_remove(geo_loc_name, "Canada: ")) %>%
  .[rownames(pca_data_filtered), , drop = FALSE]  # Ensure matching rows

# Step 7: Create PCA scores for plotting
pca_scores <- as.data.frame(pca_result_filtered$x) %>%
  mutate(Group = group_data$geo_loc_name) %>%
  mutate(Group = factor(Group, levels = region_names))

# Calculate variance explained
explained_variance <- round(100 * pca_result_filtered$sdev^2 / sum(pca_result_filtered$sdev^2), 1)

# Step 8: Make PCA plot with points colored by region
# Extend arrows for parameters and make sample dots smaller and semi-transparent
pca_plot <- ggplot() +
  # Plot PCA sample dots, making them smaller and semi-transparent
  geom_point(data = pca_scores %>%
               filter(abs(PC1) <= quantile(abs(PC1), 0.99) & abs(PC2) <= quantile(abs(PC2), 0.99)), # Remove extreme 1% of points
             aes(x = PC1, y = PC2, color = Group),
             size = 0.6) +                                                                          # Adjust size
  # Add parameter arrows with extended length
  geom_segment(data = as.data.frame(pca_result_filtered$rotation),
               aes(x = 0, y = 0, xend = PC1 * 8, yend = PC2 * 8),  # Extend arrow length
               arrow = arrow(length = unit(0.3, "cm")), color = "grey40", linewidth = 0.5) +
  # Add parameter labels, avoiding overlap
  geom_text_repel(data = as.data.frame(pca_result_filtered$rotation) %>%
                    mutate(parameter = str_remove(rownames(pca_result_filtered$rotation), "\\[.*\\]")),
                  aes(x = PC1*8, y = PC2*8, label = parameter),
                  color = "black", size = 3, fontface = 'bold',
                  max.overlaps = Inf,                # all labels are displayed
                  box.padding = 1,                 # padding around labels
                  point.padding = 0,                 # padding around arrowheads
                  force = 15,                        # repelling force for label spacing
                  force_pull = 5,                    # keep labels associated with arrows
                  segment.color = "grey70",          # connectors
                  segment.size = 0.3,                # connectors
                  seed = 7) +
  # Apply consistent color scheme for regions
  scale_color_manual(values = custom_colors, name = "Regions") +
  # Add labels and title
  labs(title = "PCA plot of significant CTD and nutrient parameters for the regions sampled",
       x = paste0("PC1 (", explained_variance[1], "%)"),
       y = paste0("PC2 (", explained_variance[2], "%)")) +
  # Enhance plot appearance
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.margin = margin(t = 40, r = 60, b = 40, l = 60),            # Add left and right margins
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold"))

print(pca_plot)

dev.off() # close the pdf

##
## Save the compiled CTD and NUT metadata to a CSV file
write.csv(ctd.nut.amun, "OUTPUT/ctd.nut.amun.csv", row.names = FALSE, fileEncoding = "UTF-8")


#############################
### EXTRA (commented out) ###

### silicate/nitrate in deeper waters
# ggplot (ctd.nut.amun %>% filter (`Depth [m]` > 300), aes(x=`Silicate [uM]`/`Nitrate [uM]`, y=`Depth [m]`)) +
#   geom_point(size=1) +
#   geom_smooth(method = "lm", color = "red", se = TRUE)

### high light transmission in deep waters
# ggplot (ctd.nut.amun, aes(x=`Light Transmission [%]`, y=`Depth [m]`)) +
#     geom_point(size=1) +
#     geom_smooth(method = "lm", color = "red", se = TRUE)


