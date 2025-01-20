library(tidyverse)
library(janitor)

source("natural_earth/get_region_name.R") # for regions
# get_region_name(48.88017, -68.16883) # test

############################################################
### Load and clean the genohub and sampling spreadsheets ###

amun_genohub <- read_csv("AMUN_DNA/meta.amun.csv") %>%
  select(1:(ncol(.) - 1)) %>%
  select("sample.id", "station", "target.depth", "box", "plate", "genohub.id", "plate.geno") %>%  
  rename ("Station" = "station")

#dna.conc in the genohub spreadsheet is not as precise as in the sample spreadsheets but they did send samples with 0.5+ ng/μl
min(read_csv("AMUN_DNA/meta.amun.csv")$dna.conc)


amun_sample_leg2 <- read_csv("AMUN_DNA/Amundsen_sampling_2023_Leg2_csv.csv") %>%
  slice(1:(n() - 2)) %>%
  mutate(pore.size = 0.22) %>%
  mutate("Sequential filtration" = 0) %>%
  mutate(Cruise = "KEBABB") %>%
  select("Cruise", "Sample label", "Station", "Sampling date", "Sampling time", "Bottle", "Target Depth", "Prefiltration mesh pore size, μm", "Sequential filtration", "pore.size", "Volume filtered, L", "DNAavg, ng/μl", "DNA Box")

amun_sample_leg4 <- read_csv("AMUN_DNA/Amundsen_sampling_2023_Leg4_csv.csv") %>%
  slice(1:(n() - 2)) %>%
  # rename_with(~ str_replace_all(., "\\?", "μ")) %>%
  pivot_longer(                                                     # getting the pore size of the filters
    cols = contains("Box"),
    names_to = "pore.size.box", 
    values_to = "DNA Box", 
    values_drop_na = TRUE) %>%
  pivot_longer(
    cols = contains("DNAavg"), 
    names_to = "pore.size.conc", 
    values_to = "DNAavg, ng/μl", 
    values_drop_na = TRUE) %>%
  mutate(pore.size = case_when(
      str_detect(pore.size.box, "0.22 μm") ~ 0.22,
      str_detect(pore.size.box, "3 μm") ~ 3,
      str_detect(pore.size.box, "20 μm") ~ 20,
      TRUE ~ NA_real_)) %>%                                         # numeric NA if none of the conditions are met
  filter(str_detect(pore.size.conc, as.character(pore.size))) %>%   # ensure proper alignment by filtering matching pore sizes
  mutate(Cruise = "FOXSIPP") %>%
  select("Cruise", "Sample label", "Station", "Sampling date", "Sampling time", "Bottle", "Target Depth", "Prefiltration mesh pore size, μm", "Sequential filtration", "pore.size", "Volume filtered, L", "DNAavg, ng/μl", "DNA Box")


# Bind the sample spreadsheets
amun_sample <- bind_rows(amun_sample_leg2, amun_sample_leg4) %>%
  rename("sample.id" = "Sample label",
         "Date" = "Sampling date",
         "target.depth" = "Target Depth",
         "prefilter.size" = "Prefiltration mesh pore size, μm",
         "seq.filt" = "Sequential filtration",
         "volume.l" = "Volume filtered, L",
         "dna.conc" = "DNAavg, ng/μl",
         "box" = "DNA Box") %>%
  # mutate(dna.conc = round(dna.conc, 2)) %>% # round to two decimal places
  mutate(volume.l = round(dna.conc, 1)) %>% # round to two decimal places
  mutate(Date = dmy(Date),                    # Convert to Date format
         Year = year(Date),                   # Extract year
         Month = month(Date)) %>%             # Extract month
                                              # do not extract any more than that b/c time/day might be different (e.g., UTC vs local)
  filter(dna.conc >= 0.5) # select DNA concentrations 0.5 ng/μl b/c only these were sent


##################################################################
### Check the genohub and sample spreadsheets for integrity ###

# check for NAs
sum(is.na(amun_genohub$sample.id)) # check for completeness
sum(is.na(amun_genohub$Station))   # check for completeness
sum(is.na(amun_genohub$target.depth))  # check for completeness
sum(is.na(amun_genohub$box))       # check for completeness

sum(is.na(amun_sample$sample.id)) # check for completeness
sum(is.na(amun_sample$Station))   # check for completeness (some station names were corrupted into dates - fixed)
sum(is.na(amun_sample$target.depth))  # check for completeness
sum(is.na(amun_sample$box))       # check for completeness
sum(is.na(amun_sample$dna.conc))  # check for completeness

# since each sample with the same label (station_depth) but different pore size was put into a different box,
# we can match the genohub and sample spreadsheets by the sample.id and box alone
length(unique(amun_genohub$sample.id)) #only 183/190
nrow(unique(amun_genohub[, c("sample.id", "box")])) # 190/190

# station_depth + box combinations are also unique for the sample spreadsheet in general
length(unique(amun_sample$sample.id)) #only 198/218
nrow(unique(amun_sample[, c("sample.id", "box")])) # 218/218

# 13% of samples (28/218) were not sent in spite of DNA being >= 0.5 ng/μl - ask EC why
not_in_genohub <- amun_sample %>%
  anti_join(amun_genohub, by = c("sample.id", "box"))


#############################
### Bind the spreadsheets ###

dna.amun <- amun_genohub %>%
  inner_join(amun_sample, by = c("sample.id", "Station", "target.depth", "box")) %>%
  select ("Cruise", "Station", "Year", "Month", "Bottle", "sample.id", "target.depth", "prefilter.size", "seq.filt", "pore.size", "volume.l", "dna.conc", "box", "plate", "genohub.id", "plate.geno") %>%
  mutate(Station = Station %>%
           str_replace("D1\"", "D1-2") %>% 
           str_replace("D1-2022", "D1_2022"))

# check if any genohub samples did not match to the ones from the sample spreadsheet, should be 0
setdiff(amun_genohub$sample.id, dna.amun$sample.id)


###############################################
### Attach CTD and NUT data to the DNA data ###

# should have ctd.nut.amun off of AMUN_CTD_NUT.R

#
## before you start, check where you had a 2 casts per station
# remove station-cast groups where you (and most other main biochemists) did not collect samples for a proper match
# checked CTD logbooks first - helpful, but, it has many mistakes and lucks useful pattern, so cannot be directly integrated and used
# then, checked individual cast spreadsheets to ensure I am filtering out the right station-casts

# in the future, record which casts you collected the samples from to avoid going back to CTD spreadsheets

ctd.logbook <- rbind(
  read.csv("AMUN_CTD/logbook/2023_04_CTD_logbook.csv"),
  read.csv("AMUN_CTD/logbook/2023_02_CTD_logbook.csv"))

cast_duplicates <- data.frame(
  Station = c("E5", "E3", "C5", "C3", "A5", "A3", "323", "163", "108", "D1", "350"),  # stations
  Cast = c(8, 5, 10, 13, 22, 19, 64, 59, 46, 30, 41))                                 # casts


#
## since ctd.nut.amun lacks river data, extract and modify river data to add to the main spreadsheet
amun_sample_river <- read_csv("AMUN_DNA/Amundsen_sampling_2023_Leg4_csv.csv") %>%
  filter(str_detect(Station, "AMD")) %>%                                              # need this to pull data for the rivers separate
  mutate(geo_loc_name = map2_chr(`Latitude`, `Longitude`, get_region_name)) %>%  # region based on the natural_earth data
  mutate(Date = dmy(`Sampling date`),                    # Convert to Date format
         Year = year(`Date`),                   # Extract year
         Month = month(`Date`),                 # Extract month
         Day = day (`Date`),
         `Sampling Time` = NA) %>%             
  mutate(`yyyy-mm-ddThh:mm:ss.sss` = paste0(sprintf("%04d", Year), "-",                              # Format Year as 4 digits
                                            sprintf("%02d", Month), "-",                             # Format Month as 2 digits
                                            sprintf("%02d", Day), "T00:00:00"))                      # Format Day as 2 digits and add "T"


#
## append ctd.nut.amun to dna.amun and add river data

intersect(names(dna.amun), names(ctd.nut.amun))                                       # check common columns

meta.amun.complete <- dna.amun %>%
  left_join(ctd.nut.amun %>%
              anti_join(cast_duplicates, by = c("Station", "Cast")),        # Exclude station-casts that were not sampled
            by = c("Cruise", "Station", "Year", "Month", "Bottle")) %>%     # Join keys
  select(all_of(c("Cruise", "Station", "Cast", "Latitude [deg]", "Longitude [deg]", "geo_loc_name",
                  "yyyy-mm-ddThh:mm:ss.sss", "Year", "Month", "Day", "Time", "Bot. Depth [m]", 
                  "Depth [m]", "Sea Pressure [db]", "Bottle")), everything()) %>%
  mutate(`Latitude [deg]` = case_when(
    Station == "AMD23-29" ~ as.numeric(amun_sample_river %>% filter(Station == "AMD23-29") %>% pull(Latitude)),
    Station == "AMD23-6" ~ as.numeric(amun_sample_river %>% filter(Station == "AMD23-6") %>% pull(Latitude)),
    TRUE ~ as.numeric(`Latitude [deg]`)),     # Retain existing values as numeric
    `Longitude [deg]` = case_when(
      Station == "AMD23-29" ~ as.numeric(amun_sample_river %>% filter(Station == "AMD23-29") %>% pull(Longitude)),
      Station == "AMD23-6" ~ as.numeric(amun_sample_river %>% filter(Station == "AMD23-6") %>% pull(Longitude)),
      TRUE ~ as.numeric(`Longitude [deg]`)),  # Retain existing values as numeric
    geo_loc_name = case_when(
      Station == "AMD23-29" ~ amun_sample_river %>% filter(Station == "AMD23-29") %>% pull(geo_loc_name),
      Station == "AMD23-6" ~ amun_sample_river %>% filter(Station == "AMD23-6") %>% pull(geo_loc_name),
      TRUE ~ geo_loc_name),                   # Retain existing values as numeric
    `yyyy-mm-ddThh:mm:ss.sss` = case_when(
      Station == "AMD23-29" ~ amun_sample_river %>% filter(Station == "AMD23-29") %>% pull(`yyyy-mm-ddThh:mm:ss.sss`),
      Station == "AMD23-6" ~ amun_sample_river %>% filter(Station == "AMD23-6") %>% pull(`yyyy-mm-ddThh:mm:ss.sss`),
      TRUE ~ `yyyy-mm-ddThh:mm:ss.sss`),      # Retain existing values
    Day = case_when(
      Station == "AMD23-29" ~ as.numeric(amun_sample_river %>% filter(Station == "AMD23-29") %>% pull(Day)),
      Station == "AMD23-6" ~ as.numeric(amun_sample_river %>% filter(Station == "AMD23-6") %>% pull(Day)),
      TRUE ~ as.numeric(Day)))                # Retain existing values as numeric


##########################
### EXPLORE AND OUTPUT ###

###
### Check bottle and DNA concentrations at the target depths for validity and potential patterns

# Open pdf
pdf("OUTPUT/amun_micb_sampling_validity.pdf", width = 12, height = 8.5) # A4 proportions

#
## Check Niskin bottles at the target depths

# Generate a gradient for depths from 20 m to 200 m
additional_depths <- c("20 m", "30 m", "40 m", "50 m", "60 m", "80 m", 
                       "120 m", "140 m", "150 m", "160 m", "180 m", "200 m")
depth_colors <- colorRampPalette(c( "orange", "magenta", "red"))(length(additional_depths))

# Combine all depths into a single color map
all_depths <- c("surf", "chl max", "bottom", additional_depths)
color_mapping <- c("Surface" = "blue", "Chl maximum" = "green", "Bottom" = "black", 
                   setNames(depth_colors, additional_depths))

bottle_check <- ggplot(meta.amun.complete %>%
                         filter(target.depth %in% all_depths) %>%
                         mutate(target.depth = factor(target.depth,
                                                      levels = all_depths,
                                                      labels = c("Surface", "Chl maximum", "Bottom", additional_depths))),
                       aes(x = target.depth, y = (Bottle))) +
  geom_boxplot(aes(color = target.depth), outlier.shape = NA) +
  geom_jitter(width = 0.3, size = 1.5, aes(color = target.depth)) +
  scale_color_manual(values = color_mapping) +
  labs(title = "Microbial sampling: bottles closed at the target depths",
       x = NULL,
       y = "Bottle number") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black", face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 15)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 30)),
        plot.margin = margin(t = 40, r = 60, b = 40, l = 80),
        legend.position = "none")

print (bottle_check)

#
## Check DNA concentrations

# Generate a gradient for depths from 20 m to 200 m
additional_depths <- c("20 m", "30 m", "40 m", "50 m", "60 m", "80 m", 
                       "120 m", "140 m", "150 m", "160 m", "180 m", "200 m")
depth_colors <- colorRampPalette(c( "orange", "magenta", "red"))(length(additional_depths))

# Combine all depths into a single color map
all_depths <- c("surf", "chl max", "bottom", additional_depths)
color_mapping <- c("Surface" = "blue", "Chl maximum" = "green", "Bottom" = "black", 
                   setNames(depth_colors, additional_depths))

dna_check <- ggplot(meta.amun.complete %>%
                      filter(target.depth %in% all_depths) %>%
                      mutate(target.depth = factor(target.depth,
                                                   levels = all_depths,
                                                   labels = c("Surface", "Chl maximum", "Bottom", additional_depths))),
                    aes(x = target.depth, y = (dna.conc))) +
  geom_boxplot(aes(color = target.depth), outlier.shape = NA) +
  geom_jitter(width = 0.3, size = 1.5, aes(color = target.depth)) +
  scale_color_manual(values = color_mapping) +
  labs(title = "Microbial sampling: DNA concentration of samples from the target depths",
       x = NULL,
       y = "DNA concentration [ng/μl]") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black", face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 15)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 30)),
        plot.margin = margin(t = 40, r = 60, b = 40, l = 80),
        legend.position = "none")

print (dna_check)

dev.off ()


###
### Save the metadata for the microbial samples to a CSV file
write.csv(meta.amun.complete, file = "OUTPUT/meta.amun.complete.csv", row.names = FALSE)


###
### Save the information on the DNA samples that could be but were not sequenced to a CSV file
write.csv(not_in_genohub, file = "OUTPUT/not_in_genohub.csv", row.names = FALSE)


###
### Create a CSV file for MIMARKS

mimarks.amun <- meta.amun.complete %>%
  mutate(organism = "marine prokaryote",
         env_broad_scale = "marine biome [ENVO:00000447]",
         env_local_scale = ifelse(grepl("AMD", Station),                  # Check if "AMD" is in Station
                                  "stream [ENVO:00000023]",               # Value if "AMD" is present
                                  "marine pelagic zone [ENVO:00000208]"), # Default value
         env_medium = ifelse(grepl("AMD", Station),                       # Check if "AMD" is in Station
                             "river water [ENVO:01000599]",               # Value if "AMD" is present
                             "sea water [ENVO:00002149]"),                # Default value
         lat_lon = paste0(sprintf("%.5f", abs(`Latitude [deg]`)),         # Format latitude to 5 decimal places
                          ifelse(`Latitude [deg]` >= 0, " N ", " S "),    # Add "N" or "S" based on the sign
                          sprintf("%.5f", abs(`Longitude [deg]`)),        # Format longitude to 5 decimal places
                          ifelse(`Longitude [deg]` >= 0, " E", " W")),    # Add "E" or "W" based on the sign
         size_frac_low = case_when(seq.filt == 0 ~ prefilter.size,  # If seq.filt is 0, use prefilter.size
                                   seq.filt == 1 & pore.size == 0.22 ~ 3,             # If seq.filt is 1 and pore.size is 0.22, set to 3
                                   seq.filt == 1 & pore.size == 3 ~ 20,               # If seq.filt is 1 and pore.size is 3, set to 20
                                   seq.filt == 1 & pore.size == 20 ~ prefilter.size,  # If seq.filt = 1 & pore.size = 20, use prefilter.size
                                   TRUE ~ NA_real_),                                  # Default to NA if none of the above conditions are met
         samp_store_temp = "-80") %>%
  select(-"Latitude [deg]", -"Longitude [deg]", -"Year", -"Month", -"Day", -"Time") %>%
  select (all_of(c("genohub.id", "organism", "yyyy-mm-ddThh:mm:ss.sss", "Depth [m]", "Bot. Depth [m]", "env_broad_scale", "env_local_scale", "env_medium", "geo_loc_name", "lat_lon", "size_frac_low", "pore.size", "volume.l", "samp_store_temp", "Temperature [deg C]", "Salinity", "Dissolved Oxygen [uM]", "PAR [uE/s/m^2]", "Fluorescence [mg/m^3]", "Nitrate [uM]", "Nitrite [uM]", "Phosphate [uM]", "Silicate [uM]")), everything()) %>%
  list(tibble(parameter = names(.) %>%
                str_remove_all(" \\[.*\\]") %>%
                janitor::make_clean_names(),
              units = names(.) %>%
                str_extract("\\[.*\\]") %>%
                str_remove_all("[\\[\\]]")) %>%
         filter(!is.na(units))) %>%
  enframe(name = NULL) %>%
  unnest(value) %>%
  rename_with(~ str_remove(., " \\[.*\\]")) %>%   # Remove text inside brackets from column names
  clean_names()

write.csv(mimarks.amun, file = "OUTPUT/mimarks.amun.csv", row.names = FALSE)

