# 1A_META.R --------------------------------------------------------------------------------------
#
# Purpose:
# Match sequenced samples to field records, bottle hydrography and laboratory nutrients.
# Derive Conservative temperature, profile features and recent wave exposure; retain all downcasts for ODV.
#
# Inputs (paths are relative to the project working directory containing CTD, DNA, EXT and OUT):
# DNA/meta.amun.dna.csv; DNA/Amundsen_sampling_2023_Leg2_csv.csv; DNA/Amundsen_sampling_2023_Leg4_csv.csv
# CTD/BTL/**/*_btl.csv; CTD/LOG/**/RosetteSheet_###.xlsx
# EXT/AM2023_Leg2_CTD_profiles.csv; EXT/AM2023_Leg4_CTD_profiles.csv
# EXT/nutrients_2023_leg_2&4.csv; EXT/2023_August_October_SWH_hourly.nc
# EXT/ne_10m/ne_10m_geography_marine_polys.shp and ne_10m_coastline.shp (+ associated files)
#
# Outputs:
# EXT/meta.amun.complete.rds -> paired-marker selection in 2A_TAX.R
# OUT/1A_META/amun_data_overview.pdf (QC); AMUN_2023_all_downcast_profiles_ODV.csv
# OUT/1A_META/amun_CTD_profiles_ODV.txt; AMUN_MIMARKS.survey.water.6.0.csv
# Optional OUT/1A_META/amun_all_CTD_TS_pycnocline_chlmax.pdf (commented profile QC)
#
# Run sections in order within this script. Later scripts load saved inputs; they do not need earlier scripts in memory.
# Heavy calculations remain optional or cached. Existing EXT and OUT folders are assumed to exist.
#---------------------------------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(sf)
library(readxl)
library(gsw)
library(patchwork)

# Results / interpretation notes:
# Shallow pycnocline depths often coincided with strong near-surface salinity gradients, consistent with freshened surface stratification.
# Temperature and salinity alone cannot identify the freshwater source. Several Leg 4 profiles had near-surface fChl-a maxima.
#
# Scientific/data references:
# Amundsen CTD dataset: https://doi.org/10.5884/12713
# Amundsen 2023 expedition report: https://amundsenscience.com/wp-content/uploads/2023/11/Expedition-Report-2023_Final.pdf
# GSW/TEOS-10: https://teos-10.github.io/GSW-R/ ; Martin et al. (2010): https://doi.org/10.3354/meps08666
# N2-maximum approach for pycnocline: Tremblay et al. (2008), https://doi.org/10.1029/2007JC004547
# SWH from ERA5: Hersbach et al. (2020), https://doi.org/10.1002/qj.3803 & https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels
# Beam attenuation calculation: https://www.seabird.com/hubfs/Resources/Application%20Notes/AN91-CStarCalCoeffs.pdf
# Natural Earth coastline: https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-coastline/


############################################
### SAMPLE ROSTER AND FILTRATION RECORDS ###
# Genohub is the authoritative list of samples that were actually sequenced, so it defines which biological samples enter the analysis.
# The two field-sampling sheets provide bottle number, target depth, filtration details and sampling date that are not all present in Genohub.
# Keep only identifiers and DNA information needed downstream; plate.geno is upper-cased so sample IDs are consistent across projects.
amun_genohub <- read_csv("DNA/meta.amun.dna.csv", show_col_types = FALSE) %>%
  select(Station = station, `Sample label` = sample.id, `DNA box` = box, plate, genohub.id, plate.geno, dna.conc) %>%
  mutate(plate.geno = toupper(plate.geno),
         plate.geno = if_else(plate.geno %in% c("7_AM1-153200", "7_AM2-BGC-ARGO200"), paste0(plate.geno, "M"), plate.geno))

# Use the same column names for both cruise-leg sampling sheets so they can be stacked safely after their different layouts are harmonized.
sample.columns <- c("Cruise", "Station", "Sample label", "Sampling date", "Bottle", "Target", "Prefiltration mesh pore size [μm]",
                    "Sequential filtration", "Pore size [μm]", "Volume filtered [L]", "DNAavg [ng/μl]", "DNA box")

# Leg 2 used one 0.22-um fraction; the last two spreadsheet rows are non-sample summary rows and are deliberately discarded.
amun_sample_leg2 <- read_csv("DNA/Amundsen_sampling_2023_Leg2_csv.csv", show_col_types = FALSE) %>% dplyr::slice(1:(n() - 2)) %>%
  mutate(Cruise = "AM2023_02", `Pore size [μm]` = 0.22, `Sequential filtration` = 0L) %>% select(all_of(sample.columns))

# Leg 4 stored the 0.22-, 3- and 20-um fractions in separate columns, so pivot_longer converts them to one row per size fraction.
amun_sample_leg4_raw <- read_csv("DNA/Amundsen_sampling_2023_Leg4_csv.csv", show_col_types = FALSE) %>% dplyr::slice(1:(n() - 2))
amun_sample_leg4 <- amun_sample_leg4_raw %>%
  pivot_longer(matches("^(0\\.22|3|20) μm (DNA box|DNAavg \\[ng/μl\\])$"), names_to = c("Pore size [μm]", ".value"),
               names_pattern = "^(0\\.22|3|20) μm (DNA box|DNAavg \\[ng/μl\\])$") %>%
  mutate(Cruise = "AM2023_04", `Pore size [μm]` = as.numeric(`Pore size [μm]`)) %>% select(all_of(sample.columns))

# Combine field records for all fractions; the sequencing roster determines which samples are retained.
amun_sample <- bind_rows(amun_sample_leg2, amun_sample_leg4) %>%
  mutate(`Volume filtered [L]` = round(`Volume filtered [L]`, 1), Date = dmy(`Sampling date`), Year = year(Date), Month = month(Date))

# Start from the sequencing roster; unmatched field records are not added, but a listed sample can still have missing metadata.
dna.amun <- amun_genohub %>% left_join(amun_sample, by = c("Station", "Sample label", "DNA box"), na_matches = "never") %>%
  mutate(`Original station` = as.character(Station), Station = replace_values(Station, "D1\"" ~ "D1.2", "D1-2022" ~ "D1_2022", "SAL-01" ~ "SAL")) %>%
  select(Cruise, Station, `Original station`, Target, Bottle, `Sample label`, Year, Month, `Sequential filtration`,
         `Prefiltration mesh pore size [μm]`, `Pore size [μm]`, `Volume filtered [L]`, `DNAavg [ng/μl]`, plate, genohub.id, plate.geno)


############################################
### BOTTLE HYDROGRAPHY AND TARGET DEPTHS ###
# These CTD/BTL files are bottle/upcast records: each row describes the water measured when a particular Niskin bottle was closed.
# Rosette logs are read separately because they contain the human-entered target labels (surface, chl max, bottom, etc.) used during sampling.
csv_files <- list.files("CTD/BTL", pattern = "_btl\\.csv$", full.names = TRUE, recursive = TRUE)
log_files <- list.files("CTD/LOG", pattern = "^RosetteSheet_\\d{3}\\.xlsx$", full.names = TRUE, recursive = TRUE)
target_depth_list <- list()
meta_list <- list()

# Read every rosette log and recover bottle position, intended target label and approximate trip pressure recorded by the shipboard team.
# Numeric target labels such as "20" can also serve as approximate closure pressures when the dedicated trip-pressure field was left blank.
for (file in log_files) {
  cruise.number <- basename(dirname(file))
  cast.number <- as.integer(str_extract(basename(file), "\\d{3}(?=\\.xlsx$)"))
  target.depth <- read_excel(file, sheet = 1, range = "G6:K29", col_types = "text",
                             col_names = c("Bottle", "% light", "Target depth", "check", "Trip pressure [db]")) %>%
    transmute(Cruise_Number = cruise.number, Cast_Number = cast.number, BOPO = as.integer(Bottle),
              `Target depth` = na_if(str_squish(`Target depth`), ""),
              `Trip pressure [db]` = round(as.numeric(str_extract(str_squish(`Trip pressure [db]`), "^-?\\d+(?:\\.\\d+)?$")), -1),
              .target_pressure = round(as.numeric(str_extract(`Target depth`, "^-?\\d+(?:\\.\\d+)?$")), -1),
              .closure_pressure = coalesce(`Trip pressure [db]`, .target_pressure)) %>% filter(!is.na(BOPO))
  target_depth_list[[file]] <- target.depth
}

# Standardize synonymous target labels so "SCM", "subSCM", "btm", etc. do not become artificial categorical differences later.
target.depth.amun <- bind_rows(target_depth_list) %>%
  mutate(`Target depth` = str_to_lower(str_squish(`Target depth`)),
         `Target depth` = replace_values(`Target depth`, "bt" ~ "bottom", "btm" ~ "bottom", "scm" ~ "chl max",
                                         "scm or 20" ~ "chl max", "subscm" ~ "chl max", "subdcm" ~ "chl max (deep)",
                                         "sub dcm" ~ "chl max (deep)", "sub dcm_01" ~ "chl max (deep, 1st)",
                                         "sub dcm_02" ~ "chl max (deep, 2nd)", "thermo" ~ "thermocline"), `Target depth` = na_if(`Target depth`, "tbd"))

# A descriptive target written once at a given cruise/cast/pressure is propagated to every bottle closed at that same pressure.
target.labels <- target.depth.amun %>% filter(is.na(.target_pressure), !is.na(`Target depth`), !is.na(`Trip pressure [db]`)) %>%
  distinct(Cruise_Number, Cast_Number, .closure_pressure, `Target label` = `Target depth`)
target.depth.amun <- target.depth.amun %>% left_join(target.labels, by = c("Cruise_Number", "Cast_Number", ".closure_pressure"),
                                                               na_matches = "never", relationship = "many-to-one") %>%
  mutate(`Target depth` = coalesce(`Target label`, `Target depth`)) %>% select(-.target_pressure, -.closure_pressure, -`Target label`)

# Read each bottle CTD file, append its cast-level header metadata to each bottle row, then attach the normalized target labels.
for (file in csv_files) {
  basic_meta_raw <- read.csv(file, nrows = 9, header = FALSE, stringsAsFactors = FALSE)
  basic_meta_keys <- gsub("^% |:.*", "", basic_meta_raw$V1)
  basic_meta_values <- gsub("^.*: ", "", basic_meta_raw$V1)
  basic_meta <- setNames(as.data.frame(t(basic_meta_values), stringsAsFactors = FALSE), basic_meta_keys)
  meta <- read.csv(file, skip = 26, check.names = FALSE) %>% select(1:15) %>% bind_cols(basic_meta[rep(1, nrow(.)), ]) %>%
    mutate(Cast_Number = as.integer(Cast_Number), BOPO = as.integer(BOPO)) %>%
    left_join(target.depth.amun, by = c("Cruise_Number", "Cast_Number", "BOPO"))
  meta_list[[file]] <- meta
}

# Natural Earth marine polygons provide broad geographic names; the coastline gives straight-line station distance to land.
sp_crs <- 3347
marine_meta <- st_read("EXT/ne_10m/ne_10m_geography_marine_polys.shp", quiet = TRUE) %>% select(marine = name) %>% st_transform(sp_crs)
coast_meta <- st_read("EXT/ne_10m/ne_10m_coastline.shp", quiet = TRUE) %>% st_transform(sp_crs)

# Combine bottle rows, standardize names/units, repair known station-name variants and derive beam attenuation from transmission.
# Amundsen Science CTD-rosette instrumentation:
# Temperature, conductivity and pressure: Sea-Bird SBE-9plus CTD system.
# Dissolved oxygen: Sea-Bird SBE-43; fChl-a: Seapoint chlorophyll fluorometer; fCDOM: Wet Labs FL(RT)D.
# Optical nitrate: Sea-Bird SUNA (sn 1546); light transmission: Wet Labs C-Star transmissometer.
# PAR and surface PAR are also included in the Amundsen Science archive.
# Dataset/instrument metadata: https://erddap.amundsenscience.com/erddap/info/amundsen12713/index.html
#
# CTD variables used here: PRES = pressure; TE90 = in-situ temperature; PSAL = Practical Salinity; OXYM = dissolved oxygen;
# NTRA = optical nitrate; FLOR = fChl-a; CDOM = fCDOM; TRAN = light transmission.
ctd.amun <- bind_rows(meta_list) %>%
  mutate(Cruise = replace_values(Cruise_Name, "2023_02_AS" ~ "AM2023_02", "2023_04_AS" ~ "AM2023_04")) %>%
  dplyr::rename("Bottle" = "BOPO", "Sigma-theta [kg/m^3]" = "SIGP", "Pressure [db]" = "PRES",
                "Depth [m]" = "DEPH", "Temperature [deg C]" = "TE90", "Oxygen solubility [ml/L]" = "OSOL",
                "Salinity" = "PSAL", "Sigma-t [kg/m^3]" = "SIGT", "fChl-a [mg/m^3]" = "FLOR",
                "fCDOM [mg/m^3]" = "CDOM", "Seawater PAR [uE/s/m^2]" = "PSAR",
                "Light transmission [%]" = "TRAN", "Dissolved oxygen [uM]" = "OXYM", "Sound velocity [m/s]" = "SVEL",
                "Nitrate (CTD) [uM]" = "NTRA", "Cast" = "Cast_Number", "Latitude [deg]" = "Initial_Latitude [deg]",
                "Longitude [deg]" = "Initial_Longitude [deg]", "Bottom depth [m]" = "Sounding [m]") %>%
  
  # Beer-Lambert attenuation: c = -ln(T)/L, where transmission T is converted from percent to fraction and the C-Star path length L is 0.25 m.
  mutate(`Beam attenuation [1/m]` = -log(`Light transmission [%]` / 100) / 0.25) %>%
  filter(Station != "test") %>%
  mutate(Cast = as.integer(Cast), `Latitude [deg]` = as.numeric(`Latitude [deg]`),
         `Longitude [deg]` = as.numeric(`Longitude [deg]`), Year = as.integer(substring(`Start_Date_Time [UTC]`, 1, 4)),
         Month = as.integer(substring(`Start_Date_Time [UTC]`, 6, 7)), Day = as.integer(substring(`Start_Date_Time [UTC]`, 9, 10)),
         Time = substring(`Start_Date_Time [UTC]`, 12), `Bottom depth [m]` = as.numeric(`Bottom depth [m]`),
         `Bottom depth [m]` = if_else(`Bottom depth [m]` < `Depth [m]`, NA_real_, `Bottom depth [m]`)) %>%
  select(-Cruise_Name, -Original_Filename, -`Start_Date_Time [UTC]`) %>%
  mutate(Station = replace_values(Station, "D1''" ~ "D1.2", "A6/195" ~ "A6", "A7/196" ~ "A7",
                                  "A8/197" ~ "A8", "A9/198" ~ "A9", "Sal-01" ~ "SAL")) %>%
  # Project station coordinates before measuring distance; the nearest marine polygon supplies the broad geographic name.
  st_as_sf(coords = c("Longitude [deg]", "Latitude [deg]"), crs = 4326, remove = FALSE) %>% st_transform(sp_crs) %>%
  mutate(`Distance to coast [km]` = {
    i <- st_nearest_feature(geometry, coast_meta)
    round(as.numeric(st_distance(geometry, st_geometry(coast_meta)[i], by_element = TRUE)) / 1000)
  }, geo_loc_name = {
    i <- st_nearest_feature(geometry, marine_meta)
    d <- as.numeric(st_distance(geometry, st_geometry(marine_meta)[i], by_element = TRUE))
    if_else(d <= 5000, paste0("Canada: ", marine_meta$marine[i]), "Canada")
  }) %>% st_drop_geometry()


###########################################################
### Identify the CTD casts containing sequenced samples ###
# Several stations were occupied more than once.  The microbial sample must inherit hydrography from the exact cast/bottle that produced that sample.
# Matching Genohub samples to cruise + station + bottle + year + month recovers the correct CTD cast without assuming one cast per station.
# These known duplicate casts were not the casts from which the biological samples were collected and therefore cannot be used for DNA matching.
cast_duplicates <- tibble(Station = c("E5", "E3", "C5", "C3", "A5", "A3", "323", "163", "108", "D1", "350"),
                          Cast = c(8, 5, 10, 13, 22, 19, 64, 59, 46, 30, 41))

# Join sequenced marine samples to the bottle table only to recover cast number; river samples remain unmatched because no CTD was used.
dna_cast_match <- dna.amun %>%
  left_join(ctd.amun %>% anti_join(cast_duplicates, by = c("Station", "Cast")) %>%
              distinct(Cruise, Station, Bottle, Year, Month, Cast), by = c("Cruise", "Station", "Bottle", "Year", "Month"), na_matches = "never", relationship = "many-to-one")

# Collapse sample-level matches to unique cruise/station/cast combinations; only these casts are used for microbial profile features and plots.
sampled_casts <- dna_cast_match %>% filter(!is.na(Cast)) %>% distinct(Cruise, Station, Cast)


################################################################################################
### Geographic station codes: source IDs remain unchanged until final metadata/display joins ###
# Station groups A-M follow the counterclockwise sampling progression around Baffin Island and adjacent waterways.
# Leg 2 uses A-J; Leg 4 continues with K-M from where the geographic progression resumed. Rivers are R1-R2.
# Station numbers follow sampling order within each geographic group.

same <- c(paste0("A", 5:9), "B1", paste0("C", 1:5), "D1.2", paste0("D", 2:5), paste0("E", 1:5))
leg2_codes <- c(setNames(same, same), "D1" = "D1.1", "D1_2022" = "D1.3", "BGC-Argo" = "F1", "153" = "F2", "158" = "F3", "163" = "F4", "Bylot-23" = "F5",
                "116" = "G1", "115" = "G2", "113" = "G3", "111" = "G4", "108" = "G5", "105" = "G6", "103" = "G7", "101" = "G8", "100" = "G9",
                "324" = "H1", "323" = "H2", "300" = "H3", "S3" = "I1", "S5" = "I2", "S8" = "I3", "S10" = "I4", "305B" = "J1", "305C" = "J2", "305D" = "J3")
leg4_codes <- c("322-GB" = "K1", "324-GB" = "K2", "SPR-005" = "K3", "SPR-004" = "K4", "SPR-003" = "K5", "FH20-001" = "L1", "330" = "L2", "332" = "L3",
                "333" = "L4", "334" = "L5", "FoxSIPP-02" = "L6", "FoxSIPP-03" = "L7", "338" = "L8", "350" = "L9", "FoxSIPP-08" = "L10", "FoxSIPP-M" = "L11", 
                "FoxSIPP-10" = "L12", "349" = "L13", "15F" = "M1", "15E" = "M2", "15C" = "M3", "HS22-001" = "M4", "HS22-003" = "M5", "HS22-005" = "M6",
                "SAL" = "M7", "14E" = "M8", "14D" = "M9", "14B" = "M10", "HS22-009" = "M11", "HS22-013" = "M12", "354" = "M13", "356" = "M14", "640" = "M15")
river_codes <- c("AMD23-29" = "R1", "AMD23-6" = "R2")
station_map <- bind_rows(tibble(Cruise = "AM2023_02", `Station source` = names(leg2_codes), `Station code` = unname(leg2_codes)),
                         tibble(Cruise = "AM2023_04", `Station source` = names(leg4_codes), `Station code` = unname(leg4_codes)),
                         tibble(Cruise = "AM2023_04", `Station source` = names(river_codes), `Station code` = unname(river_codes))) %>%
  mutate(`Station group` = str_extract(`Station code`, "^[A-Z]+"))


##############################################################
### Complete CTD downcasts: retain ALL casts for later ODV ###
# These are the calibrated continuous DOWNCAST profiles from Amundsen Science (dataset doi:10.5884/12713), not the bottle/upcast values above.
# Downcasts are used for vertical structure because they profile continuously before bottle collection and provide full water-column context.
# Every valid Leg 2 and Leg 4 profile is retained here; the microbial subset is created only after the complete archive has been saved/exported.
ctd.profile.variables <- c("filename", "cast_number", "station", "time", "latitude", "longitude", "PRES", "depth", "TE90", "PSAL",
                           "OXYM", "NTRA", "FLOR", "CDOM", "TRAN", "TURB", "PSAR", "SPAR", "SIGT")
# Optional download/update of both complete calibrated cruise downcast files.
# ctd.profile.url <- function(leg) paste0("https://erddap.amundsenscience.com/erddap/tabledap/amundsen12713.csv0?",
#                                         paste(ctd.profile.variables, collapse = "%2C"), "&filename=~%22AS_CTD_2023_", leg, "_.*%22")
# download.file(ctd.profile.url("02"), "EXT/AM2023_Leg2_CTD_profiles.csv", mode = "wb", method = "libcurl")
# download.file(ctd.profile.url("04"), "EXT/AM2023_Leg4_CTD_profiles.csv", mode = "wb", method = "libcurl")
# Read both cruise legs, standardize a few known station aliases, remove negative/invalid sea pressures, and sort each profile from shallow to deep.
# TEOS-10 converts pressure to depth (.z), Practical to Absolute Salinity (.SA), in-situ to Conservative Temperature (.CT), and calculates sigma0.
# TEOS-10 helper variables are used only for hydrographic calculations; the original measured variables remain in the table.
ctd.profile.all <- bind_rows(
  read_csv("EXT/AM2023_Leg2_CTD_profiles.csv", col_names = ctd.profile.variables, show_col_types = FALSE) %>% mutate(Cruise = "AM2023_02"),
  read_csv("EXT/AM2023_Leg4_CTD_profiles.csv", col_names = ctd.profile.variables, show_col_types = FALSE) %>% mutate(Cruise = "AM2023_04")) %>%
  mutate(Cast = as.integer(cast_number), Station = str_squish(as.character(station)),
         Station = replace_values(Station, "A6/195" ~ "A6", "A7/196" ~ "A7", "A8/197" ~ "A8", "A9/198" ~ "A9", "Sal-01" ~ "SAL"),
         Station = if_else(Cruise == "AM2023_02" & Station == "D1" & Cast == 31, "D1.2", Station)) %>%
  filter(is.finite(PRES), PRES >= 0) %>%
  mutate(.z = -gsw_z_from_p(PRES, latitude), .SA = gsw_SA_from_SP(PSAL, PRES, longitude, latitude),
         .CT = gsw_CT_from_t(.SA, TE90, PRES), .sigma0 = gsw_sigma0(.SA, .CT)) %>%
  arrange(Cruise, Cast, PRES)
# The ODV export retains original station names and adds analysis codes where available; unsampled profiles remain fully available.
ctd.profile.odv <- ctd.profile.all %>% left_join(station_map, by = c("Cruise", "Station" = "Station source"))
write_csv(ctd.profile.odv, "OUT/1A_META/AMUN_2023_all_downcast_profiles_ODV.csv")
# Restrict only the HYDROGRAPHIC ANALYSIS object to sampled casts; this does not modify or discard anything from ctd.profile.all/ctd.profile.odv.
ctd.profile.amun <- ctd.profile.all %>% semi_join(sampled_casts, by = c("Cruise", "Station", "Cast"))


################################################################################
### Pycnocline and chlorophyll-fluorescence maximum depths for sampled casts ###
# N2 is calculated from TEOS-10 Absolute Salinity and Conservative Temperature without smoothing and is not bridged across >2-m gaps.
# Only the depths of maximum positive N2 and maximum raw fChl-a within the upper 100 m are retained; boundary maxima remain valid.
n2_fun <- function(d, max_gap = 2) {
  d <- d %>% distinct(PRES, .keep_all = TRUE) %>% arrange(PRES)
  empty <- tibble(depth = double(), N2 = double(), segment = integer())
  if (nrow(d) < 2) return(empty)
  ok <- is.finite(d$PRES) & is.finite(d$.z) & is.finite(d$.SA) & is.finite(d$.CT)
  d$.seg <- cumsum(!ok | c(TRUE, !head(ok, -1) | diff(d$.z) > max_gap))
  d <- d[ok, ]
  if (nrow(d) < 2) return(empty)
  raw <- rep(NA_real_, nrow(d) - 1L)
  for (i in split(seq_len(nrow(d)), d$.seg)) if (length(i) > 1)
    raw[head(i, -1)] <- gsw_Nsquared(d$.SA[i], d$.CT[i], d$PRES[i], median(d$latitude[i], na.rm = TRUE))$N2
  tibble(depth = (head(d$.z, -1) + tail(d$.z, -1)) / 2, N2 = raw, segment = head(d$.seg, -1))
}

pycnocline_fun <- function(d, max_depth = 100) {
  n <- n2_fun(d)
  i <- which(is.finite(n$N2) & n$N2 > 0 & n$depth <= max_depth)
  if (!length(i)) NA_real_ else n$depth[i[which.max(n$N2[i])]]
}

chlmax_fun <- function(d, max_depth = 100) {
  d <- d %>% distinct(PRES, .keep_all = TRUE) %>% filter(is.finite(.z), is.finite(FLOR), .z <= max_depth)
  if (!nrow(d)) NA_real_ else d$.z[which.max(d$FLOR)]
}

profile_features.amun <- ctd.profile.amun %>% group_by(Cruise, Station, Cast) %>%
  group_modify(~ tibble(`Pycnocline depth [m]` = pycnocline_fun(.x), `Chl-a maximum depth [m]` = chlmax_fun(.x))) %>% ungroup()


# ###################################################################################
# ### QC for pycnocline and chl max: hydrographic profiles for sampled casts only ###
# # Raw upper-100-m; temperature, salinity, N2 and fChl-a are shown for visual QC; red/green lines mark pycnocline/chl-max depths
# # Inspect once and comment out
# plot_hydro_cast <- function(cruise, station, cast, display) {
#   d <- ctd.profile.amun %>% filter(Cruise == .env$cruise, Station == .env$station, Cast == .env$cast) %>%
#     distinct(PRES, .keep_all = TRUE) %>% arrange(PRES)
#   ok_h <- is.finite(d$.z) & is.finite(d$.SA) & is.finite(d$.CT)
#   d$.seg <- cumsum(!ok_h | c(TRUE, !head(ok_h, -1) | diff(d$.z) > 2))
#   dh <- d[ok_h, ]
#   ok_c <- is.finite(d$.z) & is.finite(d$FLOR)
#   d$.chl_seg <- cumsum(!ok_c | c(TRUE, !head(ok_c, -1) | diff(d$.z) > 2))
#   n <- n2_fun(d)
#   m <- profile_features.amun %>% filter(Cruise == .env$cruise, Station == .env$station, Cast == .env$cast)
#   decorate <- function(p) p +
#     geom_hline(yintercept = m$`Pycnocline depth [m]`, colour = "red3", linetype = "dotted", linewidth = .65, na.rm = TRUE) +
#     geom_hline(yintercept = m$`Chl-a maximum depth [m]`, colour = "darkgreen", linetype = "dotdash", linewidth = .65, na.rm = TRUE) +
#     scale_y_reverse(limits = c(100, 0), breaks = seq(0, 100, 20), expand = expansion(mult = 0)) +
#     labs(y = "Depth [m]") + theme_minimal(base_size = 11)
#   panel <- function(var, lab) {
#     dd <- dh %>% filter(.z <= 100, is.finite(.data[[var]])) %>% group_by(.seg) %>% filter(n() > 1) %>% ungroup()
#     decorate(ggplot(dd, aes(x = .data[[var]], y = .z, group = .seg)) + geom_path(linewidth = .65) + labs(x = lab))
#   }
# 
#   nr <- n %>% filter(depth <= 100, is.finite(N2)) %>% group_by(segment) %>% filter(n() > 1) %>% ungroup()
#   cr <- d %>% filter(.z <= 100, is.finite(FLOR)) %>% group_by(.chl_seg) %>% filter(n() > 1) %>% ungroup()
#   pn <- ggplot(nr, aes(x = N2 * 1e4, y = depth, group = segment)) + geom_path(linewidth = .65) +
#     geom_vline(xintercept = 0, colour = "grey75") + labs(x = expression(N^2~"["*10^{-4}~s^{-2}*"]"))
#   pc <- ggplot(cr, aes(x = FLOR, y = .z, group = .chl_seg)) + geom_path(linewidth = .65) +
#     labs(x = expression("fChl-a [mg m"^{-3}*"]"))
#   label <- if_else(display == station, display, paste0(display, " [", station, "]"))
#   (panel("TE90", "Temperature [deg C]") | panel("PSAL", "Salinity") | decorate(pn) | decorate(pc)) +
#     plot_annotation(title = paste(cruise, label, cast, sep = " | "),
#                     subtitle = sprintf("Pycnocline = %.1f m | Chl-a maximum = %.1f m",
#                                        m$`Pycnocline depth [m]`, m$`Chl-a maximum depth [m]`))
# }
# 
# hydro_plot_order <- profile_features.amun %>% left_join(station_map, by = c("Cruise", "Station" = "Station source")) %>%
#   mutate(.display = coalesce(`Station code`, Station), .group = str_extract(.display, "^[A-Z]+"),
#          .number = as.numeric(str_extract(.display, "\\d+(?:\\.\\d+)?"))) %>%
#   arrange(match(.group, c(LETTERS, "R")), .number, Cruise, Cast)
# pdf("OUT/1A_META/amun_all_CTD_TS_pycnocline_chlmax.pdf", width = 15, height = 8.5)
# for (i in seq_len(nrow(hydro_plot_order)))
#   print(plot_hydro_cast(hydro_plot_order$Cruise[i], hydro_plot_order$Station[i],
#                         hydro_plot_order$Cast[i], hydro_plot_order$.display[i]))
# dev.off()


#############################
### NUTRIENT MEASUREMENTS ###
# Nutrient chemistry comes from a separate laboratory table; names and identifiers are harmonized to the CTD convention before joining.
# Known station/cast corrections are applied explicitly so a naming difference cannot silently create a false "missing nutrient" sample.
nut.amun <- read_csv("EXT/nutrients_2023_leg_2&4.csv", show_col_types = FALSE) %>%
  select(-contains("QF"), -contains("yyyy"), -contains("[degrees_"), -Type, -`Depth [m]`, -`Bot.Depth [m]`) %>%
  mutate(Cruise = replace_values(Cruise, "AN2023_02" ~ "AM2023_02", "AN2023_04" ~ "AM2023_04"),
         Cast = as.integer(Cast), Bottle = as.integer(Bottle),
         Station = replace_values(Station, "D1''" ~ "D1.2", "HSS22-003" ~ "HS22-003", "HSS22-005" ~ "HS22-005",
                                  "195" ~ "A6", "196" ~ "A7", "197" ~ "A8", "198" ~ "A9")) %>%
  mutate(Station = case_when(Cruise == "AM2023_04" & Cast == 37 ~ "FoxSIPP-08",
                             Cruise == "AM2023_04" & Cast == 57 ~ "HS22-006", TRUE ~ Station)) %>%
  rename_with(~ paste0(.x, " [uM]"), all_of(c("Nitrate", "Nitrite", "Phosphate", "Silicate")))


##############################################
### MATCH BOTTLE HYDROGRAPHY AND NUTRIENTS ###
# Attach the two profile-feature depths at cast level, then exact-bottle nutrients; same-target propagation is only used when chemistry is missing.
nutrient_vars <- c("Nitrate [uM]", "Nitrite [uM]", "Phosphate [uM]", "Silicate [uM]")
nutrient_sources <- paste0(nutrient_vars, " source bottle")
ctd.nut.amun <- ctd.amun %>% anti_join(cast_duplicates, by = c("Station", "Cast")) %>%
  semi_join(sampled_casts, by = c("Cruise", "Station", "Cast")) %>%
  left_join(profile_features.amun, by = c("Cruise", "Station", "Cast"), relationship = "many-to-one") %>%
  left_join(nut.amun, by = c("Cruise", "Station", "Cast", "Bottle"), relationship = "many-to-one") %>%
  # Fill each nutrient together with its source bottle, including partially missing chemistry.
  mutate(across(all_of(nutrient_vars), ~ if_else(!is.na(.x), Bottle, NA_integer_), .names = "{.col} source bottle"),
         .target.group = if_else(is.na(`Target depth`), paste0("Missing target depth_", row_number()), `Target depth`)) %>%
  group_by(Cruise, Station, Cast, .target.group) %>% fill(all_of(c(nutrient_vars, nutrient_sources)), .direction = "downup") %>%
  mutate(`Nut bottle` = na_if(apply(pick(all_of(nutrient_sources)), 1, \(x) paste(sort(unique(na.omit(x))), collapse = ", ")), ""),
         `Nut match` = case_when(if_all(all_of(nutrient_sources), is.na) ~ NA_character_,
                                if_all(all_of(nutrient_sources), ~ is.na(.x) | .x == Bottle) ~ "Same bottle",
                                if_any(all_of(nutrient_sources), ~ .x == Bottle) ~ "Same bottle + different bottle, same target depth",
                                TRUE ~ "Different bottle, same target depth")) %>%
  ungroup() %>% select(-.target.group) %>%
  select(Cruise, Station, Cast, Bottle, `Latitude [deg]`, `Longitude [deg]`, geo_loc_name, Year, Month, Day, Time,
         `Bottom depth [m]`, `Pycnocline depth [m]`, `Chl-a maximum depth [m]`, `Depth [m]`, `Pressure [db]`,
         `Temperature [deg C]`, Salinity, `Sigma-t [kg/m^3]`, `Sigma-theta [kg/m^3]`,
         `fChl-a [mg/m^3]`, `fCDOM [mg/m^3]`, `Beam attenuation [1/m]`,
         `Nitrate (CTD) [uM]`, `Nitrate [uM]`, `Nitrite [uM]`, `Phosphate [uM]`, `Silicate [uM]`, `Nut bottle`, `Nut match`, all_of(nutrient_sources),
         `Dissolved oxygen [uM]`, `Distance to coast [km]`)

# Downcast values at bottle pressures: medians within ±1 db of each bottle pressure in the same cast (downcast code = bottle column).
qc_vars <- c(TE90 = "Temperature [deg C]", PSAL = "Salinity", OXYM = "Dissolved oxygen [uM]", NTRA = "Nitrate (CTD) [uM]",
             FLOR = "fChl-a [mg/m^3]", CDOM = "fCDOM [mg/m^3]", TRAN = "Light transmission [%]")
downcast_at_bottles <- function(btl) btl %>% distinct(Cruise, Cast, Bottle, `Pressure [db]`) %>%
  mutate(p_lo = `Pressure [db]` - 1, p_hi = `Pressure [db]` + 1) %>%
  inner_join(ctd.profile.all %>% select(Cruise, Cast, PRES, all_of(names(qc_vars))), by = join_by(Cruise, Cast, between(y$PRES, x$p_lo, x$p_hi))) %>%
  group_by(Cruise, Cast, Bottle) %>% summarise(across(all_of(names(qc_vars)), ~ median(.x, na.rm = TRUE)), .groups = "drop")

# Quality check: bottle (upcast) values against downcast values for all bottles; PAR is excluded because light differs between descent and ascent.
# Leg 4 bottle fChl-a exceeds the processed downcasts by ~1.16 mg m-3 below 100 m, whereas the other variables agree in both legs.
bottle_downcast <- ctd.amun %>% left_join(downcast_at_bottles(ctd.amun), by = c("Cruise", "Cast", "Bottle"))
imap_dfr(qc_vars, \(btl_col, dc_col) bottle_downcast %>%
           transmute(Cruise, Variable = btl_col, Below_100 = `Pressure [db]` > 100, Upcast = .data[[btl_col]], Downcast = .data[[dc_col]])) %>%
  filter(is.finite(Upcast), is.finite(Downcast)) %>% group_by(Cruise, Variable) %>%
  summarise(Bottles = n(), `Median difference` = median(Upcast - Downcast), `Median difference >100 db` = median((Upcast - Downcast)[Below_100]),
            R2 = cor(Upcast, Downcast)^2, .groups = "drop") %>% mutate(across(where(is.double), ~ signif(.x, 3))) %>% print(n = Inf)

# Sample fChl-a is taken from the downcast; Leg 2 bottles above the first downcast record keep the bottle value (Leg 2 bottle and downcast agree).
# Bottle values are retained for reference.
ctd.nut.amun <- ctd.nut.amun %>% left_join(downcast_at_bottles(ctd.nut.amun) %>% select(Cruise, Cast, Bottle, FLOR),
                                           by = c("Cruise", "Cast", "Bottle"), relationship = "one-to-one") %>%
  mutate(`fChl-a bottle [mg/m^3]` = `fChl-a [mg/m^3]`,
         `fChl-a [mg/m^3]` = coalesce(FLOR, if_else(Cruise == "AM2023_02", `fChl-a bottle [mg/m^3]`, NA_real_))) %>% select(-FLOR)


######################################################################
### Attach DNA data to CTD/NUT data and add river-station metadata ###
# This is the point where environmental measurements are finally reduced from all bottles in sampled casts to the exact Genohub microbial samples.
# Marine samples inherit metadata by cruise + station + cast + bottle + year + month; rivers are patched separately because they lack CTD casts.
# River position, time and nominal 0.5-m depth come from the reviewed field sheet because rivers were not sampled with the CTD-rosette.
# The same coast-distance/geographic-name method is applied so river and marine samples share a compatible final metadata schema.
river_metadata <- amun_sample_leg4_raw %>% filter(str_detect(Station, "^AMD")) %>%
  mutate(Date = dmy(`Sampling date`), Day = day(Date)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>% st_transform(sp_crs) %>%
  mutate(`Distance to coast [km]` = {
    i <- st_nearest_feature(geometry, coast_meta)
    round(as.numeric(st_distance(geometry, st_geometry(coast_meta)[i], by_element = TRUE)) / 1000)
  }, geo_loc_name = {
    i <- st_nearest_feature(geometry, marine_meta)
    d <- as.numeric(st_distance(geometry, st_geometry(marine_meta)[i], by_element = TRUE))
    if_else(d <= 5000, paste0("Canada: ", marine_meta$marine[i]), "Canada")
  }) %>% st_drop_geometry() %>%
  transmute(Station, `Latitude [deg]` = Latitude, `Longitude [deg]` = Longitude, geo_loc_name, Day,
            Time = as.character(`Sampling time`), `Depth [m]` = 0.5, `Distance to coast [km]`)

# Final table contains only Genohub microbial samples; CTD/nutrient information comes from their corresponding sampled bottle and cast.
# Station renaming is deliberately done HERE, after all source-based joins are finished, so the study transect codes cannot affect data integration.

meta.amun.complete <- dna_cast_match %>%
  left_join(ctd.nut.amun, by = c("Cruise", "Station", "Cast", "Bottle", "Year", "Month"),
            na_matches = "never", relationship = "many-to-one") %>%
  mutate(`Latitude [deg]` = as.numeric(`Latitude [deg]`), `Longitude [deg]` = as.numeric(`Longitude [deg]`),
         Day = as.numeric(Day)) %>%
  rows_patch(river_metadata, by = "Station", unmatched = "ignore") %>%
  left_join(station_map, by = c("Cruise", "Station" = "Station source")) %>%
  mutate(`Station source` = Station, Station = coalesce(`Station code`, Station),
         Sampling = if_else(Station %in% c("R1", "R2"), "PUMP", "ROS-B")) %>% select(-`Station code`) %>%
  # DNA concentration is numerically identical in ng/ul and mg/L; the renamed unit is used to match the cross-project metadata convention.
  dplyr::rename(`DNA [mg/L]` = `DNAavg [ng/μl]`) %>%
  separate(Time, into = c("Hour", "Minute", "Second"), sep = ":", remove = FALSE, convert = TRUE,
           extra = "merge", fill = "right") %>%
  mutate(across(c(Cruise, Station, `Original station`, Sampling, Target, `Sample label`, Time, plate, genohub.id, plate.geno), as.character),
         across(c(Cast, Bottle, Year, Month, Day, Hour, Minute, Second, `Sequential filtration`), as.integer),
         Second = replace_na(Second, 0L),
         Time = if_else(is.na(Hour) | is.na(Minute), NA_character_, sprintf("%02d:%02d:%02d", Hour, Minute, Second)),
         # DOY is fractional day-of-year, preserving time-of-day for analyses that need cruise progression rather than calendar date alone.
         DOY = yday(make_date(Year, Month, Day)) + (Hour * 3600 + Minute * 60 + Second) / 86400,
         # Recalculate a consistent TEOS-10 thermodynamic set at the actual bottle depth; AOU = oxygen saturation concentration minus measured O2.
         .SA = gsw_SA_from_SP(Salinity, `Pressure [db]`, `Longitude [deg]`, `Latitude [deg]`),
         `Conservative temperature [deg C]` = gsw_CT_from_t(.SA, `Temperature [deg C]`, `Pressure [db]`),
         .CT = `Conservative temperature [deg C]`, `Density [kg/m^3]` = gsw_rho(.SA, .CT, `Pressure [db]`),
         `Sigma-t [kg/m^3]` = coalesce(`Sigma-t [kg/m^3]`, gsw_rho(.SA, gsw_CT_from_t(.SA, `Temperature [deg C]`, 0), 0) - 1000),
         `Sigma-theta [kg/m^3]` = coalesce(`Sigma-theta [kg/m^3]`, gsw_sigma0(.SA, .CT)),
         # AOU (apparent oxygen utilization) is oxygen solubility minus measured oxygen after converting uM to umol/kg using in-situ density.
         `AOU [umol/kg]` = gsw_O2sol(.SA, .CT, `Pressure [db]`, `Longitude [deg]`, `Latitude [deg]`) -
           `Dissolved oxygen [uM]` * 1000 / `Density [kg/m^3]`) %>%
  select(any_of(c("Cruise", "Station", "Station group", "Station source", "Original station", "Sampling", "Cast", "Latitude [deg]", "Longitude [deg]",
                  "geo_loc_name", "Bottom depth [m]", "Target", "Bottle", "Sample label", "Year", "Month", "Day",
                  "Hour", "Minute", "Second", "Time", "DOY", "Depth [m]", "Pressure [db]", "Conservative temperature [deg C]",
                  "Temperature [deg C]", "Salinity", "Bottle salinity", "Density [kg/m^3]", "Sigma-t [kg/m^3]",
                  "Sigma-theta [kg/m^3]", "fChl-a [mg/m^3]", "fChl-a bottle [mg/m^3]", "fCDOM [mg/m^3]",
                  "Beam attenuation [1/m]", "Nitrate (CTD) [uM]", "Nitrate [uM]", "Nitrite [uM]", "Phosphate [uM]", "Silicate [uM]",
                  "Dissolved oxygen [uM]", "AOU [umol/kg]", "Nut bottle", "Nut match",
                  "Pycnocline depth [m]", "Chl-a maximum depth [m]",
                  nutrient_sources, "Distance to coast [km]", "Sequential filtration", "Prefiltration mesh pore size [μm]", "Pore size [μm]",
                  "Volume filtered [L]", "DNA [mg/L]", "plate", "genohub.id", "plate.geno")))


#############################################################################################################
### Recent wave exposure: maximum hourly ERA5 SWH during the 48 h preceding the recorded marine cast time ###
# The CTD header supplies cast start time; all bottles in that cast share this exposure. SWH is not a direct measurement of mixing.
# Bilinear interpolation is spatial only; hourly values are retained unchanged in the time window (cast time - 48 h, cast time].

# # One-time ERA5 preprocessing: convert the downloaded GRIB to an hourly SWH NetCDF used below
# grib_file <- "EXT/2023_August_October_waves.grib"; swh_file <- "EXT/2023_August_October_SWH.nc"
# if (!file.exists(swh_file)) {
#   swh_array <- unique(na.omit(str_extract(as.character(unlist(terra::ar_info(grib_file, what = "arrays"))), "SWH_[A-Za-z0-9_-]+")))
#   array_hour <- as.integer(str_match(swh_array, "_(\\d+)$")[, 2]); array_hour[is.na(array_hour)] <- 1L; swh_array <- swh_array[order(array_hour)]
#   sf::gdal_utils("mdimtranslate", grib_file, swh_file, options = c("-of", "NetCDF", as.vector(rbind("-array", swh_array))), quiet = FALSE)
# }
# swh0 <- terra::rast(terra::sds(swh_file, md = TRUE))
# terra::time(swh0) <- seq(as.POSIXct("2023-08-01 00:00:00", tz = "UTC"), by = "hour", length.out = terra::nlyr(swh0))
# names(swh0) <- paste0("SWH_", format(terra::time(swh0), "%Y%m%d_%H%M"))
# terra::writeCDF(swh0, "EXT/2023_August_October_SWH_hourly.nc", varname = "SWH", longname = "Significant wave height", unit = "m",
#                 split = FALSE, overwrite = TRUE, prec = "float", compression = 4)

swh <- terra::rast("EXT/2023_August_October_SWH_hourly.nc")
terra::time(swh) <- as.POSIXct(terra::time(swh), tz = "UTC")
names(swh) <- paste0("SWH_", format(terra::time(swh), "%Y%m%d_%H%M"))
wave_sites <- meta.amun.complete %>% filter(!is.na(Cast), is.finite(DOY), is.finite(`Longitude [deg]`), is.finite(`Latitude [deg]`)) %>%
  distinct(Cruise, Station, Cast, Year, DOY, `Longitude [deg]`, `Latitude [deg]`) %>%
  transmute(Cruise, Station, Cast, Sample_time = as.POSIXct(paste0(Year, "-01-01"), tz = "UTC") + (DOY - 1) * 86400,
            Longitude = `Longitude [deg]`, Latitude = `Latitude [deg]`, ID = row_number())
wave_points <- terra::vect(wave_sites, geom = c("Longitude", "Latitude"), crs = "EPSG:4326") %>% terra::project(terra::crs(swh))
sample_wave <- terra::extract(swh, wave_points, method = "bilinear") %>% as_tibble() %>%
  pivot_longer(-ID, names_to = "Layer", values_to = "SWH") %>% left_join(wave_sites, by = "ID") %>%
  left_join(tibble(Layer = names(swh), Wave_time = as.POSIXct(terra::time(swh), tz = "UTC")), by = "Layer")
sample_swh <- sample_wave %>% filter(Wave_time > Sample_time - hours(48), Wave_time <= Sample_time) %>%
  group_by(Cruise, Station, Cast) %>% summarise(`SWHmax [m]` = if (any(is.finite(SWH))) max(SWH[is.finite(SWH)]) else NA_real_, .groups = "drop")
meta.amun.complete <- meta.amun.complete %>% select(-any_of("SWHmax [m]")) %>%
  left_join(sample_swh, by = c("Cruise", "Station", "Cast"))
# Save metadata for 2A_TAX.R, which selects paired-marker samples and adds read totals.
# Diagnostic plots and submission exports below do not alter this saved analysis input.
saveRDS(meta.amun.complete, "EXT/meta.amun.complete.rds")


###############################################################################################
### Overview PDF: station map, nitrate/beam QC, environmental correlations, metadata ranges ###
stations <- meta.amun.complete %>% filter(is.finite(`Longitude [deg]`), is.finite(`Latitude [deg]`)) %>%
  distinct(Station, geo_loc_name, `Longitude [deg]`, `Latitude [deg]`) %>%
  st_as_sf(coords = c("Longitude [deg]", "Latitude [deg]"), crs = 4326) %>% st_transform(sp_crs)
station_xy <- bind_cols(st_drop_geometry(stations), as_tibble(st_coordinates(stations))) %>%
  mutate(.group = str_extract(Station, "^[A-Z]+"), .number = as.numeric(str_extract(Station, "\\d+(?:\\.\\d+)?"))) %>% arrange(.group, .number)
map_xlim <- range(station_xy$X) + diff(range(station_xy$X)) * c(-.02, .02); map_ylim <- range(station_xy$Y) + diff(range(station_xy$Y)) * c(-.02, .02)
sampling_map <- ggplot() +
  geom_sf(data = marine_meta, fill = "aliceblue", colour = "grey85", linewidth = .15) +
  geom_sf(data = coast_meta, colour = "grey45", linewidth = .25) +
  geom_point(data = station_xy, aes(X, Y, fill = geo_loc_name), shape = 21, size = 2.6, colour = "black", stroke = .35) +
  ggrepel::geom_label_repel(data = filter(station_xy, .group != "G"), aes(X, Y, label = Station), direction = "y", nudge_x = 6000,
                            hjust = 0, size = 2, fontface = "bold", fill = scales::alpha("white", .9), box.padding = .12,
                            point.padding = .08, force = .5, force_pull = 10, max.overlaps = Inf, segment.color = NA,
                            label.padding = unit(.1, "lines"), label.r = unit(.08, "lines"), seed = 1, max.time = Inf) +
  geom_label(data = filter(station_xy, .group == "G"), aes(X, Y, label = Station), nudge_x = 5000, nudge_y = -15000, hjust = 0,
             size = 2, fontface = "bold", fill = scales::alpha("white", .9), label.padding = unit(.1, "lines"), label.r = unit(.08, "lines")) +
  scale_fill_brewer(name = NULL, palette = "Dark2") + guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  coord_sf(crs = st_crs(stations), default_crs = st_crs(stations), xlim = map_xlim, ylim = map_ylim, expand = FALSE, clip = "on") +
  labs(x = NULL, y = NULL) + theme_minimal(base_size = 10) +
  theme(panel.grid.major = element_line(colour = "grey90", linewidth = .2), legend.position = "bottom",
        legend.margin = margin(t = 10), legend.box.margin = margin(t = 5), plot.margin = margin(12, 12, 8, 12))


### ### Raw environmental variables for overview correlations
check_vars <- c("Depth [m]", "Conservative temperature [deg C]", "Salinity", "fChl-a [mg/m^3]", "fCDOM [mg/m^3]",
                "Beam attenuation [1/m]", "Nitrate [uM]", "Nitrite [uM]", "Phosphate [uM]", "Silicate [uM]",
                "Distance to coast [km]", "Latitude [deg]", "SWHmax [m]", "Pycnocline depth [m]",
                "Chl-a maximum depth [m]", "Volume filtered [L]", "DNA [mg/L]", "Nitrate (CTD) [uM]")
cor_data <- meta.amun.complete %>% select(Cruise, all_of(check_vars))


### ### Seven relationships of direct interest
cor_pairs <- tribble(~Relationship, ~X, ~Y,
                     "Conservative temperature vs salinity", "Conservative temperature [deg C]", "Salinity",
                     "fChl-a vs beam attenuation", "fChl-a [mg/m^3]", "Beam attenuation [1/m]",
                     "Nitrate vs phosphate", "Nitrate [uM]", "Phosphate [uM]",
                     "Distance to coast vs max. SWH", "Distance to coast [km]", "SWHmax [m]",
                     "Pycnocline vs chl-a maximum depth", "Pycnocline depth [m]", "Chl-a maximum depth [m]",
                     "Volume filtered vs DNA", "Volume filtered [L]", "DNA [mg/L]")

make_cor_plot <- \(Relationship, X, Y) {
  d <- cor_data %>% transmute(Cruise, x = .data[[X]], y = .data[[Y]]) %>% filter(is.finite(x), is.finite(y))
  s <- d %>% group_by(Cruise) %>% summarise(n = n(), r = cor(x, y), .groups = "drop")
  ggplot(d, aes(x, y)) + geom_point(size = 1.4, alpha = .6) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = .6) + facet_wrap(~Cruise) +
    geom_text(data = s, aes(x = -Inf, y = Inf, label = sprintf("n = %d; r = %.2f", n, r)),
              inherit.aes = FALSE, hjust = -.05, vjust = 1.15, size = 3) +
    labs(title = Relationship, x = X, y = Y) + theme_bw(base_size = 10)
}

pair_plots <- pmap(cor_pairs, make_cor_plot)


### ### Exact-bottle laboratory vs CTD nitrate graph
nitrate_check <- nut.amun %>% select(Cruise, Station, Cast, Bottle, `Nitrate [uM]`) %>%
  inner_join(ctd.amun %>% select(Cruise, Station, Cast, Bottle, `Nitrate (CTD) [uM]`), by = c("Cruise", "Station", "Cast", "Bottle")) %>%
  filter(is.finite(`Nitrate [uM]`), is.finite(`Nitrate (CTD) [uM]`)) %>% distinct(Cruise, Station, Cast, Bottle, .keep_all = TRUE)

# The pooled R2 and complete-pair count support the Methods nitrate QC statement; per-leg r values remain on the plots.
nitrate_qc_overall <- nitrate_check %>% summarise(n = n(), R2 = cor(`Nitrate [uM]`, `Nitrate (CTD) [uM]`)^2)
print(nitrate_qc_overall)
nitrate_qc <- nitrate_check %>% group_by(Cruise) %>% summarise(n = n(), r = cor(`Nitrate [uM]`, `Nitrate (CTD) [uM]`), .groups = "drop")

nitrate_plot <- ggplot(nitrate_check, aes(`Nitrate [uM]`, `Nitrate (CTD) [uM]`)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") + geom_point(size = 1.4, alpha = .6) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = .6) + facet_wrap(~Cruise) +
  geom_text(data = nitrate_qc, aes(x = -Inf, y = Inf, label = sprintf("n = %d; r = %.2f", n, r)),
            inherit.aes = FALSE, hjust = -.05, vjust = 1.15, size = 3) +
  labs(title = "Laboratory vs CTD nitrate", x = "Laboratory nitrate [uM]", y = "CTD nitrate [uM]") + theme_bw(base_size = 10)

pair_plots[[7]] <- nitrate_plot


### ### Correlation matrices: ALL variables, separately by cruise, lower triangle only
make_cor_matrix <- \(cruise) {
  d <- cor_data %>% filter(Cruise == cruise) %>% select(-Cruise)
  r <- cor(d, use = "pairwise.complete.obs"); vars <- colnames(r)
  z <- as.data.frame(as.table(r)) %>% as_tibble() %>%
    transmute(X = as.character(Var1), Y = as.character(Var2), r = Freq, i = match(X, vars), j = match(Y, vars)) %>%
    filter(i > j) %>% mutate(X = factor(X, levels = vars), Y = factor(Y, levels = rev(vars)),
                             fill = if_else(abs(r) >= .5, r, NA_real_))
  
  ggplot(z, aes(X, Y, fill = fill)) + geom_tile(colour = "grey90", linewidth = .25) +
    geom_text(aes(label = if_else(abs(r) >= .5, sprintf("%.2f", r), "")), size = 2.7) +
    scale_fill_gradientn(colours = c("blue4", "white", "white", "red3"), values = scales::rescale(c(-1, -.5, .5, 1)),
                         limits = c(-1, 1), breaks = c(-1, -.5, .5, 1), na.value = "white", name = "Pearson r") +
    labs(title = cruise, x = NULL, y = NULL, subtitle = "|r| < 0.5 shown in white") + coord_equal() + theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 55, hjust = 1), plot.title = element_text(face = "bold"))
}


### ### Save overview PDF
pdf("OUT/1A_META/amun_data_overview.pdf", width = 9, height = 9)
print(sampling_map)
print(wrap_plots(pair_plots[1:4], ncol = 2))
print(wrap_plots(pair_plots[5:7], ncol = 2))
for (cr in unique(cor_data$Cruise)) print(make_cor_matrix(cr))

par(mfrow = c(3, 2), mar = c(5, 6, 4, 3))
for (col_name in names(meta.amun.complete)[sapply(meta.amun.complete, is.numeric)]) {
  x <- meta.amun.complete[[col_name]]; x <- x[is.finite(x)]; if (!length(x)) next
  hist(x, main = NULL, xlab = col_name, col = "lightblue", border = "black")
  mtext(paste("Min:", round(min(x), 2)), side = 3, line = 2, adj = 0, cex = .8)
  mtext(paste("Max:", round(max(x), 2)), side = 3, line = 1, adj = 0, cex = .8)
  mtext(paste("Mean +/- SD:", round(mean(x), 2), "+/-", round(sd(x), 2)), side = 3, line = 0, adj = 0, cex = .8)
}
dev.off()


###############################
### NCBI BIOSAMPLE METADATA ###
# Translate the internal metadata table into NCBI MIMARKS survey-water fields for direct BioSample submission.
# River samples receive freshwater ENVO terms, while CTD-rosette samples receive marine ENVO terms appropriate to the sampled environment.
bioproject_accession <- "PRJNA421293"
fmt_num <- function(x, digits = 3) if_else(is.na(x), NA_character_, sub("\\.?0+$", "", sprintf(paste0("%.", digits, "f"), x)))
# Derive the upper bound of each microbial size fraction and construct descriptive NCBI sample titles and collection timestamps.
amun <- meta.amun.complete %>%
  mutate(.river_sample = Station %in% c("R1", "R2"),
         .size_frac_up = case_when(`Sequential filtration` == 0 ~ `Prefiltration mesh pore size [μm]`,
                                   `Sequential filtration` == 1 & `Pore size [μm]` == 0.22 ~ 3,
                                   `Sequential filtration` == 1 & `Pore size [μm]` == 3 ~ 20,
                                   `Sequential filtration` == 1 & `Pore size [μm]` == 20 ~ `Prefiltration mesh pore size [μm]`,
                                   TRUE ~ NA_real_),
         .pore_id = fmt_num(`Pore size [μm]`),
         sample_title = if_else(is.na(Cruise) | is.na(Station) | is.na(Target) | is.na(.pore_id), NA_character_,
                                sprintf("Canadian Arctic microbial communities: %s, station %s, %s, %s um filter",
                                        Cruise, Station, Target, .pore_id)),
         collection_date = case_when(!.river_sample & !is.na(Year) & !is.na(Month) & !is.na(Day) & !is.na(Time) ~
                                       sprintf("%04d-%02d-%02dT%sZ", Year, Month, Day, Time),
                                     !is.na(Year) & !is.na(Month) & !is.na(Day) ~
                                       sprintf("%04d-%02d-%02d", Year, Month, Day), TRUE ~ NA_character_))
# Build the MIMARKS table; the parameter-description field records units for NCBI columns that otherwise contain unitless text.
mimarks.amun <- amun %>%
  transmute(`*sample_name` = genohub.id, sample_title, bioproject_accession = .env$bioproject_accession,
            source_material_id = `Sample label`, `*organism` = if_else(.river_sample, "freshwater metagenome", "marine metagenome"),
            `*collection_date` = collection_date, `*depth` = `Depth [m]`,
            `*env_broad_scale` = if_else(.river_sample, "freshwater river biome [ENVO:01000253]", "marine biome [ENVO:00000447]"),
            `*env_local_scale` = if_else(.river_sample, "stream [ENVO:00000023]", "marine pelagic zone [ENVO:00000208]"),
            `*env_medium` = if_else(.river_sample, "river water [ENVO:01000599]", "sea water [ENVO:00002149]"),
            `*geo_loc_name` = geo_loc_name,
            `*lat_lon` = if_else(is.na(`Latitude [deg]`) | is.na(`Longitude [deg]`), NA_character_,
                                 sprintf("%.5f %s %.5f %s", abs(`Latitude [deg]`),
                                         if_else(`Latitude [deg]` >= 0, "N", "S"), abs(`Longitude [deg]`),
                                         if_else(`Longitude [deg]` >= 0, "E", "W"))),
            tot_depth_water_col = `Bottom depth [m]`,
            size_frac_low = `Pore size [μm]`, size_frac_up = .size_frac_up, samp_size = `Volume filtered [L]`,
            samp_store_temp = -80, temp = `Temperature [deg C]`, salinity = Salinity,
            diss_oxygen = `Dissolved oxygen [uM]`, fluor = `fChl-a [mg/m^3]`,
            cdom = `fCDOM [mg/m^3]`, beam_attenuation = `Beam attenuation [1/m]`,
            nitrate = `Nitrate [uM]`, nitrite = `Nitrite [uM]`, phosphate = `Phosphate [uM]`,
            silicate = `Silicate [uM]`, cruise = Cruise, station = Station, station_group = `Station group`, cast = Cast, bottle = Bottle,
            target = Target, dna_conc = `DNA [mg/L]`, sigma_t = `Sigma-t [kg/m^3]`,
            sigma_theta = `Sigma-theta [kg/m^3]`, nitrate_ctd = `Nitrate (CTD) [uM]`,
            `parameter (units)` = paste("collection_date (UTC for marine samples; date only for river samples);",
                                        "depth and tot_depth_water_col (m); size_frac_low and size_frac_up (um);",
                                        "samp_size (L); samp_store_temp and temp (degrees C); salinity (unitless);",
                                        "diss_oxygen, nitrate, nitrite, phosphate, silicate and nitrate_ctd (umol L-1);",
                                        "fluor and cdom (mg m-3); beam_attenuation (m-1); dna_conc (mg L-1);",
                                        "sigma_t and sigma_theta (kg m-3)"))
write.csv(mimarks.amun, "OUT/1A_META/AMUN_MIMARKS.survey.water.6.0.csv", row.names = FALSE, na = "NA", fileEncoding = "UTF-8")


###########################################
### ODV EXPORT OF ALL DOWNCAST PROFILES ###
odv.profile <- ctd.profile.odv %>%
  mutate(Cruise_name = Cruise, Station_name = station, Cruise = as.integer(str_extract(Cruise, "\\d+$")),
         Station = Cast, Type = "C", `yyyy-mm-ddThh:mm:ss` = format(time, "%Y-%m-%dT%H:%M:%S", tz = "UTC"),
         `Longitude [degrees_east]` = longitude, `Latitude [degrees_north]` = latitude) %>%
  transmute(Cruise, Station, Type, `yyyy-mm-ddThh:mm:ss`, `Longitude [degrees_east]`, `Latitude [degrees_north]`,
            Cruise_name, Station_name, `Station code`, `Station group`, Cast,
            `Depth [m]` = depth, `Pressure [db]` = PRES, `Temperature [deg C]` = TE90, Salinity = PSAL,
            `Absolute salinity [g/kg]` = .SA, `Conservative temperature [deg C]` = .CT,
            `Sigma-theta [kg/m^3]` = .sigma0, `Dissolved oxygen [uM]` = OXYM,
            `Nitrate (CTD) [uM]` = NTRA, `fChl-a [mg/m^3]` = FLOR,
            `fCDOM [mg/m^3]` = CDOM, `Transmission [%]` = TRAN, `Turbidity` = TURB, PSAR, SPAR)
write.table(odv.profile, "OUT/1A_META/amun_CTD_profiles_ODV.txt", sep = "\t", quote = FALSE, row.names = FALSE, na = "")
