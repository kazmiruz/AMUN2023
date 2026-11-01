# Amundsen 2023 analysis pipeline

This repository contains seven R scripts for the Amundsen 2023 microbial community study in the eastern Canadian Arctic.

The CTD, DNA, EXT, and OUT folders contain the input data, supporting files, and analysis outputs. They are available as a 7z archive on [CanWIN](https://canwin-datahub.ad.umanitoba.ca/data/en/project/micb-amun2023). The complete collection on [Google Drive](https://drive.google.com/drive/folders/1ki3yW1-yE-W7HCAZRh4MeBzsJ4weOFR1) includes everything in the archive plus optional PR2 reference files, raw BLAST results, and higher-resolution bathymetry. Extract the four folders into the same project directory as the seven R scripts, keeping all folder names and internal structures unchanged.

## Scripts and folders

| Script | Purpose |
| --- | --- |
| `1A_META.R` | Match DNA samples to field records, CTD bottles and nutrients; assign analysis station codes while retaining source names; derive Conservative temperature, pycnocline and chlorophyll-maximum depths, SWHmax, coast distance, beam attenuation and apparent oxygen utilization. Export metadata, ODV profiles, QC plots and an NCBI MIMARKS table. |
| `2A_TAX.R` | Clean SILVA-based 16S taxonomy, load cached PR2-based 18S assignments, retain paired-marker samples, and save original counts, taxonomy, and original/non-rare relative-abundance tables. |
| `3A_CLUS.R` | Cluster samples and ASVs; retain seven broad sample clusters and one supported subcluster; identify cluster and cluster-pair bioindicators; draw heatmaps and dendrograms. |
| `4A_PLANKT.R` | Select rarefaction depths, calculate Shannon diversity, compare diversity among clusters, and summarize major taxonomic groups. |
| `5A_INDIC.R` | Plot prevalent eukaryotic and prokaryotic bioindicators with cluster associations and sample read counts. |
| `6A_ENV.R` | Draw station/cluster maps and environmental plots; compare environmental properties; test Shannon-environment correlations; export supplementary tables and taxonomic summaries. Save shared inputs for 7A. |
| `7A_dbRDA.R` | For CL3-CL6, fit environmental dbRDAs for eukaryotic and prokaryotic composition; select prevalent eukaryotic ASVs in a partial dbRDA conditioned on the selected prokaryotic environmental predictors; partition variation in prokaryotic composition between these predictor sets. |

- `CTD/`: bottle records and rosette sampling logs.
- `DNA/`: sequencing-derived ASV counts and classifier results, plus DNA sampling metadata. These scripts start from processed ASV tables, not raw FASTQ reads.
- `EXT/`: supporting data, reference/spatial files, and intermediate or cached RDS files. Several scripts update this folder.
- `OUT/`: figures, tables and submission files, organized by script. Download the complete folder with its existing subfolders; the scripts do not create output directories.

## Running the scripts

Set the working directory to the project directory containing the seven scripts and the four folders. Run **1A -> 2A -> 3A -> 4A -> 5A -> 6A -> 7A**. Each script can also run in a fresh R session when the files listed in its header are available. Existing outputs with the same filenames are overwritten.

The recorded analysis versions are R 4.6.1 and vegan 2.7-5. Use these and the package versions listed below when reproducing the deposited results. **dplyr >=1.2.0 is required** for `replace_values()` and `recode_values()`.

Required packages are `tidyverse`, `lubridate`, `sf`, `readxl`, `gsw`, `terra`, `patchwork`, `ggrepel`, `gplots`, `ComplexHeatmap`, `viridis`, `cluster`, `dendextend`, `vegan`, `indicspecies`, `permute`, `car`, `ggdendro`, `ggspatial`, `marmap`, `dunn.test`, `openxlsx`, `ggrastr`, `ggnewscale` and `svglite`. `rlang` and `scales` are installed with tidyverse dependencies. Install `ComplexHeatmap` through Bioconductor. Optional BLAST regeneration also uses Bioconductor's `Biostrings` and external BLAST+ executables. The PATH example in 2A is Windows-specific.

The following regeneration blocks are commented out; enable them when their saved inputs/results need rebuilding:

| Script | Optional block |
| --- | --- |
| 1A | Download calibrated CTD profiles from ERDDAP; convert ERA5 GRIB to hourly SWH NetCDF; draw the per-cast hydrographic QC PDF. |
| 2A | Rebuild confidence-labelled taxonomy from `taxout.rds` and `bootout.rds`; build the PR2 BLAST database and regenerate `hits_pr2.rds`; filter and annotate the hits to save `EXT/hits_pr2_filt.rds`. |
| 3A | Recompute the 9,999-permutation bioindicator analysis and save `EXT/indicator_res.rds`. |
| 4A | Reassess rarefaction depths, average Shannon diversity across 1,000 random rarefactions, and save diversity/depth-assessment RDS files. |
| 6A | Download bathymetry if the selected resolution is not already cached. `bathy_res <- 4` is used for fast bathymetry generation, but `bathy_res <- 0.5` is used for the manuscript. |

**Cached results:** 2A loads filtered and annotated PR2 results from `EXT/hits_pr2_filt.rds`. Normal runs do not require `EXT/PR2/`, `EXT/hits_pr2.rds` or `EXT/tax.18S.rds`. The deposited indicator and Shannon caches are ready to use. 3A checks abundance values and cluster assignments; 4A checks original counts. A mismatch or missing input fingerprint stops the script with instructions to rebuild the cache using its commented block. Regenerate explicitly if analysis settings change, even when inputs do not.

2A removes organelle/eukaryote sequences before treating low-confidence 16S assignments as unresolved. Rebuild the PR2 results and filtered cache if the 18S query sequences, PR2 reference data or BLAST/filtering settings change. If input data or filtering settings change, rerun affected scripts in order and regenerate the 3A/4A caches when their inputs change.

7A loads `EXT/amun_environment_setup.rds` from 6A and `EXT/proptab.amun.rds` from 2A, then recomputes model selection and all permutation tests on every run using fixed seeds.

## Analysis conventions

Relative abundance (RA) is calculated separately for eukaryotes and prokaryotes. Original RAs are retained without renormalization after rare-ASV exclusion. Composition summaries and Shannon diversity use original ASV tables; clustering, bioindicators and dbRDA use the non-rare table (mean RA >0.01% across samples).

The nine sample environmental variables are assessed for log1p transformation in 6A. Pycnocline depth and SWHmax intentionally remain untransformed. Dunn comparisons use BH correction within each environmental variable; Shannon-environment correlations use BH correction separately for eukaryotes and prokaryotes. Tukey tests use their own multiplicity adjustment.

In 7A, candidate ASVs reach RA >=2% in at least two CL3-CL6 samples; bioindicator status is not a selection criterion. ASVs are selected while conditioning on the environmental predictors retained for prokaryotes. Variation partitioning uses those two selected predictor sets. Model and partition percentages use adjusted R2; axis and marginal contributions are unadjusted percentages of total variation. Final-model and marginal p-values are post-selection tests and do not repeat selection within permutations; full-candidate model p-values are also reported. The models describe associations, not causal influence or independently validated prediction.

## Outputs

| Manuscript output | Script and output location |
| --- | --- |
| Figures 1-3 | 6A: station map, combined cluster figure and environmental properties in `OUT/6A_ENV/`. |
| Figure 4 | 4A: `OUT/4A_PLANKT/plankton_communities.pdf` and `.png`. |
| Figure 5 | 5A: `OUT/5A_INDIC/amun_bioindicators.pdf` and `.png`. |
| Figure 6 | 7A: `OUT/7A_dbRDA/amun_CL3-6_dbRDA.pdf`, `.png` and `.svg`. |
| Figure S1 | 3A: `OUT/3A_CLUS/amun_complex_hm.pdf` and `.png`. |
| Tables S1-S7 | 6A: `OUT/6A_ENV/Supplementary materials Tables S1-S7.xlsx`. |
| dbRDA statistics | 7A: `OUT/7A_dbRDA/amun_CL3-6_dbRDA_statistics.xlsx`. |

Tables S1-S7 contain metadata; original and non-rare ASV RAs; environmental summaries and comparisons; Shannon diversity comparisons; and Shannon-environment correlations. Tables S2-S3 contain **RA fractions (0-1), not percentages**, displayed to six decimals, with smaller nonzero values in scientific notation. Spearman coefficients use two decimals. Other measurements/percentages use one decimal, counts use whole numbers, and p-values/small nonzero measurements use scientific notation. Stored numeric values retain full precision; missing entries are `NA`.

6A also saves `amun_group_summaries.xlsx`, `amun_genera_summaries.xlsx` and `amun_bioindicator_summaries.xlsx`, covering displayed major groups, genus-level labels (including unclassified taxa), and prevalent bioindicator ASVs. Group and genus-level summaries include all constituent ASVs in the original tables. Means include zeros.

- `Cluster`, `Station group`, `Layer` and `Station group-Layer`: sample counts and mean RA (%).
- `Samples`: individual RA (%), cluster, station group, layer, environmental conditions and read counts.
- `Above 0% RA`, `Above 25% RA`, `Above 50% RA` and `Above 75% RA`: sample/station counts and depth, Conservative temperature and salinity ranges for samples strictly above each threshold.

Station-group sheets include transects and station groups such as L and M.

7A exports six sheets: `Models`, `Fractions`, `Predictor contributions`, `Selection steps`, `Candidate ASVs` and `Displayed axes`. `Models` includes the full-candidate sample-permutation test alongside final-model sample and whole-station tests. `Fractions` reports environment-only, ASV-only, shared and unexplained fractions, plus total explained variation. Fitted models remain in the R session and are not saved by 7A.

## Data sources

Source links are retained in script comments: [Amundsen Science CTD data](https://doi.org/10.5884/12713), [PR2](https://pr2-database.org/), [SILVA](https://www.arb-silva.de/), [Natural Earth](https://www.naturalearthdata.com/), [ERA5](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels), and NOAA bathymetry accessed through `marmap`. The 1A MIMARKS export uses BioProject `PRJNA421293`.

For questions, contact kazmiruz@myumanitoba.ca or cryomics@gmail.com.

References for packages
-----------------------

Campitelli, E. (2025). ggnewscale: Multiple fill and colour scales in 'ggplot2' (R package version 0.5.2) [Computer software]. https://CRAN.R-project.org/package=ggnewscale (https://doi.org/10.32614/CRAN.package.ggnewscale)

De Cáceres, M., & Legendre, P. (2009). Associations between species and groups of sites: indices and statistical inference. Ecology, 90, 3566-3574. https://doi.org/10.1890/08-1823.1

de Vries, A., & Ripley, B. D. (2024). ggdendro: Create dendrograms and tree diagrams using 'ggplot2' (R package version 0.2.0) [Computer software]. https://CRAN.R-project.org/package=ggdendro (https://doi.org/10.32614/CRAN.package.ggdendro)

Dinno, A. (2026). dunn.test: Dunn's test of multiple comparisons using rank sums (R package version 1.4.1) [Computer software]. https://CRAN.R-project.org/package=dunn.test (https://doi.org/10.32614/CRAN.package.dunn.test)

Dunnington, D. (2025). ggspatial: Spatial data framework for ggplot2 (R package version 1.1.10) [Computer software]. https://CRAN.R-project.org/package=ggspatial (https://doi.org/10.32614/CRAN.package.ggspatial)

Fox, J., & Weisberg, S. (2019). An R companion to applied regression (3rd ed.). Sage. https://www.john-fox.ca/Companion/

Galili, T. (2015). dendextend: An R package for visualizing, adjusting, and comparing trees of hierarchical clustering. Bioinformatics. https://doi.org/10.1093/bioinformatics/btv428

Garnier, S., Ross, N., Rudis, R., Camargo, A. P., Sciaini, M., & Scherer, C. (2024). viridis(Lite): Colorblind-friendly color maps for R (R package version 0.6.5) [Computer software]. https://sjmgarnier.github.io/viridis/ (https://doi.org/10.5281/zenodo.4679423)

Grolemund, G., & Wickham, H. (2011). Dates and times made easy with lubridate. Journal of Statistical Software, 40(3), 1-25. https://www.jstatsoft.org/v40/i03/

Gu, Z., Eils, R., & Schlesner, M. (2016). Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics. https://doi.org/10.1093/bioinformatics/btw313

Gu, Z. (2022). Complex heatmap visualization. iMeta. https://doi.org/10.1002/imt2.43

Hijmans, R., Brown, A., & Barbosa, M. (2026). terra: Spatial data analysis (R package version 1.9-34) [Computer software]. https://CRAN.R-project.org/package=terra (https://doi.org/10.32614/CRAN.package.terra)

Kelley, D., Richards, C., & WG127 SCOR/IAPSO. (2024). gsw: Gibbs Sea Water functions (R package version 1.2-0) [Computer software]. https://CRAN.R-project.org/package=gsw (https://doi.org/10.32614/CRAN.package.gsw)

Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M., & Hornik, K. (2026). cluster: Cluster analysis basics and extensions (R package version 2.1.8.2) [Computer software]. https://CRAN.R-project.org/package=cluster

Müller, K., & Wickham, H. (2026). tibble: Simple data frames (R package version 3.3.1) [Computer software]. https://CRAN.R-project.org/package=tibble (https://doi.org/10.32614/CRAN.package.tibble)

Oksanen, J., Simpson, G. L., Blanchet, F. G., Kindt, R., Legendre, P., Minchin, P. R., O'Hara, R. B., Solymos, P., Stevens, M. H. H., Szöcs, E., Wagner, H., Bedward, M., Bolker, B., Borcard, D., Carvalho, G., De Cáceres, M., Durand, S., Evangelista, H. B. A., Hannigan, G., Hill, M. O., Lahti, L., Martino, C., Ouellette, M.-H., Ribeiro Cunha, E., Smith, T., Stier, A., Ter Braak, C. J. F., & Weedon, J. (2026). vegan: Community ecology package (R package version 2.7-5) [Computer software]. https://CRAN.R-project.org/package=vegan (https://doi.org/10.32614/CRAN.package.vegan)

Pante, E., Simon-Bouhet, B., & Irisson, J.-O. (2025). marmap: Import, plot and analyze bathymetric and topographic data (R package version 1.0.12) [Computer software]. https://CRAN.R-project.org/package=marmap (https://doi.org/10.32614/CRAN.package.marmap)

Pedersen, T. L. (2025). patchwork: The composer of plots (R package version 1.3.2) [Computer software]. https://CRAN.R-project.org/package=patchwork (https://doi.org/10.32614/CRAN.package.patchwork)

Pebesma, E. (2018). Simple features for R: Standardized support for spatial vector data. The R Journal, 10(1), 439-446. https://doi.org/10.32614/RJ-2018-009

Pebesma, E., & Bivand, R. (2023). Spatial data science: With applications in R. Chapman & Hall/CRC. https://r-spatial.org/book/ (https://doi.org/10.1201/9780429459016)

Petukhov, V., van den Brand, T., & Biederstedt, E. (2023). ggrastr: Rasterize layers for 'ggplot2' (R package version 1.0.2) [Computer software]. https://CRAN.R-project.org/package=ggrastr (https://doi.org/10.32614/CRAN.package.ggrastr)

R Core Team. (2026). R: A language and environment for statistical computing. R Foundation for Statistical Computing. https://www.R-project.org/ (https://doi.org/10.32614/R.manuals)

Schauberger, P., & Walker, A. (2025). openxlsx: Read, write and edit xlsx files (R package version 4.2.8.1) [Computer software]. https://CRAN.R-project.org/package=openxlsx (https://doi.org/10.32614/CRAN.package.openxlsx)

Simpson, G. (2026). permute: Functions for generating restricted permutations of data (R package version 0.9-10) [Computer software]. https://CRAN.R-project.org/package=permute (https://doi.org/10.32614/CRAN.package.permute)

Slowikowski, K. (2026). ggrepel: Automatically position non-overlapping text labels with 'ggplot2' (R package version 0.9.8) [Computer software]. https://CRAN.R-project.org/package=ggrepel (https://doi.org/10.32614/CRAN.package.ggrepel)

Warnes, G. R., Bolker, B., Bonebakker, L., Gentleman, R., Huber, W., Liaw, A., Lumley, T., Maechler, M., Magnusson, A., Moeller, S., Schwartz, M., Venables, B., & Galili, T. (2025). gplots: Various R programming tools for plotting data (R package version 3.3.0) [Computer software]. https://CRAN.R-project.org/package=gplots (https://doi.org/10.32614/CRAN.package.gplots)

Wickham, H. (2016). ggplot2: Elegant graphics for data analysis. Springer-Verlag New York. https://ggplot2.tidyverse.org

Wickham, H. (2025). stringr: Simple, consistent wrappers for common string operations (R package version 1.6.0) [Computer software]. https://CRAN.R-project.org/package=stringr (https://doi.org/10.32614/CRAN.package.stringr)

Wickham, H., Averick, M., Bryan, J., Chang, W., McGowan, L. D., François, R., Grolemund, G., Hayes, A., Henry, L., Hester, J., Kuhn, M., Pedersen, T. L., Miller, E., Bache, S. M., Müller, K., Ooms, J., Robinson, D., Seidel, D. P., Spinu, V., Takahashi, K., Vaughan, D., Wilke, C., Woo, K., & Yutani, H. (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686. https://doi.org/10.21105/joss.01686

Wickham, H., & Bryan, J. (2026). readxl: Read Excel files (R package version 1.5.0) [Computer software]. https://CRAN.R-project.org/package=readxl (https://doi.org/10.32614/CRAN.package.readxl)

Wickham, H., François, R., Henry, L., Müller, K., & Vaughan, D. (2026). dplyr: A grammar of data manipulation (R package version 1.2.1) [Computer software]. https://CRAN.R-project.org/package=dplyr (https://doi.org/10.32614/CRAN.package.dplyr)

Wickham, H., & Henry, L. (2026). purrr: Functional programming tools (R package version 1.2.2) [Computer software]. https://CRAN.R-project.org/package=purrr (https://doi.org/10.32614/CRAN.package.purrr)

Wickham, H., Henry, L., Pedersen, T., Luciani, T., Decorde, M., & Lise, V. (2025). svglite: An 'SVG' graphics device (R package version 2.2.2) [Computer software]. https://CRAN.R-project.org/package=svglite (https://doi.org/10.32614/CRAN.package.svglite)

Wickham, H., Hester, J., & Bryan, J. (2026). readr: Read rectangular text data (R package version 2.2.0) [Computer software]. https://CRAN.R-project.org/package=readr (https://doi.org/10.32614/CRAN.package.readr)

Wickham, H., Pedersen, T., & Seidel, D. (2025). scales: Scale functions for visualization (R package version 1.4.0) [Computer software]. https://CRAN.R-project.org/package=scales (https://doi.org/10.32614/CRAN.package.scales)

Wickham, H., Vaughan, D., & Girlich, M. (2025). tidyr: Tidy messy data (R package version 1.3.2) [Computer software]. https://CRAN.R-project.org/package=tidyr (https://doi.org/10.32614/CRAN.package.tidyr)

Additional software and reference databases
------------------------------------------

Camacho, C., et al. (2009). BLAST+: architecture and applications. BMC Bioinformatics, 10, 421. https://doi.org/10.1186/1471-2105-10-421

Guillou, L., et al. (2013). The Protist Ribosomal Reference database (PR2): a catalog of unicellular eukaryote Small Sub-Unit rRNA sequences with curated taxonomy. Nucleic Acids Research, 41, D597-D604. https://doi.org/10.1093/nar/gks1160

Quast, C., et al. (2013). The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucleic Acids Research, 41, D590-D596. https://doi.org/10.1093/nar/gks1219

Biostrings: Efficient manipulation of biological strings. https://doi.org/10.18129/B9.bioc.Biostrings.