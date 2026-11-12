# 2A_TAX.R ---------------------------------------------------------------------------------------
#
# Purpose:
# Clean 16S taxonomy and assign 18S taxonomy with PR2 BLAST references, stopping at genus.
# Select non-sequentially filtered samples with both markers; preserve original ASVs for composition summaries.
# Create non-rare ASV tables for clustering, indicator analysis and dbRDA.
#
# Inputs (paths are relative to the project working directory containing EXT, OUT and DNA):
# EXT/meta.amun.complete.rds (from 1A_META.R)
# DNA/515F-806RB/seqtab_final.rds; DNA/Euk1391F-EukBr/seqtab_final.rds
# EXT/tax.16S.rds; EXT/hits_pr2_filt.rds (cached taxonomy inputs)
# Optional filtered-cache regeneration: EXT/PR2/pr2_version_5.1.1_merged.csv; EXT/hits_pr2.rds
# Optional BLAST rerun: EXT/tax.18S.rds and PR2 BLAST database in EXT/PR2/
# Optional classifier regeneration: taxout.rds and bootout.rds in each DNA marker folder
#
# Outputs:
# EXT/meta.amun.rds (paired-marker metadata and original read totals)
# EXT/proptab.amun.rds (non-rare ASVs; original within-marker RA denominators)
# EXT/proptab.amun.original.rds (all cleaned ASVs)
# EXT/amun_taxonomy.rds (full taxonomy and sequence lookup for 4A and 6A)
# EXT/amun_original_counts.rds (full integer counts for optional Shannon recalculation in 4A)
#
# Run sections in order within this script. Later scripts load saved inputs; they do not need earlier scripts in memory.
# Heavy calculations remain optional or cached. Existing EXT and OUT folders are assumed to exist.
#---------------------------------------------------------------------------------------------------

library(tidyverse)

## Windows PATH example for BLAST+; on other systems add the executable directory to PATH using that system's conventions.
# Sys.setenv(PATH=paste("C:/Program Files/NCBI/blast-2.17.0+/bin", Sys.getenv("PATH"), sep=";"))
# Sys.which("blastn") # check

# vocab: qseqid - Query Seq-id; sseqid - Subject Seq-id; length - Alignment length; mismatch - Number of mismatches;
#        gapopen - # of gap openings (gap: a space introduced into an alignment to compensate for insertions/deletions in one seq rel-ve to another);
#        qstart/qend - Start/End of alignment in query; sstart/send - Start/End of alignment in subject;
#        evalue - the smaller, the better the match - the likelihood that a given sequence match is purely by chance -
#                 if 10, up to 10 hits can be expected to be found just by chance, given the same size of a random db;
#                 if E<1e-50 (or 1×10-50), there should be an extremely high confidence that the database match is a result of homologous relat-ps;
#                 if E is btwn 0.01 and 1e-50, the match can be considered a result of homology;
#        bitscore - the higher, the better the match - the req-d size of a seq db in which the current match could be found just by chance;
#                   unlike evalue, doesn't decrease with increasing db size, doesn't depend on db size at all
#                   hence, provides better rule-of-thumb for inferring homology: for avg length proteins, 50+ is almost always significant;
#        pident - Percentage of identical matches - # of matches / # of positions - e.g., for CCCCCCCCCC (target) and CCCACCACCC (query), it is 80%.


####################################################
### OPTIONAL REGENERATION OF CLASSIFIER TAXONOMY ###

# ## 16S
# # tax.16S is sequences vs Domain/Phylum/Class/Order/Family/Genus
# # boot.16S is bootstrap confidence scores for the tax.16S table
# # bootstrapping is a resampling technique used to estimate confidence intervals / variability by repeatedly sampling with replacement
# tax.16S = data.frame(readRDS("DNA/515F-806RB/taxout.rds")$`ref_dada2_silva_v138-2_train_set_uniq.fasta.gz`)
# boot.16S = data.frame(readRDS("DNA/515F-806RB/bootout.rds")$`ref_dada2_silva_v138-2_train_set_uniq.fasta.gz`)
# tax.16S <- tax.16S %>% mutate(across(everything(), ~ paste0(.x, " [", boot.16S[[cur_column()]], "]")))
# saveRDS(tax.16S, "EXT/tax.16S.rds")
# 
# ## 18S - do the same as for 16S, notice how there are more taxonomic classifications (columns)
# tax.18S = data.frame(readRDS("DNA/Euk1391F-EukBr/taxout.rds")$ref_dada2_pr2_version_5.0.0_SSU.fasta.gz)
# boot.18S = data.frame(readRDS("DNA/Euk1391F-EukBr/bootout.rds")$ref_dada2_pr2_version_5.0.0_SSU.fasta.gz)
# tax.18S <- tax.18S %>% mutate(across(everything(), ~ paste0(.x, " [", boot.18S[[cur_column()]], "]")))
# saveRDS(tax.18S, "EXT/tax.18S.rds")


########################################
### Load raw ASV tables and taxonomy ###

# Load raw ASV tables for 16S and 18S sequences
seqtab.16S.raw <- readRDS("DNA/515F-806RB/seqtab_final.rds") %>% { rownames(.) <- gsub("716S","7-16S",rownames(.)); . }  # 716S -> 7-16S in rownames
seqtab.18S.raw <- readRDS("DNA/Euk1391F-EukBr/seqtab_final.rds")

# Load raw taxonomic identification for 16 and 18S sequences + inspect
tax.16S.raw <- readRDS("EXT/tax.16S.rds") # naive Bayesian classifier against the SILVA database
# hist(tax.16S.raw$Genus %>% str_extract("(?<=\\[)\\d+(?=\\])") %>% as.integer(), breaks = 100, main = "16S bootstraps", xlab = "Bootstrap (%)")

# Needed only for the optional BLAST rerun below.
# tax.18S.raw <- readRDS("EXT/tax.18S.rds") # naive Bayesian classifier against the PR2 database
# hist(tax.18S.raw$Class %>% str_extract("(?<=\\[)\\d+(?=\\])") %>% as.integer(), breaks = 100, main = "18S bootstraps", xlab = "Bootstrap (%)")
# for 18S, naive Bayesian classifier didn't produce good tax. identif-n; likely reason: not enough coverage for the 1391f-EukBr primer pair used
# hence, here I implement a custom pipeline that applies BLAST against the PR2 database instead


####################
### Clean up 16S ###

ranks16 <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
fill_unclassified <- function(df, ranks) {                                                # function for assigning “unclassified <lowest known rank>” to NAs
  for (k in 2:length(ranks)) {
    cur <- ranks[k]; prev <- ranks[k - 1]; x <- df[[cur]]; p <- df[[prev]]
    df[[cur]] <- ifelse(is.na(x) & !is.na(p), ifelse(grepl("^uncl. ", p), p, paste0("uncl. ", p)), x)       # checks if uncl. assig-t is there to not repeat
  }
  df
}

# Remove raw organelle/eukaryote assignments before low-confidence labels are masked.
tax.16S <- tax.16S.raw %>%
  filter(!if_any(all_of(ranks16), ~ str_detect(coalesce(.x, ""), "Eukaryota|Chloroplast|Mitochondria"))) %>%
  mutate(across(all_of(ranks16), \(x) {                                         # assign NA to bootstrap < 80 and names containing Incertae Sedis or Unknown
    b <- str_extract(x, "(?<=\\[)\\d+(?=\\])") %>% as.integer()
    x <- ifelse(b < 80, NA_character_, x)
    x <- ifelse(str_detect(x, "Incertae Sedis|Unknown|Candidatus Nitrosoarchaeum"), NA_character_, x)    # Candidatus Nitrosoarchaeum is treated as unresolved (SILVA 138.2 release notes)
    x
  })) %>%
  filter(!is.na(Domain)) %>% 
  fill_unclassified(ranks16) %>% mutate(across(all_of(ranks16), ~ str_remove(.x, "\\s*\\[\\d+\\]")))           # fill NAs w/ unclassified and rmv bootstraps

seqtab.16S <- seqtab.16S.raw[, rownames(tax.16S), drop = FALSE]


#############################
### BLAST 18S against PR2 ###

# PR2 reference data for the optional BLAST rerun below.
# pr2 <- read_csv("EXT/PR2/pr2_version_5.1.1_merged.csv", show_col_types = FALSE)

## RUN THIS ONCE IF NEED TO BLAST AGAINST PR2 AGAIN ##

# fa_file <- "EXT/PR2/PR2_5.1.1.fasta"
# fa <- Biostrings::DNAStringSet(pr2$sequence)
# names(fa) <- pr2$pr2_accession
# Biostrings::writeXStringSet(fa, filepath = fa_file, width = 80)
# db_pr2 <- "EXT/PR2/PR2_5_1_1"
# system2("makeblastdb", c("-in", fa_file, "-dbtype", "nucl", "-parse_seqids", "-out", "EXT/PR2/PR2_5_1_1"))
#
# # run all sequences from tax.18S
# xs <- Biostrings::DNAStringSet(rownames(tax.18S.raw))
# names(xs) <- paste0("Q", seq_along(xs))
#
# predict_pr2 <- function(xs, db, batch = 100, max_hits = 300) {
#   idx_list <- split(seq_along(xs), ceiling(seq_along(xs) / batch))
#   fmt_str <- "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qcovs qcovus"
#   fmt <- shQuote(fmt_str)   # quote on every platform; otherwise system2() lets the shell split the multi-word -outfmt string on macOS/Linux
#   cols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
#             "qcovhsp", "qcovs", "qcovus")
#   out <- vector("list", length(idx_list))
#   total <- length(xs)
#   done <- 0
# 
#   for (k in seq_along(idx_list)) {
#     i  <- idx_list[[k]]
#     fa <- tempfile(fileext = ".fa")
#     Biostrings::writeXStringSet(xs[i], filepath = fa)
#     raw <- system2("blastn", c("-db", db, "-max_target_seqs", as.character(max_hits), "-query", fa, "-outfmt", fmt), stdout = TRUE, stderr = TRUE)
#     raw <- raw[nzchar(raw)]
#     if (!length(raw)) next
#     df <- read_tsv(I(paste(raw, collapse = "\n")), col_names = cols, show_col_types = FALSE) %>%
#       mutate(across(c(pident, evalue, bitscore, qcovhsp, qcovs, qcovus), as.numeric),
#              across(c(length, mismatch, gapopen, qstart, qend, sstart, send), as.integer))
# 
#     # attach original sequences
#     map <- tibble(qseqid = names(xs)[i], seq = as.character(xs[i]))
#     df <- left_join(map, df, by = "qseqid") %>% relocate(seq)
# 
#     # quick progress: one line per query
#     first <- df[!duplicated(df$qseqid), ]
#     for (j in seq_len(nrow(first))) {
#       r <- first[j, ]
#       done <- done + 1
#       cat(sprintf("[%d/%d] %s -> %s | pid=%.1f len=%d e=%g bits=%.1f\n", done, total, r$qseqid, r$sseqid, r$pident, r$length, r$evalue, r$bitscore))
#     }
#     out[[k]] <- df
#   }
#   bind_rows(out)
# }
# hits_pr2 <- predict_pr2(xs, db = db_pr2)   # takes 50 minutes
# write_rds(hits_pr2, "EXT/hits_pr2.rds")


#################################
### Annotate and clean up 18S ###

# PR2 taxonomy lookup
ranks18 <- c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")
ranks18_assign <- setdiff(ranks18, "Species")                            # assignments stop at Genus; species names remain reference information

# # To regenerate taxonomic assignments (hits_pr2_filt), uncomment the block below.
#
# pr2 <- read_csv("EXT/PR2/pr2_version_5.1.1_merged.csv", show_col_types = FALSE)
# hits_pr2 <- read_rds("EXT/hits_pr2.rds")     # load raw BLAST hits only to regenerate the saved assignments
# ranks18_eff <- c("Genus", "Family", "Order", "Class", "Subdivision", "Division", "Supergroup")
# pr2_tax <- pr2 %>% transmute(sseqid = pr2_accession, across(all_of(tolower(ranks18)))) %>% rename_with(str_to_title, all_of(tolower(ranks18)))
#
# # qcovhsp is the query coverage of a single alignment, equal to (qend - qstart + 1) / nchar(seq) * 100 rounded up;
# # qcovs and qcovus pool all alignments to the same subject. No query coverage exceeds 100.
#
# hits_pr2_filt <- hits_pr2 %>%
#   filter(qcovhsp == 100, pident >= 80) %>%                              # retain alignments with 100% query coverage and ≥80% identity
#   filter(bitscore == max(bitscore), .by = qseqid) %>%                   # highest bitscore per qseqid
#   filter(pident   == max(pident), .by = qseqid) %>%                     # highest identity among those
#   filter(mismatch == min(mismatch), .by = qseqid) %>%                   # fewest mismatches
#   filter(gapopen  == min(gapopen), .by = qseqid) %>%                    # fewest gap openings
#   left_join(pr2_tax, by = "sseqid")                                     # join taxonomy only after removing lower-ranked alignments
#
# # Check whole columns instead of using rowwise()/c_across() separately within every ASV
# eff_idx <- rep.int(length(ranks18_eff) + 1L, nrow(hits_pr2_filt))
# eff_tax <- rep(NA_character_, nrow(hits_pr2_filt))
# for (k in rev(seq_along(ranks18_eff))) {                                # move towards Genus so the deepest resolved annotation wins
#   v <- hits_pr2_filt[[ranks18_eff[k]]]
#   known <- !is.na(v) & !str_detect(v, "_X")
#   eff_idx[known] <- k
#   eff_tax[known] <- v[known]
# }
#
# hits_pr2_filt <- hits_pr2_filt %>% mutate(eff_idx = eff_idx, eff_tax = eff_tax) %>%
#   filter(eff_idx == min(eff_idx), .by = qseqid) %>%                     # prefer the deepest resolved annotation among equally ranked alignments
#   mutate(tax_n = n_distinct(sseqid), .by = c(qseqid, eff_tax)) %>%      # each reference gets one vote, even if it has multiple alignments
#   filter(tax_n == max(tax_n), .by = qseqid) %>%                         # retain all equally frequent assignments at that rank
#   mutate(across(all_of(ranks18_assign), ~ ifelse(str_detect(.x, "_X"), NA_character_, .x)))
#
# # Summarise all ASVs together; reference species are pooled but never enter selection or consensus
# first_cols <- setdiff(names(hits_pr2_filt), c("qseqid", ranks18, "sseqid", "sstart", "send", "eff_idx", "eff_tax", "tax_n"))
# hits_pr2_filt <- hits_pr2_filt %>% summarise(
#   across(all_of(first_cols), dplyr::first),
#   n_best_hits = n_distinct(sseqid),
#   across(all_of(c("sseqid", "sstart", "send")), ~ paste(unique(.x), collapse = "; ")),
#   Species_hits = paste(unique(Species), collapse = "; "),
#   across(all_of(ranks18_assign), list(
#     taxon = ~ if (!anyNA(.x) && length(unique(.x)) == 1L) .x[1] else NA_character_,
#     conflict = ~ length(unique(.x[!is.na(.x)])) > 1L), .names = "{.col}__{.fn}"), .by = qseqid) %>%
#   rename_with(~ str_remove(.x, "__taxon$"), ends_with("__taxon"))
#
# # A conflict at any ancestral rank clears that rank and all lower ranks; an unresolved name alone is not a conflict
# conflict <- as.matrix(hits_pr2_filt[paste0(ranks18_assign, "__conflict")])
# blocked <- rep(FALSE, nrow(hits_pr2_filt))
# assigned <- rep(NA_character_, nrow(hits_pr2_filt))
# for (k in seq_along(ranks18_assign)) {
#   r <- ranks18_assign[k]
#   blocked <- blocked | conflict[, k]
#   hits_pr2_filt[[r]][blocked] <- NA_character_
#   assigned[!is.na(hits_pr2_filt[[r]])] <- r
# }
# hits_pr2_filt <- hits_pr2_filt %>% mutate(Assigned_rank = assigned) %>% select(-ends_with("__conflict")) %>% filter(Domain == "Eukaryota")
# write_rds(hits_pr2_filt, "EXT/hits_pr2_filt.rds", compress = "xz")

# Load the saved filtered and annotated BLAST results for normal runs.
hits_pr2_filt <- read_rds("EXT/hits_pr2_filt.rds")

# exploring hits_pr2 and the filtered version
# unique(hits_pr2 %>% left_join(pr2_tax, by = "sseqid") %>% .$Domain)  # all distinct domains in the unfiltered
# plot(as.numeric(hits_pr2_filt$gapopen), hits_pr2_filt$bitscore)      # gap opening vs bit score
# plot(hits_pr2_filt$qcovs, hits_pr2_filt$qcovus)                      # qcovs vs qcovus - same
# hist(hits_pr2_filt$n_best_hits, breaks = 1000, xlim = c(0, 10))      # distinct selected reference sequences

tax.18S <- hits_pr2_filt %>% select(seq, all_of(ranks18_assign), Species_hits) %>% column_to_rownames("seq") %>%
  fill_unclassified(ranks18_assign) %>% filter(Genus != "Homo")        # fill assigned ranks only; pooled reference species stay in Species_hits
seqtab.18S <- seqtab.18S.raw[, rownames(tax.18S), drop = FALSE]

# Check and display the distribution of non-zero human DNA (Homo) read counts per sample that were removed from the 18S dataset
# seqtab.18S.raw[, hits_pr2_filt$seq[hits_pr2_filt$Genus %in% "Homo"], drop = FALSE] %>% .[rowSums(.) > 0, , drop = FALSE] %>% as.data.frame() %>%
#   rownames_to_column("Sample_ID") %>% mutate(Total_Homo_Reads = rowSums(select(., -Sample_ID))) %>% arrange(desc(Total_Homo_Reads)) %>% View()


#########################################################################################################################
### Select the AMUN samples that were collected by non-sequential filtration and have both 16S and 18S data available ###

extract_plate_geno <- function(seqtab) {                                         ## plate.geno extractor
  rownames(seqtab) %>% tibble(value = .) %>%
    separate(value, into = c("project", "primer", "plate", "well"), sep = "\\|", fill = "right") %>%
    separate(plate, into = paste0("p", 1:6), sep = "-", fill = "right") %>%
    mutate(plate.geno = toupper(paste(p2, well, sep = "_"))) %>% pull(plate.geno)
}
plates.16S <- extract_plate_geno(seqtab.16S)
plates.18S <- extract_plate_geno(seqtab.18S)

# 115/118 have both 16S and 18S; get only the ones with non-sequential filtration
meta.amun.complete <- readRDS("EXT/meta.amun.complete.rds")
meta.amun <- meta.amun.complete %>% filter(`Sequential filtration` == 0, plate.geno %in% plates.16S, plate.geno %in% plates.18S)
sam.16S.pick <- match(meta.amun$plate.geno, plates.16S); sam.18S.pick <- match(meta.amun$plate.geno, plates.18S)

# meta.amun.complete %>% filter(`Sequential filtration` == 0, plate.geno %in% plates.16S, !plate.geno %in% plates.18S) # 2 have only 16S
# meta.amun.complete %>% filter(`Sequential filtration` == 0, plate.geno %in% plates.18S, !plate.geno %in% plates.16S) # 1 has only 18S


########################################
### ORIGINAL AND NON-RARE ASV TABLES ###

## AMUN-only count tables (still with sequences as colnames); remove ASVs absent from every AMUN sample collected using non-sequential filtration
tab.16S <- seqtab.16S[sam.16S.pick, , drop = FALSE]; tab.18S <- seqtab.18S[sam.18S.pick, , drop = FALSE]
tab.16S <- tab.16S[, colSums(tab.16S) > 0, drop = FALSE]; tab.18S <- tab.18S[, colSums(tab.18S) > 0, drop = FALSE]
# ncol(tab.16S) # how many 16S ASVs originally?  5588
# ncol(tab.18S) # how many 18S ASVs originally?   710
# much lower prokaryotic and eukaryotic diversity than in WK due to the latter being sampled across a broader environmental and salinity gradient

## Calculate relative abundances using all cleaned ASVs; keep the original table for group and genus summaries.
proptab.original <- cbind(prop.table(as.matrix(tab.16S), margin = 1), prop.table(as.matrix(tab.18S), margin = 1)) %>% as.data.frame()

## Retain non-rare ASVs for clustering, bioindicators and ASV-based dbRDA; do not renormalize after filtering.
proptab.small <- proptab.original %>% select(where(~ mean(.x) > 1e-4))

## ASV metadata: which sequence, which marker, ASV ID
n16 <- ncol(tab.16S)
all_seqs <- c(colnames(tab.16S), colnames(tab.18S))
base_asv <- tibble(idx = match(colnames(proptab.small), all_seqs), Sequence = all_seqs[idx], ASV = paste0("ASV_", seq_along(idx)))

## 16S ASVs: taxonomy only
asv16 <- base_asv %>% filter(idx <= n16) %>% select(-idx) %>% left_join(tax.16S %>% rownames_to_column("Sequence"), by = "Sequence")

## 18S ASVs: taxonomy + BLAST quality info
asv18 <- base_asv %>% filter(idx > n16) %>% select(-idx) %>% left_join(tax.18S %>% rownames_to_column("Sequence"), by = "Sequence") %>%
  left_join(hits_pr2_filt %>% transmute(Sequence = seq, qcovhsp, bitscore, evalue, pident, mismatch, gapopen), by = "Sequence")

## Build final column labels and apply to proportional table
asv_labels <- bind_rows(asv16 %>% mutate(TaxLabel = paste(Domain, Phylum, Order, Genus, sep = "; ")),
                        asv18 %>% mutate(TaxLabel = paste(Domain, Subdivision, Order, Genus, sep = "; "))) %>%
  arrange(match(Sequence, colnames(proptab.small))) %>% transmute(TaxLabel_ASV = paste(TaxLabel, ASV, sep = "; ")) # put labels in the same order as proptab

proptab.small.taxan <- proptab.small %>% `colnames<-`(asv_labels$TaxLabel_ASV)

## Keep all original ASVs for major group and genus summaries, including those excluded by the ASV rarity cutoff.
## Existing ASV numbers are preserved; rare ASVs receive additional numbers after the non-rare set.
base_asv.original <- tibble(idx = seq_along(all_seqs), Sequence = all_seqs) %>% left_join(base_asv %>% select(idx, ASV), by = "idx") %>%
  mutate(ASV = coalesce(ASV, paste0("ASV_", nrow(base_asv) + cumsum(is.na(ASV)))))

asv16.original <- base_asv.original %>% filter(idx <= n16) %>% select(-idx) %>% left_join(tax.16S %>% rownames_to_column("Sequence"), by = "Sequence")
asv18.original <- base_asv.original %>% filter(idx > n16) %>% select(-idx) %>% left_join(tax.18S %>% rownames_to_column("Sequence"), by = "Sequence")

asv_labels.original <- bind_rows(asv16.original %>% mutate(TaxLabel = paste(Domain, Phylum, Order, Genus, sep = "; ")),
                                 asv18.original %>% mutate(TaxLabel = paste(Domain, Subdivision, Order, Genus, sep = "; "))) %>%
  arrange(match(Sequence, colnames(proptab.original))) %>% transmute(TaxLabel_ASV = paste(TaxLabel, ASV, sep = "; "))

## Load metadata for the samples and use it to give more descriptive sample names
meta.amun <- meta.amun %>% mutate(updated_rowname = paste0(Station, " (", round(`Depth [m]`), " m)")) %>%
  column_to_rownames("updated_rowname") %>%
  mutate(`Eukaryotic Reads` = rowSums(tab.18S), `Prokaryotic Reads` = rowSums(tab.16S)) # total reads per sample

proptab.amun <- proptab.small.taxan %>% `rownames<-`(rownames(meta.amun))
proptab.amun.original <- proptab.original %>% `colnames<-`(asv_labels.original$TaxLabel_ASV) %>% `rownames<-`(rownames(meta.amun))

# Save paired-marker metadata and the two RA tables as explicit inputs to 3A_CLUS.R and later scripts.
# RDS preserves row names, column names, data types, and object classes exactly.
saveRDS(meta.amun, "EXT/meta.amun.rds")
saveRDS(proptab.amun, "EXT/proptab.amun.rds")
saveRDS(proptab.amun.original, "EXT/proptab.amun.original.rds")


###############################
### CHECK NON-RARE COVERAGE ###
# Report the number of retained ASVs and their share of original read abundance for Methods 2.4.1.
# Uses: Original and non-rare relative abundance tables; no renormalization after filtering.
# Produces: A printed coverage summary.
nonrare_coverage <- tibble(Domain = c("Eukaryote", "Prokaryote"), ASVs = c(nrow(asv18), nrow(asv16)),
                           Mean_RA_pct = c(mean(rowSums(proptab.amun[, startsWith(names(proptab.amun), "Eukaryota;")])),
                                           mean(rowSums(proptab.amun[, !startsWith(names(proptab.amun), "Eukaryota;")])))) %>%
  mutate(Mean_RA_pct = 100 * Mean_RA_pct)
print(nonrare_coverage)
# 260 eukaryotic and 380 prokaryotic non-rare ASVs represent approximately 98.8% and 96.9% of RA, respectively.


#########################################
### SAVE TAXONOMY AND ORIGINAL COUNTS ###
# Make taxonomy available to 4A and original integer counts available for optional Shannon recalculation.
# Uses: Taxonomy and count tables for the 115 selected non-sequential samples.
# Produces: EXT/amun_taxonomy.rds; EXT/amun_original_counts.rds.
saveRDS(list(asv16.original = asv16.original, asv18.original = asv18.original, base_asv.original = base_asv.original), "EXT/amun_taxonomy.rds")
rownames(tab.16S) <- rownames(tab.18S) <- rownames(meta.amun)
saveRDS(list(Eukaryota = tab.18S, Prokaryota = tab.16S), "EXT/amun_original_counts.rds")
