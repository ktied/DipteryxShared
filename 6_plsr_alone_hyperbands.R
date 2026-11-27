####################################################

# This script 
# PLSR on the FieldSpec scans directly  
# I also tried a separate pixel-based approach in season1_sugar_treelevel_plsr.R.
# Ground-measured spectra paired with real sugar values, no cross-sensor
# VNIR only (500-900nm), used as measured -- no
# smoothing here, but probably should to match Martins?  
# elsewhere I smoothed airborne imagery with SG smoothing .
#

# Handled instead with a tree-grouped LOO-CV loop (holds out all of one
# tree's scans at once), not pls::plsr()'s built-in LOO.
#

#Things you might do: other nutrients, smoothing. not sure Martins bands now

####################################################

library(dplyr)
library(tidyr)
library(pls)        # plsr(), loading.weights()
library(ggplot2)
library(terra)      # vect() -- reads the tree-crown shapefile
library(prospectr)  # spliceCorrection() -- Part B
library(plsVarSel)  # VIP()              -- Part B
# install once if needed:
# install.packages(c("dplyr", "tidyr", "pls", "ggplot2", "terra", "prospectr", "plsVarSel"))

# this is for me - moving bt server and you can ignore it 
on_server <- FALSE  # <<< TOGGLE: TRUE = server (~/ home-relative paths),
# FALSE = local (//10.126.19.90/... UNC paths)

# Turns a network path into a server home path when on_server is TRUE
path_root <- function(unc_path) {
  if (on_server) sub("^//[^/]+/", "~/", unc_path) else unc_path
}

# All figures and .rds outputs go here - replace with your own :) 
out_folder_hyper <- path_root("//10.126.19.90/EAS_ind/ktiedeman/Dipteryx/data/BCI_fieldseason-01/Intermediate/")

## LOAD + CLEAN FIELDSPEC DATA
# Tree crowns: used only for the crosswalk from old tree IDs (TreeID2324) - 2023 - 2024 
# to the new tree IDs (newID) newID is the ids from 2024-2025 

alldipt2 <- vect(path_root("//10.126.19.90/EAS_ind/ktiedeman/Dipteryx/data/CrownStudy2425/TreeCrowns/Fall2025/DipteryxTreeCrowns_alledit_2026_01_19(1).shp"))


newid_to_oldid <- as.data.frame(alldipt2)[, c("newID", "TreeID2324")] %>%
  dplyr::mutate(newID = as.character(newID))
cat("Crowns loaded:", nrow(alldipt2), "| crosswalk rows:", nrow(newid_to_oldid), "\n")

# FieldSpec scans (median scan per spot) with all nutrient values attached
EndMem <- readRDS(path_root("//10.126.19.90/EAS_shared/nutritional_landscapes/working/processed/2024-season-01-bci/Intermediate/Endmembers_medianscansperspot_withAllNut_20250828.rds"))

# Spots flagged as outliers in earlier checks - in the fieldspec folder 
bark_outliers  <- read.csv(path_root("//10.126.19.90/EAS_ind/ktiedeman/Dipteryx/data/BCI_fieldseason-01/Intermediate/Endmembers_bark_outliers.csv"))
fruit_outliers <- read.csv(path_root("//10.126.19.90/EAS_ind/ktiedeman/Dipteryx/data/BCI_fieldseason-01/Intermediate/Endmembers_fruit_outliers.csv"))

# Cleaning: drop scans of unknown height, fruit collected from the ground
# ("bottom"), and the outlier spots
dataCL    <- subset(EndMem, height_collected != "for fun, unknown height")
dataCLrip <- subset(dataCL, !grepl("bottom", dataCL$height_collected))
dataCLrip <- subset(dataCLrip, !dataCLrip$ID_spot %in%
                      c(bark_outliers$ID_spot, fruit_outliers$ID_spot))
dataCLrip$Date <- as.Date(dataCLrip$Date)


head(dataCLrip)
# "Crude.Protein" "NDF"           "ADF"           "protein_fiber" "sugar_cor"     "pooled_sugar" 

# FieldSpec band columns are named x350, x351, ... x2500 (wavelength in nm)
fieldspec_band_cols <- grep("^x[0-9]+$", names(dataCLrip), value = TRUE)
fieldspec_wl_all    <- as.numeric(gsub("x", "", fieldspec_band_cols))


cat("Fieldspec band columns:", length(fieldspec_band_cols),
    "| range:", min(fieldspec_wl_all), "-", max(fieldspec_wl_all), "nm\n")

## MATCH TO THE AIRBORNE HYPERSPECTRAL BANDS
# The FieldSpec measures every 1 nm; the airborne sensor has fewer, wider
# bands. To make the FieldSpec models comparable with the imagery, keep
# only the FieldSpec band nearest to each airborne band center
# (A simple subset, not a full resampling with band widths.)
#could also do a smoothing instead... 

use_hyper_bands <- TRUE  # <<< TOGGLE: TRUE = keep only bands matching the planemounted sensor

# Band-number -> wavelength lookup tables for the airborne VNIR and SWIR sensors
hyper_vnir_lookup <- read.csv(path_root("//10.126.19.90/EAS_ind/ktiedeman/Dipteryx/data/VNIR_bandno_to_wavelength.csv"))
hyper_swir_lookup <- read.csv(path_root("//10.126.19.90/EAS_ind/ktiedeman/Dipteryx/data/swir_band_lookup.csv"))

hyper_vnir_wl <- as.numeric(gsub(" nm", "", hyper_vnir_lookup$Wave_num))
hyper_swir_wl <- as.numeric(gsub(" nm", "", hyper_swir_lookup$Wavelength))
hyper_swir_wl <- hyper_swir_wl[hyper_swir_wl >= 991]  # VNIR and SWIR overlap; use VNIR below 991 nm
hyper_wl      <- c(hyper_vnir_wl, hyper_swir_wl)

# Nearest FieldSpec wavelength to each airborne band center
hyper_fs_wl <- unique(sapply(hyper_wl, function(w) fieldspec_wl_all[which.min(abs(fieldspec_wl_all - w))]))
in_hyper    <- if (use_hyper_bands) fieldspec_wl_all %in% hyper_fs_wl else rep(TRUE, length(fieldspec_wl_all))
band_tag    <- if (use_hyper_bands) "_hyperbands" else ""   # added to output file names
cat("Airborne bands:", length(hyper_vnir_wl), "VNIR +", length(hyper_swir_wl), "SWIR",
    "-> matching FieldSpec bands:", length(hyper_fs_wl), "\n")


####################################################################

## PART A: SUGAR, VNIR ONLY (500-900 nm), spectra as measured - no smoothing

## A1: SUGAR SCANS -- every scan, all dates


##########We scanned each fruit 3 - 4 times in different spots - so here I pull only spot 2
# Fruit scans from spot 2 with a sugar value (0 = not measured, so dropped),
# matched to a newID

sugar_scans <- dataCLrip[dataCLrip$type == "Fruit" & dataCLrip$spot == 2 &
                           !is.na(dataCLrip$sugar_cor) & dataCLrip$sugar_cor != 0, ]
sugar_scans <- merge(sugar_scans, newid_to_oldid, by.x = "TreeID", by.y = "TreeID2324")

#now should have ~1000 scans of fruit rather than a 3k+ 

cat("\nSugar scans (spot==2, non-NA, non-zero, matched to a newID):", nrow(sugar_scans), "\n")
cat("Distinct trees:", dplyr::n_distinct(sugar_scans$newID), "| distinct dates:", dplyr::n_distinct(sugar_scans$Date), "\n")

cat("Sugar range (all scans):", round(range(sugar_scans$sugar_cor), 2), "\n")

## A2: KEEP ONLY THE VNIR BANDS (500-900 nm)
# No smoothing: FieldSpec scans are close-range and taken under controlled
# conditions, so they are much cleaner than the airborne imagery.

vnir_keep <- fieldspec_wl_all >= 500 & fieldspec_wl_all <= 900 & in_hyper  # in_hyper: airborne bands only (see above)
fieldspec_vnir_cols <- fieldspec_band_cols[vnir_keep]
fieldspec_vnir_wl   <- fieldspec_wl_all[vnir_keep]
cat("VNIR-only (500-900nm), bands kept:", length(fieldspec_vnir_cols), "\n") 

# (now uses only the bands that match the hyperspec -- see the "Match to the airborne hyperspectral bands" section above)

## A3: TRAINING TABLE -- a subset 
sugar_scans_vnir <- sugar_scans[, c("newID", "TreeID", "Date", "sugar_cor", fieldspec_vnir_cols)]
sugar_scans_vnir <- sugar_scans_vnir[complete.cases(sugar_scans_vnir[, fieldspec_vnir_cols]), ]
cat("\nPLSR training table (every scan, all dates):", nrow(sugar_scans_vnir), "rows,",
    dplyr::n_distinct(sugar_scans_vnir$newID), "trees\n")



## A4: LEAVE-ONE-TREE-OUT CV

# For each tree: fit a fresh PLSR on all other trees, predict this tree's
# scans. Each fit tries 1 to max_ncomp components, so afterwards we have a
# held-out prediction for every scan at every ncomp.

# Formula: sugar_cor ~ x500 + x501 + ... + x900
plsr_formula_fs <- as.formula(paste("sugar_cor ~", paste(fieldspec_vnir_cols, collapse = " + ")))

# Can't use more components than (number of trees - 2)
max_ncomp <- min(10, length(fieldspec_vnir_cols) - 1, dplyr::n_distinct(sugar_scans_vnir$newID) - 2)
cat("\nMax PLS components tried in LOTO-CV:", max_ncomp, "\n")

trees_fs <- unique(sugar_scans_vnir$newID)
loto_results_fs <- vector("list", length(trees_fs)); names(loto_results_fs) <- trees_fs

cat("\nRunning tree-grouped LOTO-CV...\n")
for (i in seq_along(trees_fs)) {
  held  <- trees_fs[i]
  train <- sugar_scans_vnir[sugar_scans_vnir$newID != held, ]  # every other tree
  test  <- sugar_scans_vnir[sugar_scans_vnir$newID == held, ]  # all scans of the held-out tree
  if (nrow(train) < 10 || nrow(test) < 1) {
    cat(sprintf("  [%d/%d] Tree %s -- skipping (too few rows)\n", i, length(trees_fs), held)); next
  }
  fit <- tryCatch(plsr(plsr_formula_fs, data = train, ncomp = max_ncomp, scale = TRUE),
                  error = function(e) { message("  Tree ", held, " failed: ", e$message); NULL })
  if (is.null(fit)) next
  # predict() gives an array: scans x 1 response x ncomp
  pred_allcomp <- tryCatch(predict(fit, newdata = test), error = function(e) NULL)
  if (is.null(pred_allcomp)) next
  loto_results_fs[[i]] <- data.frame(
    newID = held, TreeID = test$TreeID, Date = test$Date, actual = test$sugar_cor,
    matrix(pred_allcomp[, 1, ], nrow = nrow(test))
  )
  cat(sprintf("  [%d/%d] Tree %s -- %d scans held out\n", i, length(trees_fs), held, nrow(test)))
}
loto_fs <- dplyr::bind_rows(loto_results_fs)
names(loto_fs)[-(1:4)] <- paste0("pred_ncomp", seq_len(max_ncomp))


## A5: CHOOSE NCOMP AND REPORT ACCURACY

# Pick the ncomp with the lowest held-out RMSE, but replace if necessary 

rmse_by_ncomp <- sapply(seq_len(max_ncomp), function(nc) {
  sqrt(mean((loto_fs$actual - loto_fs[[paste0("pred_ncomp", nc)]])^2, na.rm = TRUE))
})
best_ncomp_fs <- which.min(rmse_by_ncomp)
cat("\nRMSE by number of components (tree-grouped LOTO-CV):\n")
print(round(rmse_by_ncomp, 3))
cat("Best number of components:", best_ncomp_fs, "\n")

loto_fs$pls_pred <- loto_fs[[paste0("pred_ncomp", best_ncomp_fs)]]

r2_pls_fs <- cor(loto_fs$actual, loto_fs$pls_pred, use = "complete.obs")^2
rmse_pls_fs <- sqrt(mean((loto_fs$actual - loto_fs$pls_pred)^2, na.rm = TRUE))
baseline_rmse_fs <- sqrt(mean((loto_fs$actual - mean(loto_fs$actual))^2))  # always guess the mean

cat("\n=== Fieldspec-only PLSR, scan-level, tree-grouped LOTO-CV results ===\n")
cat("R² (scan-level):             ", round(r2_pls_fs, 3), "\n")
cat("RMSE (scan-level):           ", round(rmse_pls_fs, 2), "\n")
cat("Baseline RMSE (grand mean):  ", round(baseline_rmse_fs, 2), "\n")
cat("Improvement over baseline:   ", round(baseline_rmse_fs - rmse_pls_fs, 2),
    "(", round(100 * (1 - rmse_pls_fs / baseline_rmse_fs), 1), "% reduction)\n")

## ---- final model, fit on all scans -- used only for the loading-weights
## figure (no single CV fold's model is "the" model). Its accuracy is the
## CV accuracy above. ----
plsr_fit_fs_final <- plsr(plsr_formula_fs, data = sugar_scans_vnir, ncomp = best_ncomp_fs, scale = TRUE)

#################################################
## A6: FIGURES


## ---- predicted vs measured (dashed line = perfect prediction) ----
fig_fs_pred <- ggplot(loto_fs, aes(x = actual, y = pls_pred)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.5, color = "#e76f51", size = 1.5) +
  labs(title = "Fieldspec-based PLSR: predicted vs actual sugar, tree-grouped LOTO-CV",
       subtitle = paste0("R² = ", round(r2_pls_fs, 3), " | RMSE = ", round(rmse_pls_fs, 2),
                         " | n = ", nrow(loto_fs), " scans, ", dplyr::n_distinct(loto_fs$newID),
                         " trees | VNIR only (500-900nm)", if (use_hyper_bands) ", airborne bands" else ""),
       x = "Measured sugar", y = "LOTO-CV predicted sugar") +
  theme_bw(base_size = 12)
print(fig_fs_pred)
#ggsave(paste0(out_folder_hyper, "FIGURE_fieldspec_plsr_sugar", band_tag, "_predicted_vs_actual.png"),
#       fig_fs_pred, width = 8, height = 6, dpi = 300)

# This is boxy as all get out 

## ---- loading weights: which bands drive component 1 ----
# Bigger distance from 0 = more influence. 
loadings_df_fs <- data.frame(
  band       = fieldspec_vnir_cols,
  wavelength = fieldspec_vnir_wl,
  loading    = as.numeric(loading.weights(plsr_fit_fs_final)[, 1])
) %>% dplyr::arrange(dplyr::desc(abs(loading)))

cat("\nTop 15 VNIR bands by |PLS loading weight| (component 1, final model on all scans):\n")
print(head(loadings_df_fs, 15))

fig_fs_loadings <- ggplot(loadings_df_fs %>% dplyr::arrange(wavelength), aes(x = wavelength, y = loading)) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_line(color = "#264653", linewidth = 0.6) +
  labs(title = "Fieldspec PLSR loading weights by wavelength (VNIR)",
       subtitle = "Component 1 -- higher |weight| = more influence on the sugar prediction",
       x = "Wavelength (nm)", y = "Loading weight") +
  theme_bw(base_size = 12)

fig_fs_loadings

## ---- top bands vs measured sugar, one point per scan ----
top_n_bands  <- 10  # <<< TOGGLE
top_bands_fs <- loadings_df_fs %>% dplyr::slice_head(n = top_n_bands) %>% dplyr::pull(band)

top_band_scatter_data <- sugar_scans_vnir %>%
  dplyr::select(TreeID, sugar_cor, dplyr::all_of(top_bands_fs)) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(top_bands_fs), names_to = "band", values_to = "reflectance") %>%
  dplyr::mutate(wl = as.numeric(gsub("x", "", band)),
                wl_label = reorder(paste0(wl, " nm"), wl))

fig_fs_top_bands <- ggplot(top_band_scatter_data, aes(x = sugar_cor, y = reflectance)) +
  geom_point(alpha = 0.4, size = 1, color = "#2a9d8f") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#e76f51", linewidth = 0.7) +
  facet_wrap(~ wl_label, scales = "free_y") +
  labs(title = paste0("Individual fruit: top ", top_n_bands, " PLSR-important bands"),
       subtitle = "One point per fruit scan",
       x = "Sugar (measured, individual fruit)", y = "Reflectance") +
  theme_bw(base_size = 11)
print(fig_fs_top_bands)

#ggsave(paste0(out_folder_hyper, "FIGURE_fieldspec_individualfruit_topbands_vs_sugar.png"),
#       fig_fs_top_bands, width = 10, height = 7, dpi = 300)


## A7: SAVE

saveRDS(plsr_fit_fs_final, paste0(out_folder_hyper, "plsr_sugar_fieldspec", band_tag, "_season1.rds"))
saveRDS(loto_fs,           paste0(out_folder_hyper, "loto_sugar_fieldspec", band_tag, "_scanlevel_season1.rds"))
saveRDS(loadings_df_fs,    paste0(out_folder_hyper, "plsr_sugar_fieldspec", band_tag, "_loadings_season1.rds"))
saveRDS(sugar_scans_vnir,  paste0(out_folder_hyper, "sugar_scans_vnir_all", band_tag, "_season1.rds"))




####################################################################
## PART B: FULL RANGE (VNIR + SWIR) -- sugar or protein
## Run once with nutrient_col = "sugar_cor", then again with the
## protein column. Output files are named by nutrient, so nothing
## is overwritten.
####################################################################
####################################################################


## B0: SETTINGS

nutrient_col <- "sugar_cor"  # <<< TOGGLE: "sugar_cor" or "Crude.Protein" (with a dot -- check names(dataCLrip))
#     (check the name printed after loading)
wl_min     <- 450    # <<< TOGGLE: shortest wavelength used (nm)
wl_max     <- 2300   # <<< TOGGLE: longest wavelength used (nm)
drop_water <- TRUE   # <<< TOGGLE: drop atmospheric water absorption bands
step_nm    <- 1      # <<< TOGGLE: 1 = every band; 5 or 10 = thinned (faster). Keep 1 when use_hyper_bands = TRUE
use_snv    <- TRUE   # <<< TOGGLE: SNV-normalize each spectrum (see B3)
splice_correct <- TRUE  # <<< TOGGLE: fix detector jumps at ~1000 and ~1830 nm (see B2)

match_airborne_bands <- TRUE  # <<< TOGGLE: also drop the airborne bad bands (900-1020 nm)



ncomp_tol  <- 0.05   # <<< TOGGLE: pick the smallest ncomp within 5% of the lowest RMSE
#     (simpler models generalize better; 0 = just take the lowest)
max_ncomp_cap <- 20  # <<< TOGGLE: most components ever tried
round_col  <- "Round"  # <<< TOGGLE: name of the collection-round column (used in Part C)

nutrient_name <- sub("_cor$", "", nutrient_col)   # e.g. "sugar"

# Short label of the settings, added to output file names so different
# runs don't overwrite each other
run_tag <- paste0(nutrient_name, "_fullrange", if (use_snv) "_snv" else "_raw",
                  if (match_airborne_bands) "_airbornebands" else "", band_tag)


########## B1: NUTRIENT SCANS -- same filters as Part A, for the chosen nutrient
#could be fiber or protein 


nut_scans <- dataCLrip[dataCLrip$type == "Fruit" & dataCLrip$spot == 2 &
                         !is.na(dataCLrip[[nutrient_col]]) & dataCLrip[[nutrient_col]] != 0, ]
nut_scans <- merge(nut_scans, newid_to_oldid, by.x = "TreeID", by.y = "TreeID2324")
cat("\n", nutrient_name, "scans:", nrow(nut_scans), "| trees:", dplyr::n_distinct(nut_scans$newID), "\n")

####################
## B2: ASD DETECTOR SPLICE CORRECTION


# The FieldSpec has three detectors joined at ~1000 nm (VNIR/SWIR1) and
# ~1830 nm (SWIR1/SWIR2), which leaves small steps in the spectrum. This
# smooths them out. It needs bands on both sides of each join, so it runs
# on the full 350-2500 nm spectrum before we drop any bands.

if (splice_correct) {
  nut_scans <- nut_scans[complete.cases(nut_scans[, fieldspec_band_cols]), ]
  spec_raw  <- as.matrix(nut_scans[, fieldspec_band_cols])
  nut_scans[, fieldspec_band_cols] <- prospectr::spliceCorrection(spec_raw, wav = fieldspec_wl_all,
                                                                  splice = c(1000, 1830), interpol.bands = 10)
  cat("Splice correction applied at 1000 and 1830 nm\n")
}


## B3: BAND SELECTION + SNV


keep <- fieldspec_wl_all >= wl_min & fieldspec_wl_all <= wl_max
if (drop_water) {
  # water in the air absorbs these bands, so they are mostly noise
  keep <- keep &
    !(fieldspec_wl_all >= 1350 & fieldspec_wl_all <= 1460) &
    !(fieldspec_wl_all >= 1790 & fieldspec_wl_all <= 1960)
}
if (match_airborne_bands) {
  keep <- keep & !(fieldspec_wl_all >= 900 & fieldspec_wl_all <= 1020)
}
keep <- keep & (fieldspec_wl_all %% step_nm == 0)
keep <- keep & in_hyper   # airborne bands only (see the "Match to the airborne hyperspectral bands" section above)

fs_cols <- fieldspec_band_cols[keep]
fs_wl   <- fieldspec_wl_all[keep]
cat("Bands kept:", length(fs_cols), "| range:", min(fs_wl), "-", max(fs_wl), "nm",
    if (drop_water) "(water bands removed)" else "", "\n")

## SNV (standard normal variate): center and scale each spectrum by its own
## mean and SD across the bands used. This removes overall brightness
## differences (distance, angle, fruit surface), so the model has to use the
## shape of the spectrum. Any new data the model is applied to (e.g. airborne
## spectra resampled to these bands) must get SNV over the  band set.


snv <- function(m) {
  m <- as.matrix(m)
  t(apply(m, 1, function(x) (x - mean(x)) / sd(x)))
}

## ---- training table: the spectra go in as a matrix column called
## 'spec', so the formula is simply nutrient ~ spec.


nut_scans_full <- nut_scans[complete.cases(nut_scans[, fs_cols]), ]

pls_df_full <- data.frame(
  newID    = nut_scans_full$newID,
  TreeID   = nut_scans_full$TreeID,
  Date     = nut_scans_full$Date,
  round    = as.character(nut_scans_full[[round_col]]),  # for Part C
  nutrient = nut_scans_full[[nutrient_col]]
)
spec_mat_full <- as.matrix(nut_scans_full[, fs_cols])
if (use_snv) {
  spec_mat_full <- snv(spec_mat_full)
  colnames(spec_mat_full) <- fs_cols
  cat("SNV applied over", length(fs_cols), "bands\n")
}
pls_df_full$spec <- I(spec_mat_full)

cat("PLSR training table:", nrow(pls_df_full), "scans,",
    dplyr::n_distinct(pls_df_full$newID), "trees\n")

#####################################
## B4: LEAVE-ONE-TREE-OUT CV (same idea as before)


max_ncomp_full <- min(max_ncomp_cap, length(fs_cols) - 1, dplyr::n_distinct(pls_df_full$newID) - 2)
cat("Max PLS components tried:", max_ncomp_full, "\n")

trees_full <- unique(pls_df_full$newID)
loto_results_full <- vector("list", length(trees_full)); names(loto_results_full) <- trees_full

for (i in seq_along(trees_full)) {
  held  <- trees_full[i]
  train <- pls_df_full[pls_df_full$newID != held, ]
  test  <- pls_df_full[pls_df_full$newID == held, ]
  if (nrow(train) < 10 || nrow(test) < 1) next
  fit <- tryCatch(plsr(nutrient ~ spec, data = train, ncomp = max_ncomp_full, scale = TRUE),
                  error = function(e) { message("  Tree ", held, " failed: ", e$message); NULL })
  if (is.null(fit)) next
  pred_allcomp <- tryCatch(predict(fit, newdata = test), error = function(e) NULL)
  if (is.null(pred_allcomp)) next
  loto_results_full[[i]] <- data.frame(
    newID = held, TreeID = test$TreeID, Date = test$Date, round = test$round,
    actual = test$nutrient,
    matrix(pred_allcomp[, 1, ], nrow = nrow(test))
  )
  cat(sprintf("  [%d/%d] Tree %s -- %d scans held out\n", i, length(trees_full), held, nrow(test)))
}

loto_full <- dplyr::bind_rows(loto_results_full)
names(loto_full)[-(1:5)] <- paste0("pred_ncomp", seq_len(max_ncomp_full))


## B5: CHOOSE NCOMP AND REPORT ACCURACY

rmse_by_ncomp_full <- sapply(seq_len(max_ncomp_full), function(nc) {
  sqrt(mean((loto_full$actual - loto_full[[paste0("pred_ncomp", nc)]])^2, na.rm = TRUE))
})

rmse_by_ncomp_full

min_ncomp_full  <- which.min(rmse_by_ncomp_full)
best_ncomp_full <- which(rmse_by_ncomp_full <= min(rmse_by_ncomp_full) * (1 + ncomp_tol))[1]
cat("\nRMSE by ncomp (tree-grouped LOTO-CV):\n"); print(round(rmse_by_ncomp_full, 3))
cat("Minimum-RMSE ncomp:", min_ncomp_full, "| chosen (within", 100 * ncomp_tol, "% of min):",
    best_ncomp_full, "\n")
if (min_ncomp_full == max_ncomp_full)
  cat("NOTE: minimum is at the cap -- consider raising max_ncomp_cap.\n")

#### RMSE curve figure (dashed line = chosen ncomp) 
rmse_curve_df <- data.frame(ncomp = seq_along(rmse_by_ncomp_full), rmse = rmse_by_ncomp_full)
fig_rmse_curve <- ggplot(rmse_curve_df, aes(ncomp, rmse)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = best_ncomp_full, linetype = "dashed", color = "#e76f51") +
  labs(title = paste0("LOTO-CV RMSE by number of components (", run_tag, ")"),
       subtitle = paste0("Dashed: chosen ncomp = ", best_ncomp_full),
       x = "Number of PLS components", y = "RMSE") +
  theme_bw(base_size = 12)
print(fig_rmse_curve)
#ggsave(paste0(out_folder_hyper, "FIGURE_fieldspec_", run_tag, "_plsr_rmse_curve.png"),
#       fig_rmse_curve, width = 7, height = 4.5, dpi = 300)

loto_full$pls_pred <- loto_full[[paste0("pred_ncomp", best_ncomp_full)]]
r2_full   <- cor(loto_full$actual, loto_full$pls_pred, use = "complete.obs")^2
rmse_full <- sqrt(mean((loto_full$actual - loto_full$pls_pred)^2, na.rm = TRUE))
base_full <- sqrt(mean((loto_full$actual - mean(loto_full$actual))^2))

cat("\n=== Fieldspec PLSR,", nutrient_name, "(", min(fs_wl), "-", max(fs_wl), "nm), tree-grouped LOTO-CV ===\n")
cat("R² (scan-level):            ", round(r2_full, 3), "\n")
cat("RMSE (scan-level):          ", round(rmse_full, 2), "\n")
cat("Baseline RMSE (grand mean): ", round(base_full, 2), "\n")
cat("Reduction vs baseline:      ", round(100 * (1 - rmse_full / base_full), 1), "%\n")

## $$$final model on all scans (for loadings / figures)
plsr_fit_full_final <- plsr(nutrient ~ spec, data = pls_df_full,
                            ncomp = best_ncomp_full, scale = TRUE)


## B6: FIGURES


range_lab <- paste0(min(fs_wl), "-", max(fs_wl), "nm", if (drop_water) ", water bands removed" else "",
                    if (match_airborne_bands) ", 900-1020nm removed" else "", if (use_snv) ", SNV" else "",
                    if (use_hyper_bands) ", airborne bands" else "")

###predicted vs measured
fig_full_pred <- ggplot(loto_full, aes(x = actual, y = pls_pred)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.5, color = "#e76f51", size = 1.5) +
  labs(title = paste0("Fieldspec-based PLSR: predicted vs actual ", nutrient_name, ", tree-grouped LOTO-CV"),
       subtitle = paste0("R² = ", round(r2_full, 3), " | RMSE = ", round(rmse_full, 2),
                         " | n = ", nrow(loto_full), " scans, ",
                         dplyr::n_distinct(loto_full$newID), " trees | ", range_lab),
       x = paste("Measured", nutrient_name), y = paste("LOTO-CV predicted", nutrient_name)) +
  theme_bw(base_size = 12)
fig_full_pred # looking less boxy 


#ggsave(paste0(out_folder_hyper, "FIGURE_fieldspec_", run_tag, "_plsr_predicted_vs_actual.png"),
#       fig_full_pred, width = 8, height = 6, dpi = 300)

## loading weights (component 1) 
loadings_df_full <- data.frame(
  band       = fs_cols,
  wavelength = fs_wl,
  loading    = as.numeric(loading.weights(plsr_fit_full_final)[, 1])
)

# 'segment' gives each run of bands its own group, so the line breaks at
# the removed water bands instead of bridging them. A jump of more than
# 20 nm between neighboring bands counts as a gap (airborne bands are a
# few nm apart).
loadings_df_full$segment <- cumsum(c(1, diff(loadings_df_full$wavelength) > 20))

fig_full_loadings <- ggplot(loadings_df_full, aes(x = wavelength, y = loading, group = segment)) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_line(color = "#264653", linewidth = 0.6) +
  labs(title = paste0("Fieldspec PLSR loading weights by wavelength (VNIR + SWIR), ", nutrient_name),
       subtitle = "Component 1 -- higher |weight| = more influence on the first latent variable",
       x = "Wavelength (nm)", y = "Loading weight") +
  theme_bw(base_size = 12)
print(fig_full_loadings)
ggsave(paste0(out_folder_hyper, "FIGURE_fieldspec_", run_tag, "_plsr_loadings.png"),
       fig_full_loadings, width = 9, height = 5, dpi = 300)

loadings_df_full <- loadings_df_full %>% dplyr::arrange(dplyr::desc(abs(loading)))
cat("\nTop 15 bands by |loading weight| (component 1):\n")
print(head(loadings_df_full, 15))

## VIP scores: importance of each band across all chosen components
## (loading weights only show one component). VIP > 1 = "important". ----
vip_df <- data.frame(wavelength = fs_wl,
                     vip = plsVarSel::VIP(plsr_fit_full_final, best_ncomp_full)) %>%
  dplyr::mutate(segment = cumsum(c(1, diff(wavelength) > 20)))

fig_vip <- ggplot(vip_df, aes(wavelength, vip, group = segment)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_line(color = "#264653", linewidth = 0.6) +
  labs(title = paste0("VIP scores, ", run_tag, " (", best_ncomp_full, " components)"),
       subtitle = "Dashed line: VIP = 1 (conventional importance threshold)",
       x = "Wavelength (nm)", y = "VIP") +
  theme_bw(base_size = 12)
print(fig_vip)
ggsave(paste0(out_folder_hyper, "FIGURE_fieldspec_", run_tag, "_plsr_vip.png"),
       fig_vip, width = 9, height = 5, dpi = 300)

##  top bands vs the nutrient, one point per scan 
top_n_bands    <- 10  # <<< TOGGLE
top_bands_full <- loadings_df_full %>% dplyr::slice_head(n = top_n_bands) %>% dplyr::pull(band)

top_band_scatter_full <- data.frame(nutrient = pls_df_full$nutrient,
                                    spec_mat_full[, top_bands_full, drop = FALSE]) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(top_bands_full), names_to = "band", values_to = "reflectance") %>%
  dplyr::mutate(wl = as.numeric(gsub("x", "", band)),
                wl_label = reorder(paste0(wl, " nm"), wl))

fig_full_top_bands <- ggplot(top_band_scatter_full, aes(x = nutrient, y = reflectance)) +
  geom_point(alpha = 0.4, size = 1, color = "#2a9d8f") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#e76f51", linewidth = 0.7) +
  facet_wrap(~ wl_label, scales = "free_y") +
  labs(title = paste0("Individual fruit: top ", top_n_bands, " PLSR-important bands (", run_tag, ")"),
       subtitle = "One point per fruit scan",
       x = paste0(nutrient_name, " (measured, individual fruit)"),
       y = if (use_snv) "SNV reflectance" else "Reflectance") +
  theme_bw(base_size = 11)
print(fig_full_top_bands)
ggsave(paste0(out_folder_hyper, "FIGURE_fieldspec_", run_tag, "_individualfruit_topbands.png"),
       fig_full_top_bands, width = 10, height = 7, dpi = 300)

##################################
## B7: SAVE

saveRDS(plsr_fit_full_final, paste0(out_folder_hyper, "plsr_fieldspec_", run_tag, "_season1.rds"))
saveRDS(loto_full,           paste0(out_folder_hyper, "loto_fieldspec_", run_tag, "_scanlevel_season1.rds"))
saveRDS(loadings_df_full,    paste0(out_folder_hyper, "plsr_fieldspec_", run_tag, "_loadings_season1.rds"))
saveRDS(pls_df_full,         paste0(out_folder_hyper, "scans_", run_tag, "_season1.rds"))


####################################################################
####################################################################
## PART C: ROUND 3 ONLY -- the round closest to the airborne flight
## Run right after Part B: uses pls_df_full and the Part B settings
## (same nutrient, bands, SNV). Fits a new model on round-3 scans only.
####################################################################
####################################################################

airborne_round <- "3"  
table(pls_df_full$round)  

pls_df_r3 <- pls_df_full[pls_df_full$round == airborne_round, ]
cat("Round", airborne_round, ":", nrow(pls_df_r3), "scans,",
    dplyr::n_distinct(pls_df_r3$newID), "trees | dates:",
    format(min(pls_df_r3$Date)), "to", format(max(pls_df_r3$Date)), "\n")

## ---- how does the Part B model (trained on all rounds) do on round 3? ----
loto_full_r3 <- loto_full[loto_full$round == airborne_round, ]
cat("Part B (all-rounds) model, round-3 scans only: R² =",
    round(cor(loto_full_r3$actual, loto_full_r3$pls_pred)^2, 3), "| RMSE =",
    round(sqrt(mean((loto_full_r3$actual - loto_full_r3$pls_pred)^2)), 2), "\n")

## =====================================================================
## C1: LEAVE-ONE-TREE-OUT CV, round 3 only
## =====================================================================
# Fewer trees than Part B, so fewer components are possible

max_ncomp_r3 <- min(max_ncomp_cap, dplyr::n_distinct(pls_df_r3$newID) - 2, nrow(pls_df_r3) - 2)
cat("Max PLS components tried:", max_ncomp_r3, "\n")

trees_r3 <- unique(pls_df_r3$newID)
loto_results_r3 <- vector("list", length(trees_r3)); names(loto_results_r3) <- trees_r3

for (i in seq_along(trees_r3)) {
  held  <- trees_r3[i]
  train <- pls_df_r3[pls_df_r3$newID != held, ]
  test  <- pls_df_r3[pls_df_r3$newID == held, ]
  if (nrow(train) < 10 || nrow(test) < 1) next
  fit <- tryCatch(plsr(nutrient ~ spec, data = train, ncomp = max_ncomp_r3, scale = TRUE),
                  error = function(e) { message("  Tree ", held, " failed: ", e$message); NULL })
  if (is.null(fit)) next
  pred_allcomp <- tryCatch(predict(fit, newdata = test), error = function(e) NULL)
  if (is.null(pred_allcomp)) next
  loto_results_r3[[i]] <- data.frame(
    newID = held, TreeID = test$TreeID, Date = test$Date, actual = test$nutrient,
    matrix(pred_allcomp[, 1, ], nrow = nrow(test))
  )
}

loto_r3 <- dplyr::bind_rows(loto_results_r3)
names(loto_r3)[-(1:4)] <- paste0("pred_ncomp", seq_len(max_ncomp_r3))


## C2: CHOOSE NCOMP AND REPORT ACCURACY 


rmse_by_ncomp_r3 <- sapply(seq_len(max_ncomp_r3), function(nc) {
  sqrt(mean((loto_r3$actual - loto_r3[[paste0("pred_ncomp", nc)]])^2, na.rm = TRUE))
})
best_ncomp_r3 <- which(rmse_by_ncomp_r3 <= min(rmse_by_ncomp_r3) * (1 + ncomp_tol))[1]
cat("\nRMSE by ncomp (round", airborne_round, "LOTO-CV):\n"); print(round(rmse_by_ncomp_r3, 3))
cat("Chosen ncomp:", best_ncomp_r3, "\n")

loto_r3$pls_pred <- loto_r3[[paste0("pred_ncomp", best_ncomp_r3)]]
r2_r3   <- cor(loto_r3$actual, loto_r3$pls_pred)^2
rmse_r3 <- sqrt(mean((loto_r3$actual - loto_r3$pls_pred)^2))
base_r3 <- sqrt(mean((loto_r3$actual - mean(loto_r3$actual))^2))

cat("\n=== Round", airborne_round, "only,", nutrient_name, ", tree-grouped LOTO-CV ===\n")
cat("R² (scan-level):            ", round(r2_r3, 3), "\n")
cat("RMSE (scan-level):          ", round(rmse_r3, 2), "\n")
cat("Baseline RMSE (round mean): ", round(base_r3, 2), "\n")
cat("Reduction vs baseline:      ", round(100 * (1 - rmse_r3 / base_r3), 1), "%\n")

## final round-3 model on all round-3 scans (e.g. to apply to the
## airborne image -- remember to resample + SNV the image spectra over the same bands first)
plsr_fit_r3_final <- plsr(nutrient ~ spec, data = pls_df_r3, ncomp = best_ncomp_r3, scale = TRUE)


## C3: FIGURE + SAVE
fig_r3_pred <- ggplot(loto_r3, aes(x = actual, y = pls_pred)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.5, color = "#e76f51", size = 1.5) +
  labs(title = paste0("Round ", airborne_round, " only: predicted vs actual ", nutrient_name,
                      ", tree-grouped LOTO-CV"),
       subtitle = paste0("R² = ", round(r2_r3, 3), " | RMSE = ", round(rmse_r3, 2),
                         " | n = ", nrow(loto_r3), " scans, ", dplyr::n_distinct(loto_r3$newID),
                         " trees | ", range_lab),
       x = paste("Measured", nutrient_name), y = paste("LOTO-CV predicted", nutrient_name)) +
  theme_bw(base_size = 12)
print(fig_r3_pred)

#ggsave(paste0(out_folder_hyper, "FIGURE_fieldspec_", run_tag, "_round", airborne_round, "_plsr_predicted_vs_actual.png"),
#       fig_r3_pred, width = 8, height = 6, dpi = 300)

saveRDS(plsr_fit_r3_final, paste0(out_folder_hyper, "plsr_fieldspec_", run_tag, "_round", airborne_round, "_season1.rds"))
saveRDS(loto_r3,           paste0(out_folder_hyper, "loto_fieldspec_", run_tag, "_round", airborne_round, "_season1.rds"))