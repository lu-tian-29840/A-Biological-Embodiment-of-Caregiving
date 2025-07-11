#@@ Author: Lu Tian
#@@ Date: May 16, 2025
#@@ Project name: Caregiving and bioage
#    Section 4/4  : Part 1: Forest plot for caregiving - bioage
#                   Part 2: Forest plot for caregiving - bioage, stratified by gender 
#                   Part 3: Weighting analysis table

####################################################################################################
# 
#                           Part 1: Forest plot for caregiving - bioage
#
####################################################################################################

# ------------------------------------------------------------
# Pair‑wise ATT GBM-IPTW + Three‑Panel Forest Plot
# ------------------------------------------------------------
# Requirements ------------------------------------------------
library(dplyr)
library(purrr)     # imap_dfr()
library(broom)     # tidy()
library(survey)    # svyglm() 
library(ggplot2)

# ------------------------------------------------------------
# 0.  Inputs exist in environment -----------------
# dat_IH, dat_IM, dat_CE, dat_su_med, dat_su_hi  – analytic
#   data frames, each containing its own binary treat_* dummy.
# design_IH, design_IM, design_CE, design_su_med, design_su_hi
#   – corresponding svydesign objects (combined ATT × sample weights).
# ------------------------------------------------------------

# 1.  Model set‑up -------------------------------------------
# --- Original covariates ---
cov_w11 <- c("r11age", "female", "black", "hisp", "other",
             "r11marcat_b", "edu", "h11itot1", "r11shlt",
             "r11smoke", "r11drinkn", "r11mobila", "r11chronic",
             "r11acts_w", "r11bmi")

# --- Treatment dummies & design map ---
trt_names <- c("treat_IM", "treat_IH", "treat_CE",
               "treat_su_med", "treat_su_hi")

svy_designs <- list(
  treat_IM     = design_IM,
  treat_IH     = design_IH,
  treat_CE     = design_CE,
  treat_su_med = design_su_med,
  treat_su_hi  = design_su_hi
)

# --- Outcomes (panels A–C) ---
dv_vec <- c("grim_residual1_z", "vidalbralo_residual1_z", "zhang_residual1_z")


# --------------------------------------------------------------------------------------------------
# 1. Mapping from treatment dummy to descriptive label
# --------------------------------------------------------------------------------------------------

treatment_labels <- c(
  treat_IM     = "Transition to medium intensity care",
  treat_IH     = "Transition to high intensity care",
  treat_CE     = "Transition out of caregiving",
  treat_su_med = "Sustained at (or to) medium intensity",
  treat_su_hi  = "Sustained at (or to) high intensity"
)

# --------------------------------------------------------------------------------------------------
# 2. Check the order of labels 
# --------------------------------------------------------------------------------------------------
# Factor levels for the y‑axis so that rows appear in the
# desired order (top‑to‑bottom in the plot facet):
#   1) Transition to medium intensity care
#   2) Transition to high intensity care
#   3) Transition out of caregiving
#   4) Sustained at (or to) medium intensity
#   5) Sustained at (or to) high intensity
# ------------------------------------------------------------

treatment_order <- c(
  "Transition to medium intensity care",
  "Transition to high intensity care",
  "Transition out of caregiving",
  "Sustained at (or to) medium intensity",
  "Sustained at (or to) high intensity"
)

# -------------------------------------------------------------------------------------------------
# 3.  Forest plot for caregiving ~ bioage. 
# --------------------------------------------------------------------------------------------------

# Storage lists ------------------------------------------------
models       <- list()
formulas     <- list()
tidy_results <- list()

# 3a. Fit survey‑weighted doubly‑robust models --------------
for (treat in trt_names) {
  rhs <- paste(c(treat, cov_w11), collapse = " + ")
  des <- svy_designs[[treat]]
  models[[treat]]   <- list()
  formulas[[treat]] <- list()
  tidy_results[[treat]] <- list()
  
  for (dv in dv_vec) {
    fml <- as.formula(paste(dv, "~", rhs))
    fit <- svyglm(fml, design = des)
    
    models[[treat]][[dv]]   <- fit
    formulas[[treat]][[dv]] <- fml
    
    tidy_results[[treat]][[dv]] <- tidy(fit, conf.int = TRUE)
  }
}

# Build nested data_list matching each design --------------
data_list <- list(
  treat_IM     = dat_IM,
  treat_IH     = dat_IH,
  treat_CE     = dat_CE,
  treat_su_med = dat_su_med,
  treat_su_hi  = dat_su_hi
)

    # Extract ATT rows & assemble plotting tibble --------------
    plot_df <- imap_dfr(tidy_results, function(outcome_list, treat) {
      dat <- data_list[[treat]]
      imap_dfr(outcome_list, function(tidy_tbl, dv) {
        att <- tidy_tbl %>% filter(term == treat)
        tibble(
          outcome   = dv,
          treatment = treatment_labels[[treat]],
          estimate  = att$estimate,
          lci       = att$conf.low,
          uci       = att$conf.high,
          pval      = att$p.value,
          n_treated = sum(dat[[treat]] == 1L, na.rm = TRUE)
        )
      })
    })

    # Factor ordering & label adjustments -----------------------
    plot_df <- plot_df %>%
      mutate(
        outcome = factor(outcome,
                         levels = c("grim_residual1_z",
                                    "vidalbralo_residual1_z",
                                    "zhang_residual1_z")),
        treatment = factor(treatment, levels = rev(treatment_order)),
        label_n = ifelse(outcome == "zhang_residual1_z",
                         paste0("n = ", n_treated), NA_character_)
      )

# 3b. Create the forest plot ------------------------------------
forest_plot <- ggplot(plot_df, aes(y = treatment)) +
  geom_segment(aes(x = lci, xend = uci, yend = treatment),
               size = 7, colour = "#00BFC4", alpha = .40) +
  geom_point(aes(x = estimate), size = 3, shape = 21,
             fill = "#00BFC4", colour = "black") +
  geom_text(aes(x = uci + 0.02, label = label_n),
            hjust = 0, size = 3.5, na.rm = TRUE, show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ outcome, ncol = 3, scales = "free_x",
             strip.position = "bottom",
             labeller = as_labeller(c(
               grim_residual1_z       = "Panel A: GrimAge",
               vidalbralo_residual1_z = "Panel B: Vidal‑Bralo DNAmAge",
               zhang_residual1_z      = "Panel C: Zhang DNAmAge"))) +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.placement  = "outside",
    strip.text       = element_text(size = 12),
    strip.background = element_blank(),
    axis.text.y      = element_text(size = 12),
    axis.text.x      = element_text(size = 12),
    panel.spacing    = unit(0.8, "lines"),
    plot.margin      = margin(6, 40, 6, 6)
  )
        
        # Save -------------------------------------------------------
        
        ggsave("Figure_Bioage_caregiving.pdf",
               plot   = forest_plot,
               device = pdf, width = 11, height = 5, units = "in")
        
        ggsave("Figure_Bioage_caregiving.png", forest_plot, width = 270, height = 130, units = "mm", dpi = 300)
        

        
####################################################################################################
# 
#                     Part 2: Forest plot for caregiving - bioage stratified by gender 
#
####################################################################################################
        
        att_rows <- list()
        
        for (gender_lab in c("Female", "Male")) {
          gender_val <- ifelse(gender_lab == "Female", 1, 0)
          for (trt in trt_names) {
            rhs <- paste(c(trt, cov_w11), collapse = " + ")
            des <- subset(svy_designs[[trt]], female == gender_val)
            for (dv in dv_vec) {
              fit <- svyglm(as.formula(paste(dv, "~", rhs)), design = des)
              att <- tidy(fit, conf.int = TRUE) %>% filter(term == trt)
              att_rows[[length(att_rows) + 1]] <- tibble(
                gender    = gender_lab,
                outcome   = dv,
                treatment = treatment_labels[[trt]],  # updated label
                estimate  = att$estimate,
                lci       = att$conf.low,
                uci       = att$conf.high,
                n_treated = sum(data_list[[trt]]$female == gender_val &
                                  data_list[[trt]][[trt]] == 1, na.rm = TRUE)
              )
            }
          }
        }
        
        plot_df_g <- bind_rows(att_rows) %>%
          mutate(
            treatment = factor(treatment, levels = treatment_order),
            outcome   = factor(outcome,   levels = dv_vec),
            rank_rev  = 6 - as.numeric(treatment),   # top row = 5, bottom = 1
            offset    = if_else(gender == "Female", -0.20, 0.20),
            y_pos     = rank_rev + offset,
            label_n   = if_else(outcome == "zhang_residual1_z",
                                paste0("n = ", n_treated), NA_character_)
          )
        
        strat_plot <- ggplot(plot_df_g) +
          geom_segment(aes(x = lci, xend = uci, y = y_pos, yend = y_pos, colour = gender),
                       size = 7, alpha = .35, show.legend = FALSE) +
          geom_point(aes(x = estimate, y = y_pos, shape = gender, fill = gender),
                     size = 3, colour = "black") +
          geom_text(aes(x = uci + 0.02, y = y_pos, label = label_n, colour = gender),
                    hjust = 0, size = 3.5, na.rm = TRUE, show.legend = FALSE) +
          geom_vline(xintercept = 0, linetype = "dashed") +
          scale_colour_manual(values = c(Female = "#F8766D", Male = "#00BFC4"), guide = "none") +
          scale_fill_manual(values   = c(Female = "#F8766D", Male = "#00BFC4"), name = "Gender") +
          scale_shape_manual(values  = c(Female = 21,       Male = 22),       name = "Gender") +
          scale_x_continuous(expand = expansion(mult = c(0.05, 0.30))) +
          scale_y_continuous(breaks = 5:1, labels = treatment_order, expand = expansion(mult = c(0.05, 0.05))) +
          facet_wrap(~ outcome, ncol = 3, scales = "free_x",
                     strip.position = "bottom",
                     labeller = as_labeller(c(
                       grim_residual1_z       = "Panel A: GrimAge",
                       vidalbralo_residual1_z = "Panel B: Vidal‑Bralo DNAmAge",
                       zhang_residual1_z      = "Panel C: Zhang DNAmAge"))) +
          labs(x = NULL, y = NULL) +
          coord_cartesian(clip = "off") +
          guides(colour = "none",
                 shape  = "none",
                 fill   = guide_legend(title = "Gender",
                                       override.aes = list(shape = c(21, 22),
                                                           fill  = c("#F8766D", "#00BFC4"),
                                                           colour = "black",
                                                           size   = 4))) +
          theme_bw() +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                strip.placement  = "outside",
                strip.text       = element_text(size = 12),
                strip.background = element_blank(),
                axis.text.y      = element_text(size = 12),
                axis.text.x      = element_text(size = 10),
                panel.spacing    = unit(0.8, "lines"),
                legend.position  = "right",
                legend.text      = element_text(size = 12),
                legend.title     = element_text(size = 12),
                plot.margin      = margin(6, 40, 6, 6))
        
        # --- Save gender‑stratified figure -------------------------
        
        ggsave("Figure_Bioage_caregiving_by_gender.pdf", strat_plot,
               width = 14.5, height = 7.5, units = "in")
        
        ggsave("Figure_Bioage_caregiving_by_gender.png", strat_plot,
               width = 315, height = 160, units = "mm", dpi = 300)        

####################################################################################################
# 
#               Part 3: Weighting analysis table
#
####################################################################################################
        
        # --- Attempt to load / install qvalue --------------------------------------
        use_qvalue <- TRUE   
        if (use_qvalue) {
          if (!requireNamespace("qvalue", quietly = TRUE)) {
            message("Package 'qvalue' not found; attempting installation …")
            
            ## Step 1: ensure BiocManager is available ------------------------------
            if (!requireNamespace("BiocManager", quietly = TRUE)) {
              message("Installing 'BiocManager' from CRAN …")
              install.packages("BiocManager", repos = "https://cloud.r-project.org")
            }
            library(BiocManager)
            
            ## Try installing qvalue via Bioconductor -------------------------------
            success_bioc <- FALSE
            tryCatch({
              BiocManager::install("qvalue", ask = FALSE, update = FALSE)
              success_bioc <- TRUE
            }, error = function(e) {
              message("Bioconductor install failed: ", e$message)
            })
            
            ## Step 2: fall back to GitHub if Bioconductor fails --------------------
            if (!success_bioc) {
              message("Attempting GitHub installation via devtools …")
              if (!requireNamespace("devtools", quietly = TRUE)) {
                install.packages("devtools", repos = "https://cloud.r-project.org")
              }
              devtools::install_github("jdstorey/qvalue")
            }
          }
          ## Load qvalue if installation succeeded ----------------------------------
          if (requireNamespace("qvalue", quietly = TRUE)) {
            library(qvalue)
          } else {
            warning("Adaptive q‑values will be skipped: 'qvalue' could not be installed.")
            use_qvalue <- FALSE
          }
        }
        
        # 1. Outcome and covariate lists --------------------------------------------
        dv_list <- c("horvath_residual1_z", "hannum_residual1_z", "levine_residual1_z",
                     "horvathskin_residual1_z", "lin_residual1_z", "weidner_residual1_z",
                     "vidalbralo_residual1_z", "grim_residual1_z", "zhang_residual1_z",
                     "yang_residual1_z", "bocklandt_residual1_z", "garagnani_residual1_z",
                     "mpoa_z")
        
        cov_w11 <- c("r11age", "female", "racecat", "r11marcat_b", "edu", "h11itot1",
                     "r11shlt", "r11smoke", "r11drinkn", "r11mobila", "r11chronic",
                     "r11acts_w", "r11bmi")
        
        # 2. Treatment groups --------------------------------------------------------
        treatment_groups <- c("Initiated (High)", "Ceased", "Initiated (Medium)",
                              "Sustained (to Medium)", "Sustained (to High)")
        
        # 3. Storage object ----------------------------------------------------------
        results_list <- list()
        
        # 4. Main loop ---------------------------------------------------------------
        for (treatment in treatment_groups) {
          cat(sprintf("\n\n--- Starting Analysis for: %s vs. No Caregiving ---\n", treatment))
          
          ## 4.1 Data subset & treatment flag ---------------------------------------
          dat_temp <- subset(df3, care_trans %in% c(treatment, "No Caregiving"))
          dat_temp$treat <- ifelse(dat_temp$care_trans == treatment, 1L, 0L)
          
          ## 4.2 GBM ATT weights -----------------------------------------------------
          set.seed(123)  # reproducible GBM each loop
          w.out_temp <- weightit(
            as.formula(paste("treat ~", paste(cov_w11, collapse = " + "))),
            data              = dat_temp,
            method            = "gbm",
            estimand          = "ATT",
            n.trees           = 2500,
            interaction.depth = 2,
            shrinkage         = 0.01,
            criterion         = "smd.max",
            abs               = TRUE
          )
          
          ## 4.3 Combine with survey weight -----------------------------------------
          dat_temp$w_att  <- get.w(w.out_temp)
          dat_temp$comb_w <- dat_temp$w_att * dat_temp$vbsi16wgtra
          dat_temp <- subset(dat_temp, !is.na(comb_w) & comb_w > 0)
          
          ## 4.4 Survey design -------------------------------------------------------
          design_temp <- svydesign(ids = ~secu, strata = ~stratum,
                                   weights = ~comb_w, data = dat_temp, nest = TRUE)
          
          rhs_dr_temp <- paste("treat +", paste(cov_w11, collapse = " + "))
          
          ## 4.5 Doubly‑robust models ----------------------------------------------
          for (dv in dv_list) {
            if (!dv %in% names(dat_temp)) {
              warning(sprintf("Variable %s not found; skipping.", dv))
              next
            }
            
            dr_model <- svyglm(as.formula(paste(dv, "~", rhs_dr_temp)),
                               design = design_temp, na.action = na.omit)
            
            coef_row <- summary(dr_model)$coefficients
            if ("treat" %in% rownames(coef_row)) {
              coef_info <- coef_row["treat", ]
              results_list[[length(results_list) + 1L]] <- data.frame(
                comparison = sprintf("%s vs. No Caregiving", treatment),
                outcome    = dv,
                estimate   = coef_info["Estimate"],
                std.error  = coef_info["Std. Error"],
                t.value    = coef_info["t value"],
                p.value    = coef_info["Pr(>|t|)"]
              )
            }
          }
        }
        
        # 5. Compile master results --------------------------------------------------
        if (length(results_list) == 0L) stop("No models fitted – check variable names and data.")
        
        results_df <- do.call(rbind, results_list) |> 
          mutate(across(where(is.numeric), round, 6))
        
        # 6. False‑Discovery‑Rate corrections ---------------------------------------
        ## 6.1 Global BH (65 tests)
        results_df$q_bh_global <- round(p.adjust(results_df$p.value, method = "BH"), 4)
        
        ## 6.2 Within‑contrast BH (5 × 13 tests)
        results_df <- results_df |> 
          group_by(comparison) |> 
          mutate(q_bh_within = round(p.adjust(p.value, method = "BH"), 4)) |> 
          ungroup()
        
        ## 6.3 Benjamini‑Yekutieli (global, arbitrary dependence)
        results_df$q_by_global <- round(p.adjust(results_df$p.value, method = "BY"), 4)
        
        ## 6.4 Storey–Tibshirani adaptive q‑values (optional)
        if (use_qvalue && requireNamespace("qvalue", quietly = TRUE)) {
          results_df$q_storey <- round(qvalue::qvalue(results_df$p.value)$qvalues, 4)
        }
        
        # 7. Final rounding for presentation ----------------------------------------
        num_cols <- c("estimate", "std.error", "t.value", "p.value",
                      "q_bh_global", "q_bh_within", "q_by_global",
                      if (use_qvalue && "q_storey" %in% names(results_df)) "q_storey")
        results_df[num_cols] <- lapply(results_df[num_cols], round, 4)
        
        # 8. Export ------------------------------------------------------------------
        write.csv(results_df, "doubly_robust_att_results_FDR.csv", row.names = FALSE)
        print(results_df)
        
        
        # --------------------------------------------------------------------------------------------------
        #                  Doubly‑Robust ATT Analysis – **Gender‑Stratification
        #
        # --------------------------------------------------------------------------------------------------
        
        
        # 1. Lists of variables -----------------------------------------------------
        dv_list <- c("horvath_residual1_z","hannum_residual1_z","levine_residual1_z",
                     "horvathskin_residual1_z","lin_residual1_z","weidner_residual1_z",
                     "vidalbralo_residual1_z","grim_residual1_z","zhang_residual1_z",
                     "yang_residual1_z","bocklandt_residual1_z","garagnani_residual1_z",
                     "mpoa_z")
        
        cov_w11 <- c("r11age","female","racecat","r11marcat_b","edu","h11itot1",
                     "r11shlt","r11smoke","r11drinkn","r11mobila","r11chronic",
                     "r11acts_w","r11bmi")
        
        treatment_groups <- c("Initiated (High)","Ceased","Initiated (Medium)",
                              "Sustained (to Medium)","Sustained (to High)")
        
        # 2. Build full‑sample survey designs for each caregiving transition --------
        message("\n=== Building full‑sample GBM weights and survey designs ===")
        svy_designs <- list()
        for (treatment in treatment_groups) {
          message("  → ", treatment, " vs. No Caregiving")
          dat <- subset(df3, care_trans %in% c(treatment, "No Caregiving"))
          dat$treat <- ifelse(dat$care_trans == treatment, 1L, 0L)
          
          set.seed(123)
          w.out <- weightit(treat ~ ., data = dat[c("treat", cov_w11)],
                            method = "gbm", estimand = "ATT",
                            n.trees = 2500, interaction.depth = 2,
                            shrinkage = 0.01, criterion = "smd.max", abs = TRUE)
          
          dat$w_att  <- get.w(w.out)
          dat$comb_w <- dat$w_att * dat$vbsi16wgtra   # grand weight à la Kim 2025
          dat <- subset(dat, comb_w > 0 & !is.na(comb_w))
          
          svy_designs[[treatment]] <- svydesign(ids = ~secu, strata = ~stratum,
                                                weights = ~comb_w, data = dat, nest = TRUE)
        }
        
        # 3. Gender‑stratified outcome models using the grand weights --------------
        gender_levels <- c(Male = 0, Female = 1)
        results_list   <- list()
        
        for (g_label in names(gender_levels)) {
          g_val <- gender_levels[[g_label]]
          for (treatment in treatment_groups) {
            des_sub <- subset(svy_designs[[treatment]], female == g_val)
            rhs     <- paste("treat +", paste(cov_w11, collapse = " + "))
            
            for (dv in dv_list) {
              if (!dv %in% names(des_sub$variables)) next
              fit <- svyglm(as.formula(paste(dv, "~", rhs)), design = des_sub, na.action = na.omit)
              co  <- summary(fit)$coefficients["treat", ]
              results_list[[length(results_list)+1]] <- data.frame(
                gender     = g_label,
                comparison = sprintf("%s vs. No Caregiving", treatment),
                outcome    = dv,
                estimate   = round(co["Estimate"], 4),
                std.error  = round(co["Std. Error"], 4),
                p.value    = round(co["Pr(>|t|)"], 4)
              )
            }
          }
        }
        
        results_df <- do.call(rbind, results_list)
        
        # 4. FDR corrections --------------------------------------------------------
        results_df$q_bh_global <- round(p.adjust(results_df$p.value, method="BH"), 4)
        results_df <- results_df %>% group_by(gender, comparison) %>%
          mutate(q_bh_within = round(p.adjust(p.value, method="BH"), 4)) %>% ungroup()
        results_df$q_by_global <- round(p.adjust(results_df$p.value, method="BY"), 4)
        if (use_qvalue) results_df$q_storey <- round(qvalue::qvalue(results_df$p.value)$qvalues, 4)
        
        # 5. Export ---------------------------------------------------------------
        write.csv(results_df, "doubly_robust_att_results_gender_FDR.csv", row.names = FALSE)
        print(results_df)
