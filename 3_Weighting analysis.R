
#@@ Author: Lu Tian
#@@ Date: May 16, 2025
#@@ Project name: Caregiving and bioage
#    Section 3/4  : Run GBM-IPTW diagnostics for approach selection
#           Part I:  Pair-wise ATT diagnostics for Initiated-High vs No-Care 
#                    Step 1: Estimate GBM propensity weights for ATT
#                    Step 2: Diagnostic Pipeline
#                           1) Balance plot (SMD vs. number of trees)
#                           2) Positivity (Overlap) assumption
#                           3) Covariate Balance Diagnostics
#                           4) Variance Ratios & Distributional Diagnostics

#           Part II: Modeling in weighted dataset of Initiated-High & No-Care 
#           Part III: Do the same analysis (Part I & II) for:
#                                             1) Ceased vs No Caregiving               
#                                             2) Initiated (Medium) vs No Caregiving    
#                                             3) Sustained (to Medium) vs No Caregiving 
#                                             4) Sustained (to High) vs No Caregiving   
#          Part IV: 1) Put 5 overlap graphs together
#                   2) Put 5 weights distribution graphs together


df3 <- readRDS("df3.rds")

# GBM-IPTW Implementation (using WeightIt and cobalt)
library(WeightIt) # weighting & MSM
library(cobalt) # balance diagnostics
library(patchwork) #  combine plots

# research question answered by GBM-IPTW is “What is the average effect of being in a given 
# caregiving category (versus no caregiving) on 2016 EAA, among those who actually experienced that 
# caregiving transition?” This is a static, cross-sectional-type question focusing on final treatment assignment.

#--------------------------------------------------------------------------------------------------
### Part 1 Pair-wise ATT diagnostics for  Initiated-High vs No-Care 
#--------------------------------------------------------------------------------------------------

# 0.  Prepare analysis sample: Initiated-High vs No-Care only
table(df3$care_trans)

# No Caregiving Initiated (Medium) Initiated (High)   Ceased Sustained (to Medium)   Sustained (to High) 
# 2412                   109                    94    186                   237                   115 

dat_IH <- subset(df3, care_trans %in% c("Initiated (High)", "No Caregiving"))
dat_IH$treat_IH <- ifelse(dat_IH$care_trans == "Initiated (High)", 1, 0)

# 1.  Estimate GBM propensity weights for ATT
set.seed(123) 
w.out <- weightit(
  treat_IH ~ r11age + female + racecat +
    r11marcat_b + edu + h11itot1 + r11shlt +
    r11smoke + r11drinkn + r11mobila +
    r11chronic + r11acts_w + r11bmi,
  data               = dat_IH,
  method             = "gbm",
  estimand           = "ATT",
  focal              = 1,
  n.trees            = 2500,        
  shrinkage          = 0.01,         # learning rate
  interaction.depth  = 2,            
  criterion          = "smd.max",    # minimize max abs(SMD)
  abs = TRUE,
  verbose = FALSE)

# 2  Diagnostic Pipeline
# 2.1 Balance plot (SMD vs. number of trees)
plot(w.out, type = "balance")
optimal_trees <- w.out$info$best.tree
print(optimal_trees) # 2493

#  2.2 Positivity (Overlap: “common support") assumption
dat_IH$ps <- w.out$ps
p_overlap_IH <- ggplot(dat_IH, aes(x = ps, fill = as.factor(treat_IH))) +
  geom_density(alpha = 0.5) +
  labs(x = "Propensity Score", y = "Density", fill = "Treatment Group") +
  theme_minimal()  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

        ggsave("Figure_S1_IH_Propensity_Overlap.pdf",
               plot = p_overlap_IH,
               device = pdf,
               width = 6, height = 4, units = "in")

      # 2.2.1 . Calculate overlap and trim based on common support
      ps_t <- dat_IH$ps[dat_IH$treat_IH == 1]
      ps_c <- dat_IH$ps[dat_IH$treat_IH == 0]
      lower <- max(min(ps_t), min(ps_c))
      upper <- min(max(ps_t), max(ps_c))
      
      cat("Overlap PS before trimming:", round(lower, 3), "-", round(upper, 3), "\n")
      # 0.02 - 0.312 
      cat("Percent retained before trimming:", round(mean(dat_IH$ps >= lower & dat_IH$ps <= upper) * 100, 1), "%\n")
      # 45.7 %


summary(w.out$ps)      
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001364 0.008918 0.017471 0.037648 0.042120 0.693697 
summary(w.out$weights)
# Min.       1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001366 0.008999 0.017782 0.070015 0.044347 1.000000 

# Have strong separation—most controls look nothing like treated, so they get near‐zero weight.
# Only the controls in the overlap region (PS ≈ 0.01–0.10) meaningfully contribute to the ATT estimate.


summary(w.out)

                  # Summary of weights
                  # 
                  # - Weight ranges:
                  #   
                  #   Min                                    Max
                  # treated 1.0000                              || 1.0000
                  # control 0.0014   |-----------|                 0.4531
                  # 
                  # - Units with the 5 most extreme weights by group:
                  #   
                  #   47     44     28     11      8
                  # treated      1      1      1      1      1
                  #             1764    830   2480   2477    388
                  # control    0.3495 0.3559 0.4148 0.4167 0.4531
                  # 
                  # - Weight statistics:
                  #   
                  #   Coef of Var  MAD Entropy # Zeros
                  # treated       0.000 0.00   0.000       0
                  # control       1.353 0.85   0.576       0
                  # 
                  # - Effective Sample Sizes:
                  #   
                  #   Control Treated
                  # Unweighted 2412.        94
                  # Weighted    852.13      94
      
                  wt_ctl <- w.out$weights[dat_IH$treat_IH == 0]               # control weights
                  cut99_ctl <- quantile(wt_ctl, 0.99)                        # 99th percentile (about 0.044)
                  share_top1_ctl <- sum(wt_ctl[wt_ctl > cut99_ctl]) / sum(wt_ctl)
                  cat("Top‑1% weight share among controls:", round(share_top1_ctl * 100, 1), "%\n")
                  # ≈  9.3 %

                  # 2.2.2  Histogram of ATT weights with 99th-pct line (Initiated High)
                  dat_IH$w_att <- get.w(w.out)
                  p_weights_IH <- ggplot(dat_IH, aes(x = w_att)) +
                    geom_histogram(bins = 40, colour = "white", fill = "steelblue", alpha = .7) +
                    geom_vline(aes(xintercept = quantile(w_att, 0.99)),
                               linetype = "dashed", linewidth = .8) +
                    annotate("text",
                             x = quantile(dat_IH$w_att, 0.99),
                             y = Inf, vjust = -0.3, hjust = 1.1, angle = 90,
                             label = paste0("99th pct = ", round(quantile(dat_IH$w_att, .99), 2))) +
                    labs(x = "ATT weight",
                         y = " ",
                         title = "B Transition to high intensity care") +
                    theme_minimal() +
                    theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())
                  
                  ggsave("Figure_IH_WEIGHT_Histogram.pdf",
                         plot = p_weights_IH,
                         device = pdf, width = 6, height = 4, units = "in")

# 2.3. Covariate Balance Diagnostics
# Examine balance before and after weighting
bal <- bal.tab(w.out, un = TRUE)
# Balance Measures
# Type Diff.Un Diff.Adj
# prop.score                Distance  1.2829   0.9155
# r11age                     Contin.  0.0791   0.0008
# female                      Binary -0.0469   0.0071
# racecat_White               Binary -0.0557   0.0111
# racecat_Black               Binary  0.0164   0.0072
# racecat_Hispanic            Binary  0.0692   0.0042
# racecat_Other               Binary -0.0299  -0.0226
# r11marcat_b_Married         Binary  0.2738   0.0190
# edu_Less than High School   Binary  0.1323   0.0155
# edu_High School Graduate    Binary  0.0961   0.0155
# edu_Some College            Binary -0.0993  -0.0137
# edu_College Degree          Binary -0.1291  -0.0174
# h11itot1                   Contin. -0.0847  -0.0251
# r11shlt                    Contin.  0.1597   0.0270
# r11smoke_Never smoker       Binary -0.0418  -0.0012
# r11smoke_Past smoker        Binary  0.0086   0.0036
# r11smoke_Current smoker     Binary  0.0333  -0.0023
# r11drinkn                  Contin. -0.0785   0.0060
# r11mobila                  Contin.  0.1154   0.0199
# r11chronic                 Contin.  0.0563  -0.0202
# r11acts_w                  Contin. -0.1452  -0.0298
# r11bmi                     Contin. -0.0512   0.0010
# 
# Effective sample sizes
# Control Treated
# Unadjusted 2412.        94
# Adjusted    852.13      94

# The balance diagnostics table demonstrates that, after weighting, all covariates exhibit extremely 
# small adjusted differences (Diff.Adj), well below the conventional threshold of 0.1, indicating 
# excellent balance between treated (“Initiated (High)”) and control (“No Caregiving”) groups on all 
# measured variables. This strong covariate balance supports the plausibility of the unconfoundedness 
# (or ignorability) assumption for observed confounders, a key requirement for valid causal inference. 
# Additionally, the reduction in the control group’s effective sample size highlights the challenge of 
# limited overlap (positivity), as only controls similar to the treated group meaningfully contribute
# to the analysis. 

# 2.4. Residual analysis
# 2.4.1. Fit the Propensity Score Model
ps_model <- glm(treat_IH ~ r11age + female + black + hisp + other +
                  r11marcat_b + edu + h11itot1 + r11shlt +
                  r11smoke + r11drinkn + r11mobila +
                  r11chronic + r11acts_w + r11bmi,
                data = dat_IH,
                family = binomial)
# 2.4.2. Compute Residuals
dat_IH$residuals <- residuals(ps_model, type = "deviance")


# 2.4.3. Visualize Residuals
# Residuals vs. Fitted values
ggplot(dat_IH, aes(x = fitted(ps_model), y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  labs(x = "Fitted Propensity Scores", y = "Deviance Residuals") +
  theme_minimal()

# Histogram of residuals
ggplot(dat_IH, aes(x = residuals)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  labs(x = "Deviance Residuals", y = "Frequency") +
  theme_minimal()          
# no substantial modifications

# 3.1. SMD plot (left)
love.plot(bal, threshold = 0.1)

# Labels

var_labs <- c(
  "r11marcat_b_Married"        = "Married/Partenered",
  "r11shlt"                    = "Self-Rated Health",
  "r11acts_w"                  = "Physical Activity",
  "edu_Less than High School"  = "Less than High School",
  "edu_College Degree"         = "College Degree",
  "r11mobila"                  = "Mobility Limitation",
  "edu_Some College"           = "Some College",
  "edu_High School Graduate"   = "High School Graduate",
  "h11itot1"                   = "Household Income",
  "r11age"                     = "Age",
  "r11drinkn"                  = "Drinking Now",
  "racecat_Hispanic"           = "Hispanic",
  "r11chronic"                 = "Chronic Disease",
  "racecat_White"              = "White",
  "r11bmi"                     = "BMI",
  "female"                     = "Female",
  "r11smoke_Never smoker"      = "Never Smoker",
  "r11smoke_Current smoker"    = "Current Smoker",
  "racecat_Other"              = "Other Race",
  "racecat_Black"              = "Black",
  "r11smoke_Past smoker"       = "Past Smoker"
)
p_smd1 <- love.plot(
  w.out,
  stats = "mean.diffs",
  threshold = 0.1,
  abs = TRUE,
  line = TRUE,
  var.order = "unadjusted",
  drop.distance = TRUE,
  colors = c("#F8766D", "#00BFC4"),
  shapes = c(19, 18),
  sample.names = c("Unweighted", "Weighted"),
  var.names     = var_labs,
  title = " "
) +
  scale_x_continuous(limits = c(0, 0.5)) +
  xlab("Absolute Standardized Mean Differences")

# 3.2. KS plot (right, no y labels)
p_ks1 <- love.plot(
  x             = w.out,
  stats         = "ks.statistics",
  thresholds    = c(ks.statistics = 0.1),
  abs           = TRUE,
  line          = TRUE,
  drop.distance = TRUE,
  var.order     = p_smd1,  # Pass the love.plot object here
  colors        = c("#F8766D", "#00BFC4"),
  shapes        = c(19, 18),
  sample.names  = c("Unweighted", "Weighted"),
  disp.ks       = TRUE,
  un            = TRUE,
  var.names = var_labs,   
  title         = " ",
  xlab          = "KS Statistic"
) +
  scale_x_continuous(limits = c(0, 0.5)) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

# 3.3. Combine, collect the legend, and set legend to right
combined_plot1 <- (p_smd1 | p_ks1) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

      combined_plot1
      
      ggsave("Figure_1_IH_SMD and KS.pdf",
             plot = combined_plot1,
             device = pdf,
             width = 9, height = 5, units = "in")

# 4. Variance Ratios & Distributional Diagnostics
# 4.1 Variance ratios:
bal.tab(w.out, stats = "variance.ratios")
# Variance ratios should be close to 1 (typically 0.5 to 2 is considered acceptable).

# 4.2 Kolmogorov–Smirnov (KS) test for covariates:
bal.tab(w.out, stats = "ks.statistics")
#  Smaller KS statistics indicate better distributional balance between groups.

#--------------------------------------------------------------------------------------------------
### Part II Modeling in weighted dataset of Initiated-High & No-Care 
#--------------------------------------------------------------------------------------------------
# Step 1: Apply weights
dat_IH$w_att <- get.w(w.out)
dat_IH$comb_w <- dat_IH$w_att * dat_IH$vbsi16wgtra  # Combine IPTW with HRS survey weight

# Step 2: Define the survey design
design_IH <- svydesign(
  ids = ~secu,
  strata = ~stratum,
  weights = ~comb_w,
  data = dat_IH,
  nest = TRUE
)

# Step 3: Define outcome variables
dv_list <- c("horvath_residual1_z",      "hannum_residual1_z",    "levine_residual1_z", 
             "horvathskin_residual1_z",  "lin_residual1_z",       "weidner_residual1_z", 
             "vidalbralo_residual1_z",   "grim_residual1_z",      "zhang_residual1_z", 
             "yang_residual1_z",        "bocklandt_residual1_z", "garagnani_residual1_z", 
             "mpoa_z")

# Step 4: Fit ATT-only models (outcome ~ treat_IH)
att_models <- setNames(
  lapply(dv_list, function(dv) {
    svyglm(
      formula = as.formula(paste(dv, "~ treat_IH")),
      design = design_IH
    )
  }),
  dv_list
)

# Step 5: Define covariates for doubly robust models
cov_w11 <- c(
  "r11age", "female", "black", "hisp", "other", "r11marcat_b", "edu",
  "h11itot1", "r11shlt", "r11smoke", "r11drinkn", "r11mobila",
  "r11chronic", "r11acts_w", "r11bmi"
)

# Step 6: Fit doubly robust models (outcome ~ treat_IH + covariates)
rhs_dr <- paste("treat_IH +", paste(cov_w11, collapse = " + "))
dr_models <- setNames(
  lapply(dv_list, function(dv) {
    svyglm(
      formula = as.formula(paste(dv, "~", rhs_dr)),
      design = design_IH
    )
  }),
  dv_list
)

# Step 7: Print summaries for ATT-only models
for (dv in dv_list) {
  cat("\n===== ATT-only summary for:", dv, "=====\n")
  print(summary(att_models[[dv]]))
}

# Step 8: Print summaries for doubly robust models
for (dv in dv_list) {
  cat("\n===== Doubly Robust summary for:", dv, "=====\n")
  print(summary(dr_models[[dv]]))
}


### A) _resudual1_z      
#===== Doubly Robust summary for: zhang_dnamage_z =====
# treat_IH                 0.319701   0.121982   2.621 0.012650 * 

# ===== Doubly Robust summary for: grim_residual1_z =====
# treat_IH                 0.212059   0.080278   2.642   0.0120 * 

# ==== Doubly Robust summary for: vidalbralo_residual1_z =====
# treat_IH                 0.0620443  0.0826307   0.751  0.45748

### B) _dff_z
#===== Doubly Robust summary for: zhang_dnamage_z =====
# treat_IH                 0.0144606  0.0072511    1.994  0.05353 . 

# ===== Doubly Robust summary for: grim_diff_z =====
# treat_IH                 0.176487   0.070580   2.501   0.0170 * 

# ===== Doubly Robust summary for: vidalbralo_diff_z =====
# treat_IH                 0.022364   0.054147   0.413 0.681971   

# ===== Doubly Robust summary for: garagnani_dnamage_z =====
# treat_IH                -0.247876   0.069025  -3.591 0.000952 ***


### C) _residual_z
# ===== Doubly Robust summary for: grim_residual_z =====
# treat_IH                 0.196451   0.078575   2.500   0.0170 *  

# ===== Doubly Robust summary for: vidalbralo_residual_z =====
# treat_IH                 0.031930   0.082959   0.385 0.702527


#--------------------------------------------------------------------------------------------------
###                                     Part III
#--------------------------------------------------------------------------------------------------        

###---------------------------  Ceased vs No Caregiving ---------------------------       

# 0. Prepare analysis sample: Ceased vs No Caregiving
dat_CE <- subset(df3, care_trans %in% c("Ceased", "No Caregiving"))
dat_CE$treat_CE <- ifelse(dat_CE$care_trans == "Ceased", 1, 0)

# 1. Estimate GBM propensity weights for ATT (Ceased)
set.seed(123)
w.out_ce <- weightit(
  treat_CE ~ r11age + female + racecat +
    r11marcat_b + edu + h11itot1 + r11shlt +
    r11smoke + r11drinkn + r11mobila +
    r11chronic + r11acts_w + r11bmi,
  data               = dat_CE,
  method             = "gbm",
  estimand           = "ATT",
  focal              = 1,
  n.trees            = 2500,
  shrinkage          = 0.01,
  interaction.depth  = 2,
  criterion          = "smd.max",
  abs                = TRUE,
  verbose            = FALSE
)

# 2. Diagnostic Pipeline
# 2.1 Balance plot (SMD vs. number of trees)
plot(w.out_ce, type = "balance")
optimal_trees_ce <- w.out_ce$info$best.tree
print(optimal_trees_ce) # 2448

# 2.2 Positivity (Overlap: “Common Support") assumption
dat_CE$ps_ce <- w.out_ce$ps
p_overlap_CE <- ggplot(dat_CE, aes(x = ps_ce, fill = as.factor(treat_CE))) +
  geom_density(alpha = 0.5) +
  labs(x = "Propensity Score", y = "Density", fill = "Treatment Group") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("Figure_S2_CE_Propensity_Overlap.pdf",
       plot = p_overlap_CE,
       device = pdf,
       width = 6, height = 4, units = "in")
          
          # 2.2.1  Calculate common-support interval and retention (Ceased vs No-Care)
          
          # Separate treated and control propensity scores
          ps_t <- dat_CE$ps_ce[dat_CE$treat_CE == 1]   # Ceased caregivers
          ps_c <- dat_CE$ps_ce[dat_CE$treat_CE == 0]   # Non-caregivers
          
          # Common-support bounds (max of mins, min of maxs)
          lower <- max(min(ps_t), min(ps_c))
          upper <- min(max(ps_t), max(ps_c))
          
          cat("Overlap PS before trimming:", round(lower, 3), "-", round(upper, 3), "\n")
            # Overlap PS before trimming: 0.039 - 0.359 
          
          # Share of full sample that lies within the support band
          retained_pct <- mean(dat_CE$ps_ce >= lower & dat_CE$ps_ce <= upper) * 100
          cat("Percent retained before trimming:", round(retained_pct, 1), "%\n")
            # Percent retained before trimming: 68.8 %

summary(w.out_ce$ps)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.005884 0.034234 0.055809 0.071508 0.089117 0.512321 
summary(w.out_ce$weights)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.005919 0.035447 0.059698 0.137397 0.102524 1.000000 
 summary(w.out_ce)
          # > summary(w.out_ce)
          # Summary of weights
          # 
          # - Weight ranges:
          #   
          #   Min                                   Max
          # treated 1.0000                              || 1.000
          # control 0.0059   |--------------|              0.559
          # 
          # - Units with the 5 most extreme weights by group:
          #   
          #   117     94     92     34    19
          # treated      1      1      1      1     1
          # 1621   2226   1364   2123  2422
          # control 0.4421 0.4445 0.5114 0.5553 0.559
          # 
          # - Weight statistics:
          #   
          #   Coef of Var   MAD Entropy # Zeros
          # treated       0.000 0.000   0.000       0
          # control       0.801 0.555   0.251       0
          # 
          # - Effective Sample Sizes:
          #   
          #   Control Treated
          # Unweighted 2412.       186
          # Weighted   1469.91     186
 
           # 2.2.2 Wights distribution
           dat_CE$w_att <- get.w(w.out_ce)
           p_weights_CE <- ggplot(dat_CE, aes(x = w_att)) +
             geom_histogram(bins = 40, colour = "white", fill = "steelblue", alpha = .7) +
             geom_vline(aes(xintercept = quantile(w_att, 0.99)),
                        linetype = "dashed", linewidth = .8) +
             annotate("text",
                      x = quantile(dat_CE$w_att, 0.99),
                      y = Inf, vjust = -0.3, hjust = 1.1, angle = 90,
                      label = paste0("99th pct = ", round(quantile(dat_CE$w_att, .99), 2))) +
             labs(x = "ATT weight",
                  y = " ",
                  title = "C Transition out of care") +
             theme_minimal() +
             theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
           
           ggsave("Figure_CE_WEIGHT_Histogram.pdf",
                  plot = p_weights_CE,
                  device = pdf, width = 6, height = 4, units = "in")
 
# 2.3 Covariate Balance Diagnostics
bal_ce <- bal.tab(w.out_ce, un = TRUE)
# Balance Measures
# Type Diff.Un Diff.Adj
# prop.score                Distance  1.0584   0.7087
# r11age                     Contin. -0.3474  -0.0117
# female                      Binary  0.0556   0.0133
# racecat_White               Binary -0.0059  -0.0016
# racecat_Black               Binary  0.0182   0.0076
# racecat_Hispanic            Binary -0.0202  -0.0063
# racecat_Other               Binary  0.0078   0.0003
# r11marcat_b_Married         Binary  0.1272   0.0114
# edu_Less than High School   Binary -0.0418  -0.0052
# edu_High School Graduate    Binary -0.0124   0.0060
# edu_Some College            Binary  0.0743   0.0025
# edu_College Degree          Binary -0.0201  -0.0033
# h11itot1                   Contin.  0.1055   0.0063
# r11shlt                    Contin. -0.0304  -0.0125
# r11smoke_Never smoker       Binary  0.0163  -0.0032
# r11smoke_Past smoker        Binary -0.0511  -0.0008
# r11smoke_Current smoker     Binary  0.0349   0.0040
# r11drinkn                  Contin.  0.0773   0.0264
# r11mobila                  Contin.  0.0233   0.0236
# r11chronic                 Contin. -0.1627  -0.0276
# r11acts_w                  Contin.  0.0693  -0.0017
# r11bmi                     Contin.  0.0178   0.0134
# 
# Effective sample sizes
# Control Treated
# Unadjusted 2412.       186
# Adjusted   1469.91     186

# Extract control weights from the Ceased vs No Caregiving model
wt_ctl_ce <- w.out_ce$weights[dat_CE$treat_CE == 0]

# Compute the 99th percentile cutoff
cut99_ctl_ce <- quantile(wt_ctl_ce, 0.99)

# Proportion of total control weight contributed by the top 1%
share_top1_ctl_ce <- sum(wt_ctl_ce[wt_ctl_ce > cut99_ctl_ce]) / sum(wt_ctl_ce)

# Print the result
cat("Top‑1% weight share among controls (Ceased vs No Caregiving):",
    round(share_top1_ctl_ce * 100, 1), "%\n")
# 5.4 %

# 2.4 Residual Analysis
ps_model_ce <- glm(
  treat_CE ~ r11age + female + black + hisp + other +
    r11marcat_b + edu + h11itot1 + r11shlt +
    r11smoke + r11drinkn + r11mobila +
    r11chronic + r11acts_w + r11bmi,
  data = dat_CE,
  family = binomial
)
dat_CE$residuals_ce <- residuals(ps_model_ce, type = "deviance")
# Residuals vs Fitted
ggplot(dat_CE, aes(x = fitted(ps_model_ce), y = residuals_ce)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  labs(x = "Fitted Propensity Scores", y = "Deviance Residuals") +
  theme_minimal()
# Histogram of residuals
ggplot(dat_CE, aes(x = residuals_ce)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  labs(x = "Deviance Residuals", y = "Frequency") +
  theme_minimal()

# 3.1 SMD plot (left)
love.plot(bal_ce, threshold = 0.1)

# 3.2 SMD and KS plots for publication-style figure
p_smd_ce <- love.plot(
  w.out_ce,
  stats = "mean.diffs",
  threshold = 0.1,
  abs = TRUE,
  line = TRUE,
  var.order = "unadjusted",
  drop.distance = TRUE,
  colors = c("#F8766D", "#00BFC4"),
  shapes = c(19, 18),
  sample.names = c("Unweighted", "Weighted"),
  var.names = var_labs, 
  title = " "
) +
  scale_x_continuous(limits = c(0, 0.5)) +
  xlab("Absolute Standardized Mean Differences")

p_ks_ce <- love.plot(
  x             = w.out_ce,
  stats         = "ks.statistics",
  thresholds    = c(ks.statistics = 0.1),
  abs           = TRUE,
  line          = TRUE,
  drop.distance = TRUE,
  var.order     = p_smd_ce,  # Pass love.plot object for consistent y-order
  colors        = c("#F8766D", "#00BFC4"),
  shapes        = c(19, 18),
  sample.names  = c("Unweighted", "Weighted"),
  disp.ks       = TRUE,
  un            = TRUE,
  var.names     = var_labs, 
  title         = " ",
  xlab          = "KS Statistic"
) +
  scale_x_continuous(limits = c(0, 0.5)) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

combined_plot_ce <- (p_smd_ce | p_ks_ce) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "right")
combined_plot_ce

    ggsave("Figure_2_CE_SMD and KS.pdf",
           plot = combined_plot_ce,
           device = pdf,
           width = 9, height = 5, units = "in")


# 4. Variance Ratios & KS
bal.tab(w.out_ce, stats = "variance.ratios")
bal.tab(w.out_ce, stats = "ks.statistics")


# 5. Modeling in weighted dataset of Ceased & No-Care

# 5.1: Apply weights
dat_CE$w_att_ce <- get.w(w.out_ce)
dat_CE$comb_w_ce <- dat_CE$w_att_ce * dat_CE$vbsi16wgtra

# 5.2: Define survey design
design_CE <- svydesign(
  ids     = ~secu,
  strata  = ~stratum,
  weights = ~comb_w_ce,
  data    = dat_CE,
  nest    = TRUE
)

# 5.3: ATT-only models (outcome ~ treat_CE)
att_models_ce <- setNames(
  lapply(dv_list, function(dv) {
    svyglm(
      formula = as.formula(paste(dv, "~ treat_CE")),
      design = design_CE
    )
  }),
  dv_list
)

# 5.4: Covariates for doubly robust models
rhs_dr_ce <- paste("treat_CE +", paste(cov_w11, collapse = " + "))

# 5.5: Doubly robust models (outcome ~ treat_CE + covariates)
dr_models_ce <- setNames(
  lapply(dv_list, function(dv) {
    svyglm(
      formula = as.formula(paste(dv, "~", rhs_dr_ce)),
      design = design_CE
    )
  }),
  dv_list
)

# 5.6: Print summaries for ATT-only models
for (dv in dv_list) {
  cat("\n===== ATT-only summary for:", dv, "=====\n")
  print(summary(att_models_ce[[dv]]))
}

# 5.7: Print summaries for doubly robust models
for (dv in dv_list) {
  cat("\n===== Doubly Robust summary for:", dv, "=====\n")
  print(summary(dr_models_ce[[dv]]))
}

### A) _resudual1_z  
# ===== Doubly Robust summary for: zhang_dnamage_z =====
# treat_CE                -0.028402   0.070100  -0.405 0.687694     

# ===== Doubly Robust summary for: grim_residual1_z =====
# treat_CE                -0.0313983  0.0535433  -0.586  0.56116    

# ===== Doubly Robust summary for: vidalbralo_residual1_z =====
# treat_CE                -0.162609   0.063808  -2.548   0.0151 * 


### B) _diff_z
# ===== Doubly Robust summary for: grim_diff_z =====
# treat_CE                 0.005352   0.045274   0.118  0.90655 
# 
# ===== Doubly Robust summary for: vidalbralo_diff_z =====
# treat_CE                -0.072930   0.040765  -1.789  0.08180 . 

### C) residual_z
# # ===== Doubly Robust summary for: vidalbralo_residual_z =====
# treat_CE                -0.1124637  0.0636180  -1.768  0.08534 .
# 
# ===== Doubly Robust summary for: grim_residual_z =====
# treat_CE                 0.006938   0.050365   0.138  0.89118  


###--------------------------- Initiated (Medium) vs No Caregiving ---------------------------
# 0. Prepare Analysis Sample: Initiated (Medium) vs No Caregiving
dat_IM <- subset(df3, care_trans %in% c("Initiated (Medium)", "No Caregiving"))
dat_IM$treat_IM <- ifelse(dat_IM$care_trans == "Initiated (Medium)", 1, 0)

# 1. Estimate GBM Propensity Weights for ATT
set.seed(123)
w.out_im <- weightit(
  treat_IM ~ r11age + female + racecat +
    r11marcat_b + edu + h11itot1 + r11shlt +
    r11smoke + r11drinkn + r11mobila +
    r11chronic + r11acts_w + r11bmi,
  data               = dat_IM,
  method             = "gbm",
  estimand           = "ATT",
  focal              = 1,
  n.trees            = 2500,
  shrinkage          = 0.01,
  interaction.depth  = 2,
  criterion          = "smd.max",
  abs                = TRUE,
  verbose            = FALSE
)

# 2. Diagnostic Pipeline
# 2.1 Balance Plot (SMD vs. Number of Trees)
plot(w.out_im, type = "balance")
optimal_trees_im <- w.out_im$info$best.tree
print(optimal_trees_im) # 2500

# 2.2 Positivity (Overlap: “Common Support”) Assumption
dat_IM$ps_im <- w.out_im$ps
p_overlap_IM <-  ggplot(dat_IM, aes(x = ps_im, fill = as.factor(treat_IM))) +
  geom_density(alpha = 0.5) +
  labs(x = "Propensity Score", y = "Density", fill = "Treatment Group") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

      ggsave("Figure_S3_IM_Propensity_Overlap.pdf",
             plot = p_overlap_IM,
             device = pdf,
             width = 6, height = 4, units = "in")
      
      # 2.2.1 Calculate overlap and retention for Initiated (Medium) vs No-Care
      
      # Separate treated and control propensity scores
      ps_t_im <- dat_IM$ps_im[dat_IM$treat_IM == 1]  # Initiated (Medium)
      ps_c_im <- dat_IM$ps_im[dat_IM$treat_IM == 0]  # Non-caregivers
      
      # Compute common-support bounds
      lower_im <- max(min(ps_t_im), min(ps_c_im))
      upper_im <- min(max(ps_t_im), max(ps_c_im))
      
      cat("Overlap PS before trimming:", round(lower_im, 3), "-", round(upper_im, 3), "\n")
      # Overlap PS before trimming: 0.011 - 0.431
      
      # Calculate retention percentage within support
      retained_pct_im <- mean(dat_IM$ps_im >= lower_im & dat_IM$ps_im <= upper_im) * 100
      cat("Percent retained before trimming:", round(retained_pct_im, 1), "%\n")
      # Percent retained before trimming: 53.2 %
      

summary(w.out_im$ps)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001358 0.005320 0.014413 0.043336 0.055060 0.608370 
summary(w.out_im$weights)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001360 0.005348 0.014971 0.081453 0.059071 1.000000 
summary(w.out_im)

        # Summary of weights
        # 
        # - Weight ranges:
        #   
        #   Min                                    Max
        # treated 1.0000                              || 1.0000
        # control 0.0014   |--------------------|        0.7584
        # 
        # - Units with the 5 most extreme weights by group:
        #   
        #   725    704    516    215     41
        # treated      1      1      1      1      1
        # 1908   1521   2364   2009   2283
        # control 0.4134 0.4459 0.5331 0.5956 0.7584
        # 
        # - Weight statistics:
        #   
        #   Coef of Var   MAD Entropy # Zeros
        # treated       0.000 0.000   0.000       0
        # control       1.551 1.007   0.767       0
        # 
        # - Effective Sample Sizes:
        #   
        #   Control Treated
        # Unweighted 2412.       109
        # Weighted    708.53     109

      # 2.2.2 Histogram of ATT weights with 99th-pct line for Initiated (Medium) vs No Caregiving
      dat_IM$w_att <- get.w(w.out_im)
      p_weights_IM <- ggplot(dat_IM, aes(x = w_att)) +
        geom_histogram(bins = 40, colour = "white", fill = "steelblue", alpha = .7) +
        geom_vline(aes(xintercept = quantile(w_att, 0.99)),
                   linetype = "dashed", linewidth = .8) +
        annotate("text",
                 x = quantile(dat_IM$w_att, 0.99),
                 y = Inf, vjust = -0.3, hjust = 1.1, angle = 90,
                 label = paste0("99th pct = ", round(quantile(dat_IM$w_att, .99), 2))) +
        labs(x = "ATT weight",
             y = "Count",
             title = "A Transition to medium intensity care") +
        theme_minimal() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      
      ggsave("Figure_IM_WEIGHT_Histogram.pdf",
             plot = p_weights_IM,
             device = pdf, width = 6, height = 4, units = "in")

# 2.3 Covariate Balance Diagnostics
bal_im <- bal.tab(w.out_im, un = TRUE)

# Balance Measures
# Type                                Diff.Un Diff.Adj
# prop.score                Distance  1.3447   0.7861
# r11age                     Contin. -1.9001  -0.1141
# female                      Binary  0.0267   0.0200
# racecat_White               Binary -0.0151   0.0263
# racecat_Black               Binary  0.0388  -0.0011
# racecat_Hispanic            Binary -0.0122  -0.0113
# racecat_Other               Binary -0.0115  -0.0139
# r11marcat_b_Married         Binary  0.0360  -0.0069
# edu_Less than High School   Binary -0.0753  -0.0079
# edu_High School Graduate    Binary -0.0711   0.0131
# edu_Some College            Binary  0.1111   0.0043
# edu_College Degree          Binary  0.0353  -0.0095
# h11itot1                   Contin.  0.0402  -0.0044
# r11shlt                    Contin. -0.1478   0.0043
# r11smoke_Never smoker       Binary -0.0241  -0.0020
# r11smoke_Past smoker        Binary  0.0388   0.0127
# r11smoke_Current smoker     Binary -0.0148  -0.0107
# r11drinkn                  Contin.  0.1039   0.0036
# r11mobila                  Contin. -0.4116  -0.0506
# r11chronic                 Contin. -0.4195  -0.0464
# r11acts_w                  Contin.  0.3900  -0.0237
# r11bmi                     Contin.  0.0269  -0.0190
# 
# Effective sample sizes
# Control Treated
# Unadjusted 2412.       109
# Adjusted    708.53     109

    # Calculate top-1% weight share among controls for Initiated (Medium)
    wt_ctl_im <- w.out_im$weights[dat_IM$treat_IM == 0]    # control weights
    cut99_ctl_im <- quantile(wt_ctl_im, 0.99)              # 99th percentile for controls
    share_top1_ctl_im <- sum(wt_ctl_im[wt_ctl_im > cut99_ctl_im]) / sum(wt_ctl_im)
    cat("Top‑1% weight share among controls:", round(share_top1_ctl_im * 100, 1), "%\n")
    # Top‑1% weight share among controls: 10.4 %


# 2.4 Residual Analysis
ps_model_im <- glm(
  treat_IM ~ r11age + female + black + hisp + other +
    r11marcat_b + edu + h11itot1 + r11shlt +
    r11smoke + r11drinkn + r11mobila +
    r11chronic + r11acts_w + r11bmi,
  data = dat_IM,
  family = binomial
)
dat_IM$residuals_im <- residuals(ps_model_im, type = "deviance")

# Residuals vs. Fitted Values
ggplot(dat_IM, aes(x = fitted(ps_model_im), y = residuals_im)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  labs(x = "Fitted Propensity Scores", y = "Deviance Residuals") +
  theme_minimal()

# Histogram of Residuals
ggplot(dat_IM, aes(x = residuals_im)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  labs(x = "Deviance Residuals", y = "Frequency") +
  theme_minimal()

# 3. Covariate Balance Plots
# 3.1 SMD Plot
love.plot(bal_im, threshold = 0.1)

# 3.2 SMD and KS Plots for  Figure
p_smd_im <- love.plot(
  w.out_im,
  stats = "mean.diffs",
  threshold = 0.1,
  abs = TRUE,
  line = TRUE,
  var.order = "unadjusted",
  drop.distance = TRUE,
  colors = c("#F8766D", "#00BFC4"),
  shapes = c(19, 18),
  sample.names = c("Unweighted", "Weighted"),
  var.names = var_labs,
  title = " "
) +
  scale_x_continuous(limits = c(0, 0.5)) +
  xlab("Absolute Standardized Mean Differences")

p_ks_im <- love.plot(
  x             = w.out_im,
  stats         = "ks.statistics",
  thresholds    = c(ks.statistics = 0.1),
  abs           = TRUE,
  line          = TRUE,
  drop.distance = TRUE,
  var.order     = p_smd_im,
  colors        = c("#F8766D", "#00BFC4"),
  shapes        = c(19, 18),
  sample.names  = c("Unweighted", "Weighted"),
  disp.ks       = TRUE,
  un            = TRUE,
  var.names     = var_labs,
  title         = " ",
  xlab          = "KS Statistic"
) +
  scale_x_continuous(limits = c(0, 0.5)) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

combined_plot_im <- (p_smd_im | p_ks_im) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "right")
combined_plot_im
    ggsave("Figure_3_IM_SMD and KS.pdf",
           plot = combined_plot_im,
           device = pdf,
           width = 9, height = 5, units = "in")

# 4. Variance Ratios & Distributional Diagnostics
bal.tab(w.out_im, stats = "variance.ratios")
bal.tab(w.out_im, stats = "ks.statistics")

# Part II: Modeling in weighted Dataset of Initiated (Medium) & No Care
# Step 1: Apply Weights
dat_IM$w_att_im <- get.w(w.out_im)
dat_IM$comb_w_im <- dat_IM$w_att_im * dat_IM$vbsi16wgtra

# Step 2: Define the Survey Design
design_IM <- svydesign(
  ids     = ~secu,
  strata  = ~stratum,
  weights = ~comb_w_im,
  data    = dat_IM,
  nest    = TRUE
)

# Step 3: Fit ATT-Only Models (Outcome ~ treat_IM)
att_models_im <- setNames(
  lapply(dv_list, function(dv) {
    svyglm(
      formula = as.formula(paste(dv, "~ treat_IM")),
      design = design_IM
    )
  }),
  dv_list
)

# Step 4: Fit Doubly Robust Models (Outcome ~ treat_IM + Covariates)
rhs_dr_im <- paste("treat_IM +", paste(cov_w11, collapse = " + "))
dr_models_im <- setNames(
  lapply(dv_list, function(dv) {
    svyglm(
      formula = as.formula(paste(dv, "~", rhs_dr_im)),
      design = design_IM
    )
  }),
  dv_list
)
for (dv in dv_list) {
  cat("\n===== ATT-only summary for:", dv, "=====\n")
  print(summary(att_models_im[[dv]]))
}
for (dv in dv_list) {
  cat("\n===== Doubly Robust summary for:", dv, "=====\n")
  print(summary(dr_models_im[[dv]]))
}

### A) _residual1_z
# ===== Doubly Robust summary for: zhang_dnamage_z =====
# treat_IM                -0.013979   0.103121  -0.136  0.89290  

# ===== Doubly Robust summary for: grim_residual1_z =====
# treat_IM                -0.045847   0.114326  -0.401  0.69071 

# ===== Doubly Robust summary for: vidalbralo_residual1_z =====
# treat_IM                -0.109915   0.111134  -0.989  0.32907 

### B) _diff_z
# ===== Doubly Robust summary for: grim_diff_z =====
# treat_IM                 0.0006908  0.0952099   0.007  0.99425  
# ===== Doubly Robust summary for: vidalbralo_diff_z =====
# treat_IM                -0.057692   0.072183  -0.799  0.42925  

### C) _residual_z
# ===== Doubly Robust summary for: grim_residual_z =====
#   treat_IM                 0.002235   0.105789   0.021   0.9833 
# 
# ===== Doubly Robust summary for: vidalbralo_residual_z =====
#   treat_IM                -0.087244   0.112185  -0.778  0.44170  


###---------------------------  Sustained (to Medium) vs No Caregiving ---------------------------


# 0. Prepare analysis sample: Sustained (to Medium) vs No Caregiving
dat_su_med <- subset(df3, care_trans %in% c("Sustained (to Medium)", "No Caregiving"))
dat_su_med$treat_su_med <- ifelse(dat_su_med$care_trans == "Sustained (to Medium)", 1, 0)

# 1. Estimate GBM propensity weights for ATT
set.seed(123)
w.out_su_med <- weightit(
  treat_su_med ~ r11age + female + racecat +
    r11marcat_b + edu + h11itot1 + r11shlt +
    r11smoke + r11drinkn + r11mobila +
    r11chronic + r11acts_w + r11bmi,
  data               = dat_su_med,
  method             = "gbm",
  estimand           = "ATT",
  focal              = 1,
  n.trees            = 2500,
  shrinkage          = 0.01,
  interaction.depth  = 2,
  criterion          = "smd.max",
  abs                = TRUE,
  verbose            = FALSE
)

# 2. Diagnostic Pipeline
# 2.1 Balance plot (SMD vs. number of trees)
plot(w.out_su_med, type = "balance")
optimal_trees_su_med <- w.out_su_med$info$best.tree
print(optimal_trees_su_med)

  

# 2.2 Positivity (Overlap: “common support") assumption
dat_su_med$ps <- w.out_su_med$ps
p_overlap_su_med <-ggplot(dat_su_med, aes(x = ps, fill = as.factor(treat_su_med))) +
  geom_density(alpha = 0.5) +
  labs(x = "Propensity Score", y = "Density", fill = "Treatment Group") +
  theme_minimal()  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

        ggsave("Figure_S4_su_med_Propensity_Overlap.pdf",
               plot = p_overlap_su_med,
               device = pdf,
               width = 6, height = 4, units = "in")
        
        # 2.2.1 Calculate overlap and trim based on common support for Sustained (to Medium)
        ps_t_su_med <- dat_su_med$ps[dat_su_med$treat_su_med == 1]   # Treated (Sustained to Medium)
        ps_c_su_med <- dat_su_med$ps[dat_su_med$treat_su_med == 0]   # Controls (No Caregiving)
        
        lower_su_med <- max(min(ps_t_su_med), min(ps_c_su_med))
        upper_su_med <- min(max(ps_t_su_med), max(ps_c_su_med))
        
        cat("Overlap PS before trimming:", round(lower_su_med, 3), "-", round(upper_su_med, 3), "\n")
        # Overlap PS before trimming: 0.01 - 0.481 
        retained_pct_su_med <- mean(w.out_su_med$ps >= lower_su_med & w.out_su_med$ps <= upper_su_med) * 100
        cat("Percent retained before trimming:", round(retained_pct_su_med, 1), "%\n")
        # Percent retained before trimming: 84.1 %

summary(w.out_su_med$ps)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00265 0.01260 0.04545 0.08977 0.13998 0.64557 
summary(w.out_su_med$weights)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.002657 0.012759 0.049559 0.172387 0.175381 1.000000

summary(w.out_su_med)
# Summary of weights
# 
# - Weight ranges:
#   
#   Min                                    Max
# treated 1.0000                              || 1.0000
# control 0.0027   |-------------------------|   0.9265
# 
# - Units with the 5 most extreme weights by group:
#   
#   166    120     79     38     23
# treated      1      1      1      1      1
# 2366   2498   2097   2338   2305
# control 0.7867 0.8034 0.8651 0.9095 0.9265
# 
# - Weight statistics:
#   
#   Coef of Var   MAD Entropy # Zeros
# treated        0.00 0.000   0.000       0
# control        1.27 0.946   0.636       0
# 
# - Effective Sample Sizes:
#   
#   Control Treated
# Unweighted 2412.       237
# Weighted    923.22     237


          # 2.2.2  Histogram of ATT weights with 99th-pct line -------------
          dat_su_med$w_att <- get.w(w.out_su_med)           # already computed earlier
          p_weights_su_med <- ggplot(dat_su_med,
                                     aes(x = w_att)) +
            geom_histogram(bins = 40, colour = "white",
                           fill = "steelblue", alpha = .7) +
            geom_vline(aes(xintercept = quantile(w_att, 0.99)),
                       linetype = "dashed", linewidth = .8) +
            annotate("text",
                     x      = quantile(dat_su_med$w_att, 0.99),
                     y      = Inf, vjust = -0.3, hjust = 1.1,
                     angle  = 90,
                     label  = paste0("99th pct = ",
                                     round(quantile(dat_su_med$w_att, .99), 2))) +
            labs(x = "ATT weight",
                 y = "Count",
                 title = "D Sustained at (or to) medium") +
            theme_minimal() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())
          
          ggsave("Figure_su_med_WEIGHT_Histogram.pdf",
                 plot   = p_weights_su_med,
                 device = pdf, width = 6, height = 4, units = "in")


wt_ctl_su_med <- w.out_su_med$weights[dat_su_med$treat_su_med == 0]
cut99_ctl_su_med <- quantile(wt_ctl_su_med, 0.99)
share_top1_ctl_su_med <- sum(wt_ctl_su_med[wt_ctl_su_med > cut99_ctl_su_med]) / sum(wt_ctl_su_med)
cat("Top‑1% weight share among controls (Sustained‑to‑Medium):",
    round(share_top1_ctl_su_med * 100, 1), "%\n")
# op‑1% weight share among controls (Sustained‑to‑Medium): 7.2 %

# 2.3. Covariate Balance Diagnostics
bal_su_med <- bal.tab(w.out_su_med, un = TRUE)
      # Balance Measures
      # Type Diff.Un Diff.Adj
      # prop.score                Distance  1.2749   0.4840
      # r11age                     Contin. -1.6898  -0.0697
      # female                      Binary  0.1596   0.0220
      # racecat_White               Binary -0.0949   0.0080
      # racecat_Black               Binary  0.0783   0.0067
      # racecat_Hispanic            Binary  0.0127  -0.0129
      # racecat_Other               Binary  0.0039  -0.0017
      # r11marcat_b_Married         Binary -0.0098   0.0008
      # edu_Less than High School   Binary -0.0665  -0.0133
      # edu_High School Graduate    Binary -0.0825  -0.0173
      # edu_Some College            Binary  0.1000   0.0073
      # edu_College Degree          Binary  0.0490   0.0233
      # h11itot1                   Contin.  0.0765   0.0298
      # r11shlt                    Contin. -0.1559  -0.0451
      # r11smoke_Never smoker       Binary  0.0054  -0.0056
      # r11smoke_Past smoker        Binary -0.0458   0.0004
      # r11smoke_Current smoker     Binary  0.0404   0.0051
      # r11drinkn                  Contin.  0.0778   0.0225
      # r11mobila                  Contin. -0.1997  -0.0095
      # r11chronic                 Contin. -0.3610  -0.0131
      # r11acts_w                  Contin.  0.2466   0.0116
      # r11bmi                     Contin.  0.1203   0.0043
      # 
      # Effective sample sizes
      # Control Treated
      # Unadjusted 2412.       237
      # Adjusted    923.22     237

# 2.4. Residual analysis
ps_model_su_med <- glm(treat_su_med ~ r11age + female + black + hisp + other +
                         r11marcat_b + edu + h11itot1 + r11shlt +
                         r11smoke + r11drinkn + r11mobila +
                         r11chronic + r11acts_w + r11bmi,
                       data = dat_su_med,
                       family = binomial)
dat_su_med$residuals <- residuals(ps_model_su_med, type = "deviance")

ggplot(dat_su_med, aes(x = fitted(ps_model_su_med), y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  labs(x = "Fitted Propensity Scores", y = "Deviance Residuals") +
  theme_minimal()

ggplot(dat_su_med, aes(x = residuals)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  labs(x = "Deviance Residuals", y = "Frequency") +
  theme_minimal()

# 3.1. SMD plot (left)
love.plot(bal_su_med, threshold = 0.1)

# 3.2. KS plot (right, no y labels)
p_smd_su_med <- love.plot(
  w.out_su_med,
  stats = "mean.diffs",
  threshold = 0.1,
  abs = TRUE,
  line = TRUE,
  var.order = "unadjusted",
  drop.distance = TRUE,
  colors = c("#F8766D", "#00BFC4"),
  shapes = c(19, 18),
  sample.names = c("Unweighted", "Weighted"),
  var.names     = var_labs,
  title = " "
) +
  scale_x_continuous(limits = c(0, 2.0)) +
  xlab("Absolute Standardized Mean Differences")

p_ks_su_med <- love.plot(
  x             = w.out_su_med,
  stats         = "ks.statistics",
  thresholds    = c(ks.statistics = 0.1),
  abs           = TRUE,
  line          = TRUE,
  drop.distance = TRUE,
  var.order     = "unadjusted",
  colors        = c("#F8766D", "#00BFC4"),
  shapes        = c(19, 18),
  sample.names  = c("Unweighted", "Weighted"),
  disp.ks       = TRUE,
  un            = TRUE,
  var.names     = var_labs,
  title         = " ",
  xlab          = "KS Statistic"
) +
  scale_x_continuous(limits = c(0, 0.5)) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

combined_plot_su_med <- (p_smd_su_med | p_ks_su_med) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

combined_plot_su_med
    ggsave("Figure_4_su_med_SMD and KS.pdf",
           plot = combined_plot_su_med,
           device = pdf,
           width = 9, height = 5, units = "in")

# 4. Variance Ratios & Distributional Diagnostics
bal.tab(w.out_su_med, stats = "variance.ratios")
bal.tab(w.out_su_med, stats = "ks.statistics")

# Part II Modeling in weighted dataset of Sustained (to Medium) & No-Care
dat_su_med$w_att <- get.w(w.out_su_med)
dat_su_med$comb_w <- dat_su_med$w_att * dat_su_med$vbsi16wgtra

design_su_med <- svydesign(
  ids = ~secu,
  strata = ~stratum,
  weights = ~comb_w,
  data = dat_su_med,
  nest = TRUE
)
# retrieve dv_list
att_models_su_med <- setNames(
  lapply(dv_list, function(dv) {
    svyglm(
      formula = as.formula(paste(dv, "~ treat_su_med")),
      design = design_su_med
    )
  }),
  dv_list
)

rhs_dr <- paste("treat_su_med +", paste(cov_w11, collapse = " + "))
dr_models_su_med <- setNames(
  lapply(dv_list, function(dv) {
    svyglm(
      formula = as.formula(paste(dv, "~", rhs_dr)),
      design = design_su_med
    )
  }),
  dv_list
)

for (dv in dv_list) {
  cat("\n===== ATT-only summary for:", dv, "=====\n")
  print(summary(att_models_su_med[[dv]]))
}

for (dv in dv_list) {
  cat("\n===== Doubly Robust summary for:", dv, "=====\n")
  print(summary(dr_models_su_med[[dv]]))
}

### A) _residual1_z
# ===== Doubly Robust summary for: zhang_dnamage_z =====
# treat_su_med            -0.082221   0.073171  -1.124  0.26838 

# ===== Doubly Robust summary for: grim_residual1_z =====
# treat_su_med            -0.116801   0.051761  -2.257 0.030031 * 

# ===== Doubly Robust summary for: vidalbralo_residual1_z =====
# treat_su_med            -0.068205   0.071234  -0.957   0.3445      

### B) _diff_z
# ===== Doubly Robust summary for: grim_diff_z =====
# treat_su_med            -0.093222   0.045866  -2.032 0.049328 * 

#===== Doubly Robust summary for: vidalbralo_diff_z =====
# treat_su_med            -0.038164   0.044092  -0.866 0.392319 

### C) _residual_z
# ===== Doubly Robust summary for: grim_residual_z =====
# treat_su_med            -0.106277   0.050938  -2.086 0.043890 * 

# ===== Doubly Robust summary for: vidalbralo_residual_z =====
# treat_su_med            -0.0643352  0.0687700  -0.936 0.355590

###--------------------------- Sustained (to High) vs No Caregiving ---------------------------        
# 0. Prepare analysis sample: Sustained (to High) vs No Caregiving
dat_su_hi <- subset(df3, care_trans %in% c("Sustained (to High)", "No Caregiving"))
dat_su_hi$treat_su_hi <- ifelse(dat_su_hi$care_trans == "Sustained (to High)", 1, 0)

# 1. Estimate GBM propensity weights for ATT
set.seed(123)
w.out_su_hi <- weightit(
  treat_su_hi ~ r11age + female + racecat +
    r11marcat_b + edu + h11itot1 + r11shlt +
    r11smoke + r11drinkn + r11mobila +
    r11chronic + r11acts_w + r11bmi,
  data               = dat_su_hi,
  method             = "gbm",
  estimand           = "ATT",
  focal              = 1,
  n.trees            = 2500,
  shrinkage          = 0.01,
  interaction.depth  = 2,
  criterion          = "smd.max",
  abs                = TRUE,
  verbose            = FALSE
)

# 2. Diagnostic Pipeline
# 2.1 Balance plot (SMD vs. number of trees)
plot(w.out_su_hi, type = "balance")
optimal_trees_su_hi <- w.out_su_hi$info$best.tree
print(optimal_trees_su_hi) # 2206

# 2.2 Positivity (Overlap: “common support") assumption
dat_su_hi$ps <- w.out_su_hi$ps
p_overlap_su_hi <-  ggplot(dat_su_hi, aes(x = ps, fill = as.factor(treat_su_hi))) +
  geom_density(alpha = 0.5) +
  labs(x = "Propensity Score", y = "Density", fill = "Treatment Group") +
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

  ggsave("Figure_S5_su_hi_Propensity_Overlap.pdf",
         plot = p_overlap_su_hi,
         device = pdf,
         width = 6, height = 4, units = "in")
  
        # 2.2.1 . Calculate overlap and trim based on common support
  # 2.2.1 Calculate overlap and trim based on common support for Sustained (to High)
  ps_t_su_hi <- dat_su_hi$ps[dat_su_hi$treat_su_hi == 1]   # Treated (Sustained to High)
  ps_c_su_hi <- dat_su_hi$ps[dat_su_hi$treat_su_hi == 0]   # Controls (No Caregiving)
  
  lower_su_hi <- max(min(ps_t_su_hi), min(ps_c_su_hi))
  upper_su_hi <- min(max(ps_t_su_hi), max(ps_c_su_hi))
  
  cat("Overlap PS before trimming:", round(lower_su_hi, 3), "-", round(upper_su_hi, 3), "\n")
  # Overlap PS before trimming: 0.009 - 0.431 
  cat("Percent retained before trimming:",
      round(mean(dat_su_hi$ps >= lower_su_hi & dat_su_hi$ps <= upper_su_hi) * 100, 1), "%\n")
  # Percent retained before trimming: 71.7 %
  
        
summary(w.out_su_hi$ps)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001155 0.008652 0.019225 0.045504 0.051471 0.675495 
summary(w.out_su_hi$weights)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001156 0.008727 0.019655 0.084796 0.055446 1.000000
summary(w.out_su_hi)
      # Summary of weights
      # 
      # - Weight ranges:
      #   
      #   Min                                    Max
      # treated 1.0000                              || 1.0000
      # control 0.0012   |--------------------|        0.7587
      # 
      # - Units with the 5 most extreme weights by group:
      #   
      #   189    140    118   107     41
      # treated      1      1      1     1      1
      # 433    924    670  1578   2469
      # control 0.4805 0.5331 0.5826 0.739 0.7587
      # 
      # - Weight statistics:
      #   
      #   Coef of Var   MAD Entropy # Zeros
      # treated       0.000 0.000   0.000       0
      # control       1.503 0.926   0.673       0
      # 
      # - Effective Sample Sizes:
      #   
      #   Control Treated
      # Unweighted 2412.       115
      # Weighted    740.25     115

      # 2.2.2  Histogram of ATT weights with 99th-pct line (Sustained to High)
      dat_su_hi$w_att <- get.w(w.out_su_hi)
      p_weights_su_hi <- ggplot(dat_su_hi, aes(x = w_att)) +
        geom_histogram(bins = 40, colour = "white", fill = "steelblue", alpha = .7) +
        geom_vline(aes(xintercept = quantile(w_att, 0.99)),
                   linetype = "dashed", size = .8) +
        annotate("text",
                 x = quantile(dat_su_hi$w_att, 0.99),
                 y = Inf, vjust = -0.3, hjust = 1.1, angle = 90,
                 label = paste0("99th pct = ", round(quantile(dat_su_hi$w_att, .99), 2))) +
        labs(x = "ATT weight",
             y = " ",
             title = "E Sustained at (or to) high") +
        theme_minimal() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      
      ggsave("Figure_su_hi_WEIGHT_Histogram.pdf",
             plot = p_weights_su_hi,
             device = pdf, width = 6, height = 4, units = "in")


# 2.3. Covariate Balance Diagnostics
bal_su_hi <- bal.tab(w.out_su_hi, un = TRUE)
    # Balance Measures
    # Type Diff.Un Diff.Adj
    # prop.score                Distance  1.2035   0.7606
    # r11age                     Contin.  0.1141   0.0110
    # female                      Binary -0.0136   0.0086
    # racecat_White               Binary -0.0245  -0.0129
    # racecat_Black               Binary -0.0060   0.0137
    # racecat_Hispanic            Binary  0.0429   0.0045
    # racecat_Other               Binary -0.0125  -0.0054
    # r11marcat_b_Married         Binary  0.3165   0.0146
    # edu_Less than High School   Binary  0.0847   0.0150
    # edu_High School Graduate    Binary  0.0290   0.0044
    # edu_Some College            Binary  0.0059  -0.0062
    # edu_College Degree          Binary -0.1196  -0.0132
    # h11itot1                   Contin. -0.1120  -0.0238
    # r11shlt                    Contin.  0.0681   0.0052
    # r11smoke_Never smoker       Binary -0.0287  -0.0266
    # r11smoke_Past smoker        Binary -0.0382   0.0084
    # r11smoke_Current smoker     Binary  0.0669   0.0183
    # r11drinkn                  Contin. -0.3700  -0.0576
    # r11mobila                  Contin.  0.1270   0.0570
    # r11chronic                 Contin.  0.0953   0.0203
    # r11acts_w                  Contin. -0.1877  -0.0192
    # r11bmi                     Contin.  0.1408  -0.0001
    # 
    # Effective sample sizes
    # Control Treated
    # Unadjusted 2412.       115
    # Adjusted    740.25     115

# Top‑1% weight share among controls (Sustained to High vs No Caregiving)
wt_ctl_su_hi <- w.out_su_hi$weights[dat_su_hi$treat_su_hi == 0]

# Compute the 99th percentile cutoff
cut99_ctl_su_hi <- quantile(wt_ctl_su_hi, 0.99, na.rm = TRUE)

# Proportion of total control weight contributed by the top 1%
share_top1_ctl_su_hi <- sum(wt_ctl_su_hi[wt_ctl_su_hi > cut99_ctl_su_hi], na.rm = TRUE) / sum(wt_ctl_su_hi, na.rm = TRUE)

# Print the result
cat("Top‑1% weight share among controls (Sustained to High vs No Caregiving):",
    round(share_top1_ctl_su_hi * 100, 1), "%\n")
# Top‑1% weight share among controls: 10.3 %

# 2.4. Residual analysis
ps_model_su_hi <- glm(treat_su_hi ~ r11age + female + black + hisp + other +
                        r11marcat_b + edu + h11itot1 + r11shlt +
                        r11smoke + r11drinkn + r11mobila +
                        r11chronic + r11acts_w + r11bmi,
                      data = dat_su_hi,
                      family = binomial)
dat_su_hi$residuals <- residuals(ps_model_su_hi, type = "deviance")

ggplot(dat_su_hi, aes(x = fitted(ps_model_su_hi), y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  labs(x = "Fitted Propensity Scores", y = "Deviance Residuals") +
  theme_minimal()

ggplot(dat_su_hi, aes(x = residuals)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  labs(x = "Deviance Residuals", y = "Frequency") +
  theme_minimal()

# 3.1. SMD plot (left)
love.plot(bal_su_hi, threshold = 0.1)

# 3.2. KS plot (right, no y labels)
p_smd_su_hi <- love.plot(
  w.out_su_hi,
  stats = "mean.diffs",
  threshold = 0.1,
  abs = TRUE,
  line = TRUE,
  var.order = "unadjusted",
  drop.distance = TRUE,
  colors = c("#F8766D", "#00BFC4"),
  shapes = c(19, 18),
  sample.names = c("Unweighted", "Weighted"),
  var.names     = var_labs,
  title = " "
) +
  scale_x_continuous(limits = c(0, 0.5)) +
  xlab("Absolute Standardized Mean Differences")

p_ks_su_hi <- love.plot(
  x             = w.out_su_hi,
  stats         = "ks.statistics",
  thresholds    = c(ks.statistics = 0.1),
  abs           = TRUE,
  line          = TRUE,
  drop.distance = TRUE,
  var.order     = "unadjusted",
  colors        = c("#F8766D", "#00BFC4"),
  shapes        = c(19, 18),
  sample.names  = c("Unweighted", "Weighted"),
  disp.ks       = TRUE,
  un            = TRUE,
  var.names     = var_labs,
  title         = " ",
  xlab          = "KS Statistic"
) +
  scale_x_continuous(limits = c(0, 0.5)) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

combined_plot_su_hi <- (p_smd_su_hi | p_ks_su_hi) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

combined_plot_su_hi
      ggsave("Figure_5_su_hi_SMD and KS.pdf",
             plot = combined_plot_su_hi,
             device = pdf,
             width = 9, height = 5, units = "in")

# 4. Variance Ratios & Distributional Diagnostics
bal.tab(w.out_su_hi, stats = "variance.ratios")
bal.tab(w.out_su_hi, stats = "ks.statistics")

# Part II Modeling in weighted dataset of Sustained (to High) & No-Care
dat_su_hi$w_att <- get.w(w.out_su_hi)
dat_su_hi$comb_w <- dat_su_hi$w_att * dat_su_hi$vbsi16wgtra

design_su_hi <- svydesign(
  ids = ~secu,
  strata = ~stratum,
  weights = ~comb_w,
  data = dat_su_hi,
  nest = TRUE
)

att_models_su_hi <- setNames(
  lapply(dv_list, function(dv) {
    svyglm(
      formula = as.formula(paste(dv, "~ treat_su_hi")),
      design = design_su_hi
    )
  }),
  dv_list
)


rhs_dr <- paste("treat_su_hi +", paste(cov_w11, collapse = " + "))
dr_models_su_hi <- setNames(
  lapply(dv_list, function(dv) {
    svyglm(
      formula = as.formula(paste(dv, "~", rhs_dr)),
      design = design_su_hi
    )
  }),
  dv_list
)

for (dv in dv_list) {
  cat("\n===== ATT-only summary for:", dv, "=====\n")
  print(summary(att_models_su_hi[[dv]]))
}

for (dv in dv_list) {
  cat("\n===== Doubly Robust summary for:", dv, "=====\n")
  print(summary(dr_models_su_hi[[dv]]))
}

### A) residual1_z
# ===== Doubly Robust summary for: zhang_dnamage_z =====
# treat_su_hi             -0.011571   0.074188  -0.156  0.87690   

# ===== Doubly Robust summary for: grim_residual1_z =====
# treat_su_hi             -0.113002   0.125868  -0.898  0.37511 

# ===== Doubly Robust summary for: vidalbralo_residual1_z =====
# treat_su_hi             0.028779   0.162276   0.177  0.86020        

# ===== Doubly Robust summary for: horvath_residual1_z =====        
# treat_su_hi             -0.226886   0.104884  -2.163   0.0371 *       

### B) _diff_z
# ===== Doubly Robust summary for: grim_diff_z =====       
# treat_su_hi             -0.1015981  0.1077614  -0.943 0.351896          

# ===== Doubly Robust summary for: vidalbralo_diff_z =====  
# treat_su_hi              0.025361   0.102560   0.247  0.80606   

# ===== Doubly Robust summary for: horvath_diff_z 
# treat_su_hi             -0.2323724  0.0962433  -2.414  0.02082 * 

### C) _residual_z
# ===== Doubly Robust summary for: grim_residual_z =====
# treat_su_hi             -0.115078   0.120206  -0.957 0.344612   
# 
# 
# ===== Doubly Robust summary for: vidalbralo_residual_z =====
# treat_su_hi              0.036997   0.161719   0.229  0.82030  
# 
# ===== Doubly Robust summary for: horvath_residual_z =====
# treat_su_hi             -2.502e-01  1.033e-01  -2.423   0.0204 *


#--------------------------------------------------------------------------------------------------
### Part IV Pair-wise ATT diagnostics for  Initiated-High vs No-Care 
#--------------------------------------------------------------------------------------------------

### 1 Combined plot
### Put 5 overlap graphs together
theme_overlap <- theme_minimal(base_size = 10) +
  theme(
    panel.grid   = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 8),
    plot.title   = element_text(face = "bold", hjust = 0),
    legend.title = element_blank()
  )

fill_scale <- scale_fill_manual(
  values = c("0" = "#F8766D", "1" = "#00BFC4"),
  labels = c("0" = "Control", "1" = "Treated"),
  name   = "Group"
)  


p1 <- p_overlap_IM     +
  ggtitle("A  Transition to medium intensity care") +
  theme_overlap +
  fill_scale +
  theme(legend.position = "none")

p2 <- p_overlap_IH     +
  ggtitle("B  Transition to high intensity care") +
  theme_overlap +
  fill_scale +
  theme(legend.position = "none")

p3 <- p_overlap_CE     +
  ggtitle("C  Transition out of care") +
  theme_overlap +
  fill_scale +
  theme(legend.position = "none")

p4 <- p_overlap_su_med +
  ggtitle("D  Sustained at (or to) medium") +
  theme_overlap +
  fill_scale +
  theme(legend.position = "none")

p5 <- p_overlap_su_hi  +
  ggtitle("E  Sustained at (or to) high") +
  theme_overlap +
  fill_scale +
  theme(legend.position = "none")

              #  Assemble final plot with legend in panel 
              final_plot <- 
                p1 + p2 + p3 +
                p4 + p5 + guide_area() +       
                plot_layout(ncol = 3, guides = "collect") &
                theme(legend.position = "bottom") &  
                theme(plot.margin = margin(5,5,5,5)) &
                labs(x = "Propensity score", y = "Density")


                ## Export at journal-ready resolution 
                ggsave("FigureS_overlap.png", final_plot, width = 235, height = 130, units = "mm", dpi = 300)
                
                ggsave("FigureS_overlap_plots.pdf",
                       plot = final_plot,
                       device = pdf,
                       width = 10, height = 5, units = "in")

### 2) Combined plot 2
### Put 5 weights distribution graphs together
                dummy_plot <- ggplot() +
                  theme_void() +
                  annotate("text", x = 0.5, y = 0.5,
                           label = "Dashed line: 99th percentile",
                           size = 5, fontface = "italic", color = "steelblue",
                           hjust = 0.5, vjust = 0.5)
                
                combined_weights <- (
                  p_weights_IM  +
                    p_weights_IH  +
                    p_weights_CE  +
                    p_weights_su_med +
                    p_weights_su_hi +
                    dummy_plot
                ) +
                  plot_layout(ncol = 3, nrow = 2, byrow = TRUE, guides = 'collect') +
                  plot_annotation(
                    title = 'Distribution of ATT Weights Across Caregiving Transition Groups',
                    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
                  )
                
                ggsave("Figure_ALL_WEIGHT_Histograms.png", combined_weights, width = 235, height = 130, units = "mm", dpi = 300)
                
                ggsave("Figure_ALL_WEIGHT_Histograms_ANNOTATED.pdf",
                       plot = combined_weights, 
                       device = pdf, 
                       width = 16, height = 8, units = "in")
