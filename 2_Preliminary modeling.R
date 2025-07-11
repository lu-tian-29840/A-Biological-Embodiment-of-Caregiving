
#@@ Author: Lu Tian
#@@ Date: May 16, 2025
#@@ Project name: Caregiving and bioage
#    Section 2/4  : OLS Modeling -- preliminary analysis
#           Part I:  1) DNAmAge_residual1_z ~ care_trans1 + wave 11 controls
#                    2) DNAmAge_residual_z ~ care_trans1 + wave 11 controls
#                    3) DNAmAge_diff_z~ care_trans1 + wave 11 controls
#         
# 
#           Part II: 1) DNAmAge_residual1_z ~ care_trans + wave 11 controls
#                    2) DNAmAge_residual_z ~ care_trans + wave 11 controls
#                    3) DNAmAge_diff_z~ care_trans + wave 11 controls

### Part I
### care_trans1 as the IV of interest

# table(df3$care_trans1)
# No Caregiving Initiated (Medium)   Initiated (High)             Ceased          Sustained 
# 2412                109                 94                186                352 
# table(df3$care_trans)
# No Caregiving   Initiated (Medium)  Initiated (High)  Ceased Sustained (to Medium) Sustained (to High) 
# 2412            109                 94                186    237                   115 

# Step 0
# DV list
names(df3)
dv_list <- c("horvath_dnamage",     "hannum_dnamage",   "levine_dnamage",
             "horvathskin_dnamage", "lin_dnamage",      "weidner_dnamage",
             "vidalbralo_dnamage",  "grim_dnamage",     "zhang_dnamage", 
             "yang_dnamage",       "bocklandt_dnamage", "garagnani_dnamage", 
             "mpoa")

dv_list <- c("horvath_residual_z",     "hannum_residual_z",    "levine_residual_z", 
             "horvathskin_residual_z", "lin_residual_z",       "weidner_residual_z", 
             "vidalbralo_residual_z",  "grim_residual_z",      "zhang_residual_z", 
             "yang_residual_z",       "bocklandt_residual_z", "garagnani_residual_z", 
             "mpoa_z")

dv_list <- c("horvath_diff_z",     "hannum_diff_z",    "levine_diff_z", 
             "horvathskin_diff_z", "lin_diff_z",       "weidner_diff_z", 
             "vidalbralo_diff_z",  "grim_diff_z",      "zhang_diff_z", 
             "yang_diff_z",        "bocklandt_diff_z", "garagnani_diff_z", 
             "mpoa_z")

dv_list <- c("horvath_residual1_z",      "hannum_residual1_z",    "levine_residual1_z", 
             "horvathskin_residual1_z",  "lin_residual1_z",       "weidner_residual1_z", 
             "vidalbralo_residual1_z",   "grim_residual1_z",      "zhang_residual1_z", 
             "yang_residual1_z",        "bocklandt_residual1_z", "garagnani_residual1_z", 
             "mpoa_z")

# Define rhs1 
rhs1 <- "r11age + female + black + hisp + other + r11marcat_b + edu + h11itot1 + r11shlt + r11smoke + r11drinkn + r11mobila + r11chronic + r11acts_w +r11bmi"

# 1) Define the survey design:
d.vbs <- svydesign(
  ids     = ~secu,            # clustering 
  strata  = ~stratum,         # stratification 
  weights = ~vbsi16wgtra,     # biomarker weight
  data    = df3,
  nest    = TRUE
)

# 2) Fit each model with svyglm():
model_list <- setNames(
  lapply(dv_list, function(dv) {
    svyglm(
      formula = as.formula(paste(dv, "~ care_trans1 +", rhs1)),
      design  = d.vbs
    )
  }),
  dv_list
)

# 3) Print design‐corrected summaries:
for (dv in dv_list) {
  cat("\n===== Summary for:", dv, "=====\n")
  print(summary(model_list[[dv]]))
}


### Part II
### care_trans as the IV of interest
# 1) Fit each model with svyglm():
model_list <- setNames(
  lapply(dv_list, function(dv) {
    svyglm(
      formula = as.formula(paste(dv, "~ care_trans +", rhs1)),
      design  = d.vbs
    )
  }),
  dv_list
)

# 2) Print design‐corrected summaries:
for (dv in dv_list) {
  cat("\n===== Summary for:", dv, "=====\n")
  print(summary(model_list[[dv]]))
}

### A） residual1_z
# ===== Summary for: zhang_dnamage_z =====
# care_transInitiated (High)       0.3215599  0.1275639   2.521 0.016724 *
#   
# ===== Summary for: grim_residual1_z =====
# care_transInitiated (High)       0.214121   0.074694   2.867 0.007171 ** 
# care_transSustained (to Medium) -0.103169   0.051891  -1.988 0.055136 .  

# ===== Summary for: vidalbralo_residual1_z ===== 
# care_transCeased                -0.1887346  0.0718918  -2.625 0.013020 * 

# ===== Summary for: horvath_residual1_z =====
# care_transSustained (to High)   -0.246077   0.126343  -1.948 0.060002 . 

### B) _diff_z
# ===== Summary for: zhang_dnamage_z =====
# care_transInitiated (High)       0.0150230  0.0070465    2.132   0.0405 * 

# ===== Summary for: grim_diff_z =====
# care_transInitiated (High)       0.179681   0.068328   2.630 0.012881 *  

# ===== Summary for: vidalbralo_diff_z =====
# care_transCeased                -0.0857520  0.0455206  -1.884  0.06843 .  

# ===== Summary for: horvath_diff_z =====
# care_transSustained (to High)   -0.245444   0.118004  -2.080 0.045372 *   

### C) residual_z
# ===== Summary for: zhang_residual_z =====
# care_transInitiated (High)       0.271172   0.114545   2.367  0.02393 *  

# ===== Summary for: grim_residual_z =====
# care_transInitiated (High)       0.200114   0.076267   2.624 0.013064 * 

# ===== Summary for: vidalbralo_residual_z =====
# care_transCeased                -0.131199   0.071216  -1.842 0.074437 . 

# ===== Summary for: horvath_residual_z =====
# care_transSustained (to High)   -0.263287   0.126543  -2.081 0.045308 *  




# rhs2 (independent variables)
rhs2 <- "r12carecat + r11carecat + r11age + female + black + hisp + other + r11marcat + edu + h11itot1 + r11shlt + r11smoke + r11drinkn + r11mobila + r11chronic + r11acts_w +r11bmi"

# Store models in a list
model_list <- list()

# Loop through each DV and run lm()
for (dv in dv_list) {
  formula <- as.formula(paste(dv, "~", rhs2))
  model_list[[dv]] <- lm(formula, data = df3)
}

# Print summaries for all models
for (dv in dv_list) {
  cat("\n===== Summary for:", dv, "=====\n")
  print(summary(model_list[[dv]]))
}
