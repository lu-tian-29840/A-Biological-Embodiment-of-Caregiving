
#@@ Author: Lu Tian
#@@ Date: May 16, 2025
#@@ Project name: Caregiving and bioage
#    Section 1/4  : Data prepare and variable creation
#           Part I: Dependent variable construction
#           Part II: Independent variable construction
#           Part III: Controls creation
#           Parr IV: Descriptive table 1 before and after survey weighting
library(knitr)        # rendering .Rmd, tables 
library(tidyverse)    # loads dplyr, ggplot2, tidyr, tibble, etc. see tidyverse.org/packages
library(Hmisc)     
library(haven)     # read .dta file
library(kableExtra) # Make nice table
library(psych)    # to describe variables
library(purrr)   # use list
library(labelled) 
library(dplyr)
library(survey)
library(jtools)  
library(naniar) # missing pattern
library(tableone) # make table 1
library(haven)


##########################################################################################      
### Part 0: 
### Load and prepare data sets
##########################################################################################   

### Step 0: 
### Load data

df1 <- read_dta("full_2.dta") %>% 
  as.data.frame() 

df2 <- df1 %>% 
  select(rahhidpn, hhid, pn, opn, secu, stratum, vbsi16wgtra, r11age, r12age, r13age, 
         r11marst, r12marst, ragender, rahispan, raracem, raedyrs,
         r11caregiving_flag, study_trk, r12caregiving_flag, 
         r12pcare_flag, r11pcare_flag, r11f001, r11f011, r12f001, r12f011,
         r12is_helper, r12is_spouse_helper, r12is_ADL_helper, r12is_IADL_helper, r12pPersonal_helper, r12pErrand_helper, 
         r11is_helper, r11is_spouse_helper, r11is_ADL_helper, r11is_IADL_helper, r11pPersonal_helper, r11pErrand_helper, 
         r11shlt, r11cesd, r11vgactx, r11bmi, r11smoken, r11smokev, r11drinkn, 
         r11hibpe, r11diabe, r11cancre, r11hearte, r11stroke, h11itot, r11mobila, r11adl6a, 
         r12shlt, r12cesd, r12vgactx, r12bmi, r12smoken, r12smokev, r12drinkn, 
         r12hibpe, r12diabe, r12cancre, r12hearte, r12stroke, h12itot, r12mobila, r12adl6a, 
         r13bcell_pct:mpoa, r11acts_w, r12acts_w)


names(df2)[names(df2) == "dnamgrimage"] <- "grim_dnamage"


##########################################################################################      
### Part I 
### Dependent variable construction
##########################################################################################      

## Method 1: 
## Difference between DNAmAge and calendar age  
dnamage_vars <- c(
  "horvath_dnamage", "hannum_dnamage", "levine_dnamage",
  "horvathskin_dnamage", "lin_dnamage", "weidner_dnamage",
  "vidalbralo_dnamage", "grim_dnamage", "yang_dnamage",
  "zhang_dnamage", "bocklandt_dnamage", "garagnani_dnamage"
)

# Loop through each DNAmAge variable to create diffs
for (var in dnamage_vars) {
  accelDiff <- sub("_dnamage$", "_diff", var)
  zDiff     <- paste0(accelDiff, "_z")
  
  # 1) raw difference
  df2[[accelDiff]] <- df2[[var]] - df2$r13age
  # 2) z-score of the difference
  df2[[zDiff]]     <- as.numeric(scale(df2[[accelDiff]]))
}

## Method 2:
## Take residuals (Lu et al. 2019) e.g.DNAmAge ~ r13age


for (var in dnamage_vars) {
  # build and fit the model
  f        <- as.formula(paste(var, "~ r13age"))
  fit      <- lm(f, data = df2)
  
  # name the residual and its z-score
  res_name <- sub("_dnamage$", "_residual", var)
  zRes     <- paste0(res_name, "_z")
  
  # 1) raw residual
  df2[[res_name]] <- resid(fit)
  # 2) z-score of the residual
  df2[[zRes]]     <- as.numeric(scale(df2[[res_name]]))
}

## Method 3:
## Take residuals (Kim et al. 2025) e.g. e.g.DNAmAge ~ r13age + 
#                                                       r13bcell_pct + r13cd8t_pct +  r13cd4t_pct + 
#                                                       r13cd8n_pct + r13nk_pct + r13mono_pc

# Step 1: drop missing in bio age weighting data AND 6 cell types
df2 <- df2[!is.na(df2$vbsi16wgtra), ]
vars <- c("r13bcell_pct", "r13cd8t_pct", "r13cd4t_pct",
          "r13cd8n_pct", "r13nk_pct",   "r13mono_pct")

df2 <- df2[complete.cases(df2[, vars]), ]

for (var in dnamage_vars) {
  # 1) build formula
  f1       <- as.formula(paste(
    var,
    "~ r13age + r13bcell_pct + r13cd8t_pct +",
    "r13cd4t_pct + r13cd8n_pct + r13nk_pct + r13mono_pct"
  ))
  fit1     <- lm(f1, data = df2)
  
  # 2) name the raw residual
  res_name1 <- sub("_dnamage$", "_residual1", var)
  df2[[res_name1]] <- resid(fit1)
  
  # 3) name and compute the z-score of that residual
  z_name1   <- paste0(res_name1, "_z")
  df2[[z_name1]] <- as.numeric(scale(df2[[res_name1]]))
}

### Part II
### Independent variable construction

## Approach 1 
df2 <- df2 %>%
  mutate(
    care_entered    = if_else(r11caregiving_flag == 0 & r12caregiving_flag == 1, 1, 0),
    care_exited     = if_else(r11caregiving_flag == 1 & r12caregiving_flag == 0, 1, 0),
    care_consistent = if_else(r11caregiving_flag == 1 & r12caregiving_flag == 1, 1, 0),
    care_never      = if_else(r11caregiving_flag == 0 & r12caregiving_flag == 0, 1, 0),
    # optional: a single categorical variable
    care_traj = case_when(
      care_entered    == 1 ~ "care_entered",
      care_exited     == 1 ~ "care_exited",
      care_consistent == 1 ~ "care_consistent",
      care_never      == 1 ~ "care_never",
      TRUE            ~ NA_character_
    )
  )

# Quick check
df2 %>% count(care_traj)
df2 %>% summarise(
  n_entered    = sum(care_entered),
  n_exited     = sum(care_exited),
  n_consistent = sum(care_consistent),
  n_never      = sum(care_never)
)
table(df2$r12caregiving_flag)
table(df2$r11caregiving_flag)


## 2. 
## Approach 2 (in model analysis)

df2 <- df2 %>%
  mutate(
    # 3‐level caregiving category in wave 12:
    r12carecat = case_when(
      r12caregiving_flag == 0                   ~ 1,  # No caregiving
      r12pcare_flag       == 1                   ~ 2,  # Parental errands/personal needs
      r12is_helper        == 1                   ~ 3,  # ADL or IADL help
      TRUE                                        ~ NA_real_
    ),
    # (Optionally) label them as a factor:
    r12carecat = factor(
      r12carecat,
      levels = 1:3,
      labels = c("None", "Parental", "ADL/IADL")
    )
  )


df2 <- df2 %>%
  mutate(
    # 3‐level caregiving category in wave 12:
    r11carecat = case_when(
      r11caregiving_flag == 0                   ~ 1,  # No caregiving
      r11pcare_flag       == 1                   ~ 2,  # Parental errands/personal needs
      r11is_helper        == 1                   ~ 3,  # ADL or IADL help
      TRUE                                        ~ NA_real_
    ),
    r11carecat = factor(
      r11carecat,
      levels = 1:3,
      labels = c("None", "Parental", "ADL/IADL")
    )
  )

table(df2$r11carecat)
table(df2$r12carecat)

# a) transition 1 
df2 <- df2 %>%
  mutate(
    care_trans = case_when(
      # No caregiving in both waves
      r11carecat == "None" & r12carecat == "None" ~ "No Caregiving",
      
      # Initiation of caregiving
      r11carecat == "None" & r12carecat == "Parental" ~ "Initiated (Medium)",
      r11carecat == "None" & r12carecat == "ADL/IADL" ~ "Initiated (High)",
      
      # Cessation of caregiving
      r11carecat %in% c("Parental", "ADL/IADL") & r12carecat == "None" ~ "Ceased",
      
      # Sustained caregiving
      r11carecat == "Parental" & r12carecat == "Parental" ~ "Sustained (to Medium)",
      r11carecat == "ADL/IADL" & r12carecat == "ADL/IADL" ~ "Sustained (to High)",
      
      # Changed effort transitions
      r11carecat == "Parental" & r12carecat == "ADL/IADL" ~ "Sustained (to High)",
      r11carecat == "ADL/IADL" & r12carecat == "Parental" ~ "Sustained (to Medium)",
      
      # Catch-all for unexpected combinations
      TRUE ~ NA_character_
    ),
    care_trans = factor(
      care_trans,
      levels = c(
        "No Caregiving",
        "Initiated (Medium)",
        "Initiated (High)",
        "Ceased",
        "Sustained (to Medium)",
        "Sustained (to High)"
      )
    )
  )


table(df2$care_trans)

# b) transition 2     (in model)

df2 <- df2 %>%
  mutate(
    care_trans1 = case_when(
      # No caregiving in both waves
      r11carecat == "None" & r12carecat == "None" ~ "No Caregiving",
      
      # Initiation of caregiving
      r11carecat == "None" & r12carecat == "Parental" ~ "Initiated (Medium)",
      r11carecat == "None" & r12carecat == "ADL/IADL" ~ "Initiated (High)",
      
      # Cessation of caregiving
      r11carecat %in% c("Parental", "ADL/IADL") & r12carecat == "None" ~ "Ceased",
      
      # Sustained caregiving (including intensity changes)
      r11carecat %in% c("Parental", "ADL/IADL") & r12carecat %in% c("Parental", "ADL/IADL") ~ "Sustained",
      
      # Catch-all for unexpected combinations
      TRUE ~ NA_character_
    ),
    care_trans1 = factor(
      care_trans1,
      levels = c(
        "No Caregiving",
        "Initiated (Medium)",
        "Initiated (High)",
        "Ceased",
        "Sustained"
      )
    )
  )
table(df2$care_trans1)


##########################################################################################      
##@ Part 3
##@ Controls
##########################################################################################    '

# 1). income
df2 <- df2 %>%
  mutate(
    h11itot1 = log(h11itot + sqrt(h11itot^2 + 1))
  )

df2 <- df2 %>%
  mutate(
    h12itot1 = log(h12itot + sqrt(h12itot^2 + 1))
  )

# 2) Race
df2$racecat <- ifelse(df2$raracem == 1 & df2$rahispan == 0, 1,     # non-hisp white
                      ifelse(df2$raracem == 2 & df2$rahispan == 0, 2,     # non-hisp black
                             ifelse(df2$rahispan == 1,  3,                
                                    ifelse(df2$raracem == 3 & df2$rahispan == 0, 4, NA))))    

df2$white <- ifelse(df2$racecat == 1, 1, 0)
df2$black <- ifelse(df2$racecat == 2, 1, 0)
df2$hisp  <- ifelse(df2$racecat == 3, 1, 0)
df2$other <- ifelse(df2$racecat == 4, 1, 0)

df2$racecat <- factor(df2$racecat,
                      levels = c(1, 2, 3, 4),
                      labels = c("White", "Black", "Hispanic", "Other"))

# 3) Gender
df2$female <- ifelse(df2$ragender == 2, 1, 0)


# 4) Maritccal status 
df2 <- df2 %>%
  mutate(
    r11marcat = case_when(
      r11marst == 1 ~ "Married",
      r11marst == 2 ~ "Div/Sep",
      r11marst == 3 ~ "Widow",
      r11marst == 4 ~ "NevMar",
      TRUE          ~ NA_character_
    ),
    # now turn into a plain factor with your desired reference order
    r11marcat = factor(
      r11marcat,
      levels = c("Married", "Div/Sep", "Widow", "NevMar")
    )
  )
table(df2$r11marcat)

df2 <- df2 %>%
  mutate(
    r12marcat = case_when(
      r12marst == 1 ~ "Married",
      r12marst == 2 ~ "Div/Sep",
      r12marst == 3 ~ "Widow",
      r12marst == 4 ~ "NevMar",
      TRUE          ~ NA_character_
    ),
    # now turn into a plain factor with your desired reference order
    r12marcat = factor(
      r12marcat,
      levels = c("Married", "Div/Sep", "Widow", "NevMar")
    )
  )
table(df2$r11marcat)



# 5). Drink
df2 <- df2 %>%
  mutate(
    # recode numeric values
    r11drinkn = case_when(
      r11drinkn == 0           ~ 0,
      r11drinkn %in% c(1, 2, 3) ~ 1,
      r11drinkn == 4           ~ 2,
      TRUE                      ~ NA_real_
    )
  ) %>%
  # set value labels
  set_value_labels(r11drinkn = c(
    "Doesn't drink" = 0,
    "1–4/day"       = 1,
    "5+/day"        = 2
  ))

df2 <- df2 %>%
  mutate(
    # recode numeric values
    r12drinkn = case_when(
      r12drinkn == 0           ~ 0,
      r12drinkn %in% c(1, 2, 3) ~ 1,
      r12drinkn == 4           ~ 2,
      TRUE                      ~ NA_real_
    )
  ) %>%
  # set value labels
  set_value_labels(r12drinkn = c(
    "Doesn't drink" = 0,
    "1–4/day"       = 1,
    "5+/day"        = 2
  ))


# 6). smoke
df2$r11smoke <- NA_integer_

# Current smoker
df2$r11smoke[df2$r11smoken == 1] <- 2

# Past smoker
df2$r11smoke[df2$r11smoken == 0 & df2$r11smokev == 1] <- 1

# Never smoked
df2$r11smoke[df2$r11smoken == 0 & df2$r11smokev == 0] <- 0

# Add labels
attr(df2$r11smoke, "labels") <- c(
  "Never smoked" = 0,
  "Past smoker" = 1,
  "Current smoker" = 2
)

df2$r12smoke <- NA_integer_

# Current smoker
df2$r12smoke[df2$r12smoken == 1] <- 2

# Past smoker
df2$r12smoke[df2$r12smoken == 0 & df2$r12smokev == 1] <- 1

# Never smoked
df2$r12smoke[df2$r12smoken == 0 & df2$r12smokev == 0] <- 0

# Add labels
attr(df2$r12smoke, "labels") <- c(
  "Never smoked" = 0,
  "Past smoker" = 1,
  "Current smoker" = 2
)

# 7). Number of health conditions. 
df2 <- df2 %>%
  mutate(across(c(r11hibpe, r11diabe, r11cancre, r11hearte, r11stroke),
                ~ ifelse(is.na(.), 0, .), .names = "clean_{col}")) %>%
  rowwise() %>%
  mutate(r11chronic = sum(c_across(c(r11hibpe, r11diabe, r11cancre, r11hearte, r11stroke)) == 1, na.rm = TRUE)) %>%
  ungroup()

df2 <- df2 %>%
  mutate(across(c(r12hibpe, r12diabe, r12cancre, r12hearte, r12stroke),
                ~ ifelse(is.na(.), 0, .), .names = "clean_{col}")) %>%
  rowwise() %>%
  mutate(r12chronic = sum(c_across(c(r12hibpe, r12diabe, r12cancre, r12hearte, r12stroke)) == 1, na.rm = TRUE)) %>%
  ungroup()

# 8. education 
df2 <- df2 %>%
  mutate(
    edu = case_when(
      raedyrs <= 11 ~ "Less than High School",
      raedyrs == 12 ~ "High School Graduate",
      raedyrs >= 13 & raedyrs <= 15 ~ "Some College",
      raedyrs >= 16 ~ "College Degree",
      TRUE ~ NA_character_  
    ),
    edu = factor(edu, levels = c("Less than High School", "High School Graduate", 
                                 "Some College", "College Degree"))
  ) %>%
  mutate(
    lshs = if_else(edu == "Less than High School", 1, 0),
    hs = if_else(edu == "High School Graduate", 1, 0),
    somecoll = if_else(edu == "Some College", 1, 0),
    coll = if_else(edu == "College Degree", 1, 0)
  )


### Part 4
### Choose variables
# df3 <- df2 %>% 
#   select(rahhidpn:r13age, raedyrs,  care_trans, r11caregiving_flag:r12caregiving_flag, 
#          r11f001, r11f011, r12f001, r12f011, 
#          r11shlt, r11cesd, r11bmi, r11drinkn, r11mobila, r11smoken, 
#          r12shlt, r12cesd, r12bmi, r12drinkn, r12mobila,
#          horvath_dnamage:garagnani_residual1_z, r12carecat, r11carecat, care_trans1, h11itot1, h12itot1, 
#          racecat:r12smoke, r12chronic, r11chronic:coll) %>% 
#   as.data.frame() 

df3 <- df2 %>% 
  select(rahhidpn:r13age, raedyrs,  care_trans, r11caregiving_flag:r12caregiving_flag, 
         r11shlt, r11cesd, r11bmi, r11drinkn, r11mobila, r11smoken, 
         r12shlt, r12cesd, r12bmi, r12drinkn, r12mobila,
         horvath_dnamage:garagnani_residual1_z, r12carecat, r11carecat, care_trans1, h11itot1, h12itot1, 
         racecat:r12smoke, r12chronic, r11chronic:coll) %>% 
  as.data.frame() 


# Calculate the percentage of missing values per column
percentage_missing_per_column <- colMeans(is.na(df3)) * 100

# Display the result
print(round(percentage_missing_per_column, 2))

# vis_miss(df3)
# gg_miss_var(df3)

colSums(is.na(df3)) # 0.49% ~ 2.56%

df3 <- df3 %>%
  mutate(
    momAlive = case_when(
      r11f001 == 1 & r12f001 == 1 ~ 1L,
      r11f001 == 1 & r12f001 == 0 ~ 0L,
      TRUE ~ NA_integer_
    ),
    dadAlive = case_when(
      r11f011 == 1 & r12f011 == 1 ~ 1L,
      r11f011 == 1 & r12f011 == 0 ~ 0L,
      TRUE ~ NA_integer_
    )
  )

table(df3$dadAlive)
table(df3$momAlive)

df3 <- df3 %>%
  mutate(
    pnewDie = case_when(
      r11f011 == 1 & r11f001 == 1 & (r12f011 == 0 | r12f001 == 0) ~ 1L,
      r11f011 == 1 & r11f001 == 1 & (r12f011 == 1 & r12f001 == 1) ~ 0L,
      TRUE ~ NA_integer_
    )
  )

table(df3$pnewDie, useNA = "ifany")
# g no individuals lost a parent between those waves 

df3 <- df3[complete.cases(df3), ]
table(df3$care_trans)    

table(df3$r11carecat)
table(df3$r12carecat)
table(df3$r11carecat, df3$r11f001)
# Binary marital status 
df3 <- df3 %>%
  mutate(
    r11marcat_b = case_when(
      r11marcat == "Married" ~ "Married",
      TRUE                   ~ "noMarried"
    ),
    r11marcat_b = factor(r11marcat_b, levels = c("noMarried", "Married"))
  )

# edu
df3$edu <- factor(df3$edu,
                  levels = c("Less than High School", 
                             "High School Graduate", 
                             "Some College", 
                             "College Degree"))

# Set r11smoke as a factor and assign value labels
df3$r11smoke <- factor(df3$r11smoke,
                       levels = c(0, 1, 2),
                       labels = c("Never smoker", "Past smoker", "Current smoker"))

# Standardize other 4 DNA methylation age
vars_to_standardize <- c("mpoa")
df3[paste0(vars_to_standardize, "_z")] <- scale(df3[vars_to_standardize])


saveRDS(df3, file = "df3.rds")
# save as RData after this step


##########################################################################################      
#           Parr IV: 
#           Descriptive table 1
##########################################################################################   
# Table 1 
# Before weighting
df_table1 <- df3 %>%
  select(
    horvath_dnamage, hannum_dnamage, levine_dnamage,
    horvathskin_dnamage, lin_dnamage, weidner_dnamage,
    vidalbralo_dnamage, grim_dnamage, zhang_dnamage,
    bocklandt_dnamage, garagnani_dnamage, mpoa,
    care_trans, r11age, female, racecat, r11marcat_b,
    edu, h11itot1, r11shlt, r11smoke, r11drinkn,
    r11mobila, r11chronic, r11acts_w, r11bmi
  )


categorical_vars <- c("care_trans", "female", "racecat", "r11marcat_b", "edu", "r11smoke", "r11drinkn")

df_table1 <- df_table1 %>%
  mutate(across(c(care_trans, female, racecat, r11marcat_b, edu, r11smoke, r11drinkn), as.factor))

tab1 <- CreateTableOne(
  vars = names(df_table1),
  factorVars = categorical_vars,
  data = df_table1
)
# Extract as data frame
tab1_df <- print(tab1, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
tab1_df <- as.data.frame(tab1_df)

# If the result is a matrix, convert to data frame
tab1_df <- as.data.frame(tab1_df)

# 2. Tidy Up and Reformat Table
# Example: Split tab1_df row names to get variable/level
tab1_df$Variable <- rownames(tab1_df)
tab1_df <- tab1_df %>%
  separate(Variable, into = c("Variable", "Level"), sep = " ", extra = "merge", fill = "right") %>%
  mutate(Level = ifelse(is.na(Level), "", Level))

# Clean up for publication style
tab1_tidy <- tab1_df %>%
  select(Variable, Level, Overall) %>%
  rename("Mean/Proportion" = Overall)

# Optionally add SD/Range/Count columns, if available
# For continuous variables, you might extract SD from the parenthesis
tab1_tidy <- tab1_tidy %>%
  mutate(
    SD = str_extract(`Mean/Proportion`, "(?<=\\().+?(?=\\))"),
    `Mean/Proportion` = str_remove(`Mean/Proportion`, " \\(.+\\)")
  )

kable(
  tab1_tidy,
  booktabs = TRUE,
  col.names = c("Variable", "Level", "Mean/Proportion", "SD")
) %>%
  kable_styling(full_width = FALSE, position = "left")

print(tab1, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE)
write.csv(print(tab1, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE), "Table1.csv")



###
### after weighting
df_table1 <- df_table1 %>%
  mutate(across(where(is.labelled), ~ as.numeric(.)))

# Or, if you prefer to convert to factors
df_table1 <- df_table1 %>%
  mutate(across(where(is.labelled), ~ as_factor(.)))


# 1. Define the Survey Design
d.vbs <- svydesign(
  ids     = ~secu,            # Primary Sampling Units (clusters)
  strata  = ~stratum,         # Stratification variable
  weights = ~vbsi16wgtra,     # Survey weights
  data    = df_table1,        # Subsetted dataset
  nest    = TRUE              # Ensures unique cluster IDs within strata
)



# 2. Summarize Continuous Variables with Weights
# List of continuous variables
continuous_vars <- setdiff(names(df_table1), categorical_vars)

# Initialize an empty list to store summaries
cont_summary_list <- list()

# Loop through each continuous variable
for (var in continuous_vars) {
  # Create a formula for the variable
  formula <- as.formula(paste0("~", var))
  
  # Calculate weighted mean
  mean_val <- svymean(formula, d.vbs, na.rm = TRUE)
  
  # Calculate weighted variance and derive standard deviation
  var_val <- svyvar(formula, d.vbs, na.rm = TRUE)
  sd_val <- sqrt(var_val)
  
  # Calculate min and max (unweighted)
  min_val <- min(df_table1[[var]], na.rm = TRUE)
  max_val <- max(df_table1[[var]], na.rm = TRUE)
  
  # Count of non-missing observations (unweighted)
  count_val <- sum(!is.na(df_table1[[var]]))
  
  # Store the results in a data frame
  cont_summary_list[[var]] <- data.frame(
    Variable = var,
    Level = "",
    Mean = round(coef(mean_val), 2),
    SD = round(coef(sd_val), 2),
    Min = round(min_val, 2),
    Max = round(max_val, 2),
    Count = count_val,
    stringsAsFactors = FALSE
  )
}

# Combine all summaries into one data frame
cont_summary <- do.call(rbind, cont_summary_list)




# 3. Summarize Categorical Variables (Unweighted Counts)
# Initialize list to store summaries
cat_summary_list <- list()

# Loop through each categorical variable
for (var in categorical_vars) {
  formula <- as.formula(paste0("~", var))
  
  # Calculate weighted proportions
  prop <- svymean(formula, design = d.vbs, na.rm = TRUE)
  
  # Extract levels and proportions
  levels <- names(prop)
  proportions <- as.numeric(prop)
  
  # Calculate unweighted counts
  counts <- table(df_table1[[var]])
  
  # Create data frame for current variable
  df <- data.frame(
    Variable = var,
    Level = gsub(paste0(var), "", levels),
    Mean = round(proportions, 4),
    SD = "",
    Min = "",
    Max = "",
    Count = as.numeric(counts),
    stringsAsFactors = FALSE
  )
  
  # Append to list
  cat_summary_list[[var]] <- df
}

# Combine all summaries into one data frame
cat_summary <- do.call(rbind, cat_summary_list)




# 4. Combine and Format the Table
# Combine summaries
table1 <- rbind(cont_summary, cat_summary)

# Convert all columns to character type for consistent formatting
table1 <- table1 %>%
  mutate(across(c(Mean, SD, Min, Max, Count), as.character))



kable(
  table1,
  booktabs = TRUE,
  col.names = c("Variable", "Level", "Weighted Mean/Proportion", "SD", "Min", "Max", "Unweighted Count"),
  caption = "Descriptive Statistics (Weighted)"
) %>%
  kable_styling(full_width = FALSE, position = "left")





min(df3$horvath_dnamage)
min(df3$r11age)


