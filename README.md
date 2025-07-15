# Code for A Biological Embodiment of Caregiving: A Longitudinal Study of Transitions and Gender Differences in Epigenetic Aging
This repository contains R scripts that reproduce analyses for the project *A Biological Embodiment of Caregiving*. The code is organized into four scripts that should be run in sequence:

1. **1_Data prepare.R** – prepares the analysis dataset and derives variables.
2. **2_Preliminary modeling.R** – fits initial GLM models exploring associations.
3. **3_Weighting analysis.R** – performs inverse probability weighting diagnostics and models.
4. **4_Forest plot.R** – creates the forest plots summarizing weighted results.

Each script relies on a variety of packages from the R ecosystem (e.g., `tidyverse`, `survey`, and `broom`). Ensure these packages are installed before running the code.

The original data are not included in this repository. Users should supply their own data in the expected format before executing the scripts.
