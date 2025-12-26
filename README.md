A R-code to run a Multiscale physiologically-based model of age-dependent CD4+ T-lymphocyte homeostasis.

R version: 4.2.3

R packages:
tidyverse (1.3.0)
scales (1.3.0)
rxode2 (2.1.2)


Models:
- CD4-hmst-mod-1 - files for model 1 (CD4+ T-Lymphocyte Cellular Kinetics Model)
  - CD4-hmst-mod-1-rxode2.txt - rxode-based model code
  - populationParameters.txt - estimated population parameter values from Monolix project
  - correlationEstimatesLin.txt - Correlation Fisher infromation matrix from Monolix project (for simylations with uncertainty)

- CD4-hmst-mod-2 - files for model 2 (Age-dependent CD4+ T-Lymphocyte Homeostatic Model (with age effect incorporated))
  - CD4-hmst-mod-2-rxode2.txt - rxode-based model code
  - populationParameters.txt - estimated population parameter values from Monolix project
  - correlationEstimatesLin.txt - Correlation Fisher infromation matrix from Monolix project (for simylations with uncertainty)
 
- CD4-hmst-mod-3 - files for model 3 (Age-dependent CD4+ T-Lymphocyte Homeostatic Model (with age and cell count effects incorporated))
  - CD4-hmst-mod-3-rxode2.txt - rxode-based model code
  - populationParameters.txt - estimated population parameter values from Monolix project
  - correlationEstimatesLin.txt - Correlation Fisher infromation matrix from Monolix project (for simylations with uncertainty)
 
Scripts:
- s01_model_simulations.R - R script to evaluate age dynamics of specific CD4+ T-lymphocyte subpopulations in blood, lymphoid tissue, gastro-intestinal tract and lungs (simulations with typical values and with uncertainty in parameter values)
- functions.R - R script with visualization functions

