## Author: Victoria Kulesh
## Description: CD4+ T-lymphocyte homeostatic model simulation

rm(list=ls())
# path to R script with functions
source("./Scripts/functions.R")


library(tidyverse)
library(scales)
library(rxode2)

theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())


#### ---- Model 1: CD4+ T-Lymphocyte Cellular Kinetics Model ------ ####

# path to folder with results
mod1_folder_path <- "./Models/CD4-hmst-mod-1/"
# name of file with model code (rxode)
mod1_rxode <- "CD4-hmst-mod-1-rxode2"

# typical value simulation
mod1_results <- mod_eval_func(folder_path = mod1_folder_path,
                              mod_p = mod1_rxode)

# visualization of results
mod1_results$base_manual_bl_cells
mod1_results$base_manual_bl
mod1_results$base_manual_lt
mod1_results$base_manual_lt_perc
mod1_results$base_manual_git_perc
mod1_results$base_manual_lung_perc



# simulation with uncertainty
mod1_results_unc100 <- mod_eval_func(folder_path = mod1_folder_path,
                                     mod_p = mod1_rxode,
                                     unc_fl_i = T, npop_i = 100, age_vec_i = c(seq(0, 2, 0.1), seq(3, 120, 1)))

# visualization of results
mod1_results_unc100$base_manual_bl_cells
mod1_results_unc100$base_manual_bl
mod1_results_unc100$base_manual_lt
mod1_results_unc100$base_manual_lt_perc
mod1_results_unc100$base_manual_git_perc
mod1_results_unc100$base_manual_lung_perc


#### ---- Model 2: Age-dependent CD4+ T-Lymphocyte Homeostatic Model (with age effect incorporated)  ------ ####

# path to folder with results
mod2_folder_path <- "./Models/CD4-hmst-mod-2/"
# name of file with model code (rxode)
mod2_rxode <- "CD4-hmst-mod-2-rxode2"

# typical value simulation
mod2_results <- mod_eval_func(folder_path = mod2_folder_path,
                              mod_p = mod2_rxode)

# visualization of results
mod2_results$base_manual_bl_cells
mod2_results$base_manual_bl
mod2_results$base_manual_lt
mod2_results$base_manual_lt_perc
mod2_results$base_manual_git_perc
mod2_results$base_manual_lung_perc



# simulation with uncertainty
mod2_results_unc100 <- mod_eval_func(folder_path = mod2_folder_path,
                                     mod_p = mod2_rxode,
                                     unc_fl_i = T, npop_i = 100, age_vec_i = c(seq(0, 2, 0.1), seq(3, 120, 1)))

# visualization of results
mod2_results_unc100$base_manual_bl_cells
mod2_results_unc100$base_manual_bl
mod2_results_unc100$base_manual_lt
mod2_results_unc100$base_manual_lt_perc
mod2_results_unc100$base_manual_git_perc
mod2_results_unc100$base_manual_lung_perc




#### ---- Model 3: Age-dependent CD4+ T-Lymphocyte Homeostatic Model (with age and cell count effects incorporated)  ------ ####

# path to folder with results
mod3_folder_path <- "./Models/CD4-hmst-mod-3/"
# name of file with model code (rxode)
mod3_rxode <- "CD4-hmst-mod-3-rxode2"

# typical value simulation
mod3_results <- mod_eval_func(folder_path = mod3_folder_path,
                              mod_p = mod3_rxode)

# visualization of results
mod3_results$base_manual_bl_cells
mod3_results$base_manual_bl
mod3_results$base_manual_lt
mod3_results$base_manual_lt_perc
mod3_results$base_manual_git_perc
mod3_results$base_manual_lung_perc



# simulation with uncertainty
mod3_results_unc100 <- mod_eval_func(folder_path = mod3_folder_path,
                                     mod_p = mod3_rxode,
                                     unc_fl_i = T, npop_i = 100, age_vec_i = c(seq(0, 2, 0.1), seq(3, 120, 1)))

# visualization of results
mod3_results_unc100$base_manual_bl_cells
mod3_results_unc100$base_manual_bl
mod3_results_unc100$base_manual_lt
mod3_results_unc100$base_manual_lt_perc
mod3_results_unc100$base_manual_git_perc
mod3_results_unc100$base_manual_lung_perc

