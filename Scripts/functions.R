## Author: Victoria Kulesh
## Functions for model evaluation

MSDcol <- c("#1a1866", "#f2b93b", "#b73b58", "#a2d620", "#5839bb", "#9c4ec7", "#3a6eba", "#efdd3c", "#69686d", "red", "darkmagenta")


funSum_sim <- list(mean   = ~mean(.),
                   median = ~median(.),
                   min    = ~min(.),
                   max    = ~max(.),
                   sd     = ~sd(.), 
                   P025   = ~quantile(., 0.025, na.rm = T),
                   P05    = ~quantile(., 0.05, na.rm = T),
                   P10    = ~quantile(., 0.10, na.rm = T),
                   P25    = ~quantile(., 0.25, na.rm = T),
                   P75    = ~quantile(., 0.75, na.rm = T),
                   P90    = ~quantile(., 0.90, na.rm = T),
                   P95    = ~quantile(., 0.95, na.rm = T),
                   P975   = ~quantile(., 0.975, na.rm = T),
                   CV     = ~sd(., na.rm = T)/mean(., na.rm = T)*100)



mod_eval_func <- function(folder_path, 
                          mod_p, 
                          age_vec_i = NULL, par_def = NULL,
                          npop_i = 1, unc_fl_i = F,
                          obs_data = NULL, carry_out = NULL) {
  
  base_proj_0 <- folder_path
  base_par_fin <- read_csv(str_c(folder_path, "populationParameters.txt"), col_types = cols())
  base_par_pop <- base_par_fin %>% filter(str_detect(parameter, "_pop")) %>% select(parameter, value) %>% deframe()
  names(base_par_pop) <- names(base_par_pop) %>% str_replace("_pop", "")
  
  mod_p_rxode <- rxode2::rxode(str_c(folder_path, mod_p, ".txt"))
  mod_p_rxode$params[!mod_p_rxode$params %in% names(base_par_pop)]
  
  par_all <- base_par_pop
  
  if (!is.null(par_def)){
    par_all[names(par_def)] <- par_def
  }
  
  if (!unc_fl_i){
    res_base <- mod_diag_func(path = NULL, par_test = par_all, mod_path = mod_p_rxode, age_vec = age_vec_i)
  } else {
    res_base <- mod_diag_func(path = base_proj_0, mod_path = mod_p_rxode, age_vec = age_vec_i, npop = npop_i, unc_fl = unc_fl_i, aggr_fl = unc_fl_i)
  }
  
  res_base_mod <- res_base %>% 
    mutate(NAME = case_when(str_detect(VAR, "_RTE4") ~ "RTE",
                            str_detect(VAR, "_N4") ~ "Naive",
                            str_detect(VAR, "_A4") ~ "Activated",
                            str_detect(VAR, "_CM4") ~ "Central-memory",
                            str_detect(VAR, "_EM4") ~ "Effector-memory",
                            str_detect(VAR, "_EFF4") ~ "Effector",
                            str_detect(VAR, "T4_") ~ "Total CD4+",
                            str_detect(VAR, "_MEM4") ~ "Total memory CD4+"),
           ORGAN = case_when(str_detect(VAR, "cells_uL|bl") ~ "Blood",
                             str_detect(VAR, "slo_cfb|slo|lt_perc") ~ "Lymphoid tissue",
                             str_detect(VAR, "git_cfb|git") ~ "GI tract",
                             str_detect(VAR, "lung_cfb|lung") ~ "Lungs"),
           TYPE = ifelse(str_detect(VAR, "_cfb|_perc"), "cfb", "abs")) %>% 
    mutate(NAME = factor(NAME, levels = c("RTE", "Naive", "Activated", "Central-memory", "Effector-memory", "Effector", "Total CD4+", "Total memory CD4+")),
           ORGAN = factor(ORGAN, levels = c("Blood", "Lymphoid tissue", "GI tract", "Lungs")))
  
  
  if (!is.null(obs_data)){
    
    obs_data_mod <- obs_data %>% 
      mutate(NAME = case_when(str_detect(DVNAME, "_RTE4") ~ "RTE",
                              str_detect(DVNAME, "_N4") ~ "Naive",
                              str_detect(DVNAME, "_A4") ~ "Activated",
                              str_detect(DVNAME, "_CM4") ~ "Central-memory",
                              str_detect(DVNAME, "_EM4") ~ "Effector-memory",
                              str_detect(DVNAME, "_EFF4") ~ "Effector",
                              str_detect(DVNAME, "T4_") ~ "Total CD4+",
                              str_detect(DVNAME, "_MEM4") ~ "Total memory CD4+"),
             ORGAN = case_when(str_detect(DVNAME, "cells_uL|bl") ~ "Blood",
                               str_detect(DVNAME, "slo_cfb|slo|lt_perc") ~ "Lymphoid tissue",
                               str_detect(DVNAME, "git_cfb|git") ~ "GI tract",
                               str_detect(DVNAME, "lung_cfb|lung") ~ "Lungs"),
             TYPE = ifelse(str_detect(DVNAME, "_cfb|_perc"), "cfb", "abs")) %>% 
      mutate(NAME = factor(NAME, levels = c("RTE", "Naive", "Activated", "Central-memory", "Effector-memory", "Effector", "Total CD4+", "Total memory CD4+")),
             ORGAN = factor(ORGAN, levels = c("Blood", "Lymphoid tissue", "GI tract", "Lungs")))
    
  } else {obs_data_mod <- NULL}
  
  

  # explore age dynamics
  base_manual_bl_cells <- age_dyn_func(res_data_i = res_base, obs_data = obs_data, aggr_fl = unc_fl_i,
                                       var_vec = c("T_RTE4_bl_cells_uL", "T_N4_bl_cells_uL",
                                                   "T_A4_bl_cells_uL", "T_CM4_bl_cells_uL",
                                                   "T_EM4_bl_cells_uL", "T_EFF4_bl_cells_uL", 
                                                   "T4_bl_cells_uL", "T4_MEM_bl_cells_uL"))
  
  base_manual_bl <- age_dyn_func(res_data_i = res_base, obs_data = obs_data, aggr_fl = unc_fl_i,
                                 var_vec = c("T_RTE4_bl", "T_N4_bl", "T_A4_bl", "T_CM4_bl", "T_EM4_bl", "T_EFF4_bl"),
                                 y_axis_name = "Number of cells")
  
  base_manual_slo <- age_dyn_func(res_data_i = res_base, obs_data = obs_data, aggr_fl = unc_fl_i,
                                  var_vec = c("T_RTE4_lt", "T_N4_lt", "T_A4_lt", "T_CM4_lt", "T_EM4_lt", "T_EFF4_lt"),
                                  y_axis_name = "Number of cells")
  
  
  base_manual_slo_cfb <- age_dyn_func(res_data_i = res_base, obs_data = obs_data, aggr_fl = unc_fl_i,
                                      var_vec = c("T_RTE4_lt_perc", "T_N4_lt_perc", "T_A4_lt_perc", "T_CM4_lt_perc", "T_EM4_lt_perc", "T_EFF4_lt_perc"),
                                      y_axis_name = "Cells, %")
  
  base_manual_git_cfb <- age_dyn_func(res_data_i = res_base, obs_data = obs_data, aggr_fl = unc_fl_i,
                                      var_vec = c("T_N4_git_perc", "T_CM4_git_perc", "T_EM4_git_perc", "T_EFF4_git_perc"),
                                      y_axis_name = "Cells, %")
  
  base_manual_lung_cfb <- age_dyn_func(res_data_i = res_base, obs_data = obs_data, aggr_fl = unc_fl_i,
                                       var_vec = c("T_N4_lung_perc", "T_CM4_lung_perc", "T_EM4_lung_perc", "T_EFF4_lung_perc"),
                                       y_axis_name = "Cells, %")
  
  
  
  if (!is.null(carry_out)){
    
    base_par_ruv <- base_par_fin %>% filter(!str_detect(parameter, "_pop|omega_|beta_")) %>% select(parameter, value) %>% deframe()
    ruv_data <- data.frame(NAME = names(base_par_ruv),
                           SIGMA = base_par_ruv) %>% 
      mutate(DVID = str_extract_all(NAME, "\\d+") %>% as.numeric())
    
    
  } else {
    
    ruv_data <- NULL
    
  }
  
  return(list(par_all = par_all,
              res_base = res_base_mod,
              obs_data = obs_data_mod,
              ruv_data = ruv_data,
              base_manual_bl = base_manual_bl,
              base_manual_bl_cells = base_manual_bl_cells,
              base_manual_lt = base_manual_slo,
              base_manual_lt_perc = base_manual_slo_cfb,
              base_manual_git_perc = base_manual_git_cfb,
              base_manual_lung_perc = base_manual_lung_cfb))
  
}




fancy_scientific <- function(l) {
  if (any(l>10^4, na.rm = T)){
    
    # turn in to character string in scientific notation
    l <- format(l, scientific = TRUE)
    # quote the part before the exponent to keep all the digits
    l <- gsub("^(.*)e", "'\\1'e", l)
    # turn the 'e+' into plotmath format
    
    l <- gsub("e", "%*%10^", l)
    
    l <- gsub('[:+:]', "", l)
    
    l <- ifelse(grepl("'0'", l) | grepl("'0.0'", l), "0", l)
    
    # return this as an expression
  } 
  
  parse(text=l)
}


age_dyn_func <- function(res_data_i = res_age_all, obs_data = NULL, var_column = "VAR",organ = "Blood", type = "abs",
                         var_vec = NULL,
                         aggr_fl = F, color_col = "LABEL", y_axis_name = "Concentration, cells/uL", 
                         cells_fl = F, facet_fl = T, insert_fl = F, ylim = NULL, age_h = 5, ylim_ins = NULL,
                         x_ins = 0.45, y_ins = 0.57, width_ins = 0.5, height_ins = 0.4,
                         title_nm = NULL){

  
  
  if (var_column == "VAR"){
    
    if (!is.null(obs_data)) {
      obs_data <- obs_data %>% mutate(VAR_pl = DVNAME)
      }
    
    res_data <- res_data_i %>% mutate(VAR_pl = VAR)
    
  } else {
    
    if (!is.null(obs_data)) {obs_data <- obs_data %>% rename(VAR_pl = var_column) %>% filter(ORGAN == organ, TYPE == type)}
    res_data <- res_data_i %>% rename(VAR_pl = var_column) %>% filter(ORGAN == organ, TYPE == type)
    
  }

  res_data <- res_data %>% rename(COLOR_pl = color_col)
  
  if (!aggr_fl){
    
    if (facet_fl){
    
    age_plot <- ggplot(data = res_data %>% filter(VAR_pl %in% var_vec),
                       aes(x = age, y = VAL, color = COLOR_pl, group = COLOR_pl))+
      geom_line(size = 1, color =  MSDcol[1])+
      scale_color_manual(values = MSDcol)+
      facet_wrap(~VAR_pl, scales = "free")
    
    
    if (!is.null(obs_data)){
      
      age_plot <- age_plot + 
        geom_point(data = obs_data %>% filter(VAR_pl %in% var_vec),
                                        aes(x = age, y = DV, group = VAR_pl), color = MSDcol[3])+
        geom_errorbar(data = obs_data %>% filter(VAR_pl %in% var_vec),
                      aes(x = age,y = DV, ymin = DVLL, ymax = DVUL, group = VAR_pl), width = 1, color =  MSDcol[3], alpha = 0.5)
      
    }
    
    } else {
      
      age_plot <- ggplot()+
        geom_line(data = res_data %>% filter(VAR_pl %in% var_vec),
                  aes(x = age, y = VAL, color = COLOR_pl, group = COLOR_pl), size = 1)+
        scale_color_manual(values = MSDcol)
      
      if (!is.null(obs_data)){
        
      
        age_plot <- age_plot + 
          geom_point(data = obs_data %>% filter(VAR_pl %in% var_vec),
                     aes(x = age, y = DV, group = VAR_pl, color = VAR_pl))+
          geom_errorbar(data = obs_data %>% filter(VAR_pl %in% var_vec),
                        aes(x = age,y = DV, ymin = DVLL, ymax = DVUL, group = VAR_pl, color = VAR_pl), width = 1, alpha = 0.5)
      }
      
      }
    
  } else {
    
    if (facet_fl){
      
      #browser()
      
      age_plot <- ggplot()
      
      if (!is.null(obs_data)){
        
        age_plot <- age_plot +
          geom_point(data = obs_data %>% filter(VAR_pl %in% var_vec),
                     aes(x = age, y = DV, group = VAR_pl), color = MSDcol[3])+
          geom_errorbar(data = obs_data %>% filter(VAR_pl %in% var_vec),
                        aes(x = age,y = DV, ymin = DVLL, ymax = DVUL, group = VAR_pl), width = 1, color =  MSDcol[3], alpha = 0.5)
      }
      
      age_plot <- age_plot +
        geom_line(data = res_data %>% filter(VAR_pl %in% var_vec),
                  aes(x = age, y = mean, color = COLOR_pl,  group = COLOR_pl), size = 1, color =  MSDcol[1])+
        geom_ribbon(data = res_data %>% filter(VAR_pl %in% var_vec),
                    aes(x = age, y = mean, ymin = P025, ymax = P975, fill = COLOR_pl,  color = COLOR_pl,  group = COLOR_pl), alpha = 0.2, color =  MSDcol[1], fill =  MSDcol[1]) +
        scale_color_manual(values = MSDcol)+
        facet_wrap(~VAR_pl, scales = "free")
      
      
    } else {
      
      
      age_plot <- ggplot()
      
      
      if (!is.null(obs_data)){
        
        age_plot <- age_plot+
          geom_point(data = obs_data %>% filter(VAR_pl %in% var_vec),
                     aes(x = age, y = DV, group = VAR_pl), color = MSDcol[3])+
          geom_errorbar(data = obs_data %>% filter(VAR_pl %in% var_vec),
                        aes(x = age,y = DV, ymin = DVLL, ymax = DVUL, group = VAR_pl), width = 1, color =  MSDcol[3], alpha = 0.5)
        
        
      }
      
      
      age_plot <- age_plot +
        geom_line(data = res_data %>% filter(VAR_pl %in% var_vec),
                  aes(x = age, y = mean, color = COLOR_pl,  group = COLOR_pl), size = 1, color =  MSDcol[1])+
        geom_ribbon(data = res_data %>% filter(VAR_pl %in% var_vec),
                    aes(x = age, y = mean, ymin = P025, ymax = P975, fill = COLOR_pl,  color = COLOR_pl,  group = COLOR_pl), alpha = 0.2, color =  MSDcol[1], fill =  MSDcol[1]) +
        scale_color_manual(values = MSDcol)
      
    }
    
  }
  
  age_plot <- age_plot+
    scale_y_continuous(name = y_axis_name, labels = fancy_scientific, breaks = pretty_breaks(), limits = ylim)+
    scale_x_continuous(name = "Age, y.o.", breaks = seq(0, 120, 20), limits = c(-1, 120))+
    ggtitle(title_nm)
  
  if (length(unique(res_data$COLOR_pl)) == 1) {
    
    age_plot <- age_plot + theme(legend.position = "none")
    
  }
  
  if (insert_fl){
     if (facet_fl) {stop("cannot do insert plot for facets")} else {
       
       main.plot <- age_plot + 
         theme(legend.position = "none")+
         geom_rect(data = data.frame(ID = 1, AGE = 1), aes(ymin = 0, ymax = ylim[2], xmin = 0, xmax = age_h), fill = MSDcol[6], alpha = 0.05)
      
       if (is.null(ylim_ins)){
         ylim_ins <- ylim
         
       }
       inset.plot <- age_plot + 
         coord_cartesian(xlim = c(-0.5, age_h))+
         scale_x_continuous(name = "", breaks = seq(0, age_h, 1)) + 
         scale_y_continuous(name = "",breaks = pretty_breaks(), labels = fancy_scientific, limits = ylim_ins) + 
         geom_rect(data = data.frame(ID = 1, AGE = 1), aes(ymin = -Inf, ymax = Inf, xmin = -Inf, xmax = Inf), fill = MSDcol[6], alpha = 0.05)+
         ggtitle(NULL)+
         theme(legend.position = "none",
               #panel.background = element_rect(fill = "transparent"),
               plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
               #panel.grid.major = element_blank(), #remove major gridlines
               panel.grid.minor = element_blank()
         ) #remove minor gridlines)
       
       age_plot <- ggdraw() +
         draw_plot(main.plot) +
         draw_plot(inset.plot, x = x_ins, y = y_ins, width = width_ins, height = height_ins)
       
        
     }
    
    
    
  }
  
  
  return(age_plot)
}


sampl_unc <- function(path_i, nsim_i = 1, mat = F, type = "SA"){
  poppar <- read_csv(str_c(path_i, "/populationParameters.txt"), col_types = cols())
  if (type == "SA") {
    se <- poppar$se_sa
  } else {
    se <- poppar$se_lin
  }
  
  if (any(is.na(se))){
    # fixed parameters
    fixed_fl <- T
    mu_f <- poppar$value[is.na(se)]
    names(mu_f) <- poppar$parameter[is.na(se)]
  } else {fixed_fl <- F} 
  
  mu <- poppar$value[!is.na(se)]
  
  pop_names <- poppar$parameter[!is.na(se)]
  se <- se[!is.na(se)]
  ccor <- read_csv(str_c(path_i, "correlationEstimates", type, ".txt"),
                   col_names = F, col_types = cols()) %>% select(where(is.numeric)) %>% as.matrix()
  sigmas <- diag(se) %*% ccor %*% diag(se)
  colnames(sigmas) <- pop_names; rownames(sigmas) <- pop_names
  if(mat){ return(sigmas) }
  popparams <- MASS::mvrnorm(n = nsim_i, mu = mu, Sigma = sigmas) %>% as.data.frame() %>% as_tibble()
  colnames(popparams) <- pop_names
  
  if (fixed_fl) {popparams <- popparams %>% merge(as_tibble(enframe(mu_f) %>% spread(name, value)))}
  
  return(popparams)
}

mod_diag_func <- function(path, par_test = NULL, 
                          cov_tbl = NULL, mod_path, 
                          type_alg = "Lin", npop = 1, unc_fl = F, age_step = 1, 
                          age_vec = NULL, aggr_fl = F,
                          var_vec = c("T_RTE4_bl_cells_uL", "T_N4_bl_cells_uL",
                                      "T_A4_bl_cells_uL", "T_CM4_bl_cells_uL",
                                      "T_EM4_bl_cells_uL", "T_EFF4_bl_cells_uL", 
                                      "T_RTE4_bl", "T_N4_bl", "T_A4_bl", "T_CM4_bl", "T_EM4_bl", "T_EFF4_bl",
                                      "T_RTE4_lt", "T_N4_lt", "T_A4_lt", "T_CM4_lt", "T_EM4_lt", "T_EFF4_lt",
                                      "T_RTE4_lt_perc", "T_N4_lt_perc", "T_A4_lt_perc", "T_CM4_lt_perc", "T_EM4_lt_perc", "T_EFF4_lt_perc",
                                      "T_N4_git_perc", "T_CM4_git_perc", "T_EM4_git_perc", "T_EFF4_git_perc",
                                      "T_N4_lung_perc", "T_CM4_lung_perc", "T_EM4_lung_perc", "T_EFF4_lung_perc",
                                      "T4_bl_cells_uL", "T4_MEM_bl_cells_uL"
                                      )){
  
  if (!is.null(path)) {
    
    par_fin <- read_csv(str_c(path, "/populationParameters.txt"), col_types = cols())
    
    ### Population parameters
    par_pop <- par_fin %>% filter(str_detect(parameter, "_pop")) %>% select(parameter, value) %>% deframe()
    names(par_pop) <- names(par_pop) %>% str_replace("_pop", "")
    
  } else {
    
    par_pop <- par_test
  }
  
  if (!is.null(age_vec)){
    age_range <- age_vec
  } else {
    age_range <- seq(0, 120, age_step)
  }
  
  stimes <- seq(0, 10^6, 10^4)
  
  if (unc_fl & !is.null(path)) {
    
    
    m_theta <- sampl_unc(path_i = path, nsim_i = npop, mat = T, type = type_alg)
    m_theta_pop <- m_theta[str_detect(rownames(m_theta), "_pop"), str_detect(colnames(m_theta), "_pop")]
    pop_par_set <- sampl_unc(path_i = path, nsim_i = npop, mat = F, type = type_alg)
    colnames(pop_par_set) <- str_replace(colnames(pop_par_set), "_pop", "")
    
    res_sim_dt_unc <- data.frame()
    for (i in 1:nrow(pop_par_set)){
      
      print(i)
      
      
      params_sim_i <- pop_par_set[i, ] %>% merge(tibble(age = age_range))
      
      ev_tbl_i <- params_sim_i %>% 
        mutate(id = 1:nrow(.), time = 0) %>% select(id, time) %>% as.data.frame() %>% 
        rxode2::et() %>% rxode2::add.sampling(stimes)
      
      res_sim_i <- rxode2::rxSolve(object = rxode2::rxode2(mod_path),
                                   params = params_sim_i,
                                   events = ev_tbl_i)
      
      res_sim_i_dt <- res_sim_i %>% 
        filter(time == 10^6) %>% 
        gather("VAR", "VAL", -"id") %>% 
        left_join(params_sim_i %>% mutate(id = 1:nrow(.)), "id") %>% 
        mutate(LABEL = "Monolix fit",
               sim.id = i)
      
      res_sim_dt_unc <- rbind(res_sim_dt_unc, res_sim_i_dt)
      
    }
    
    if (aggr_fl){
      
      res_sim_dt_unc_aggr <- res_sim_dt_unc %>% filter(VAR %in% var_vec) %>% 
        group_by(age, VAR, LABEL) %>% summarise_at(vars(VAL), funSum_sim) %>% ungroup()
      
      result_dt <- res_sim_dt_unc_aggr
    } else {
      
      result_dt <- res_sim_dt_unc
      
    }
    
    
  } else {
    
    if (is.null(cov_tbl)){
      
      params_sim <- data.frame(par_pop) %>% t() %>% merge(tibble(age = age_range))
      
    } else {
      
      params_sim <- data.frame(par_pop) %>% t() %>% merge(tibble(age = age_range)) %>% merge(cov_tbl)
      
    }
    
    ev_tbl <- params_sim %>% 
      mutate(id = 1:nrow(.), time = 0) %>% select(id, time) %>% as.data.frame() %>% 
      rxode2::et() %>% rxode2::add.sampling(stimes)
    
    mod <- rxode2::rxode2(mod_path)
    vars <- c(mod$lhs, mod$state)
    vars_new <- str_c(vars, "_pl")
      
    res_sim <- rxode2::rxSolve(object = rxode2::rxode2(mod_path),
                               params = params_sim,
                               events = ev_tbl)
    #browser()
    res_sim_dt <- res_sim %>% 
      filter(time == 10^6) %>%
      mutate(across(all_of(vars), 
                    ~ ., 
                    .names = "{.col}_pl")) %>% 
      gather("VAR", "VAL", -c("id", vars_new)) %>% 
      left_join(params_sim %>% mutate(id = 1:nrow(.)), "id") %>% 
      mutate(LABEL = "Monolix fit")
    
    result_dt <- res_sim_dt
    
  }
  
  return(result_dt)
  
}


