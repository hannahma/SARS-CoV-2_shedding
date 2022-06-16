library(tidyverse)
library(patchwork)
library(lubridate)

load("./data.RData")
source("./code/00_functionS2.R")

dat_nest <- dat_shed %>% nest(data = everything()) %>% 
  mutate(agecat = "all ages")
dat_nest 


dat_nest_mod_ev <-  dat_nest %>% # run models
  mutate(
    model_elisa_vax  = map(data,  model_pred_ev),
    etr_ev           = map_df(data, model_ETR_ev), 
    pred_ev          = map2(model_elisa_vax, "elisa_vax", predict_shed)  # predict pct & iqr
    
  ) %>% select(-c(data, contains("model_")))
dat_nest_mod_ev

dat_nest_mod_ev %>% select(agecat, "pred_ev") %>% unnest("pred_ev")

# plot
shed_max <- 70
x_iqr    <- 45

p_serovax <- plot_add_ETRS_ev_shed(dat_nest_mod_ev, shed_max, x_iqr) +
  theme(legend.position = "bottom")  

p_serovax

save(p_serovax, file = "./results/ggplot_p_serovax.RData")
