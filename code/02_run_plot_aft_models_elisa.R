

# load functions
source("./code/00_functions.R")
library(patchwork)
library(lubridate)

# load data----
load("../2022_covid_shedding/data/dat_shed_2022-04-07.RData")

# exclude vaxx'd >= 2 weeks prior to elisa ----
# nv = no vaccinated 
dat_shed_nv <- dat_shed %>% filter(vax==0)
dat_shed_nv


# group data by age----
nv_nest <- dat_shed_nv %>% 
  group_by(agecat) %>% 
  nest() %>% 
  full_join(dat_shed_nv %>% 
              mutate(agecat = "all ages") %>% 
              group_by(agecat) %>% 
              nest()) %>%
  arrange(agecat)
nv_nest 

# run models on nested data----
nv_nest_mod <-  nv_nest %>% 
  mutate(
    model_elisa = map(data,  model_pred_elisa),
    model_titer = map(data,  model_pred_titer),
    
    pred_elisa  = map2(model_elisa, "prior_elisa", predict_shed),
    pred_titer  = map2(model_titer, "prior_titer", predict_shed),
    
    
    etr_elisa   = map_df(data, model_ETR_elisa),
    etr_titer   = map_df(data, model_ETR_titer)
  ) %>% 
  select(-c(data, contains("model_")))
nv_nest_mod


# figure 2----
shed_max <- 70
x_iqr <- 45
p_shed_nv <- (plot_add_ETRS_elisa_shed(nv_nest_mod %>% 
                                      filter(agecat == "all ages"), shed_max, x_iqr) /
             plot_add_ETRS_titers_shed(nv_nest_mod %>% 
                                         filter(agecat == "all ages"), shed_max, x_iqr) ) +
  plot_annotation(tag_levels = c("A"),
                  #title = "excluding people with any vaccination >=14 days prior to elisa"
                  ) 
p_shed_nv

load("./results/ggplot_p_serovax.RData")
p_serovax <- p_serovax +
  guides(
    fill=guide_legend(nrow=3,byrow=TRUE),
    fill=guide_legend(nrow=3,byrow=TRUE),
    linetype=guide_legend(nrow=3,byrow=TRUE)
  )

library(cowplot)
col1 <- cowplot::plot_grid(p_shed_nv)
col2 <- cowplot::plot_grid(NULL, (p_serovax +
                                    plot_annotation(tag_levels = list(c("C")))), NULL,
                           rel_heights = c(0.5,2,0.5), ncol = 1)

p2_v2 <- cowplot::plot_grid(col1, col2, ncol = 2, rel_widths = c(2,1.5))
p2_v2

ggsave(p2_v2,
       filename = "./results/fig_2_shed.pdf",
       width = 10, height = 10)


# figure 3----
x_iqr_age <- 36
p_shed_age_nv <- (plot_add_ETRS_elisa_shed(nv_nest_mod %>% 
                                          filter(agecat != "all ages"), shed_max, x_iqr_age) /
                 plot_add_ETRS_titers_shed(nv_nest_mod %>% 
                                             filter(agecat != "all ages"), shed_max, x_iqr_age) )  +
  plot_annotation(tag_levels = c("A"),
                  #title = "excluding people with any vaccination >=14 days prior to infection"
                  ) 
p_shed_age_nv

ggsave(p_shed_age_nv,
       filename = "./results/fig_3_shed_age.pdf",
       width = 11, #14, #11,
       height = 8)




