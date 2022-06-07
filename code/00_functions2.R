
# shedding model functions

library(tidyverse)
library(survival)
library(SurvRegCensCov)
library(patchwork)

# exposures:
# 1) any vax
# 2) fully vax
# 3) prior elisa
# 4) prior infection (yes/no)
# 5) prior infection (1,2,3)
# 6) elisa_vax (1,2,3 = seroneg/novax, fully vax, seropos)


theme_hm  <- function(plot){
  theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          
          strip.text.x = element_text(size = 12),
          strip.background = element_blank()) 
}
plot_shed <- function(data, exposure, shed_max){
  
  data %>%
    # coord_cartesian isn't working with coord_flip,
    # so, manually cutting off the data:
    mutate(shed_time = if_else(shed_time > {{shed_max}}, 
                               NA_real_, shed_time)) %>% 
    ggplot(aes(
      x = 1 - pct,
      y = shed_time,
    )
    ) +
    coord_flip() +
    geom_line(aes(color = {{exposure}},
                  linetype = {{exposure}})) +
    geom_ribbon(alpha = .2, inherit.aes = F,
                aes(ymin = pmax(rep(0.0001,length(l)), l),
                    ymax = pmin(rep({{shed_max}},    length(u)), u),
                    x = 1 - pct,
                    fill = {{exposure}})
    ) +
    facet_wrap(~agecat, #scales = "free", 
               ncol = 1) +
    labs(
      x = "proportion RT-PCR positive",
      y = "shedding duration (days)"
    ) +
    scale_x_continuous(breaks = seq(0, 1, .5))  # (flipped coords)
  
  
}


# * model objects
# name: [model_pred_] [shortened var name]:
model_pred_ev <- function(data){   
  survival::survreg(
    shed_interval ~ factor(elisa_vax),
    data = {{data}},
    dist = "weibull"
  )
}


# * ETRs
# name: [model_ETR_] [shortened var name]:
model_ETR_ev    <- function(data){
  SurvRegCensCov::WeibullReg(
    formula = shed_interval ~ factor(elisa_vax), 
    data    = {{data}}
  )$ETR %>% 
    as_tibble() %>% 
    mutate(across(c(ETR, LB, UB), round, digits = 2),
           etr_ci = paste0(ETR, " (", LB, ", ", UB, ")")) %>% 
    filter(row_number()==1)
} # selecting only 1 ETR (fully vax vs seroneg & not fully vax)



predict_shed_curve <- function(shed_model, exposure){ #  nest function within function=
  
  pred_she_nested <- function(shed_model, exp_value){
    stats::predict(
      {{shed_model}}, 
      newdata = tibble({{exposure}}:= {{exp_value}}),
      type    = 'quantile',
      p       = seq(0.01, 0.99, by = 0.01),
      se      = T
    ) %>%
      as_tibble() %>%
      rename(shed_time = fit) %>%
      
      mutate(pct = seq(0.01, 0.99, by = 0.01),
             exp_value = {{exp_value}},
             
             l = shed_time - 2*se.fit,
             u = shed_time + 2*se.fit) %>%
      relocate(shed_time, l, u) 
  }
  
   if (exposure=="elisa_vax")
    pred_she_nested({{shed_model}}, 1) %>% 
    full_join(pred_she_nested({{shed_model}}, 2)) %>% 
    full_join(pred_she_nested({{shed_model}}, 3)) %>% 
    rename(elisa_vax = exp_value)
  
}
predict_shed_iqr   <- function(shed_model, exposure){
   if (exposure=="elisa_vax")
    
    stats::predict(
      {{shed_model}}, 
      newdata = tibble(elisa_vax = c(1,2,3)),
      type    = 'quantile',
      p       = c(0.5, 0.25, 0.75),
      se      = F
    ) %>%     
    as_tibble() %>% 
    mutate(
      elisa_vax =  c(1,2,3)
    ) %>% 
    rename_with(~paste0("pct_", c(0.5, 0.25, 0.75)), 
                c("V1","V2","V3")) %>% 
    mutate(
      across(contains("pct"), round, digits = 1),
      across(contains("pct"), format, nsmall = 1),
      
      elisa_vax_iqr = paste0(   #### update name for var
        elisa_vax, ": ",
        pct_0.5, " (", pct_0.25, ", ", pct_0.75, ")"
      )
    ) 

  
}
predict_shed       <- function(shed_model, exposure){
  predict_shed_curve({{shed_model}}, {{exposure}}) %>% 
    full_join(predict_shed_iqr({{shed_model}}, {{exposure}}))
}





#plots
plot_ev_shed          <- function(data, shed_max){
  data %>% 
    unnest(pred_ev) %>%               # exp-SPECIFIC
    mutate(elisa_vax = factor(elisa_vax, 
                              levels = c(1,2,3),
                              labels = c("seroneg & not fully vax",
                                         "fully vax >= 14 days before infection",
                                         "seropositive")
                              )) %>% 
    
    plot_shed(elisa_vax, {{shed_max}}) + # GENERAL SHED PLOT FUNCTION
    
    # add etrs:
    geom_text(data = 
                data %>% 
                unnest(etr_ev) %>%    # exp-SPECIFIC
                mutate( # define where to put the ETRs on the plot
                  pct = 0.2,
                  shed_time = 40,
                ) , 
              #size = 5,
              aes(label = paste0("ETR fully vax vs \nseroneg & not fully vax: \n", etr_ci)) ) + # update
    
    # general plot 
    #ylim(0,70) +
    theme_hm() +
    labs(color     = "prior exposure", 
         fill      = "prior exposure",
         linetype  = "prior exposure")  +
    scale_linetype_manual(values = c("dotted","longdash", "solid"))
}
plot_add_ETRS_ev_shed <- function(data, shed_max, iqr_x_position){
  
  how_far_below <- 1.2#1.4
  
  ev_iqr <- data %>% 
    unnest(pred_ev) %>%               # exp-SPECIFIC
    mutate(elisa_vax = factor(elisa_vax, 
                              levels = c(1,2,3),
                              labels = c("seroneg & not fully vax",
                                         "fully vax >= 14 days before infection",
                                         "seropositive")
                              )) %>% 
    distinct(agecat, elisa_vax, pct_0.5, pct_0.25, pct_0.75, elisa_vax_iqr) %>% 
    mutate(across(c(pct_0.5, pct_0.25, pct_0.75), as.numeric)) %>% 
    mutate( elisa_vax_iqr = str_sub(elisa_vax_iqr, 4, -1))
  
  plot_ev_shed(data, {{shed_max}}) + # exp-specific
    
    geom_pointrange(
      data = ev_iqr %>% 
        mutate(pct = how_far_below, # place below y-axis (y = 1-pct)
               shed_time = pct_0.5),
      position = position_dodge(0.2), show.legend = F,
      aes(ymin = pct_0.25,
          ymax = pct_0.75,
          color = elisa_vax,
          linetype = elisa_vax)
    ) +
    geom_vline(xintercept = 0, color = "grey40") +
    geom_text(data = ev_iqr %>% 
                mutate(pct = how_far_below, # place below y-axis (y = 1-pct)
                       shed_time = pct_0.5),
              position = position_dodge(0.2), show.legend = F,
              aes(label = elisa_vax_iqr, 
                  color = elisa_vax,
                  y = {{iqr_x_position}}), 
              #size = 5,
              hjust = 0) +
    annotate("text", 
             y = 25, 
             x = -0.05, 
             size = 3,
             label = "mean (IQR)", hjust = 0.5) 
  
}


