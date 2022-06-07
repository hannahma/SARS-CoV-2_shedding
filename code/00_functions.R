
# shedding model functions

library(tidyverse)
library(survival)
library(SurvRegCensCov)
library(patchwork)

theme_hm <- function(plot){
  theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          
          strip.text.x = element_text(size = 12),
          
          #strip.background = element_rect(fill = "grey90")) 
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
               nrow = 1) +
    labs(
      x = "proportion RT-PCR positive",
      y = "shedding duration (days)"
    ) +
    scale_x_continuous(breaks = seq(0, 1, .5))  # (flipped coords)
  
  
}

model_pred_elisa <- function(data){
  survival::survreg(
    shed_interval ~ factor(prior_elisa),
    data = {{data}},
    dist = "weibull"
  )
}
model_pred_titer <- function(data){   
  survival::survreg(
    shed_interval ~ prior_titer,
    data = {{data}},
    dist = "weibull"
  )
}


model_ETR_elisa <- function(data){
  SurvRegCensCov::WeibullReg(
    formula = shed_interval ~ factor(prior_elisa), 
    data    = {{data}}
  )$ETR %>% 
    as_tibble() %>% 
    mutate(across(c(ETR, LB, UB), round, digits = 2),
           etr_ci = paste0(ETR, " (", LB, ", ", UB, ")"))
}
model_ETR_titer <- function(data){
  SurvRegCensCov::WeibullReg(
    formula = shed_interval ~ prior_titer, 
    data    = {{data}}
  )$ETR %>% 
    as_tibble() %>% 
    mutate(across(c(ETR, LB, UB), round, digits = 2),
           etr_ci = paste0(ETR, " (", LB, ", ", UB, ")"))
}



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
  
  if (exposure=="prior_elisa")
    
    pred_she_nested({{shed_model}}, 1) %>% 
    full_join(
      pred_she_nested({{shed_model}}, 0)
    ) %>% 
    rename(prior_elisa = exp_value)
  
  else if (exposure=="prior_titer")
    
    pred_she_nested({{shed_model}}, 0) %>% 
    full_join(pred_she_nested({{shed_model}}, 2)) %>% 
    full_join(pred_she_nested({{shed_model}}, 3)) %>% 
    full_join(pred_she_nested({{shed_model}}, 4)) %>% 
    full_join(pred_she_nested({{shed_model}}, 5)) %>% 
    rename(prior_titer = exp_value)
  
}
predict_shed_iqr <- function(shed_model, exposure){
  if (exposure=="prior_elisa") 
    stats::predict(
      {{shed_model}}, 
      newdata = tibble(prior_elisa = c(0,1)),
      type    = 'quantile',
      p       = c(0.5, 0.25, 0.75),
      se      = F
    ) %>%     
    as_tibble() %>% 
    mutate(
      prior_elisa =  c(0,1)
    ) %>% 
    rename_with(~paste0("pct_", c(0.5, 0.25, 0.75)), 
                starts_with("V")) %>% 
    mutate(
      across(contains("pct"), round, digits = 1),
      across(contains("pct"), format, nsmall = 1),
      
      prior_elisa_iqr = paste0(
        prior_elisa, ": ",
        pct_0.5, " (", pct_0.25, ", ", pct_0.75, ")"
      )
    ) #%>% select(-contains("pct"))
  
  else if (exposure=="prior_titer")
    
    stats::predict(
      {{shed_model}}, 
      newdata = tibble(prior_titer = c(0,2,3,4,5)),
      type    = 'quantile',
      p       = c(0.5, 0.25, 0.75),
      se      = F
    ) %>%     
    as_tibble() %>% 
    mutate(
      prior_titer =  c(0,2,3,4,5)
    ) %>% 
    rename_with(~paste0("pct_", c(0.5, 0.25, 0.75)), 
                starts_with("V")) %>% 
    mutate(
      across(contains("pct"), round, digits = 1),
      across(contains("pct"), format, nsmall = 1),
      
      prior_titer_iqr = paste0(
        prior_titer, ": ",
        pct_0.5, " (", pct_0.25, ", ", pct_0.75, ")"
      )
    ) #%>% select(-contains("pct"))
  
  
  
}
predict_shed <- function(shed_model, exposure){
  predict_shed_curve({{shed_model}}, {{exposure}}) %>% 
    full_join(predict_shed_iqr({{shed_model}}, {{exposure}}))
}




# * ELISA
plot_elisa_shed <- function(data, shed_max){
  data %>% 
    unnest(pred_elisa) %>%               # ELISA-SPECIFIC
    mutate(prior_elisa = factor(prior_elisa, 
                                levels = c(0,1),
                                labels = c("negative","positive"))) %>% 
    
    plot_shed(prior_elisa, {{shed_max}}) + # GENERAL SHED PLOT FUNCTION
    
    # add etrs:
    geom_text(data = 
                data %>% 
                unnest(etr_elisa) %>%    # ELISA-SPECIFIC
                mutate( # define where to put the ETRs on the plot
                  pct = 0.2,
                  shed_time = 40,
                ) , 
              #size = 5,
              aes(label = paste0("ETR: \n", etr_ci)) ) +
    
    # general plot 
    #ylim(0,70) +
    theme_hm() +
    labs(color = "prior \nserostatus", 
         fill  = "prior \nserostatus",
         linetype  = "prior \nserostatus")  +
    scale_linetype_manual(values = c("dotted", "solid"))
}

# add IQRs below elisa shedding plots
plot_add_ETRS_elisa_shed <- function(data, shed_max, iqr_x_position){
  
  how_far_below <- 1.2#1.4
  
  elisa_mean_iqr <- data %>% 
    unnest(pred_elisa) %>% 
    mutate(prior_elisa = factor(prior_elisa, 
                                levels = c(0,1),
                                labels = c("negative","positive"))) %>% 
    distinct(agecat, prior_elisa, pct_0.5, pct_0.25, pct_0.75, prior_elisa_iqr) %>% 
    mutate(across(c(pct_0.5, pct_0.25, pct_0.75), as.numeric)) %>% 
    mutate( prior_elisa_iqr = str_sub(prior_elisa_iqr, 4, -1))
  
  plot_elisa_shed(data, {{shed_max}}) +
    
    geom_pointrange(
      data = elisa_mean_iqr %>% 
        mutate(pct = how_far_below, # place below y-axis (y = 1-pct)
               shed_time = pct_0.5),
      position = position_dodge(0.2), show.legend = F,
      aes(ymin = pct_0.25,
          ymax = pct_0.75,
          color = prior_elisa,
          linetype = prior_elisa)
    ) +
    geom_vline(xintercept = 0, color = "grey40") +
    geom_text(data = elisa_mean_iqr %>% 
                mutate(pct = how_far_below, # place below y-axis (y = 1-pct)
                       shed_time = pct_0.5),
              position = position_dodge(0.4),show.legend = F,
              aes(label = prior_elisa_iqr, 
                  color = prior_elisa,
                  y = {{iqr_x_position}}), 
              #size = 5,
              hjust = 0) +
    annotate("text", 
             y = 25, 
             x = -0.05, 
             size = 3,
             label = "mean (IQR)", hjust = 0.5) 
  
}



# * TITERS
plot_titer_shed <- function(data, shed_max){
  data %>% 
    unnest(pred_titer) %>%               # TITERS-SPECIFIC
    #mutate(prior_titer_bt = 4^prior_titer*5) %>% 
    mutate(prior_titer = 
             factor(prior_titer,  
                    levels = c(0,2,3,4,5),
                    labels = c("negative","80","320","1280","5120"))) %>% 
    
    
    
    plot_shed(prior_titer, {{shed_max}}) + # GENERAL SHED PLOT FUNCTION
    
    # add etr:
    geom_text(data = 
                data %>% 
                unnest(etr_titer) %>%    # TITERS-SPECIFIC
                mutate( # define where to put the ETRs on the plot
                  pct = 0.2,
                  shed_time = 40,
                ) , 
              #size = 5,
              aes(label = paste0("ETR 4-fold rise: \n", etr_ci)) ) +
    
    # general plot 
    #ylim(0,60) +
    theme_hm() +
    labs(color = "prior \ntiter", 
         fill  = "prior \ntiter",
         linetype  = "prior \ntiter") +
    scale_linetype_manual(values = c("dotted",
                                     "dotdash",
                                     "longdash",
                                     "twodash",
                                     "solid"))
}

# add IQRs below TITERS shedding plots
plot_add_ETRS_titers_shed <- function(data, shed_max, iqr_x_position){
  
  how_far_below <- 1.4

  titer_mean_iqr <- data %>% 
    unnest(pred_titer) %>% 
    mutate(prior_titer =
             factor(prior_titer,
                    levels = c(0,2,3,4,5),
                    labels = c("negative","80","320","1280","5120"))) %>% 
    distinct(agecat, prior_titer, pct_0.5, pct_0.25, pct_0.75, prior_titer_iqr) %>% 
    mutate(across(c(pct_0.5, pct_0.25, pct_0.75), as.numeric)) %>% 
    mutate(prior_titer_iqr = str_sub(prior_titer_iqr, 4, -1)) 

  
  plot_titer_shed(data, {{shed_max}}) +
    
    
    geom_pointrange(
      data = titer_mean_iqr %>%
        mutate(pct = how_far_below, # place below y-axis (y = 1-pct)
               shed_time = pct_0.5),
      position = position_dodge(0.5), show.legend = F,
      size = 0.5,
      aes(ymin = pct_0.25,
          ymax = pct_0.75,
          color = prior_titer,
          linetype = prior_titer)
    ) +
    geom_vline(xintercept = 0, color = "grey40") +


    geom_text(data = titer_mean_iqr %>%
                mutate(pct = how_far_below), # place below y-axis (y = 1-pct)
              position = position_dodge(0.5), show.legend = F,
              aes(label = prior_titer_iqr,
                  color = prior_titer,
                  y = {{iqr_x_position}}),
              #size = 5,
              hjust = 0) +
    annotate("text",
             y = 25,
             x = -0.05,
             size = 3,
             label = "mean (IQR)", hjust = 0.5)
  
}


