library(here)
library(tidyverse)
library(cowplot)
library(broom)
library(scales)

log10_negative <- function(x){
  if(x>0) return(log10(x))
  if(x==0) return(0)
  return(log10(-x))
}

# read data
d.plaques <- read_csv(here("data/adsorption/adsorption_plaques.csv"))
d.cfu <- read_csv(here("data/adsorption/adsorption_cfu.csv"))

# select best count (highest) 
d.plot <- 
  d.plaques %>% 
  filter(t.sample!=0) %>% 
  group_by(exp.date, flask, t.sample) %>% 
  slice_max(plaques) %>% 
  select(exp.date, flask, cells, t.sample, pfu.ml) %>% 
  mutate(t.sample = paste0("t", t.sample)) %>%
  filter(!is.na(pfu.ml)) %>% 
  pivot_wider(names_from = t.sample, values_from = pfu.ml)

#add T0
d.plot <-
  d.plaques %>% 
  filter(t.sample==0) %>% 
  filter(volPlated.ml==0.1) %>% 
  group_by(exp.date, t.sample) %>% 
  summarise(m=mean(pfu.ml),v =sd(pfu.ml)) %>% 
  select(exp.date, t0=m) %>% 
  left_join(d.plot, .) %>% 
  relocate("t0", .before = "t5")

# add cfu
d.plot <- 
  left_join(d.plot,d.cfu, by = c("exp.date", "flask", "cells"))

# veg and spores only (stationary as veg)
d.plot <- d.plot %>% 
  filter(cells !=  "Log-phase") %>% 
  mutate(cells = str_replace_all(cells,"Stationary-phase", "Vegetative"))

# Percent adsorption
d.plot <- d.plot%>% 
  mutate(perc.ads = t5/t0) %>% 
  # mutate(perc.ads = if_else(perc.ads>100, 100, perc.ads)) %>% 
  mutate(perc.ads = 1 - perc.ads)  

# one-sample t-test
t.result <- bind_rows(
  d.plot %>% 
    filter(cells == "Spores") %>% 
    pull(perc.ads) %>%
    t.test(., alternative = "g") %>% 
    broom::tidy(),
  
  d.plot %>% 
    filter(cells == "Vegetative") %>% 
    pull(perc.ads) %>%
    t.test(., alternative = "g") %>% 
    broom::tidy()) %>% 
  mutate(cells = c("Spores", "Vegetative"))

p <- d.plot %>% 
  group_by(cells) %>% 
  summarise(m = mean(perc.ads), v = sd(perc.ads)) %>% 

 ggplot(aes(cells))+
  
  # CI
  geom_linerange(data = t.result,
                 aes(y = estimate, ymin=conf.low, ymax=conf.high), 
                 color = "grey90", linewidth = 10)+
  
  geom_hline(yintercept = 0, color = "black")+

  geom_pointrange(aes(y = m, ymax = m+v, ymin = m-v), 
                  stroke = 1,  size=1, shape = 21, fill = "white")+
  geom_jitter(data = d.plot, aes(y= perc.ads),
              color = "grey40", width = 0.05, height = 0, size = 2)+
  scale_y_continuous(labels = scales::percent, 
                     breaks = c(0,.25,.50,.75,1),
                     limits = c(NA, 1))+
  # scale_x_discrete(limits = c( "Vegetative","Spores"))+
  
  theme_classic(base_size = 15)+
  panel_border(color = "black")+
  
  labs(y = "Phage adsorption",
       x='Host cell type')+
  theme(axis.text = element_text(face = "bold"))
p
ggsave(here("plots/SPO1_adsorptionPERC_delta6.png"), p,
       width =4, height =4)
