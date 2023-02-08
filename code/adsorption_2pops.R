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



# adsorption rate ---------------------------------------------------------

k.rate <- function(cfu.ml, pfu.t0, pfu.t, t){
  k <- (2.3/(cfu.ml*t)) * log (pfu.t0 / pfu.t)
  return(k)
}

d.rate <- 
  d.plot %>% 
  pivot_longer(cols = c(t5, t10), names_to = "t.end", values_to = "pfu.ml") %>% 
  mutate(k = k.rate(CFU.ml, t0, pfu.ml, parse_number(t.end) )) %>% 
  mutate (cells = if_else(cells == "Spores", "spore", "vegetative")) 




t.result <- bind_rows(
  d.rate %>% 
  filter(t.end == "t5") %>% 
  filter(cells == "spore") %>% 
  pull(k) %>%
  t.test(., alternative = "g") %>% 
  broom::tidy(),

d.rate %>% 
  filter(t.end == "t5") %>% 
  filter(cells == "vegetative") %>% 
  pull(k) %>% 
  t.test(., alternative = "g") %>% 
  broom::tidy()
) %>% 
  mutate(cells = c("spore", "vegetative"))

d.sum <- d.rate %>% 
  filter(t.end == "t5") %>% 
  # average
  group_by(flask, cells) %>% 
  summarise(m=mean(k), v = sd(k), n= n()) %>% 
  mutate(cells = fct_reorder(cells, flask)) 

p <- d.rate %>%   
  filter(t.end == "t5") %>% 
  mutate(cells = fct_reorder(cells, flask)) %>% 
  ggplot(aes(cells, k))+
  
  # CI
  geom_linerange(data = t.result,
                  aes(y = estimate, ymin=conf.low, ymax=conf.high), 
                  color = "grey90", size = 10)+
  
  geom_hline(yintercept = 0, color = "black")+
  geom_pointrange(data = d.sum,
                  aes(y =m, ymax = m+v, ymin = m-v),
                  stroke = 1,  size=1, shape = 21, fill = "white")+
  geom_jitter(color = "grey40", width = 0.1, height = 0, size = 2)+
  theme_classic(base_size = 15)+
  panel_border(color = "black")+
  
  scale_y_continuous(label= function(x) {ifelse(x<=0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))} )+
  labs(y = bquote('phage adsorption rate ' (mL^-1~min^-1)),
       x='host cell type')+
  theme(axis.text = element_text(face = "bold"))

p  
ggsave(here("plots/SPO1_adsorptionRATE_delta6.png"), p,
       width =4, height =4)
  
  # ggplot(aes(cells,k, color = exp.date))+
  # geom_hline(yintercept = 0, color = "black")+
  # geom_point()+
  # theme_classic()+
  # # scale_y_log10() +
  # panel_border(color = "black")



# Plot titer --------------------------------------------------------------

d.plot %>% 
  pivot_longer(cols = starts_with("t"), names_to = "t.sample", values_to = "pfu.ml") %>% 
  mutate(t.sample = parse_number(t.sample)) %>% 
  filter(t.sample < 8) %>% 
  ggplot(aes(t.sample, pfu.ml, color = exp.date))+
  geom_point()+
  geom_line()+
  facet_wrap(~cells)+
  theme_classic()+
  scale_y_log10()+
  panel_border(color = "black")

# Plot relative titer ------------------------------------------------------

d.plot %>% 
  pivot_longer(cols = starts_with("t"), names_to = "t.sample", values_to = "pfu.ml") %>% 
  mutate(t.sample = parse_number(t.sample)) %>% 
  
  # normalize to T0
  arrange(t.sample) %>% 
  group_by(exp.date, flask, cells) %>% 
  mutate(rel.titer =  pfu.ml/first(pfu.ml)) %>% 
  filter(t.sample < 8) %>% 
  
  ggplot(aes(t.sample, rel.titer, color = exp.date))+
  geom_point()+
  geom_line()+
  facet_wrap(~cells)+
  theme_classic()+
  panel_border(color = "black")


# Plot relative titer average ---------------------------------------------

p <- d.plot %>% 
  pivot_longer(cols = starts_with("t"), names_to = "t.sample", values_to = "pfu.ml") %>% 
  mutate(t.sample = parse_number(t.sample)) %>% 
  
  # normalize to T0
  arrange(t.sample) %>% 
  group_by(exp.date, flask, cells) %>% 
  mutate(rel.titer =  100*pfu.ml/first(pfu.ml)) %>% 
  
  # average
  group_by(flask, cells, t.sample) %>% 
  summarise(m=mean(rel.titer), v = sd(rel.titer), n= n()) %>% 
  mutate (cells = if_else(cells == "Spores", "spore", "vegetative")) %>%
  mutate(cells = fct_reorder(cells, flask)) %>% 
  filter(t.sample < 8) %>% 
  
  ggplot(aes(t.sample, m, group = cells))+
  geom_line(size = 1)+
  geom_pointrange(aes(ymax = m+v, ymin = m-v, shape= cells, fill = cells),
                  stroke = 1,  size=1)+

  theme_classic(base_size = 15)+
  panel_border(color = "black")+
  ylab('% unadsorbed phage')+ xlab('Time (min)')+
  scale_x_continuous(breaks=c(0,5))+
  ylim(0,NA)+
  # labs(caption = "mean?SD")+
  # scale_color_manual(values = c("black","black"))+
  scale_fill_manual(values = c("white", "black"))+
  scale_shape_manual(values=c(21, 22))+
  theme(axis.text = element_text(face = "bold"),
        legend.position = c(0.3,0.2),
        legend.title = element_blank())


p

ggsave(here("plots/SPO1_adsorption_delta6.png"), p,
       width =4, height =4)
