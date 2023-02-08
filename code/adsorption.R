library(here)
library(tidyverse)
library(cowplot)

# read data
d.plaques <- read_csv(here("data/adsorption/adsorption_plaques.csv"))

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
  left_join(d.plot, .)



# Plot titer --------------------------------------------------------------

d.plot %>% 
  pivot_longer(cols = starts_with("t"), names_to = "t.sample", values_to = "pfu.ml") %>% 
  mutate(t.sample = parse_number(t.sample)) %>% 
  ggplot(aes(t.sample, pfu.ml, color = exp.date))+
  geom_point()+
  geom_line()+
  facet_wrap(~cells)+
  theme_classic()+
  panel_border(color = "black")

# Plot relative titer ------------------------------------------------------

d.plot %>% 
  pivot_longer(cols = starts_with("t"), names_to = "t.sample", values_to = "pfu.ml") %>% 
  mutate(t.sample = parse_number(t.sample)) %>% 
  
  # normalize to T0
  arrange(t.sample) %>% 
  group_by(exp.date, flask, cells) %>% 
  mutate(rel.titer =  pfu.ml/first(pfu.ml)) %>% 
  
  
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
  mutate(cells = paste0(cells, " (n=", n, ")")) %>% 
  mutate(cells = fct_reorder(cells, flask)) %>% 
  
  ggplot(aes(t.sample, m, color = cells))+
  geom_line(size = 1)+
  geom_pointrange(aes(ymax = m+v, ymin = m-v, shape= cells),
                  stroke = 1, fill = "white", size=1)+

  theme_classic(base_size = 15)+
  # panel_border(color = "black")+
  ylab('% unadsorbed phage')+ xlab('Time (min)')+
  scale_x_continuous(breaks=c(0,5,10))+
  ylim(0,NA)+
  labs(caption = "mean?SD")+
  scale_color_manual(values=c("#FFC000", "#87a6d9", "#87cfd9"))+
  scale_shape_manual(values=c(21, 22,24))+
  theme(axis.text = element_text(face = "bold"))

ggsave(here("plots/SPO1_adsorption_delta6.png"), p,
       width =5, height =3)
