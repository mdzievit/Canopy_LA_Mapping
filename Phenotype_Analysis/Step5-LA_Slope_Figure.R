library(tidyverse)
library(viridis)
library(boot)
library(cowplot)

par_blues <- read_tsv("Phenotype_Analysis/BLUEs/Comb_Analysis_BLUEs_Pop_Parents.txt") %>% 
  mutate(lin_rsq = inv.logit(lin_rsq_t)) %>% 
  select(-lin_rsq_t) %>% 
  gather(Trait,Value,-Gen2,-Pedigree)

par_la_traits <- par_blues %>% 
  filter(Trait %in% c("ULA","TUELA","FUELA","MLA")) %>% 
  mutate(Coded_Trait = case_when(
    Trait == "ULA" ~ 5,
    Trait == "TUELA" ~ 3,
    Trait == "FUELA" ~ 1,
    Trait == "MLA" ~ -2
    )) %>% 
  left_join(par_blues %>%
              filter(Trait %in% c("lin_rsq","lin_a","lin_b1","Avg_LA")) %>% 
              spread(Trait,Value))

par_la_traits_format <- par_la_traits %>% 
  filter(Gen2 %in% c(1,3)) %>% 
  mutate(Pop = 4) %>% 
  bind_rows(par_la_traits %>% 
              filter(Gen2 %in% c(2,3)) %>% 
              mutate(Pop = 5)) %>% 
  mutate(Y_Start = ((-2 * lin_b1) + lin_a),
         Y_End = ((5 * lin_b1) + lin_a))

dh_blups <- read_tsv("Phenotype_Analysis/BLUPs/Comb_Analysis_BLUPs_Pop.txt")%>%
  mutate(lin_rsq = inv.logit(lin_rsq_t)) %>% 
  select(-lin_rsq_t) %>% 
  gather(Trait,Value,-Gen2,-Pedigree,-Pop)


la_traits <- dh_blups %>% 
  filter(Trait %in% c("ULA","TUELA","FUELA","MLA")) %>% 
  mutate(Coded_Trait = case_when(
    Trait == "ULA" ~ 5,
    Trait == "TUELA" ~ 3,
    Trait == "FUELA" ~ 1,
    Trait == "MLA" ~ -2
    )) %>% 
  left_join(dh_blups %>%
              filter(Trait %in% c("lin_rsq","lin_a","lin_b1","Avg_LA")) %>% 
              spread(Trait,Value)) %>% 
  mutate(Y_Start = ((-2 * lin_b1) + lin_a),
         Y_End = ((5 * lin_b1) + lin_a))

(la_figure2 <- la_traits %>%
    # filter(Gen2 %in% c(4:10)) %>% 
    ggplot(aes(group = factor(Gen2),
               color = Avg_LA)) +
    geom_point(aes(x = Coded_Trait,
                   y = Value),
               size = 4,
               show.legend = TRUE) +
    geom_segment(aes(x = -2,y = Y_Start,
                     xend = 5, yend = Y_End,
                     color = Avg_LA),
                 size = .5,
                 show.legend = FALSE) +
    geom_point(data = par_la_traits_format,
               aes(x = Coded_Trait,
                   y = Value),
               size = 6,
               alpha = .75,
               color = "Black",
               show.legend = FALSE) +
    geom_segment(data = par_la_traits_format,
                 aes(x = -2,y = Y_Start,
                     xend = 5, yend = Y_End,
                     color = Avg_LA),
                 size = 2,
                 color = "Black",
                 show.legend = FALSE) +
    theme_bw() +
    coord_flip() +
    scale_y_continuous(breaks = seq(from = 30, to = 90, by = 5), 
                       minor_breaks = NULL,
                       limits  = c(30,90)) +
    facet_wrap(~Pop,
               labeller = labeller(Pop = c('4' = "B73 Population",
                                           '5' = "Mo17 Population"))) +
    scale_color_viridis(guide = "colorbar",
                        name = "Canopy\nAverage") +
    theme(aspect.ratio = 1.25,
          panel.grid = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 15),
          strip.text.x = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 15, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          legend.position = "left",
          legend.box.margin = margin(-100,-100,-100,-100)))
grad_legend <- get_legend(la_figure2)

###need to plot as continuous, steal legend, and then plot it. get_legend()
colors <- ggplot_build(la_figure2)
dh_colors <- colors$data[1] %>%
  as.data.frame() %>%
  select(x,group,PANEL,colour) %>%
  rename(Coded_Trait = x,
         Gen2 = group,
         Pop = PANEL) %>%
  mutate(Pop = (as.numeric(as.character(Pop)) + 3),
         Gen2 = Gen2 + 3,
         size_pt = 2.5,
         size_ln = .6) %>% 
  right_join(la_traits, by = c("Coded_Trait","Gen2","Pop"))

par_colors <- colors$data[3] %>% 
  as.data.frame() %>% 
  select(x,group,PANEL,colour) %>%
  rename(Coded_Trait = x,
         Gen2 = group,
         Pop = PANEL) %>% 
  mutate(Pop = (as.numeric(as.character(Pop)) + 3),
         colour = case_when(
           Gen2 == 1 ~ "Blue",
           Gen2 == 2 ~ "Black",
           Gen2 == 3 ~ "Red"),
         size_pt = 4.5,
         size_ln = 1.75) %>% 
  right_join(par_la_traits, by = c("Coded_Trait","Gen2")) %>% 
  mutate(Y_Start = ((-2 * lin_b1) + lin_a),
         Y_End = ((5 * lin_b1) + lin_a))

col_list <- dh_colors %>% 
  select(colour,Avg_LA) %>% 
  unique() %>% 
  arrange(Avg_LA) %>% 
  pull(colour)

##Need to add a negative for the slope and intercept in order to reverse the plot, this
##is only for plotting purposes. Does Not change the data.
(la_figure3 <- dh_colors %>%
    ggplot(aes(group = factor(Gen2))) +
    geom_point(aes(x = Coded_Trait,
                   y = Value,
                   color = colour,
                   size = size_pt),
               alpha = .75,
               show.legend = TRUE) +
    geom_segment(aes(x = -2,y = Y_Start,
                     xend = 5, yend = Y_End,
                     color = colour,
                     size = size_ln),
                 alpha = .75,
                 show.legend = FALSE) +
    geom_point(data = par_colors,
               aes(x = Coded_Trait,
                   y = Value,
                   color = colour,
                   size = size_pt),
               show.legend = TRUE) +
    geom_segment(data = par_colors,
                 aes(x = -2,y = Y_Start,
                     xend = 5, yend = Y_End,
                     color = colour,
                     size = size_ln),
                 show.legend = FALSE) +
    ylab("Leaf Angle") +
    scale_x_continuous(name = "",
                       breaks = c(6,5,3,1,0,-2), limits = c(-4,6),
                       labels = c("Flag\nLeaf",
                                  "One below\nflag leaf",
                                  "Three above\near leaf",
                                  "One above\near leaf",
                                  "Ear\nleaf",
                                  "Two below\near leaf"),
                       position = "top") +
    theme_bw() +
    coord_flip() +
    scale_y_reverse(breaks = seq(from = 30, to = 90, by = 10),
                    limits = c(90,30),
                    minor_breaks = NULL,
                    name = "Leaf angle (y-axis)") +
    facet_wrap(~Pop,
               labeller = labeller(Pop = c('4' = "B73 Population",
                                           '5' = "Mo17 Population"))) +
    scale_size_identity() +
    scale_color_identity(guide = "legend",
                         labels = c("B73","PHW30","Mo17"),
                         breaks = c("Blue","Red","Black"),
                         name = "Parents") +
    theme(aspect.ratio = 2,
          panel.grid = element_blank(),
          plot.margin = margin(0,-100,0,-100),
          axis.text = element_text(size = 15, face = "bold"),
          axis.title = element_text(size = 15, hjust = 0),
          strip.text.x = element_text(size = 18, face = "bold"),
          strip.background = element_blank(),
          legend.text = element_text(size = 24, face = "bold"),
          legend.title = element_text(size = 24, face = "bold"),
          legend.position = "bottom") +
    guides(colour = guide_legend(override.aes = list(size = 6))))


(final_la_plot <- plot_grid(grad_legend,la_figure3,
                            nrow = 1,
                            rel_widths = c(.3,2)))


ggsave("Phenotype_Analysis/Figures/Figure3-LA_Slope_Plot.png",
       plot = final_la_plot,
       width = 12.75,
       height = 11,
       dpi = 1000)
