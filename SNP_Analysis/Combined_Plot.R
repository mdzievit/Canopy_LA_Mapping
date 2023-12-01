library(tidyverse)
library(cowplot)

traits <- read_tsv("QTL_Analysis/Trait_Key.txt") %>% 
  rename(TraitName = Trait,
         Trait = Coded_Trait)

lod_cutoffs <- read_tsv("QTL_Analysis/B73/B73_LOD_Cutoffs.txt") %>% 
  mutate(Population = "B73") %>% 
  bind_rows(read_tsv("QTL_Analysis/Mo17/Mo17_LOD_Cutoffs.txt") %>% 
              mutate(Population = "Mo17")) %>% 
  select(Trait,Population,LOD_Cutoff = LOD)

data <- read_csv(file = "QTL_Analysis/B73/DH_B73-QTL_Results.csv") %>%
  mutate(Population = "B73") %>% 
  bind_rows(read_csv(file = "QTL_Analysis/Mo17/DH_Mo17-QTL_Results.csv") %>% 
              mutate(Population = "Mo17")) %>% 
  mutate(Color = Chromosome %% 2) %>% 
  rename(LOD = `LOD(H0:H1)`,
         Position = `Position(cM)`,
         Additive = `Additive(H1)`) %>% 
  select(Population,Trait,Chromosome,Position,LOD,Color) %>% 
  left_join(lod_cutoffs,
            by = c("Trait","Population")) %>% 
  left_join(traits,
            by = "Trait") %>% 
  select(-Trait) %>% 
  mutate(TraitID = case_when(
    TraitName == "ULA" ~ 1,
    TraitName == "TUELA" ~ 2,
    TraitName == "FUELA" ~ 3,
    TraitName == "MLA" ~ 4,
    TraitName == "Avg_LA" ~ 5,
    TraitName == "lin_a" ~ 6,
    TraitName == "lin_b1" ~ 7
  )) %>% 
  select(-TraitName)

##############################Supplemental Figure #######################
##This is going to create supplemental figure 1 which compares the single marker scan
##with the LOD mapping

##Created LOD map with all chromosomes
##Figures out the midpoints now for the LOD data with all the chromosomes
midpoints_all <- data %>%
  select(Population,Chromosome,Position) %>%
  unique() %>%
  arrange(Population,Chromosome,Position) %>%
  group_by(Population) %>%
  mutate(Order = row_number()) %>%
  ungroup() %>%
  group_by(Population, Chromosome) %>%
  summarise(Max = max(Order),
            Min = min(Order),
            Break = (Max + Min) / 2) %>% 
  ungroup()

max_all <- midpoints_all %>% 
  select(Population,Chromosome,Max,Min)
##This creates a plotting function to plot a genetic map given a population name and traitID

LOD_plotter <- function(data,Pop,Trait, midpoints) {
  midpoints <- midpoints %>%
    filter(Population == Pop)
  
  data_filt <- data %>%
    filter(Population == Pop,
           TraitID == Trait)
  
  cur_lod_cutoff <- data_filt$LOD_Cutoff %>% unique()
  
  data_filt %>% 
    arrange(Chromosome,Position) %>%
    mutate(Order = row_number()) %>%
    ggplot() +
    geom_vline(data = midpoints %>%
                 filter(Chromosome != 10),
               aes(xintercept = Max),
               linetype = "dotted",
               size = 1,
               color = "light gray") +
    geom_line(aes(x = Order, 
                  y = LOD, 
                  group = Chromosome, 
                  color = factor(Color)),
              linewidth = .75,
              show.legend = FALSE) +
    geom_hline(aes(yintercept = cur_lod_cutoff),color = "red") +
    xlab("Chromosome") +
    scale_y_continuous(limits = c(0,12),
                       breaks = seq(from = 0, to = 12, by = 2))  +
    scale_x_continuous(labels = c(1:10),
                       breaks = midpoints %>%
                         pull(Break),
                       expand = c(0.025,0)) +
    scale_color_manual(values = c("black","blue")) +
    theme(text = element_text(size = 8),
          axis.text = element_text(size = 8)) +
}

##Creating the supplementary comparison file of the single marker scan
##with the lnikage mapping results
##This imports the single marker scan data
man_data <- read_tsv("SNP_Analysis/Filtered_SNPs/B73_Population/B73_single_marker_scan_pvalues_BLUPs.txt") %>% 
  mutate(Population = "B73") %>% 
  bind_rows(read_tsv("SNP_Analysis/Filtered_SNPs/Mo17_Population/Mo17_single_marker_scan_pvalues_BLUPs.txt") %>% 
              mutate(Population = "Mo17")) %>%
  filter(!(Trait %in% c("Ear_Shank","lin_rsq_t"))) %>% 
  mutate(TraitID = case_when(
    Trait == "ULA" ~ 1,
    Trait == "TUELA" ~ 2,
    Trait == "FUELA" ~ 3,
    Trait == "MLA" ~ 4,
    Trait == "Avg_LA" ~ 5,
    Trait == "lin_a" ~ 6,
    Trait == "lin_b1" ~ 7
  )) %>% 
  select(TraitID,CHR,BP,Population,p.value) %>% 
  rename(P = p.value)

renumber_man <- man_data %>%
  select(Population,CHR,BP) %>%
  unique() %>%
  arrange(Population,CHR,BP) %>% 
  group_by(Population) %>%
  mutate(Order = row_number()) %>% 
  group_by(Population, CHR) %>%
  summarise(Max2 = max(Order),
            Min2 = min(Order)) %>% 
  ungroup() %>% 
  mutate(Chromosome = CHR) %>% 
  left_join(max_all,by = c("Population","Chromosome")) %>%
  ungroup() %>%
  mutate(Ratio = (Max - Min)/(Max2 - Min2)) %>% 
  select(Population,CHR,Min,Ratio)

midpoints_man <- renumber_man %>% 
  left_join(man_data %>% 
              select(Population,CHR,BP) %>% 
              unique(),
            by = c("Population","CHR")) %>% 
  arrange(Population,CHR,BP) %>% 
  group_by(Population,CHR) %>%
  mutate(Order = (((row_number() - 1) * Ratio) + Min)) %>% 
  summarise(Max = max(Order),
            Min = min(Order),
            Break = (Max + Min)/2) %>% 
  ungroup()

manhattan_plotter <- function(data,Pop,Trait, midpoints,renumber) {
  midpoints <- midpoints %>%
    filter(Population == Pop)
  
  data_filt <- data %>%
    filter(Population == Pop,
           TraitID == Trait)
  
  pval_cut <- -log10(0.05/nrow(data_filt))
  
  data_filt %>% 
    arrange(CHR,BP) %>%
    left_join(renumber,by = c("Population","CHR")) %>%
    arrange(Population,CHR,BP) %>% 
    group_by(Population,CHR) %>%
    mutate(Order = (((row_number() - 1) * Ratio) + Min),
           Color = CHR %% 2) %>% 
    ungroup() %>% 
    ggplot() +
    geom_vline(data = midpoints %>%
                 filter(CHR != 10),
               aes(xintercept = Max),
               linetype = "dotted",
               linewidth = .75,
               color = "light gray") +
    geom_hline(aes(yintercept = pval_cut),color = "red") +
    geom_point(aes(x = Order, y = -log10(P),
                   color = factor(Color)),
               show.legend = FALSE,
               size = .5) +
    xlab("Chromosome") +
    scale_x_continuous(labels = c(1:10),
                       breaks = midpoints %>%
                         pull(Break),
                       expand = c(0.025,0)) +
    scale_y_continuous(limits = c(0,12),
                       breaks = seq(from = 0, to = 12, by = 2))  +
    scale_color_manual(values = c("black","blue")) +
    theme(text = element_text(size = 8),
          axis.text = element_text(size = 8)) +
    theme_bw()
}

plots <- list()

for (i in 1:14) {
  if (i < 8) {
    plot <- plot_grid(LOD_plotter(data = data,
                                  Pop = "B73",
                                  Trait = i,
                                  midpoints = midpoints_all),
                      manhattan_plotter(data = man_data,
                                        Pop = "B73",
                                        Trait = i,
                                        midpoints = midpoints_man,
                                        renumber = renumber_man),
                      ncol = 1)
  } else {
    plot <- plot_grid(LOD_plotter(data = data,
                                  Pop = "Mo17",
                                  Trait = i - 7,
                                  midpoints = midpoints_all),
                      manhattan_plotter(data = man_data,
                                        Pop = "Mo17",
                                        Trait = i - 7,
                                        midpoints = midpoints_man,
                                        renumber = renumber_man),
                      ncol = 1)
  }
  plots[[i]] <- plot
}


popLabel_Size <- 6
locDiff <- .01

xpos = .2

plot_grid(plots[[1]],
          plots[[2]],
          ncol = 1,
          labels = c("A","B"),
          label_size = 10) +
  draw_label(label = "B73 population: One below flag",
             size = popLabel_Size,
             fontface = "bold",
             x = xpos,
             y = 1 - locDiff) +
  draw_label(label = "B73 population: Three above ear",
             size = popLabel_Size,
             fontface = "bold",
             x = xpos,
             y = .5 - locDiff)
ggsave("SNP_Analysis/Combined_Comparison_Figure_Supplement_color_pg1.png",
       width = 7,
       height = 9,
       dpi = 600)

plot_grid(plots[[3]],
          plots[[4]],
          ncol = 1,
          labels = c("C","D"),
          label_size = 10) +
  draw_label(label = "B73 population: One above ear",
             size = popLabel_Size,
             fontface = "bold",
             x = xpos,
             y = 1 - locDiff) +
  draw_label(label = "B73 population: Two below ear",
             size = popLabel_Size,
             fontface = "bold",
             x = xpos,
             y = .5 - locDiff)

ggsave("SNP_Analysis/Combined_Comparison_Figure_Supplement_color_pg2.png",
       width = 7,
       height = 9,
       dpi = 600)

plot_grid(plots[[5]],
          plots[[6]],
          ncol = 1,
          labels = c("E","F"),
          label_size = 10) +
  draw_label(label = "B73 population: Canopy average",
             size = popLabel_Size,
             fontface = "bold",
             x = xpos,
             y = 1 - locDiff) +
  draw_label(label = "B73 population: Canopy intercept",
             size = popLabel_Size,
             fontface = "bold",
             x = xpos,
             y = .5 - locDiff)
ggsave("SNP_Analysis/Combined_Comparison_Figure_Supplement_color_pg3.png",
       width = 7,
       height = 9,
       dpi = 600)

plot_grid(plots[[7]],
          plots[[8]],
          ncol = 1,
          labels = c("G","H"),
          label_size = 10) +
  draw_label(label = "B73 population: Canopy slope",
             size = popLabel_Size,
             fontface = "bold",
             x = xpos,
             y = 1 - locDiff) +
  draw_label(label = "Mo17 population: One below flag",
             size = popLabel_Size,
             fontface = "bold",
             x = xpos,
             y = .5 - locDiff)

ggsave("SNP_Analysis/Combined_Comparison_Figure_Supplement_color_pg4.png",
       width = 7,
       height = 9,
       dpi = 600)

plot_grid(plots[[9]],
          plots[[10]],
          ncol = 1,
          labels = c("I","J"),
          label_size = 10) +
  draw_label(label = "Mo17 population: Three above ear",
             size = popLabel_Size,
             fontface = "bold",
             x = xpos,
             y = 1 - locDiff) +
  draw_label(label = "Mo17 population: One above ear",
             size = popLabel_Size,
             fontface = "bold",
             x = xpos,
             y = .5 - locDiff)
ggsave("SNP_Analysis/Combined_Comparison_Figure_Supplement_color_pg5.png",
       width = 7,
       height = 9,
       dpi = 600)

plot_grid(plots[[11]],
          plots[[12]],
          ncol = 1,
          labels = c("K","L"),
          label_size = 10) +
  draw_label(label = "Mo17 population: Two below ear",
             size = popLabel_Size,
             fontface = "bold",
             x = xpos,
             y = 1 - locDiff) +
  draw_label(label = "Mo17 population: Canopy average",
             size = popLabel_Size,
             fontface = "bold",
             x = xpos,
             y = .5 - locDiff)
ggsave("SNP_Analysis/Combined_Comparison_Figure_Supplement_color_pg6.png",
       width = 7,
       height = 9,
       dpi = 600)
plot_grid(plots[[13]],
          plots[[14]],
          ncol = 1,
          labels = c("M","N"),
          label_size = 10) +
  draw_label(label = "Mo17 population: Canopy intercept",
             size = popLabel_Size,
             fontface = "bold",
             x = xpos,
             y = 1 - locDiff) +
  draw_label(label = "Mo17 population: Canopy slope",
             size = popLabel_Size,
             fontface = "bold",
             x = xpos,
             y = .5 - locDiff)
ggsave("SNP_Analysis/Combined_Comparison_Figure_Supplement_color_pg7.png",
       width = 7,
       height = 9,
       dpi = 600)
