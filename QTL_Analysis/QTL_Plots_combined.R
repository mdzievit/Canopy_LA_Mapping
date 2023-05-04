library(tidyverse)
library(broom)
library(cowplot)
library(grid)
library(gridExtra)

traits <- read_tsv("QTL_Analysis/Trait_Key.txt") %>% 
  rename(TraitName = Trait,
         Trait = Coded_Trait) %>% 
  mutate(Fancy_Trait = case_when(
    Trait == 9 ~ "One below flag",
    Trait == 8 ~ "Three above ear",
    Trait == 3 ~ "One above ear",
    Trait == 7 ~ "Two below ear",
    Trait == 6 ~ "R-squared",
    Trait == 4 ~ "Canopy intercept",
    Trait == 5 ~ "Canopy slope",
    Trait == 1 ~ "Canopy average"))

###Compiled the LOD score outputs from IcIM software into one file to create 
###the QTL maps (for plotting figure 2)
##Input the data
raw_data <- read_csv(file = "QTL_Analysis/B73/DH_B73-QTL_Results.csv") %>%
  mutate(Pop = "B73") %>% 
  bind_rows(read_csv(file = "QTL_Analysis/Mo17/DH_Mo17-QTL_Results.csv") %>% 
              mutate(Pop = "Mo17")) %>% 
  mutate(Color = Chromosome %% 2) %>% 
  rename(LOD = `LOD(H0:H1)`,
         Position = `Position(cM)`,
         Additive = `Additive(H1)`)

sig_lod <- read_tsv("QTL_Analysis/B73/B73_LOD_Cutoffs.txt") %>% 
  mutate(Pop = "B73") %>% 
  bind_rows(read_tsv("QTL_Analysis/Mo17/Mo17_LOD_Cutoffs.txt") %>% 
              mutate(Pop = "Mo17")) %>% 
  select(Trait,Pop,LOD) %>% 
  rename(Sig_LOD = LOD)

###summarise the significant QTL and output
sig_qtl <- raw_data %>%
  left_join(sig_lod) %>% 
  mutate(Sig = ifelse(LOD > Sig_LOD,TRUE,FALSE)) %>% 
  filter(Sig)

# write_tsv(x = sig_qtl %>%
#             filter(Pop == "B73"),
#           path = "B73/B73_Sig_QTL.txt")

# write_tsv(x = sig_qtl %>% 
#             filter(Pop == "Mo17"),
#           path = "Mo17/Mo17_Sig_QTL.txt")

##Removing chr 9 and 10 since we don't have any significant QTL there. It makes plotting easier
data <- raw_data %>% 
  filter(!(Chromosome %in% c(9,10)))

##Since we have a consensus map, we can pull the markers from all the populations and 
##use that as the order (for plotting). We want the genetic positions to align across
##the populations as we plot
unique_order <- data %>%
  select(Chromosome,Position) %>%
  unique() %>%
  arrange(Chromosome,Position) %>%
  mutate(Order = row_number())


##This pulls the order of the markers into the data frame so we can plot
data_formatted <- data %>%
  left_join(unique_order,by = c("Chromosome","Position")) %>% 
  select(Trait,Pop,Chromosome,Marker,Position,LOD,Additive,Color,Order) %>% 
  left_join(traits) %>% 
  left_join(sig_lod, by = c("Trait","Pop"))

##This calculates the chromosome midpoints and max. It is for plotting purposes, so we have
##the chromsome name in the middle of the chr and a break there
midpoints <- data_formatted %>%
  select(Chromosome,Order,LOD) %>%
  group_by(Chromosome) %>%
  summarise(Max = max(Order),
            Min = min(Order),
            Max_LOD = max(LOD)) %>%
  mutate(Break = round(((Max-Min)/2 + Min),0)) %>%
  ungroup()


##This is going to plot each LOD map (4) to a data.frame. We are going to add stuff specifically
##to each map   
cols <- c("B73" = "blue",
          "Mo17" = "black")
plot_maps <- function(chr, Pop, lod_pos) {
  
  maps <- data_formatted %>%
    filter(Chromosome %in% chr) %>% 
    group_by(TraitName) %>%
    do(
      plots = ggplot(data = .) +
        geom_hline(aes(yintercept = Sig_LOD, linetype = factor(Pop)), color = "red", size = .25, show.legend = FALSE) +
        geom_line(aes(x = Order, y = LOD, group = interaction(Chromosome,Pop), color = factor(Pop), linetype = factor(Pop)),
                  size = .35, show.legend = TRUE) +
        theme_bw() +
        scale_color_manual(values = cols,
                           name = "Population") +
        scale_linetype_discrete(name = "Population") +
        ggtitle(paste(unique(.$Fancy_Trait), sep = "")) +
        scale_y_continuous(name = " \n",
                           #limits = c(0,max(.$LOD) + 1),
                           limits = c(0,12),
                           breaks = seq(from = 0, to = 12, by = 2))  +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text.x = element_text(face = "bold"),
              plot.title = element_text(size = 10, face = "bold"),
              legend.direction = "horizontal",
              legend.position = c(.99,1.35),
              legend.justification = c(1,1),
              legend.text = element_text(size = 6, face = "bold"),
              legend.title = element_text(size = 7, face = "bold"),
              legend.margin = margin(c(0,0,0,0)))
    )
  
  
  ##This is going to add the chromsome breaks and chromsome lines
  for (i in 1:dim(maps)[1]) {
    holder <- midpoints %>% filter(Chromosome %in% chr)
    maps$plots[i][[1]] <- maps$plots[i][[1]] +
      scale_x_continuous(breaks = holder$Break,
                         labels = holder$Chromosome,
                         expand = c(.005,.005),
                         limits = c(min(unique_order %>% filter(Chromosome %in% chr) %>%  pull(Order)),
                                    max(unique_order %>% filter(Chromosome %in% chr) %>% pull(Order))),
                         name = element_blank()) +
      geom_vline(data = holder,
                 aes(xintercept = Max),color = "lightGray", 
                 size = .5, linetype = "dotted")
  }
  
  ##This is going to actually plot the maps and add some labels
  popLabel_Size <- 7
  locDiff <- .04
  
  allPlots <- plot_grid(maps$plots[5][[1]],
                        maps$plots[4][[1]],
                        maps$plots[2][[1]],
                        maps$plots[3][[1]],
                        maps$plots[1][[1]],
                        maps$plots[6][[1]],
                        maps$plots[7][[1]],
                        ncol = 1) + 
    draw_plot_label(LETTERS[1:7],
                    rep(0,7),
                    c(1.0,0.855,0.715,0.57,0.425,0.285,0.14),
                    size = 12) +
    draw_label(label = "LOD",
               x = lod_pos, y = .5,
               size = 11,
               hjust = 0,vjust = 0,
               angle = 90,
               fontface = "bold")
  
  allPlots <- add_sub(allPlots,"Selected chromosomes genetic position (cM)", 
                      vpadding = grid::unit(0, "lines"),
                      fontface = "bold",
                      size = 11,
                      y = 1.3, x = .3, hjust = 0)
  return(allPlots)
}

ggsave(filename = "QTL_Analysis/Fig4-QTL_Maps_consensus_combined.png",
       plot = plot_grid(plot_maps(chr = c(1:10), 
                                  Pop = "Both",
                                  lod_pos = 0.05)),
       dpi = 1200,
       width = 7,
       height = 10)
