library(tidyverse)
library(PerformanceAnalytics)
library(RColorBrewer)
library(GGally)
library(lattice)
library(ggpubr)
library(cowplot)

my_lims <- c(30,90)
my_size <- 7
cor_fun <- cor_fun <- function(data, mapping, method="pearson", ndp=3, sz=5, stars=TRUE, ...){
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  corr <- cor.test(x, y, method=method)
  est <- corr$estimate
  lb.size <- 2.5 
  
  if(stars){
    stars <- c("***", "**", "*", "")[findInterval(corr$p.value, c(0, 0.001, 0.01, 0.05, 1))]
    lbl <- paste0("Corr:","\n",
                  round(est, ndp), stars)
  }else{
    lbl <- round(est, ndp)
  }
  
  ggplot(data=data, mapping=mapping) + 
    annotate("text", x=mean(x, na.rm=TRUE), y=mean(y, na.rm=TRUE), label=lbl, size=lb.size,...)+
    theme_bw(base_size = 7) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
}

limitRange <- function(data, mapping,my_color,...) { 
  ggplot(data = data, mapping = mapping, ...) + 
    geom_point(...,color = my_color,
               size = 0.5,alpha = 0.50) + 
    geom_smooth(method = "lm", se = FALSE, color = 'red', size = 0.5) +
    scale_y_continuous(limits = my_lims,breaks = c(40,60,80)) +
    scale_x_continuous(limits = my_lims,breaks = c(40,60,80)) +
    theme_bw(base_size = 7)
}

barDiagDens <- function(data, mapping,my_color,...) { 
  if(my_color == 'black') {
    my_color <- "gray"
  }
  my_line_color <- "black"
  ggplot(data = data, mapping = mapping, ...) + 
    geom_histogram(aes(y=..density..), fill = my_color,color = my_line_color) +
    geom_density(col = "red") +
    scale_x_continuous(limits = my_lims,breaks = c(40,60,80)) +
    theme_bw(base_size = 7) 
}
gen_key <- read_tsv("Phenotype_Analysis/Genotype_Coder_Key.txt") %>% 
  rename(Gen2 = Gen)

trait_key <- read_tsv("Phenotype_Analysis/SAS_Output_Files/LA_Results/Coded_Traits.txt") %>% 
  mutate(Fancy_Trait = case_when(
    Trait == "ULA" ~ "One below flag",
    Trait == "TUELA" ~ "Three above ear",
    Trait == "FUELA" ~ "One above ear",
    Trait == "MLA" ~ "Two below ear",
    Coded_Trait == 5 ~ "R-squared",
    Coded_Trait == 6 ~ "Canopy intercept",
    Coded_Trait == 7 ~ "Canopy slope",
    Coded_Trait == 8 ~ "Canopy average"))


comb_data <- read_tsv("Phenotype_Analysis/SAS_Output_Files/LA_Results/BLUPs_Results.txt",
                      na = c(".","_")) %>% 
  bind_rows(read_tsv("Phenotype_Analysis/SAS_Output_Files/LA_Results/BLUEs_Results_Parents.txt",
                     na = c(".","_")))


calc_gen_eff <- function(data, 
                         G_effect = NULL,
                         PG_effect = NULL,
                         no_parents = TRUE) {
  if(no_parents) {
    data <- data %>% 
      filter(!Gen2 %in% c(1:3))
    
    mu <-  data %>% 
      filter(Effect == "Intercept") %>% 
      select(Coded_Trait,Estimate) %>% 
      rename(mu = Estimate)
    
    g <- data %>%
      filter(Effect == G_effect) %>% 
      select(Coded_Trait,Pop,Gen2,Effect,Estimate) %>% 
      rename(G = Estimate) %>% 
      select(-Effect)
    
    pg <- data %>% 
      filter(Effect == PG_effect,
             Pop %in% c(4,5)) %>% 
      select(Coded_Trait,Pop,Effect,Estimate) %>% 
      rename(Pop_G = Estimate) %>% 
      select(-Effect)
    
    all <- mu %>% 
      right_join(g) %>% 
      left_join(pg)
    
    return(all)
  } else {
    data2 <- data %>% 
      filter(Grp %in% c(1,2,3))
    
    mu <-  data %>% 
      filter(Effect == "Intercept") %>% 
      select(Coded_Trait,Estimate) %>% 
      rename(mu = Estimate)
    
    g <- data2 %>%
      filter(Effect == "Grp") %>%
      select(Coded_Trait,Grp,Estimate) %>% 
      rename(G = Estimate,
             Gen2 = Grp)
    
    all <- mu %>% 
      right_join(g)
    
    return(all)
  }
  
}

dh_blups <- calc_gen_eff(data = comb_data,
                         G_effect = "Switch*Gen2*Grp(Pop)",
                         PG_effect = "Switch*Pop*Grp",
                         no_parents = TRUE) %>% 
  mutate(P_G = mu + G + Pop_G) %>% 
  select(-Pop_G,-G,-mu) %>% 
  left_join(gen_key) %>% 
  left_join(trait_key) %>% 
  rename(P = P_G)

write_tsv(x = dh_blups %>% 
            select(-Coded_Trait,-Fancy_Trait) %>% 
            pivot_wider(names_from = "Trait",
                   values_from = "P"),
          path = "Phenotype_Analysis/BLUPs/Comb_Analysis_BLUPs_Pop.txt")

traits <- dh_blups %>%
  select(Gen2,Pop,Fancy_Trait,P) %>% 
  spread(Fancy_Trait,P) %>% 
  select(Gen2,
         Pop,
         `One below flag`,
         `Three above ear`,
         `One above ear`,
         `Two below ear`,
         `Canopy average`,
         `Canopy intercept`,
         `Canopy slope`) %>% 
  mutate(Pop = ifelse(Pop == 4,"B73","Mo17"))

dh_summary <- traits %>% 
  gather(Trait,Value,-Gen2,-Pop) %>% 
  group_by(Pop,Trait) %>% 
  summarise(Avg = mean(Value),
            Max = max(Value),
            Min = min(Value),
            SD = sd(Value),
            n = n(),
            sdErr = SD/sqrt(n),
            .groups = 'drop')

slope_lims <- range(traits$`Canopy slope`)
slope_lims_text <- round(scales::rescale(c(30,40,60,80,90),to = slope_lims),1)[2:4]


traits$`Canopy slope` <- scales::rescale(traits$`Canopy slope`,to = my_lims)

b73_traits <- traits %>% 
  filter(Pop == 'B73')

b73_pairs <- GGally::ggpairs(data = b73_traits,columns = 3:9,
                             upper = list(continuous = wrap(cor_fun,
                                                            color = "blue")),
                             diag=list(continuous=wrap(barDiagDens, my_color = "blue")),
                             labeller = label_wrap_gen(10),
                             lower = list(continuous = wrap(limitRange,
                                                            my_color = "blue"))) +
  theme(strip.text = element_text(size = 6),
        strip.background = element_rect(fill = "white",
                                        color = 'black',
                                        linewidth = 0.2),
        panel.spacing=grid::unit(0,"lines")) +
  ggtitle(label = "B73 Population")


for(i in 1:7) {
  if(i == 7) {
    b73_pairs[7,i] <- b73_pairs[7,i] +
      scale_x_continuous(labels = slope_lims_text,breaks = c(40,60,80),limits = my_lims)    
  } else {
    b73_pairs[7,i] <- b73_pairs[7,i] +
      scale_y_continuous(labels = slope_lims_text,breaks = c(40,60,80),limits = my_lims)
  }
}

mo17_traits <- traits %>% 
  filter(Pop == 'Mo17')

mo17_pairs <- GGally::ggpairs(data = mo17_traits,columns = 3:9,
                             upper = list(continuous = wrap(cor_fun,
                                                            color = "purple")),
                             diag=list(continuous=wrap(barDiagDens, my_color = "purple")),
                             labeller = label_wrap_gen(10),
                             lower = list(continuous = wrap(limitRange,
                                                            my_color = "purple"))) +
  theme(strip.text = element_text(size = 5),
        strip.background = element_rect(fill = "white",
                                        color = 'black',
                                        linewidth = 0.2),
        panel.spacing=grid::unit(0,"lines")) +
  ggtitle(label = "Mo17 Population")

for(i in 1:7) {
  if(i == 7) {
    mo17_pairs[7,i] <- mo17_pairs[7,i] +
      scale_x_continuous(labels = slope_lims_text,breaks = c(40,60,80),limits = my_lims)    
  } else {
    mo17_pairs[7,i] <- mo17_pairs[7,i] +
      scale_y_continuous(labels = slope_lims_text,breaks = c(40,60,80),limits = my_lims)
  }
}

comb_traits <- traits

comb_pairs <- GGally::ggpairs(data = comb_traits,columns = 3:9,
                              upper = list(continuous = wrap(cor_fun,
                                                             color = "black")),
                              diag=list(continuous=wrap(barDiagDens, my_color = "black")),
                              labeller = label_wrap_gen(10),
                              lower = list(continuous = wrap(limitRange,
                                                             my_color = "black"))) +
  theme(strip.text = element_text(size = 5),
        strip.background = element_rect(fill = "white",
                                        color = 'black',
                                        linewidth  = 0.2),
        panel.spacing=grid::unit(0,"lines")) +
  ggtitle(label = "Combined Population")

for(i in 1:7) {
  if(i == 7) {
    comb_pairs[7,i] <- comb_pairs[7,i] +
      scale_x_continuous(labels = slope_lims_text,breaks = c(40,60,80),limits = my_lims)    
  } else {
    comb_pairs[7,i] <- comb_pairs[7,i] +
      scale_y_continuous(labels = slope_lims_text,breaks = c(40,60,80),limits = my_lims)
  }
}

all_plots <- plot_grid(
  ggmatrix_gtable(b73_pairs),
  ggmatrix_gtable(mo17_pairs),
  ggmatrix_gtable(comb_pairs),
  nrow = 1
)

ggsave(filename = "Phenotype_Analysis/Figures/Figure2-Correlation_BLUPs_all.png",
       plot = all_plots,
       width = 18,
       height = 6,
       dpi = 1200)


par_blues <- calc_gen_eff(data = comb_data,
                          no_parents = FALSE) %>% 
  mutate(P_G = mu + G) %>% 
  select(-G,-mu) %>% 
  left_join(trait_key) %>% 
  rename(P = P_G) %>% 
  mutate(Pedigree = case_when(
    Gen2 == 1 ~ "B73",
    Gen2 == 2 ~ "Mo17",
    Gen2 == 3 ~ "PWH30"))

write_tsv(x = par_blues %>% 
            select(-Coded_Trait,-Fancy_Trait) %>% 
            spread(Trait,P),
          path = "Phenotype_Analysis/BLUEs/Comb_Analysis_BLUEs_Pop_Parents.txt")

####summary
####

##f1_summary

f1 <- read_tsv("Phenotype_Analysis/F1-2017_2018.txt") %>% 
  mutate(Coded_Trait = case_when(
    Trait == "ULA" ~ 1,
    Trait == "TUELA" ~ 2,
    Trait == "FUELA" ~ 3,
    Trait == "MLA" ~ 4
  )) %>% 
  mutate(Stem_Pos = ifelse(Coded_Trait == 1,5,
                           ifelse(Coded_Trait == 2,3,
                                  ifelse(Coded_Trait == 3,1,
                                         ifelse(Coded_Trait == 4,-2,Coded_Trait))))) %>%
  group_by(Sample,Pedigree) %>%
  do(.,as_tibble(cbind(lin_a = (summary(lm(.$Value ~ .$Stem_Pos))$coefficients[1]),
                       lin_b1 = (summary(lm(.$Value ~ .$Stem_Pos))$coefficients[2]),
                       Avg_LA = mean(.$Value)))) %>% 
  ungroup() %>% 
  gather(Trait,Value,-Pedigree,-Sample) %>%
  bind_rows(read_tsv("Phenotype_Analysis/F1-2017_2018.txt") %>% 
              select(-Year))

f1_summary <- f1 %>%
  mutate(F1 = case_when(
    str_detect(Pedigree,"B73") ~ "B73",
    TRUE ~ "Mo17")) |> 
  group_by(Trait,F1) %>% 
  summarise(Avg = mean(Value),
            Max = max(Value),
            Min = min(Value),
            SD = sd(Value),
            n = n(),
            sdErr = SD/sqrt(n),
            .groups = 'drop')

write_tsv(x = dh_summary,
          path = "Phenotype_Analysis/Table-Summary/DH_Summary.txt")

write_tsv(x = f1_summary,
          path = "Phenotype_Analysis/Table-Summary/f1_Summary.txt")
