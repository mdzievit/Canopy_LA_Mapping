library(tidyverse)
library(PerformanceAnalytics)
library(RColorBrewer)

gen_key <- read_tsv("Phenotype_Analysis/Genotype_Coder_Key.txt") %>% 
  rename(Gen2 = Gen) %>% 
  mutate(Pop = ifelse(Pop == 4,1,ifelse(Pop == 5,2,Pop)))

trait_key <- read_tsv("Phenotype_Analysis/SAS_Output_Files/LA_Results/Coded_Traits.txt")


comb_data <- read_tsv("Phenotype_Analysis/SAS_Output_Files/LA_Results/BLUEs_Results_Full.txt",
                      na = c(".","_"))

calc_blues <- function(data, 
                       G_effect = NULL,
                       PG_effect = NULL,
                       no_parents = TRUE) {
  if(no_parents) {
  data <- data %>% 
    filter(!Gen2 %in% c(1:3))
  
  mu <-  data %>% 
    filter(Effect == "Intercept") %>% 
    select(Coded_Trait,Env,Estimate) %>% 
    rename(mu = Estimate)
  
  g <- data %>%
    filter(Effect == G_effect) %>% 
    select(Coded_Trait,Env,Pop,Gen2,Effect,Estimate) %>% 
    rename(G = Estimate) %>% 
    select(-Effect)
  
  pg <- data %>% 
    filter(Effect == PG_effect) %>% 
    select(Coded_Trait,Env,Pop,Effect,Estimate) %>% 
    rename(Pop_G = Estimate) %>% 
    select(-Effect)
  
  all <- mu %>% 
    right_join(g,
               by = c("Coded_Trait","Env")) %>% 
    left_join(pg,
              by = c("Coded_Trait","Env","Pop"))
  
  return(all)
  } else {
    data2 <- data %>% 
      filter(Grp %in% c(1,2,3))
    
    mu <-  data %>% 
      filter(Effect == "Intercept") %>% 
      select(Coded_Trait,Env,Estimate) %>% 
      rename(mu = Estimate)
    
    g <- data2 %>%
      filter(Effect == "Grp") %>%
      select(Coded_Trait,Grp,Estimate) %>% 
      rename(G = Estimate,
             Gen2 = Grp)
    
    all <- mu %>% 
      right_join(g,
                 by = 'Coded_Trait')
    
    return(all)
  }
  
}

blues <- calc_blues(data = comb_data,
                    G_effect = "Switch*Gen2*Grp(Pop)",
                    PG_effect = "Switch*Pop*Grp",
                    no_parents = TRUE) %>% 
  mutate(P_G = mu + G + Pop_G) %>% 
  select(-Pop_G,-G,-mu,-Env) %>% 
  left_join(gen_key,
            by = c("Pop","Gen2")) %>% 
  left_join(trait_key,
            by = c("Coded_Trait")) %>% 
  rename(P = P_G)

write_tsv(x = blues %>% 
            select(-Coded_Trait) %>% 
            pivot_wider(names_from = 'Trait',
                        values_from = 'P'),
          path = "Phenotype_Analysis/BLUEs/Comb_Analysis_BLUEs_Pop.txt")

traits <- blues %>%
  select(Gen2,Trait,P) %>% 
  spread(Trait,P) %>% 
  select(ULA,TUELA,FUELA,MLA,lin_b1,lin_a,lin_rsq_t)

par_blues <- calc_blues(data = comb_data,
                        no_parents = FALSE) %>% 
  mutate(P_G = mu + G) %>% 
  select(-G,-mu,-Env) %>% 
  left_join(trait_key,
            by = 'Coded_Trait') %>% 
  rename(P = P_G) %>% 
  mutate(Pedigree = case_when(
    Gen2 == 1 ~ "B73",
    Gen2 == 2 ~ "Mo17",
    Gen2 == 3 ~ "PWH30"))

write_tsv(x = par_blues %>% 
            select(-Coded_Trait) %>% 
            pivot_wider(names_from = 'Trait',
                        values_from = 'P'),
          path = "Phenotype_Analysis/BLUEs/Comb_Analysis_BLUEs_Pop_Parents.txt")
