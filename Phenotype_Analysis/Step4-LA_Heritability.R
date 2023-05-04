library(tidyverse)

trait_key <- read_tsv("Phenotype_Analysis/SAS_Output_Files/LA_Results/Coded_Traits.txt")

field <- read_tsv("Phenotype_Analysis/SAS_Output_Files/LA_Results/Field_Layouts.txt")

field_sum <- field %>%
  select(Coded_Trait,Env,Range) %>% 
  unique() %>% 
  group_by(Coded_Trait,Env) %>% 
  summarise(n = n())

env_err_var <- field %>% 
  select(Coded_Trait,Trait) %>% 
  unique() %>% 
  left_join(field %>%   
              select(Coded_Trait,Env) %>% 
              unique() %>% 
              group_by(Coded_Trait) %>% 
              summarise(Num_Env = n()),
            by = "Coded_Trait")

varComp <- read_tsv("Phenotype_Analysis/SAS_Output_Files/LA_Results/Variance_Components.txt",
                    na = ".")

her.line <- varComp %>%
  select(Coded_Trait,CovParm,Estimate) %>% 
  spread(CovParm,Estimate) %>%
  rename(G = 'Switch*Gen2*Grp(Pop)',
         GE = 'Swi*Env*Gen*Grp(Pop)',
         Pop_G = 'Switch*Pop*Grp',
         Pop_GE = 'Switch*Env*Pop*Grp') %>% 
  left_join(env_err_var, by = "Coded_Trait") %>% 
  mutate(h = ((Pop_G + G)/(Pop_G + G + (GE / Num_Env) + (Pop_GE / Num_Env) + (Residual/Num_Env)))) %>% 
  left_join(trait_key)

var_spd <- varComp %>% 
  select(Coded_Trait,CovParm,Estimate,ProbZ) %>% 
  left_join(trait_key) %>% 
  mutate(New_Estimate = case_when(
    ProbZ < .01 ~ paste(Estimate,"**",sep = ""),
    ProbZ < .05 ~ paste(Estimate,"*",sep = ""),
    TRUE ~ as.character(Estimate))
  ) %>% 
  select(-Estimate,-ProbZ) %>% 
  spread(CovParm,New_Estimate) %>% 
  left_join(env_err_var) %>% 
  left_join(her.line %>%
              select(Trait,h))

write_tsv(x = var_spd,
          path = "Phenotype_Analysis/Table-Summary/Table1-Model_Results.txt")
