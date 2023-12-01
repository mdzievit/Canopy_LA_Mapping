library(tidyverse)
library(boot)


data_2016 <- read_tsv("Phenotype_Analysis/Trait_Files/LA-DH-2016-Tidy.txt")

data_2017 <- read_tsv("Phenotype_Analysis/Trait_Files/LA-DH-2017-Tidy.txt")

##In 2018 a DH line was out of seed so we used B73 as a placeholder. We
##need to remove that extra one because it causes problems downstream in the model
##Removing Range 6 Pass 27 recid == 4085951
data_2018 <- read_tsv("Phenotype_Analysis/Trait_Files/LA-DH-2018-Tidy.txt") %>% 
  mutate(Year = 2018) %>% 
  rename(Value = LA) %>% 
  select(RecId,Pedigree,Year,Loc,Range, Pass,Trait,Value) %>%
  filter(RecId != 4085951)


combined_data <- data_2016 %>%
  bind_rows(data_2017) %>%
  bind_rows(data_2018) %>% 
  select(-RecId) %>%
  mutate(Coded_Trait = case_when(
    Trait == "ULA" ~ 1,
    Trait == "TUELA" ~ 2,
    Trait == "FUELA" ~ 3,
    Trait == "MLA" ~ 4))

models <- combined_data %>%
  filter(Year != 2016) %>% 
  mutate(Stem_Pos = ifelse(Coded_Trait == 1,5,
                           ifelse(Coded_Trait == 2,3,
                                  ifelse(Coded_Trait == 3,1,
                                         ifelse(Coded_Trait == 4,-2,Coded_Trait))))) %>%
  group_by(Pedigree,Year,Loc,Range,Pass) %>%
  do(.,as_tibble(cbind(lin_rsq = (summary(lm(.$Value ~ .$Stem_Pos))$r.squared),
                       lin_a = (summary(lm(.$Value ~ .$Stem_Pos))$coefficients[1]),
                       lin_b1 = (summary(lm(.$Value ~ .$Stem_Pos))$coefficients[2]),
                       Avg_LA = mean(.$Value)))) %>% 
  ungroup() %>% 
  mutate(lin_rsq_t = logit(lin_rsq)) %>%
  select(-lin_rsq) %>% 
  gather(Trait,Value,-Pedigree,-Year,-Loc,-Range,-Pass) %>% 
  mutate(Coded_Trait = as.numeric(case_when(
    Trait == "lin_rsq_t" ~ "5",
    Trait == "lin_a" ~ "6",
    Trait == "lin_b1" ~ "7",
    Trait == "Avg_LA" ~ "8")))


combined_data2 <- combined_data %>%
  bind_rows(models)

##Converting genotypes (Pedigree) into a number
##and adding which pop the genotype came from. Converting to a number.
##Group by pop, then sort by genotype. Ungroup. Give row_number + 3 to 
##each genotype name. This allows parents to be 1,2,3

genotypes <- combined_data2 %>%
  select(Pedigree) %>%
  unique() %>%
  filter(!Pedigree %in% c("B73","Mo17","PHW30")) %>%
  mutate(Pop = str_match(Pedigree, "B73|Mo17")[,1],
         Pop = ifelse(Pop == "B73",4,5)) %>%
  group_by(Pop) %>%
  arrange(Pedigree) %>%
  mutate() %>% 
  ungroup() %>% 
  mutate(Gen = row_number() + 3)


# #Write this key
write_tsv(x = genotypes,
          path = "Phenotype_Analysis/Genotype_Coder_Key.txt")

new_ranges <- combined_data2 %>%
  select(Loc,Year,Range) %>%
  unique() %>%
  group_by(Year,Loc) %>%
  mutate(New_Range = row_number())

coded_data <- combined_data2 %>%
  left_join(genotypes, by = "Pedigree") %>%
  filter(!Pedigree %in% c("B73","Mo17","PHW30")) %>%
  bind_rows(combined_data2 %>%
              filter(Pedigree %in% c("B73","Mo17","PHW30")) %>%
              mutate(Gen = case_when(
                Pedigree == "B73" ~ 1,
                Pedigree == "Mo17" ~ 2,
                Pedigree == "PHW30" ~ 3),
                Pop = 3)) %>%
  left_join(combined_data2 %>%
              select(Loc,Year,Range) %>%
              unique() %>%
              group_by(Year,Loc) %>%
              mutate(New_Range = row_number()) %>%
              ungroup(),
            by = c("Year","Loc","Range")) %>%
  select(-Range) %>%
  rename(Range = New_Range) %>% 
  mutate(Env = paste(Year,"_Loc",Loc,sep = ""))

write_tsv(coded_data,path = "Phenotype_Analysis/Full_LA_data-2016_2017_2018.txt",
          na = ".")

traits <- coded_data %>% 
  select(Trait,Coded_Trait) %>% 
  unique()

write_tsv(traits,
          path = "Phenotype_Analysis/SAS_Output_Files/LA_Results/Coded_Traits.txt")

write_tsv(path = "Phenotype_Analysis/SAS_Output_Files/LA_Results/Field_Layouts.txt",
          x = coded_data %>%
            na.omit() %>%
            select(Coded_Trait,Trait,Year,Loc,Env,Range,Pass,Pedigree,Pop,Gen) %>%
            unique())

parent_avg <-  coded_data %>% 
  filter(Gen %in% c(1,2,3)) %>% 
  group_by(Gen,Trait) %>% 
  summarise(Avg_Value = mean(Value,
                             na.rm = TRUE),
            StdErr = sd(Value,na.rm = TRUE)/sqrt(n()),
            .groups = 'drop')

write_tsv(parent_avg,path = "Phenotype_Analysis/Parent_Avg_All_Traits.txt",
          na = ".")
