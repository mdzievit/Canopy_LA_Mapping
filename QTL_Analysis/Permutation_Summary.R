library(tidyverse)


##Mo17 permutation file reading
files <- tibble(filename = list.files("QTL_Analysis/Mo17/Permutation_Results/",
                    full.names = TRUE))

permutations <- files %>%
  group_by(filename) %>% 
  mutate(file_contents = map(filename, ~ read_tsv(.,skip = 41,
                                                  col_names = FALSE))) %>% 
  unnest() %>% 
  separate(`X1`, into = c("Rep","Max"), sep = "        ") %>% 
  separate(filename, sep = "-C-", into = c("Delete","Trait")) %>% 
  separate(Trait, sep = "-", into = c("Trait","Delete2")) %>% 
  select(-(starts_with("Delete"))) %>% 
  mutate(Trait = str_replace(Trait,"t","")) %>% 
  mutate_if(is.character,as.numeric)


perm_sum <- permutations %>%
  arrange(Trait,desc(Max)) %>% 
  group_by(Trait) %>%
  mutate(Rank = row_number()) %>% 
  ungroup() %>% 
  filter(Rank == 50) %>% 
  mutate(LOD = (Max / (log(10) * 2)))

write_tsv(path = "QTL_Analysis/Mo17/Mo17_LOD_Cutoffs.txt",
          x = perm_sum)


##B73 permutation file reading
files <- tibble(filename = list.files("QTL_Analysis/B73/Permutation_Results/",
                                      full.names = TRUE))

permutations <- files %>%
  group_by(filename) %>% 
  mutate(file_contents = map(filename, ~ read_tsv(.,skip = 41,
                                                  col_names = FALSE))) %>% 
  unnest() %>% 
  separate(`X1`, into = c("Rep","Max"), sep = "        ") %>% 
  separate(filename, sep = "-C-", into = c("Delete","Trait")) %>% 
  separate(Trait, sep = "-", into = c("Trait","Delete2")) %>% 
  select(-(starts_with("Delete"))) %>% 
  mutate(Trait = str_replace(Trait,"t","")) %>% 
  mutate_if(is.character,as.numeric)


perm_sum <- permutations %>%
  arrange(Trait,desc(Max)) %>% 
  group_by(Trait) %>%
  mutate(Rank = row_number()) %>% 
  ungroup() %>% 
  filter(Rank == 50) %>% 
  mutate(LOD = (Max / (log(10) * 2)))

write_tsv(path = "QTL_Analysis/B73/B73_LOD_Cutoffs.txt",
          x = perm_sum)
