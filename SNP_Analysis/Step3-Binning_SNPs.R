library(tidyverse)

b73map <- read_tsv("SNP_Analysis/Filtered_SNPs/B73_Population/Genotype_Corrector/B73_Progeny_Corrected.map")

b73map_bin <- read_tsv("SNP_Analysis/Filtered_SNPs/B73_Population/Genotype_Corrector/B73_Progeny_Corrected_Bin_1mismatch.map")


binSNPs <- b73map_bin %>% 
  select(locus_name) %>% 
  separate(locus_name,into = c("locus_name","bins"),
           sep = '\\(') %>% 
  select(-bins) %>%
  mutate(Bin = row_number())

allSNPs_bins <- b73map %>% 
  select(locus_name) %>% 
  left_join(binSNPs) %>% 
  fill(Bin,.direction = "down")

write_tsv(path = "SNP_Analysis/Filtered_SNPs/B73_Population/Genotype_Corrector/B73_SNP_Bins.txt",
          x = allSNPs_bins)

###Mo17 Population #####

mo17map <- read_tsv("SNP_Analysis/Filtered_SNPs/Mo17_Population/Genotype_Corrector/Mo17_Progeny_Corrected.map")

mo17map_bin <- read_tsv("SNP_Analysis/Filtered_SNPs/Mo17_Population/Genotype_Corrector/Mo17_Progeny_Corrected_Bin_1mistmatch.map")


binSNPs <- mo17map_bin %>% 
  select(locus_name) %>% 
  separate(locus_name,into = c("locus_name","bins"),
           sep = '\\(') %>% 
  select(-bins) %>%
  mutate(Bin = row_number())

allSNPs_bins <- mo17map %>% 
  select(locus_name) %>% 
  left_join(binSNPs) %>% 
  fill(Bin,.direction = "down")

write_tsv(path = "SNP_Analysis/Filtered_SNPs/Mo17_Population/Genotype_Corrector/Mo17_SNP_Bins.txt",
          x = allSNPs_bins)
