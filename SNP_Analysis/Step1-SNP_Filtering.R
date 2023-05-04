library(tidyverse)

data <- read_tsv("SNP_Analysis/Filtered_SNPs/Parental_SNPs/Parent_SNPs_filtered.recode.vcf") %>%
  select(-QUAL,-FILTER,-INFO,-FORMAT) %>% 
  gather(Genotype, Call, -'#CHROM',-POS,-ID,-REF,-ALT) %>% 
  rename(CHROM = '#CHROM') %>%
  mutate(Genotype = str_replace(Genotype,"-","_")) %>% 
  separate(col = "Genotype",
           sep = "_",
           into = c("Genotype","Type"))

raw_parent_file <- read_tsv("SNP_Analysis/Filtered_SNPs/Parent_File_nostat.vcf",
                            comment = "##") %>% 
  select('#CHROM':FORMAT)

unique_calls <- data %>% 
  filter(Call != "./.") %>% 
  select(CHROM,POS,Genotype,Call) %>%
  arrange(CHROM,POS,Genotype,Call) %>% 
  group_by(CHROM,POS,Genotype,Call) %>%
  mutate(n = n()) %>% 
  unique() %>%
  ungroup() %>% 
  group_by(CHROM,POS,Genotype) %>% 
  mutate(Total = sum(n),
         Per = n/Total) %>% 
  ungroup() %>% 
  filter(Per > .50)


final_snps <- unique_calls %>% 
  select(-n,-Total,-Per) %>% 
  spread(Genotype,Call) %>% 
  gather(Parent_B,Call,-CHROM,-POS,-PHW30) %>% 
  na.omit() %>% 
  group_by(CHROM,POS,Parent_B) %>% 
  mutate(Remove = ifelse(PHW30 == Call,FALSE,TRUE)) %>% 
  filter(Remove) %>% 
  ungroup()

final_snps %>% 
  group_by(CHROM,POS) %>% 
  summarise(n = n()) %>% 
  filter(n > 1) %>%
  ungroup() %>% 
  summarise(overlap = n())

final_snps %>% 
  group_by(Parent_B) %>% 
  summarise(Total_SNPs = n())


B73_SNPs <- data %>%
  select(CHROM:ALT) %>% 
  unique() %>% 
  right_join(final_snps %>%
               select(CHROM:Call) %>% 
               filter(Parent_B == "B73") %>% 
               spread(Parent_B,Call)%>% 
               unique())

Mo17_SNPs <- data %>%
  select(CHROM:ALT) %>% 
  unique() %>% 
  right_join(final_snps %>%
               select(CHROM:Call) %>% 
               filter(Parent_B == "Mo17") %>% 
               spread(Parent_B,Call)%>% 
               unique())


b73_prog_data <- read_delim("SNP_Analysis/Filtered_SNPs/B73_Population/Proj1_B73_Pop_Raw_No-Het.frq",
                            col_names = paste0("Col",1:6),
                            delim = "\t",
                            skip = 1) %>% 
  separate(Col5,
           into = (c("Allele1","Freq1")),
           sep = ":") %>% 
  separate(Col6,
           into = c("Allele2","Freq2"),
           sep  = ":") %>% 
  rename(CHROM = Col1,
         POS = Col2,
         N_ALLELES = Col3,
         N_CHR = Col4) %>% 
  na.omit() %>% 
  mutate_at(vars(matches("CHROM|POS|N_CHR|N_ALLELES|Freq1|Freq2")),.funs = as.numeric)

snps1 <- unique_calls %>% 
  select(-n,-Total,-Per) %>% 
  spread(Genotype,Call) %>%
  left_join(b73_prog_data %>% 
              select(CHROM,POS,Allele1,Allele2,Freq1,Freq2)) %>% 
  filter(between(Freq1, 0.1,0.9),
         between(Freq2,0.1,0.9)) %>% 
  left_join(data %>% 
              select(CHROM,POS,REF,ALT) %>%
              unique()) %>% 
  select(CHROM:PHW30) %>% 
  mutate(B73 = case_when(
    is.na(B73) & !is.na(PHW30) & PHW30 == "1/1" ~ "0/0",
    is.na(B73) & !is.na(PHW30) & PHW30 == "0/0" ~ "1/1",
    TRUE ~ B73),
    PHW30 = case_when(
      is.na(PHW30) & !is.na(B73) & B73 == "1/1" ~ "0/0",
      is.na(PHW30) & !is.na(B73) & B73 == "0/0" ~ "1/1",
      TRUE ~ PHW30)) %>% 
  group_by(CHROM,POS) %>% 
  mutate(Remove = ifelse(PHW30 == B73,FALSE,TRUE)) %>% 
  filter(Remove) %>% 
  ungroup()

B73_SNPs <- data %>%
  select(CHROM:ALT) %>% 
  unique() %>% 
  right_join(snps1 %>%
               select(-Mo17) %>%
               unique())

write_tsv(file = "SNP_Analysis/Filtered_SNPs/Parental_SNPs/B73_Polymorphic_SNPs.txt", 
          x = B73_SNPs %>% 
            select(CHROM,POS) %>% 
            arrange(CHROM,POS))

B73_vcf <- raw_parent_file %>% 
  right_join(B73_SNPs %>% 
               rename('#CHROM' = CHROM) %>% 
               select(-Remove))

write_tsv(x = B73_vcf,
          path = "SNP_Analysis/Filtered_SNPs/B73_Population/B73_Parent_file.vcf")
############################################
mo17_prog_data <- read_delim("SNP_Analysis/Filtered_SNPs/Mo17_Population/Proj1_Mo17_Pop_Raw_No-Het.frq",
                           col_names = paste0("Col",1:6),
                           delim = "\t",
                           skip = 1) %>% 
  separate(Col5,
           into = (c("Allele1","Freq1")),
           sep = ":") %>% 
  separate(Col6,
           into = c("Allele2","Freq2"),
           sep  = ":") %>% 
  rename(CHROM = Col1,
         POS = Col2,
         N_ALLELES = Col3,
         N_CHR = Col4) %>% 
  na.omit() %>% 
  mutate_at(vars(matches("CHROM|POS|N_CHR|N_ALLELES|Freq1|Freq2")),funs(as.numeric))

snps2 <- unique_calls %>% 
  select(-n,-Total,-Per) %>% 
  spread(Genotype,Call) %>%
  left_join(mo17_prog_data %>% 
              select(CHROM,POS,Allele1,Allele2,Freq1,Freq2)) %>% 
  filter(between(Freq1, 0.3,0.7),
         between(Freq2,0.3,0.7)) %>% 
  left_join(data %>% 
              select(CHROM,POS,REF,ALT) %>%
              unique()) %>% 
  select(CHROM:PHW30) %>% 
  mutate(Mo17 = case_when(
    is.na(Mo17) & !is.na(PHW30) & PHW30 == "1/1" ~ "0/0",
    is.na(Mo17) & !is.na(PHW30) & PHW30 == "0/0" ~ "1/1",
    TRUE ~ Mo17),
    PHW30 = case_when(
      is.na(PHW30) & !is.na(Mo17) & Mo17 == "1/1" ~ "0/0",
      is.na(PHW30) & !is.na(Mo17) & Mo17 == "0/0" ~ "1/1",
      TRUE ~ PHW30)) %>% 
  group_by(CHROM,POS) %>% 
  mutate(Remove = ifelse(PHW30 == Mo17,FALSE,TRUE)) %>% 
  filter(Remove) %>% 
  ungroup()

Mo17_SNPs <- data %>%
  select(CHROM:ALT) %>% 
  unique() %>% 
  right_join(snps2 %>%
               select(-B73) %>%
               unique())

write_tsv(path = "SNP_Analysis/Filtered_SNPs/Parental_SNPs/Mo17_Polymorphic_SNPs.txt", 
          x = Mo17_SNPs %>% 
            select(CHROM,POS) %>% 
            arrange)

Mo17_vcf <- raw_parent_file %>% 
  right_join(Mo17_SNPs %>% 
               rename('#CHROM' = CHROM) %>% 
               select(-Remove))

write_tsv(x = Mo17_vcf,
          path = "SNP_Analysis/Filtered_SNPs/Mo17_Population/Mo17_Parent_file.vcf")

###Use this section to estimate genotyping errors

(error <- unique_calls %>% 
    rename(True_Call = Call) %>% 
    select(CHROM:True_Call) %>% 
    left_join(data %>% 
                filter(Call != "./.")) %>% 
    mutate(Incorrect = ifelse(Call == True_Call,0,1)) %>% 
    group_by(Genotype) %>% 
    summarise(Total_Incorrect = sum(Incorrect),
              n = n()) %>% 
    mutate(Error = Total_Incorrect/n))

write_tsv(x = error,
          path = "SNP_Analysis/Filtered_SNPs/Parental_SNPs/error_rate.txt")
