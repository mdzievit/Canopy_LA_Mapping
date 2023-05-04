library(tidyverse)
library(LPmerge)

##I editted the code to the LPmerge function.
##The only edits I did was had LPmerge function report back the 
##statistics for each interval so I could chose which interval to use for the chromosome
##Also it looks like some packages weren't working, so fixed those.
source(file = "SNP_Analysis/LPMerge_MD_Edit.R")



#### These bin files were the output from the Genotype Corrector Pipeline from Miao et al 2018
#### Unfortunately, the exact calls to run the algorithm were not saved, however, the explicit
#### options and criteria used were dictated in my manuscript. Additionally, an estimate of
#### genotyping errors was provided in the parental snps folder

bins <- read_tsv("SNP_Analysis/Filtered_SNPs/B73_Population/Genotype_Corrector/B73_SNP_Bins.txt") %>%
  mutate(Population = "B73") %>% 
  bind_rows(read_tsv("SNP_Analysis/Filtered_SNPs/Mo17_Population/Genotype_Corrector/Mo17_SNP_Bins.txt") %>% 
              mutate(Population  = "Mo17")) %>% 
  rename(Marker = locus_name)

##Output the total number of bins
bins %>% group_by(Population) %>% summarise(Max = max(Bin))

b73 <- read_csv("SNP_Analysis/Filtered_SNPs/B73_Population/RQTL_Analysis/B73_Progeny_Corrected.map_Bin_1mismatch_clean_dup.rqtl.csv")
b73_map <- b73[(1:2),-1] %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = "Marker") %>% 
  rename(Chromosome = 'V1',
         Position = 'V2') %>% 
  mutate(Chromosome = as.numeric(as.character(Chromosome)),
         Position = as.numeric(as.character(Position))) %>% 
  mutate(Marker = str_replace(Marker,"\\(.*",""),
         Population = "B73")

b73 <- b73[-(1:2),]

mo17 <- read_csv("SNP_Analysis/Filtered_SNPs/Mo17_Population/RQTL_Analysis/Mo17_Progeny_Corrected.map_Bin_1mismatch_clean_dup.rqtl.csv")
mo17_map <- mo17[(1:2),-1] %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = "Marker") %>% 
  rename(Chromosome = 'V1',
         Position = 'V2') %>% 
  mutate(Chromosome = as.numeric(as.character(Chromosome)),
         Position = as.numeric(as.character(Position))) %>% 
  mutate(Marker = str_replace(Marker,"\\(.*",""),
         Population = "Mo17")

mo17 <- mo17[-(1:2),]


##This reads in the genetic maps that are ordered long format.
final_geneticMap <- b73_map %>%
  bind_rows(mo17_map) %>% 
  group_by(Population,Chromosome) %>%
  mutate(Max = max(Position)) %>%
  ungroup()

##This pulls all the original SNPs from the unbinned maps, and then uses it to fill in the genetic 
##position based on the bin it is in. Since binned SNPs are almost identical, we wanted to pull
##in more data to find overlapping SNPs between the two maps
geneticMap <- read_tsv("SNP_Analysis/Filtered_SNPs/Consensus_Map/All_SNPs.txt") %>% 
  rename(Marker = locus_name) %>% 
  left_join(final_geneticMap) %>% 
  separate(Marker,
           into = c("Chr","BP"),
           sep = "-",
           remove = FALSE) %>%
  mutate(Chr = as.numeric(str_replace(Chr,"chr","")),
         BP = as.numeric(BP)) %>%
  left_join(bins) %>% 
  arrange(Population,Chr,Bin,Position) %>%
  group_by(Population,Bin) %>% 
  fill(Position,
       .direction = "down") %>% 
  fill(Chromosome,
       .direction = "down") %>% 
  select(-Chr,-BP) %>% 
  ungroup() %>% 
  group_by(Population,Chromosome) %>%
  mutate(Max = max(Position)) %>%
  ungroup() %>% 
  na.omit()


##This builds the consensus map for each chromosome. It uses the LPmerge function that 
##I editted above to construct it
##We wanted to minimize RMSE and pulled the first ranked one as the map of choice
consensusMap <- NULL
for (i in 1:10) {
  maxK <- 6
  mapList <- list("Mo17" = geneticMap %>%
                    filter(Population == "Mo17",
                           Chromosome == i) %>%
                    select(Marker,Position) %>%
                    unique() %>% 
                    as.data.frame(),
                  "B73" = geneticMap %>%
                    filter(Population == "B73",
                           Chromosome == i) %>%
                    select(Marker,Position) %>%
                    unique() %>% 
                    as.data.frame())
  makeConMap <- LPmerge(mapList,max.interval = 1:maxK)
  mapTally <- NULL
  for (i in 1:maxK){
    mapTally <- bind_rows(mapTally,
                          as.data.frame(makeConMap[i * 2]))
  }
  map_to_select <- mapTally %>% 
    filter(map == "sd") %>%
    arrange(RMSE) %>% 
    mutate(Rnk = rank(RMSE,ties.method = "first")) %>% 
    filter(Rnk == 1) %>% 
    pull(k)
  consensusMap <- bind_rows(consensusMap,
                            as.data.frame(makeConMap[(map_to_select * 2) - 1]))
}


write_tsv(consensusMap,
          "SNP_Analysis/Filtered_SNPs/Consensus_Map/Full_ConsensusMap.txt")
# consensusMap <- read_tsv("SNP_Analysis/Filtered_SNPs/Consensus_Map/Full_ConsensusMap.txt")


#####Consensus Summary######
##This is used to summarize the genetic maps for each one to report back in the manuscript

##Total genetic distance
consensusMap %>% 
  separate(marker,into = c("Chr","BP"), sep  = "-",
           remove = FALSE) %>%
  mutate(Chr = as.numeric(str_replace(Chr,"chr","")),
         BP = as.numeric(BP)) %>% 
  group_by(Chr) %>% 
  summarise(Length = max(consensus)) %>% 
  summarise(Total = sum(Length)) %>% 
  pull(Total)

##Genetic distance of B73 Population
##Summarizes the two maps:
geneticMap %>%
  select(Population,Max) %>%
  unique() %>%
  group_by(Population) %>% 
  summarise(Total = sum(Max)) %>% 
  print.data.frame()

##Number of bins used for each population
final_geneticMap %>% 
  group_by(Population) %>% 
  summarise(n = n())

##I think that this finds all the SNPs that are in common to both maps, and then pulls them. Should be the
##SNPs to use for the joint linkage map

overlap <- consensusMap %>% 
  gather(Pop,Present,-consensus,-marker) %>% 
  na.omit() %>% 
  group_by(marker) %>% 
  summarise(n = n()) %>%
  ungroup() %>% 
  filter(n > 1) %>% 
  left_join(consensusMap) %>% 
  separate(marker,into = c("Chr","BP"), sep  = "-",
           remove = FALSE) %>%
  mutate(Chr = as.numeric(str_replace(Chr,"chr","")),
         BP = as.numeric(BP)) %>%
  arrange(Chr,consensus) %>%
  gather(Population,Pop_Pos,-(marker:consensus)) %>% 
  rename(Marker = marker) %>% 
  left_join(geneticMap %>% 
              select(Marker,Population,Bin))
overlap_bins <- overlap %>% 
  select(Chr,consensus) %>% 
  unique() %>% 
  arrange(Chr,consensus) %>%
  mutate(Joint_Bin = row_number()) %>% 
  right_join(overlap)
  

##Pulls the b73_markers from the map, and then fills in the cM within each bin for markers in LD
b73_markers <- geneticMap %>% 
  select(Marker,Population,Bin) %>%
  filter(Population == "B73") %>%
  left_join(overlap_bins %>% 
              select(Chr,Marker,consensus,Population,Joint_Bin) %>% 
              filter(Population == "B73") %>% 
              unique() %>% 
              arrange(Chr,consensus)) %>% 
  arrange(Chr,consensus,Bin) %>% 
  group_by(Bin) %>% 
  fill(consensus,
       .direction = "down") %>% 
  fill(Joint_Bin,
       .direction = "down") %>% 
  select(-Chr) %>% 
  na.omit()




##Pulls the mo17_markers from the map, and then fills in the cM within each bin for markers in LD
mo17_markers <- geneticMap %>% 
  select(Marker,Population,Bin) %>%
  filter(Population == "Mo17") %>%
  left_join(overlap_bins %>% 
              select(Chr,Marker,consensus,Population,Joint_Bin) %>% 
              filter(Population == "Mo17") %>% 
              unique() %>% 
              arrange(Chr,consensus)) %>% 
  arrange(Chr,consensus,Bin) %>% 
  group_by(Bin) %>% 
  fill(consensus,
       .direction = "down") %>% 
  fill(Joint_Bin,
       .direction = "down") %>% 
  select(-Chr) %>% 
  na.omit()

b73_gen_markers <- tibble(Marker = colnames(b73)[-1]) %>%
  mutate(Marker = str_replace(Marker,"\\(.*","")) %>% 
  left_join(b73_markers) %>% 
  na.omit()

mo17_gen_markers <- tibble(Marker = colnames(mo17)[-1]) %>%
  mutate(Marker = str_replace(Marker,"\\(.*","")) %>% 
  left_join(mo17_markers) %>% 
  na.omit()


final_markers <- Reduce(intersect,list(b73_gen_markers %>% 
                                         pull(Joint_Bin),
                                       mo17_gen_markers %>% 
                                         pull(Joint_Bin)))

b73_final_markers <- b73_gen_markers %>% 
  filter(Joint_Bin %in% final_markers) %>% 
  select(-Bin)

mo17_final_markers <- mo17_gen_markers %>% 
  filter(Joint_Bin %in% final_markers) %>% 
  select(-Bin)

##Summary of the markers for missing data. Will have some markers with the same genetic distance, so we 
##want to remove markers with most missing data. Keep 1 marker per physical marker
b73_sum <- b73 %>%
  gather(Marker,Call,-id) %>%
  mutate(Marker = str_replace(Marker,"\\(.*","")) %>% 
  separate(Marker,into = c("Chr","BP"), sep  = "-",
           remove = FALSE) %>%
  mutate(Chr = as.numeric(str_replace(Chr,"chr","")),
         BP = as.numeric(BP)) %>% 
  group_by(Chr,Marker,Call) %>% 
  summarise(num = n()) %>% 
  filter(Call == "-")

##Filtering the raw markers and pulling out the ones that have consensus cM
b73_filt <- b73 %>%
  gather(Marker,Call,-id) %>%
  mutate(Marker = str_replace(Marker,"\\(.*","")) %>% 
  separate(Marker,into = c("Chr","BP"), sep  = "-",
           remove = FALSE) %>%
  mutate(Chr = as.numeric(str_replace(Chr,"chr","")),
         BP = as.numeric(BP)) %>% 
  left_join(b73_sum %>% 
              select(Marker,Chr,num)) %>% 
  mutate(num = replace_na(num,0)) %>% 
  filter(Marker %in% (b73_final_markers %>% 
           pull(Marker))) %>% 
  left_join(b73_final_markers)

##Summary of the remaining markrs, useful to have markers that are in LD, then we are ranking the markers based on missing data
b73_mark_sum <- b73_filt %>% 
  select(Marker,Chr,Joint_Bin,num) %>% 
  unique() %>% 
  arrange(Chr,Joint_Bin,num) %>% 
  group_by(Chr,Joint_Bin) %>% 
  mutate(Rank = row_number()) %>% 
  ungroup()
  

b73_final <- b73_filt %>% 
  filter(Marker %in% (b73_mark_sum %>%
                        filter(Rank == 1) %>% 
                        pull(Marker))) %>% 
  select(-(Marker:BP),-(num:consensus))


#################Create the consensus map for the invidiausl populations
b73_ind_cons <- tibble(marker = colnames(b73)[-1]) %>%
  mutate(marker = str_replace(marker,"\\(.*","")) %>% 
  left_join(consensusMap %>% 
              select(marker,consensus)) %>% 
  na.omit() %>% 
  left_join(b73 %>%
              gather(marker,Call,-id) %>% 
              mutate(marker = str_replace(marker,"\\(.*",""))) %>% 
  separate(col = marker,
           into = c("Chr","Pos"),
           sep = "-",
           remove = FALSE) %>%
  mutate(Chr2 = as.numeric(str_replace(Chr,"chr","")),
         Call = ifelse(Call == "-","-1",Call)) %>% 
  arrange(Chr2,consensus) %>% 
  spread(id,Call)

b73_chr_sum <- tibble(marker = colnames(b73)[-1]) %>% 
  mutate(marker = str_replace(marker,"\\(.*","")) %>% 
  left_join(consensusMap %>% 
              select(marker,consensus)) %>% 
  na.omit() %>% 
  separate(col = marker,
           into = c("Chr","Pos"),
           sep = "-",
           remove = FALSE) %>%
  mutate(Chr = as.numeric(str_replace(Chr,"chr",""))) %>% 
  arrange(Chr) %>% 
  select(marker,Chr) %>% 
  group_by(Chr) %>% 
  summarise(n = n())

write_tsv("SNP_Analysis/Filtered_SNPs/B73_Population/Consensus_Map/B73_Consensus_Map-All_Traits_Genotype_Spread.txt",
          x = b73_ind_cons)
write_tsv("SNP_Analysis/Filtered_SNPs/B73_Population/Consensus_Map/B73_Consensus_Map-All_Traits_Chr_Sum.txt",
          x = b73_chr_sum)

##Summary of the markers for missing data. Will have some markers with the same genetic distance, so we 
##want to remove markers with most missing data. Keep 1 marker per physical marker
mo17_sum <- mo17 %>%
  gather(Marker,Call,-id) %>%
  mutate(Marker = str_replace(Marker,"\\(.*","")) %>%
  separate(Marker,into = c("Chr","BP"), sep  = "-",
           remove = FALSE) %>%
  mutate(Chr = as.numeric(str_replace(Chr,"chr","")),
         BP = as.numeric(BP)) %>% 
  group_by(Chr,Marker,Call) %>% 
  summarise(num = n()) %>% 
  filter(Call == "-")

##Filtering the raw markers and pulling out the ones that have consensus cM
mo17_filt <- mo17 %>%
  gather(Marker,Call,-id) %>%
  mutate(Marker = str_replace(Marker,"\\(.*","")) %>%
  separate(Marker,into = c("Chr","BP"), sep  = "-",
           remove = FALSE) %>%
  mutate(Chr = as.numeric(str_replace(Chr,"chr","")),
         BP = as.numeric(BP)) %>% 
  left_join(mo17_sum %>% 
              select(Marker,Chr,num)) %>% 
  mutate(num = replace_na(num,0)) %>% 
  filter(Marker %in% (mo17_final_markers %>% 
                        pull(Marker))) %>% 
  left_join(mo17_final_markers)

##Summary of the remaining markrs, useful to have markers that are in LD, then we are ranking the markers based on missing data
mo17_mark_sum <- mo17_filt %>% 
  select(Marker,Chr,Joint_Bin,num) %>% 
  unique() %>% 
  arrange(Chr,Joint_Bin,num) %>% 
  group_by(Chr,Joint_Bin) %>% 
  mutate(Rank = row_number()) %>% 
  ungroup()


mo17_final <- mo17_filt %>% 
  filter(Marker %in% (mo17_mark_sum %>%
                        filter(Rank == 1) %>% 
                        pull(Marker))) %>% 
  select(-(Marker:BP),-(num:consensus))

#################Create the consensus map for the invidiausl populations
mo17_ind_cons <- tibble(marker = colnames(mo17)[-1]) %>%
  mutate(marker = str_replace(marker,"\\(.*","")) %>%
  mutate(marker = str_replace(marker,"\\(.*","")) %>% 
  left_join(consensusMap %>% 
              select(marker,consensus)) %>% 
  na.omit() %>% 
  left_join(mo17 %>%
              gather(marker,Call,-id) %>% 
              mutate(marker = str_replace(marker,"\\(.*",""))) %>% 
  separate(col = marker,
           into = c("Chr","Pos"),
           sep = "-",
           remove = FALSE) %>%
  mutate(Chr2 = as.numeric(str_replace(Chr,"chr","")),
         Call = ifelse(Call == "-","-1",Call)) %>% 
  arrange(Chr2,consensus) %>% 
  spread(id,Call)

mo17_chr_sum <- tibble(marker = colnames(mo17)[-1]) %>%
  mutate(marker = str_replace(marker,"\\(.*","")) %>%
  left_join(consensusMap %>% 
              select(marker,consensus)) %>% 
  na.omit() %>% 
  separate(col = marker,
           into = c("Chr","Pos"),
           sep = "-",
           remove = FALSE) %>% 
  mutate(Chr = as.numeric(str_replace(Chr,"chr",""))) %>% 
  arrange(Chr) %>%
  select(marker,Chr) %>% 
  group_by(Chr) %>% 
  summarise(n = n())

write_tsv("SNP_Analysis/Filtered_SNPs/mo17_Population/Consensus_Map/mo17_Consensus_Map-All_Traits_Genotype_Spread.txt",
          x = mo17_ind_cons)
write_tsv("SNP_Analysis/Filtered_SNPs/mo17_Population/Consensus_Map/mo17_Consensus_Map-All_Traits_Chr_Sum.txt",
          x = mo17_chr_sum)

#####Final markers format#####
#####
final_gen_calls <- b73_final %>% 
  bind_rows(mo17_final) %>%
  arrange(Joint_Bin) %>% 
  mutate(Call = case_when(
    Call == "AA" ~ 2,
    Call == "BB" ~ 0,
    Call == "-" ~ -1
    ),
  Joint_Bin_Name = paste("JB-",Joint_Bin, sep = "")) %>% 
  spread(id,Call)

write_tsv("SNP_Analysis/Filtered_SNPs/Consensus_Map/Consensus_Genotype_Calls.txt",
          x = final_gen_calls)


final_map <- b73_final_markers %>% 
  filter(Joint_Bin %in% final_markers) %>% 
  separate(Marker,into = c("Chr","BP"), sep  = "-",
           remove = FALSE) %>%
  arrange(Joint_Bin) %>%
  mutate(Chr = as.numeric(str_replace(Chr,"chr","")),
         BP = as.numeric(BP),
         Joint_Bin_Name = paste("JB-",Joint_Bin, sep = "")) %>% 
  select(Joint_Bin,Joint_Bin_Name,Chr,consensus) %>% 
  unique()
chr_summary <- final_map %>%
  group_by(Chr) %>% 
  summarise(n = n())

write_tsv("SNP_Analysis/Filtered_SNPs/Consensus_Map/Final_Consensus_Map.txt",
          x = final_map)

write_tsv("SNP_Analysis/Filtered_SNPs/Consensus_Map/Final_Consensus_chr_summary.txt",
          x = chr_summary)


##Total markers in LD with the final bins
final_map_all <- b73_markers %>% 
  filter(Joint_Bin %in% final_markers) %>% 
  bind_rows(mo17_markers %>% 
              filter(Joint_Bin %in% final_markers))

write_tsv("SNP_Analysis/Filtered_SNPs/Consensus_Map/Final_Consensus_Map_All.txt",
          x = final_map_all)
