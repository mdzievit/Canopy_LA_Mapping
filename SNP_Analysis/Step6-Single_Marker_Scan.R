library(tidyverse)
library(broom)
library(cowplot)

print(sessionInfo())

#######Single marker scan for B73 populations ############

##Read phenotype data
pheno <- read_tsv("Phenotype_Analysis/BLUPs/Comb_Analysis_BLUPs_Pop.txt") %>%
  select(-Gen2,-Pop) %>% 
  rename(Gen = Pedigree) %>% 
  gather(Trait,BLUP,-Gen) %>% 
  mutate(Gen = str_replace_all(string = str_remove_all(Gen," |[(]|[)]"),
                               pattern = "/|-",
                               replacement = "."))
  

##Read coded SNP data in before it sliding window
data <- read_tsv("SNP_Analysis/Filtered_SNPs/B73_Population/B73_Progeny_Filter1.map",
                 na = "-") %>%
  rename(SNP = locus_name) %>%
  gather(Gen,Call,-SNP) %>%
  mutate(Gen = str_replace_all(Gen,c("Proj1." = "",
                                     ".R1.fastq.gz" = "")),
         Call = as.numeric(as.factor(Call))) %>%
  left_join(pheno,
            by = "Gen",
            relationship = "many-to-many")

b73_snps <- data %>%
  select(SNP) %>%
  unique() %>%
  separate(SNP,into = c("CHR","BP"), sep = "-",
           remove = FALSE) %>%
  mutate(CHR = as.numeric(str_replace(CHR,"chr","")),
         BP = as.numeric(BP),
         SNP2 = row_number())
write_tsv(path = "SNP_Analysis/Filtered_SNPs/B73_Population/B73_single_marker_scan_snps.txt",
          x = b73_snps)

model_out <-  data %>%
  group_by(Trait,SNP) %>%
  do(tidy(aov(BLUP ~ Call,
              data = .,
              na.action = na.omit))) %>%
  ungroup()

b73_pvalues <- model_out %>%
  filter(term == "Call") %>%
  left_join(b73_snps) %>%
  arrange(p.value)

write_tsv(path = "SNP_Analysis/Filtered_SNPs/B73_Population/B73_single_marker_scan_pvalues_BLUPs.txt",
          x = b73_pvalues)

rm(data,model_out)
gc()

##Read this in if we already calculated p-values
b73_pvalues <- read_tsv(file = "SNP_Analysis/Filtered_SNPs/B73_Population/B73_single_marker_scan_pvalues_BLUPs.txt")
b73_snps <- read_tsv(file = "SNP_Analysis/Filtered_SNPs/B73_Population/B73_single_marker_scan_snps.txt",
                     col_types = cols(SNP = "c",
                                      CHR = "i",
                                      BP = "d",
                                      SNP2 = "i"))
##This creates the mid points for this new ordered data
b73_midpoints_all <- b73_snps %>%
  select(CHR,BP) %>%
  unique() %>%
  mutate(Order = row_number()) %>%
  ungroup() %>%
  group_by(CHR) %>%
  summarise(Max = max(Order),
            Min = min(Order),
            Break = (Max + Min) / 2) %>% 
  ungroup()

b73_renumber_man <- b73_snps %>%
  select(CHR,BP) %>%
  unique() %>%
  mutate(Order = row_number()) %>% 
  group_by(CHR) %>%
  summarise(Max2 = max(Order),
            Min2 = min(Order)) %>% 
  ungroup() %>% 
  mutate(Chromosome = CHR) %>% 
  left_join(b73_midpoints_all) %>%
  ungroup() %>%
  mutate(Ratio = (Max - Min)/(Max2 - Min2)) %>% 
  select(CHR,Min,Ratio)


b73_midpoints_man <- b73_renumber_man %>% 
  left_join(b73_pvalues %>% 
              select(CHR,BP) %>% 
              unique()) %>% 
  arrange(CHR,BP) %>% 
  group_by(CHR) %>%
  mutate(Order = (((row_number() - 1) * Ratio) + Min)) %>% 
  summarise(Max = max(Order),
            Min = min(Order),
            Break = (Max + Min)/2) %>% 
  ungroup()

manhattan_plotter <- function(data,cur_Trait, midpoints,renumber) {
  data %>%
    filter(Trait == cur_Trait) %>%
    arrange(CHR,BP) %>%
    left_join(renumber,by = c("CHR")) %>%
    arrange(CHR,BP) %>% 
    group_by(CHR) %>%
    mutate(Order = (((row_number() - 1) * Ratio) + Min),
           Color = CHR %% 2) %>% 
    ungroup() %>% 
    ggplot() +
    geom_vline(data = midpoints %>%
                 filter(CHR != 10),
               aes(xintercept = Max),
               linetype = "dotted",
               size = .75,
               color = "light gray") +
    geom_point(aes(x = Order, y = -log10(p.value),
                   color = factor(Color)),
               show.legend = FALSE,
               size = .5) +
    xlab("Chromosome") +
    scale_x_continuous(labels = c(1:10),
                       breaks = midpoints %>%
                         pull(Break),
                       expand = c(0.025,0)) +
    scale_color_manual(values = c("black","blue")) +
    theme_bw() +
    theme(text = element_text(size = 8),
          axis.text = element_text(size = 8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}

(b73_ULA <- manhattan_plotter(data = b73_pvalues,
                              cur_Trait = "ULA",
                              midpoints = b73_midpoints_man,
                              renumber = b73_renumber_man))
(b73_TUELA <- manhattan_plotter(data = b73_pvalues,
                                cur_Trait = "TUELA",
                                midpoints = b73_midpoints_man,
                                renumber = b73_renumber_man))
(b73_FUELA <- manhattan_plotter(data = b73_pvalues,
                                cur_Trait = "FUELA",
                                midpoints = b73_midpoints_man,
                                renumber = b73_renumber_man))
(b73_MLA <- manhattan_plotter(data = b73_pvalues,
                              cur_Trait = "MLA",
                              midpoints = b73_midpoints_man,
                              renumber = b73_renumber_man))
(b73_slope <- manhattan_plotter(data = b73_pvalues,
                              cur_Trait = "lin_b1",
                              midpoints = b73_midpoints_man,
                              renumber = b73_renumber_man))
(b73_inter <- manhattan_plotter(data = b73_pvalues,
                                cur_Trait = "lin_a",
                                midpoints = b73_midpoints_man,
                                renumber = b73_renumber_man))
(b73_avg <- manhattan_plotter(data = b73_pvalues,
                                cur_Trait = "Avg_LA",
                                midpoints = b73_midpoints_man,
                                renumber = b73_renumber_man))



plot_grid(b73_ULA,b73_TUELA,b73_FUELA,b73_MLA,ncol = 1)

manhattan_plotter(data = b73_pvalues,
                  cur_Trait = "Ear_Shank",
                  midpoints = b73_midpoints_man,
                  renumber = b73_renumber_man)

ggsave(filename = "SNP_Analysis/Filtered_SNPs/B73_Population/Combined_ComparisonFigure_Supplemental_color_BLUPs.png",
       plot = plot_grid(b73_ULA,b73_TUELA,b73_FUELA,b73_MLA,ncol = 1),
       dpi = 600,
       width = 7,
       height = 9)

rm(b73_pvalues,b73_snps)
#######Single marker scan for Mo17 populations ############
#######
##Read phenotype data in already above

##Read coded SNP data in before it sliding window
data <- read_tsv("SNP_Analysis/Filtered_SNPs/Mo17_Population/Mo17_Progeny_Filter1.map",
                 na = "-") %>%
  rename(SNP = locus_name) %>%
  gather(Gen,Call,-SNP) %>%
  mutate(Gen = str_replace_all(Gen,c("Proj1." = "",
                                     ".R1.fastq.gz" = "")),
         Call = as.numeric(as.factor(Call))) %>%
  left_join(pheno,
            by = "Gen",
            relationship = "many-to-many")



mo17_snps <- data %>%
  select(SNP) %>%
  unique() %>%
  separate(SNP,into = c("CHR","BP"), sep = "-",
           remove = FALSE) %>%
  mutate(CHR = as.numeric(str_replace(CHR,"chr","")),
         BP = as.numeric(BP),
         SNP2 = row_number())
write_tsv(path = "SNP_Analysis/Filtered_SNPs/Mo17_Population/Mo17_single_marker_scan_snps.txt",
          x = mo17_snps)

model_out <-  data %>%
  group_by(Trait,SNP) %>%
  do(tidy(aov(BLUP ~ Call,
              data = .,
              na.action = na.omit))) %>%
  ungroup()

mo17_pvalues <- model_out %>%
  filter(term == "Call") %>%
  left_join(mo17_snps) %>%
  arrange(p.value)

write_tsv(path = "SNP_Analysis/Filtered_SNPs/Mo17_Population/Mo17_single_marker_scan_pvalues_BLUPs.txt",
          x = mo17_pvalues)

rm(data,model_out)
gc()


##Run this, already ran the single marker scan
mo17_pvalues <- read_tsv(file = "SNP_Analysis/Filtered_SNPs/Mo17_Population/Mo17_single_marker_scan_pvalues_BLUPs.txt")
mo17_snps <- read_tsv(file = "SNP_Analysis/Filtered_SNPs/Mo17_Population/Mo17_single_marker_scan_snps.txt",
                      col_types = cols(SNP = "c",
                                       CHR = "i",
                                       BP = "d",
                                       SNP2 = "i"))
##This creates the mid points for this new ordered data
mo17_midpoints_all <- mo17_snps %>%
  select(CHR,BP) %>%
  unique() %>%
  mutate(Order = row_number()) %>%
  ungroup() %>%
  group_by(CHR) %>%
  summarise(Max = max(Order),
            Min = min(Order),
            Break = (Max + Min) / 2) %>% 
  ungroup()

mo17_renumber_man <- mo17_snps %>%
  select(CHR,BP) %>%
  unique() %>%
  mutate(Order = row_number()) %>% 
  group_by(CHR) %>%
  summarise(Max2 = max(Order),
            Min2 = min(Order)) %>% 
  ungroup() %>% 
  mutate(Chromosome = CHR) %>% 
  left_join(mo17_midpoints_all) %>%
  ungroup() %>%
  mutate(Ratio = (Max - Min)/(Max2 - Min2)) %>% 
  select(CHR,Min,Ratio)


mo17_midpoints_man <- mo17_renumber_man %>% 
  left_join(mo17_pvalues %>% 
              select(CHR,BP) %>% 
              unique()) %>% 
  arrange(CHR,BP) %>% 
  group_by(CHR) %>%
  mutate(Order = (((row_number() - 1) * Ratio) + Min)) %>% 
  summarise(Max = max(Order),
            Min = min(Order),
            Break = (Max + Min)/2) %>% 
  ungroup()

(mo17_ULA <- manhattan_plotter(data = mo17_pvalues,
                               cur_Trait = "ULA",
                               midpoints = mo17_midpoints_man,
                               renumber = mo17_renumber_man))
(mo17_TUELA <- manhattan_plotter(data = mo17_pvalues,
                                 cur_Trait = "TUELA",
                                 midpoints = mo17_midpoints_man,
                                 renumber = mo17_renumber_man))
(mo17_FUELA <- manhattan_plotter(data = mo17_pvalues,
                                 cur_Trait = "FUELA",
                                 midpoints = mo17_midpoints_man,
                                 renumber = mo17_renumber_man))
(mo17_MLA <- manhattan_plotter(data = mo17_pvalues,
                               cur_Trait = "MLA",
                               midpoints = mo17_midpoints_man,
                               renumber = mo17_renumber_man))

(mo17_MLA <- manhattan_plotter(data = mo17_pvalues,
                               cur_Trait = "Ear_Shank",
                               midpoints = mo17_midpoints_man,
                               renumber = mo17_renumber_man))


plot_grid(mo17_ULA,mo17_TUELA,mo17_FUELA,mo17_MLA,ncol = 1)


ggsave(filename = "SNP_Analysis/Filtered_SNPs/Mo17_Population/Combined_ComparisonFigure_Supplemental_color_BLUPs.png",
       plot = plot_grid(mo17_ULA,mo17_TUELA,mo17_FUELA,mo17_MLA,ncol = 1),
       dpi = 600,
       width = 7,
       height = 9)




comb <- plot_grid(plot_grid(b73_ULA,b73_TUELA,b73_FUELA,b73_MLA,ncol = 1),
                  plot_grid(mo17_ULA,mo17_TUELA,mo17_FUELA,mo17_MLA,ncol = 1),
                  ncol = 2)

ggsave(filename = "SNP_Analysis/Filtered_SNPs/Single_Marker_Both_Pops_BLUPs.png",
       plot = comb,
       dpi = 600)