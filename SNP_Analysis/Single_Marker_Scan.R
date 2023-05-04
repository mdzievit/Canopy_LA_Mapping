library(tidyverse)
library(qqman)
library(doParallel)

print(sessionInfo())

#######Single marker scan for B73 populations ############

##Read coded SNP data in before it sliding window
data <- read_tsv("Filtered_SNPs/Proj1_B73_Pop_Rough-Filter_No-Het_new-filter.txt",
                 na = "N",
                 col_names = FALSE)
snps <- as.tibble(t(data[1,][-1])) %>% 
  rename(BP = V1) %>% 
  mutate(BP = as.numeric(BP))

snps$CHR <- 0

snps <- as.matrix(snps)
chr <- 1
for(i in 1:(dim(snps)[1])) {
  if (i == 1) {
    snps[i,2] <- chr
  } else if ((snps[i,1] > snps[i - 1,1])) {
    snps[i,2] <- chr
  } else {
    chr <- chr + 1
    snps[i,2] <- chr
  }
}
##Create a SNP matrix
data_m <- data[-1,]
data_m <- data_m %>%
  arrange(X1) %>% 
  rename(Lines = X1)

##Format snps for manhattan plot
snps <- as.tibble(snps) %>% 
  mutate(SNP = row_number()) %>% 
  select(SNP,CHR,BP)


##Read phenotype data in. Contains the F2 and F2:3 data
pheno <- read_tsv("Filtered_SNPs/Phenotype_Data/Comb_Analysis_BLUEs.txt")

pheno_format <- pheno %>% 
  mutate(value = str_replace_all(Pedigree,c(" " = "",
                                            "/" = ".",
                                            " " = ".",
                                            "\\(" = "",
                                            "\\)" = "",
                                            "-" = ".")))
##Format the Lines
lines <- data_m$Lines %>%
  as.tibble() %>% 
  mutate(value = str_replace_all(value,c("Proj1." = "",
                                         ".R1.fastq.gz" = "")))

spd_pheno <- lines %>% 
  left_join(pheno_format %>%
              filter(Pop == 1) %>% 
              select(Trait,value,P),
            by = "value") %>% 
  select(Trait,P,value) %>% 
  spread(Trait,P,sep = "Trait_")

##Bind the genotype and phenotype data together
data_m$Lines <- lines %>% pull(value)
data_total <- cbind(spd_pheno[-1],data_m)
rm(data)
gc()

##Run loop that performs an t-test on each SNP
cl <- makeCluster(3)
registerDoParallel(cl)
pvalues_out <- foreach(trait=(1:ncol(spd_pheno) - 1),
                       .combine = 'rbind') %:%
  foreach(snp=(ncol(spd_pheno) + 1):(ncol(data_total)),
          .combine='rbind',
          .packages='tidyverse') %dopar% {
            #test = aov(data_total[,trait] ~ as.numeric(as.factor(unlist(data_total[,i]))))
            #holder = summary(test)[[1]][["Pr(>F)"]]
            #test <- lm(data_total[,trait] ~ as.numeric(as.factor(unlist(data_total[,i]))), na.action = na.omit)
            ##pvalues[i - 4, trait] <- summary(test)$coefficients[7]
            
            test_tibble <- tibble(x = as.numeric(as.factor(unlist(data_total[,snp]))),
                                  y = data_total[,trait]) %>%
              na.omit()
            test <- t.test(x = test_tibble %>% 
                             filter(x == 1) %>% 
                             pull(y),
                           y = test_tibble %>% 
                             filter(x == 2) %>% 
                             pull(y),
                           var.equal = FALSE)
            #pvalues[i - 5, trait] <- test$p.value
            return(c(trait,snp,test$p.value))
          }
stopCluster(cl)


pvalues_b73 <- pvalues_out %>% 
  as.tibble() %>% 
  rename(Trait = V1,
         SNP = V2,
         P = V3) %>%
  mutate(SNP = SNP - ncol(spd_pheno)) %>% 
  left_join(snps,
            by = "SNP") %>% 
  spread(Trait,P,sep = "_") %>% 
  rename(P = Trait_1)
##Binds SNP information to the pvalues
colnames(pvalues_b73)[4:7] = rep("P",4)

manhattan(pvalues_b73[,c(1:4)])

#######Single marker scan for Mo17 population ############

##Read coded SNP data in before it sliding window
data <- read_tsv("Mo17-Genotypes.txt",
                 col_names = FALSE) %>%
  t()
colnames(data) <- data[1,]
data <- data[-1,]

##These are the snp positions formatted for manhattan plot
snps = read_tsv("Mo17-SNPs.txt")

##Read phenotype data in. Contains the F2 and F2:3 data
pheno <- read_tsv("Mo17-Phenotypes.txt")

##Bind the genotype and phenotype data together
data <- cbind(pheno[,2:3],data)

##Set up the f2 and f3 pvalue vector for loop
len <- dim(data)[2]
f2_pvalues = c(1:(len - 3))
f3_pvalues = c(1:(len - 3))

##Run loop that performs an f-test on each SNP
for (i in 4:ncol(data) )
{
  test = aov(data$F2_MLA ~ data[,i], )
  holder = summary(test)[[1]][["Pr(>F)"]]
  f2_pvalues[i - 3] = holder[1]
  
  test = aov(data$F3_MLA ~ data[,i])
  holder = summary(test)[[1]][["Pr(>F)"]]
  f3_pvalues[i-3] = holder[1]
  
}

##Binds SNP information to the pvalues
snps_f2_mo17 = cbind(snps,f2_pvalues)
colnames(snps_f2_mo17)[4] = "P"

snps_f3_mo17 = cbind(snps,f3_pvalues)
colnames(snps_f3_mo17)[4] = "P"

manhattan(snps_f2_mo17)
manhattan(snps_f3_mo17)
##Output the datafiles

write_tsv(snps_f2_b73 %>%
            rename(F2 = P) %>%
            left_join(snps_f3_b73 %>%
                        rename(F3 = P),
                      by = c("SNP","CHR","BP")),
          path = "pvalues_B73.txt")
write_tsv(snps_f2_mo17 %>%
            rename(F2 = P) %>%
            left_join(snps_f3_mo17 %>%
                        rename(F3 = P),
                      by = c("SNP","CHR","BP")),
          path = "pvalues_Mo17.txt")

