library(qtl)
library(tidyverse)
library(ASMap)


##Read in all of the data from the binned markers (After genotype corrector)
##Use this to remove individuals and markers
mapthis_b <- read.cross(format = "csv",
                      file = "SNP_Analysis/Filtered_SNPs/B73_Population/RQTL_Analysis/B73_Progeny_Corrected.map_Bin_1mismatch.rqtl.csv",
                      na.strings = c("-","X"),
                      crosstype = "riself",
                      alleles = c("A","B"),
                      genotype = c("A","B","AB"),
                      map.function = "kosambi",
                      estimate.map = FALSE)

##Remove markers with more than 90% missing data
pos_cutoff <- summary(mapthis_b)$n.ind *.90
nt.bymar <- ntyped(mapthis_b, "mar")
todrop1 <- names(nt.bymar[nt.bymar < pos_cutoff])
mapthis_b <- drop.markers(mapthis_b, todrop1)

##Removes individuals with more than x missing data
indv_cutoff <- sum(summary(mapthis_b)$n.mar) *0.90
inds_removed <- ntyped(mapthis_b)>indv_cutoff

write_tsv(x = as_tibble(as.list(inds_removed)),
          path = "SNP_Analysis/Filtered_SNPs/B73_Population/RQTL_Analysis/B73_Inds_Removed.txt")

mapthis_b <- subset(mapthis_b, ind=(ntyped(mapthis_b)>indv_cutoff))

mapthis_b_gm <- mstmap.cross(mapthis_b,
                           bychr=TRUE,
                           anchor = TRUE,
                           dist.fun = "kosambi",
                           objective.fun = "COUNT",
                           id = "id",
                           mvest.bc = TRUE,
                           detectBadData = TRUE,
                           trace = TRUE,
                           noMap.dist = 15,
                           miss.thresh = .25,
                           noMap.size = 2)

mapthis_b_gm <- calc.errorlod(mapthis_b_gm, error.prob=0.005,
                         map.function = "kosambi",
                         version = "new")
print(toperr_b <- top.errorlod(mapthis_b_gm, cutoff=4))

mapthis_b_gm.clean <- mapthis_b_gm
for(i in 1:nrow(toperr_b)) {
  chr <- toperr_b$chr[i]
  id <- toperr_b$id[i]
  mar <- toperr_b$marker[i]
  mapthis_b_gm.clean$geno[[chr]]$data[mapthis_b_gm$pheno$id==id, mar] <- NA
}             

mapthis_b_gm.clean <- mstmap.cross(mapthis_b_gm.clean,
                           bychr=TRUE,
                           anchor = TRUE,
                           dist.fun = "kosambi",
                           objective.fun = "COUNT",
                           id = "id",
                           mvest.bc = TRUE,
                           detectBadData = TRUE,
                           trace = TRUE,
                           noMap.dist = 15,
                           miss.thresh = .25,
                           noMap.size = 2)
summary.map(mapthis_b_gm.clean)

write.cross(cross = mapthis_b_gm.clean,
            format = "csv",
            chr = c(1:10),
            filestem = "SNP_Analysis/Filtered_SNPs/B73_Population/RQTL_Analysis/B73_Progeny_Corrected.map_Bin_1mismatch_clean_dup.rqtl"
)



#######Precleaning the Mo17 data###########

##Read in all of the data from the binned markers (After genotype corrector)
##Use this to remove individuals and markers
mapthis_m <- read.cross(format = "csv",
                      file = "SNP_Analysis/Filtered_SNPs/Mo17_Population/RQTL_Analysis/Mo17_Progeny_Corrected.map_1mismatch.rqtl.csv",
                      na.strings = c("-","X"),
                      crosstype = "riself",
                      alleles = c("A","B"),
                      genotype = c("A","B"),
                      map.function = "kosambi",
                      estimate.map = FALSE)

##Remove markers with more than 90% missing data
pos_cutoff <- summary(mapthis_m)$n.ind *.90
nt.bymar <- ntyped(mapthis_m, "mar")
todrop1 <- names(nt.bymar[nt.bymar < pos_cutoff])
mapthis_m <- drop.markers(mapthis_m, todrop1)

write_tsv(x = as.tibble(todrop1),
          path = "SNP_Analysis/Filtered_SNPs/Mo17_Population/RQTL_Analysis/Mo17_Markers_Removed_Round1.txt")

##Removes individuals with more than x missing data
indv_cutoff <- sum(summary(mapthis_m)$n.mar) * 0.90
inds_removed <- ntyped(mapthis_m)>indv_cutoff
write_tsv(x = as.tibble(as.list(inds_removed)),
          path = "SNP_Analysis/Filtered_SNPs/Mo17_Population/RQTL_Analysis/Mo17_Inds_Removed.txt")

mapthis_m <- subset(mapthis_m, ind=(ntyped(mapthis_m)>indv_cutoff))

mapthis_m_gm <- mstmap.cross(mapthis_m,
                           bychr=TRUE,
                           anchor = TRUE,
                           dist.fun = "kosambi",
                           objective.fun = "COUNT",
                           id = "id",
                           mvest.bc = TRUE,
                           detectBadData = TRUE,
                           trace = TRUE,
                           noMap.dist = 15,
                           miss.thresh = .25,
                           noMap.size = 2)

mapthis_m_gm <- calc.errorlod(mapthis_m_gm, error.prob = 0.005,
                            map.function = "kosambi",
                            version = "new")

mapthis_m_gm.clean <- mapthis_m_gm
print(toperr_m <- top.errorlod(mapthis_m_gm, cutoff=4))

for(i in 1:nrow(toperr_m)) {
  chr <- toperr_m$chr[i]
  id <- toperr_m$id[i]
  mar <- toperr_m$marker[i]
  mapthis_m_gm.clean$geno[[chr]]$data[mapthis_m_gm$pheno$id==id, mar] <- NA
}             

mapthis_m_gm.clean <- mstmap.cross(mapthis_m_gm.clean,
                                 bychr = TRUE,
                                 anchor = TRUE,
                                 dist.fun = "kosambi",
                                 objective.fun = "COUNT",
                                 id = "id",
                                 mvest.bc = TRUE,
                                 detectBadData = TRUE,
                                 trace = TRUE,
                                 noMap.dist = 15,
                                 miss.thresh = .25,
                                 noMap.size = 2)


summary.map(mapthis_m_gm.clean)

write.cross(cross = mapthis_m_gm.clean,
            format = "csv",
            chr = c(1:10),
            filestem = "SNP_Analysis/Filtered_SNPs/Mo17_Population/RQTL_Analysis/Mo17_Progeny_Corrected.map_Bin_1mismatch_clean_dup.rqtl"
)
