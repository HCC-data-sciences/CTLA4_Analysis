library(dplyr)
library(tidyverse)

mela <- read.csv(file='results.txt', sep='\t', header=TRUE, row.names=1)
colo <- read.csv(file='results_colorectal.txt', sep='\t', header=TRUE, row.names=1)
brca_ctla <- read.csv(file='results_brca.txt', sep='\t', header=TRUE, row.names=1)
brca_gitr <- read.csv(file='results_brca_gitr.txt', sep='\t', header=TRUE, row.names=1)

mela <- mela %>% select(Coef.Group2vsGroup1, Coef.Group2vsGroup3, Coef.Group1vsGroup3, 
                        P.value.Group2vsGroup1, P.value.Group2vsGroup3, P.value.Group1vsGroup3)

colnames(mela) <- c('LogFC_Group2_vs_Group1', 'LogFC_Group2_vs_Group3', 'LogFC_Group1_vs_Group3',
                    'P_Group2_vs_Group1', 'P_Group2_vs_Group3', 'P_Group1_vs_Group3')
mela$gene <- rownames(mela)
rownames(mela) <- NULL



colo <- colo %>% select(Coef.Group2vsGroup1, Coef.Group2vsGroup3, Coef.Group1vsGroup3, 
                        P.value.Group2vsGroup1, P.value.Group2vsGroup3, P.value.Group1vsGroup3)
colnames(colo) <- c('LogFC_Group2_vs_Group1', 'LogFC_Group2_vs_Group3', 'LogFC_Group1_vs_Group3',
                    'P_Group2_vs_Group1', 'P_Group2_vs_Group3', 'P_Group1_vs_Group3')
colo$gene <- rownames(colo)
rownames(colo) <- NULL



brca_ctla <- brca_ctla %>% select(Coef.Group2vsGroup1, Coef.Group2vsGroup3, Coef.Group1vsGroup3, 
                                  P.value.Group2vsGroup1, P.value.Group2vsGroup3, 
                                  P.value.Group1vsGroup3)
colnames(brca_ctla) <- c('LogFC_Group2_vs_Group1', 'LogFC_Group2_vs_Group3', 
                         'LogFC_Group1_vs_Group3',
                    'P_Group2_vs_Group1', 'P_Group2_vs_Group3', 'P_Group1_vs_Group3')
brca_ctla$gene <- rownames(brca_ctla)
rownames(brca_ctla) <- NULL




brca_gitr <- brca_gitr %>% select(Coef.Group2vsGroup1, Coef.Group2vsGroup3, Coef.Group1vsGroup3, 
                                  P.value.Group2vsGroup1, P.value.Group2vsGroup3, 
                                  P.value.Group1vsGroup3)
colnames(brca_gitr) <- c('LogFC_Group2_vs_Group1', 'LogFC_Group2_vs_Group3', 
                         'LogFC_Group1_vs_Group3',
                    'P_Group2_vs_Group1', 'P_Group2_vs_Group3', 'P_Group1_vs_Group3')
brca_gitr$gene <- rownames(brca_gitr)
rownames(brca_gitr) <- NULL




mela_up <- mela %>% filter(
                           (LogFC_Group2_vs_Group1 > 0.5 & P_Group2_vs_Group1 < 0.005) |
                           (LogFC_Group2_vs_Group3 > 0.5 & P_Group2_vs_Group3 < 0.005) |
                           (LogFC_Group1_vs_Group3 > 0.5 & P_Group1_vs_Group3 < 0.005)
                       )


mela_down <- mela %>% filter(
                           (LogFC_Group2_vs_Group1 < -0.5 & P_Group2_vs_Group1 < 0.005) |
                           (LogFC_Group2_vs_Group3 < -0.5 & P_Group2_vs_Group3 < 0.005) |
                           (LogFC_Group1_vs_Group3 < -0.5 & P_Group1_vs_Group3 < 0.005)
                       )


mela_both <- mela %>% filter(
                           (abs(LogFC_Group2_vs_Group1) > 0.5 & P_Group2_vs_Group1 < 0.005) |
                           (abs(LogFC_Group2_vs_Group3) > 0.5 & P_Group2_vs_Group3 < 0.005) |
                           (abs(LogFC_Group1_vs_Group3) > 0.5 & P_Group1_vs_Group3 < 0.005)
                       )

write.table(mela_up, file='ipa_tables/melanoma_ctla4_upregulated.tsv', row.names=FALSE, quote=FALSE, sep='\t')
write.table(mela_down, file='ipa_tables/melanoma_ctla4_downregulated.tsv', row.names=FALSE, quote=FALSE, sep='\t')
write.table(mela_both, file='ipa_tables/melanoma_ctla4_both.tsv', row.names=FALSE, quote=FALSE, sep='\t')
###############################

colo_up <- colo %>% filter(
                           (LogFC_Group2_vs_Group1 > 0.5 & P_Group2_vs_Group1 < 0.005) |
                           (LogFC_Group2_vs_Group3 > 0.5 & P_Group2_vs_Group3 < 0.005) |
                           (LogFC_Group1_vs_Group3 > 0.5 & P_Group1_vs_Group3 < 0.005)
                       )


colo_down <- colo %>% filter(
                           (LogFC_Group2_vs_Group1 < -0.5 & P_Group2_vs_Group1 < 0.005) |
                           (LogFC_Group2_vs_Group3 < -0.5 & P_Group2_vs_Group3 < 0.005) |
                           (LogFC_Group1_vs_Group3 < -0.5 & P_Group1_vs_Group3 < 0.005)
                       )


colo_both <- colo %>% filter(
                           (abs(LogFC_Group2_vs_Group1) > 0.5 & P_Group2_vs_Group1 < 0.005) |
                           (abs(LogFC_Group2_vs_Group3) > 0.5 & P_Group2_vs_Group3 < 0.005) |
                           (abs(LogFC_Group1_vs_Group3) > 0.5 & P_Group1_vs_Group3 < 0.005)
                       )

write.table(colo_up, file='ipa_tables/colorectal_ctla4_upregulated.tsv', row.names=FALSE, quote=FALSE, sep='\t')
write.table(colo_down, file='ipa_tables/colorectal_ctla4_downregulated.tsv', row.names=FALSE, quote=FALSE, sep='\t')
write.table(colo_both, file='ipa_tables/colorectal_ctla4_both.tsv', row.names=FALSE, quote=FALSE, sep='\t')
##############################
brca_ctla_up <- brca_ctla %>% filter(
                           (LogFC_Group2_vs_Group1 > 0.5 & P_Group2_vs_Group1 < 0.005) |
                           (LogFC_Group2_vs_Group3 > 0.5 & P_Group2_vs_Group3 < 0.005) |
                           (LogFC_Group1_vs_Group3 > 0.5 & P_Group1_vs_Group3 < 0.005)
                       )


brca_ctla_down <- brca_ctla %>% filter(
                           (LogFC_Group2_vs_Group1 < -0.5 & P_Group2_vs_Group1 < 0.005) |
                           (LogFC_Group2_vs_Group3 < -0.5 & P_Group2_vs_Group3 < 0.005) |
                           (LogFC_Group1_vs_Group3 < -0.5 & P_Group1_vs_Group3 < 0.005)
                       )


brca_ctla_both <- brca_ctla %>% filter(
                           (abs(LogFC_Group2_vs_Group1) > 0.5 & P_Group2_vs_Group1 < 0.005) |
                           (abs(LogFC_Group2_vs_Group3) > 0.5 & P_Group2_vs_Group3 < 0.005) |
                           (abs(LogFC_Group1_vs_Group3) > 0.5 & P_Group1_vs_Group3 < 0.005)
                       )
write.table(brca_ctla_up, file='ipa_tables/brca_ctla4_upregulated.tsv', row.names=FALSE, quote=FALSE, sep='\t')
write.table(brca_ctla_down, file='ipa_tables/brca_ctla4_downregulated.tsv', row.names=FALSE, quote=FALSE, sep='\t')
write.table(brca_ctla_both, file='ipa_tables/brca_ctla4_both.tsv', row.names=FALSE, quote=FALSE, sep='\t')
##############################
brca_gitr_up <- brca_gitr %>% filter(
                           (LogFC_Group2_vs_Group1 > 0.5 & P_Group2_vs_Group1 < 0.005) |
                           (LogFC_Group2_vs_Group3 > 0.5 & P_Group2_vs_Group3 < 0.005) |
                           (LogFC_Group1_vs_Group3 > 0.5 & P_Group1_vs_Group3 < 0.005)
                       )


brca_gitr_down <- brca_gitr %>% filter(
                           (LogFC_Group2_vs_Group1 < -0.5 & P_Group2_vs_Group1 < 0.005) |
                           (LogFC_Group2_vs_Group3 < -0.5 & P_Group2_vs_Group3 < 0.005) |
                           (LogFC_Group1_vs_Group3 < -0.5 & P_Group1_vs_Group3 < 0.005)
                       )


brca_gitr_both <- brca_gitr %>% filter(
                           (abs(LogFC_Group2_vs_Group1) > 0.5 & P_Group2_vs_Group1 < 0.005) |
                           (abs(LogFC_Group2_vs_Group3) > 0.5 & P_Group2_vs_Group3 < 0.005) |
                           (abs(LogFC_Group1_vs_Group3) > 0.5 & P_Group1_vs_Group3 < 0.005)
                       )

write.table(brca_gitr_up, file='ipa_tables/brca_gitr_upregulated.tsv', row.names=FALSE, quote=FALSE, sep='\t')
write.table(brca_gitr_down, file='ipa_tables/brca_gitr_downregulated.tsv', row.names=FALSE, quote=FALSE, sep='\t')
write.table(brca_gitr_both, file='ipa_tables/brca_gitr_both.tsv', row.names=FALSE, quote=FALSE, sep='\t')

