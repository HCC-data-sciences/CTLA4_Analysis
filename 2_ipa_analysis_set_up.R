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




mela_up_c1 <- mela %>% filter(
                           (LogFC_Group2_vs_Group1 > 0.5 & P_Group2_vs_Group1 < 0.005)
                       )

mela_up_c2 <- mela %>% filter(
                           (LogFC_Group2_vs_Group3 > 0.5 & P_Group2_vs_Group3 < 0.005)
                       )

mela_up_c3 <- mela %>% filter(
                           (LogFC_Group1_vs_Group3 > 0.5 & P_Group1_vs_Group3 < 0.005)
                       )

mela_up_c1 <- mela_up_c1 %>% select(gene, LogFC_Group2_vs_Group1, P_Group2_vs_Group1)
mela_up_c2 <- mela_up_c2 %>% select(gene, LogFC_Group2_vs_Group3, P_Group2_vs_Group3)
mela_up_c3 <- mela_up_c3 %>% select(gene, LogFC_Group1_vs_Group3, P_Group1_vs_Group3)
mela_all_up <- bind_rows(mela_up_c1, mela_up_c2, mela_up_c3)

write.table(mela_all_up, file='new_ipa_table/melanoma_all_upregulated.tsv', row.names=FALSE, quote=FALSE, sep='\t')

mela_down_c1 <- mela %>% filter(
                           (LogFC_Group2_vs_Group1 < -0.5 & P_Group2_vs_Group1 < 0.005)
                       )

mela_down_c2 <- mela %>% filter(
                           (LogFC_Group2_vs_Group3 < -0.5 & P_Group2_vs_Group3 < 0.005)
                       )

mela_down_c3 <- mela %>% filter(
                           (LogFC_Group1_vs_Group3 < -0.5 & P_Group1_vs_Group3 < 0.005)
                       )

mela_down_c1 <- mela_down_c1 %>% select(gene, LogFC_Group2_vs_Group1, P_Group2_vs_Group1)
mela_down_c2 <- mela_down_c2 %>% select(gene, LogFC_Group2_vs_Group3, P_Group2_vs_Group3)
mela_down_c3 <- mela_down_c3 %>% select(gene, LogFC_Group1_vs_Group3, P_Group1_vs_Group3)
mela_all_down <- bind_rows(mela_down_c1, mela_down_c2, mela_down_c3)

write.table(mela_all_down, file='new_ipa_table/melanoma_all_downregulated.tsv', row.names=FALSE, quote=FALSE, sep='\t')

mela_both_c1 <- mela %>% filter(
                           (abs(LogFC_Group2_vs_Group1) > 0.5 & P_Group2_vs_Group1 < 0.005)
                       )


mela_both_c2 <- mela %>% filter(
                           (abs(LogFC_Group2_vs_Group3) > 0.5 & P_Group2_vs_Group3 < 0.005)
                       )


mela_both_c3 <- mela %>% filter(
                           (abs(LogFC_Group1_vs_Group3) > 0.5 & P_Group1_vs_Group3 < 0.005)
                       )


mela_both_c1 <- mela_both_c1 %>% select(gene, LogFC_Group2_vs_Group1, P_Group2_vs_Group1)
mela_both_c2 <- mela_both_c2 %>% select(gene, LogFC_Group2_vs_Group3, P_Group2_vs_Group3)
mela_both_c3 <- mela_both_c3 %>% select(gene, LogFC_Group1_vs_Group3, P_Group1_vs_Group3)
mela_all_both <- bind_rows(mela_both_c1, mela_both_c2, mela_both_c3)
write.table(mela_all_both, file='new_ipa_table/melanoma_all_up_and_down_regulated.tsv', row.names=FALSE, quote=FALSE, sep='\t')



###############################


colo_up_c1 <- colo %>% filter(
                           (LogFC_Group2_vs_Group1 > 0.5 & P_Group2_vs_Group1 < 0.005)
                       )

colo_up_c2 <- colo %>% filter(
                           (LogFC_Group2_vs_Group3 > 0.5 & P_Group2_vs_Group3 < 0.005)
                       )

colo_up_c3 <- colo %>% filter(
                           (LogFC_Group1_vs_Group3 > 0.5 & P_Group1_vs_Group3 < 0.005)
                       )

colo_up_c1 <- colo_up_c1 %>% select(gene, LogFC_Group2_vs_Group1, P_Group2_vs_Group1)
colo_up_c2 <- colo_up_c2 %>% select(gene, LogFC_Group2_vs_Group3, P_Group2_vs_Group3)
colo_up_c3 <- colo_up_c3 %>% select(gene, LogFC_Group1_vs_Group3, P_Group1_vs_Group3)
colo_all_up <- bind_rows(colo_up_c1, colo_up_c2, colo_up_c3)
write.table(colo_all_up, file='new_ipa_table/colorectal_all_upregulated.tsv', row.names=FALSE, quote=FALSE, sep='\t')

colo_down_c1 <- colo %>% filter(
                           (LogFC_Group2_vs_Group1 < -0.5 & P_Group2_vs_Group1 < 0.005)
                       )

colo_down_c2 <- colo %>% filter(
                           (LogFC_Group2_vs_Group3 < -0.5 & P_Group2_vs_Group3 < 0.005)
                       )

colo_down_c3 <- colo %>% filter(
                           (LogFC_Group1_vs_Group3 < -0.5 & P_Group1_vs_Group3 < 0.005)
                       )

colo_down_c1 <- colo_down_c1 %>% select(gene, LogFC_Group2_vs_Group1, P_Group2_vs_Group1)
colo_down_c2 <- colo_down_c2 %>% select(gene, LogFC_Group2_vs_Group3, P_Group2_vs_Group3)
colo_down_c3 <- colo_down_c3 %>% select(gene, LogFC_Group1_vs_Group3, P_Group1_vs_Group3)
colo_all_down <- bind_rows(colo_down_c1, colo_down_c2, colo_down_c3)
write.table(colo_all_down, file='new_ipa_table/colorectal_all_downregulated.tsv', row.names=FALSE, quote=FALSE, sep='\t')


colo_both_c1 <- colo %>% filter(
                           (abs(LogFC_Group2_vs_Group1) > 0.5 & P_Group2_vs_Group1 < 0.005)
                       )


colo_both_c2 <- colo %>% filter(
                           (abs(LogFC_Group2_vs_Group3) > 0.5 & P_Group2_vs_Group3 < 0.005)
                       )


colo_both_c3 <- colo %>% filter(
                           (abs(LogFC_Group1_vs_Group3) > 0.5 & P_Group1_vs_Group3 < 0.005)
                       )

colo_both_c1 <- colo_both_c1 %>% select(gene, LogFC_Group2_vs_Group1, P_Group2_vs_Group1)
colo_both_c2 <- colo_both_c2 %>% select(gene, LogFC_Group2_vs_Group3, P_Group2_vs_Group3)
colo_both_c3 <- colo_both_c3 %>% select(gene, LogFC_Group1_vs_Group3, P_Group1_vs_Group3)
colo_all_both <- bind_rows(colo_both_c1, colo_both_c2, colo_both_c3)
write.table(colo_all_both, file='new_ipa_table/colorectal_all_up_and_downregulated.tsv', row.names=FALSE, quote=FALSE, sep='\t')



##############################


brca_up_c1 <- brca_ctla %>% filter(
                           (LogFC_Group2_vs_Group1 > 0.5 & P_Group2_vs_Group1 < 0.005)
                       )

brca_up_c2 <- brca_ctla %>% filter(
                           (LogFC_Group2_vs_Group3 > 0.5 & P_Group2_vs_Group3 < 0.005)
                       )

brca_up_c3 <- brca_ctla %>% filter(
                           (LogFC_Group1_vs_Group3 > 0.5 & P_Group1_vs_Group3 < 0.005)
                       )

brca_up_c1 <- brca_up_c1 %>% select(gene, LogFC_Group2_vs_Group1, P_Group2_vs_Group1)
brca_up_c2 <- brca_up_c2 %>% select(gene, LogFC_Group2_vs_Group3, P_Group2_vs_Group3)
brca_up_c3 <- brca_up_c3 %>% select(gene, LogFC_Group1_vs_Group3, P_Group1_vs_Group3)
brca_all_up <- bind_rows(brca_up_c1, brca_up_c2, brca_up_c3)
write.table(brca_all_up, file='new_ipa_table/brca_ctla_all_upregulated.tsv', row.names=FALSE, quote=FALSE, sep='\t')

brca_down_c1 <- brca_ctla %>% filter(
                           (LogFC_Group2_vs_Group1 < -0.5 & P_Group2_vs_Group1 < 0.005)
                       )

brca_down_c2 <- brca_ctla %>% filter(
                           (LogFC_Group2_vs_Group3 < -0.5 & P_Group2_vs_Group3 < 0.005)
                       )

brca_down_c3 <- brca_ctla %>% filter(
                           (LogFC_Group1_vs_Group3 < -0.5 & P_Group1_vs_Group3 < 0.005)
                       )

brca_down_c1 <- brca_down_c1 %>% select(gene, LogFC_Group2_vs_Group1, P_Group2_vs_Group1)
brca_down_c2 <- brca_down_c2 %>% select(gene, LogFC_Group2_vs_Group3, P_Group2_vs_Group3)
brca_down_c3 <- brca_down_c3 %>% select(gene, LogFC_Group1_vs_Group3, P_Group1_vs_Group3)
brca_all_down <- bind_rows(brca_down_c1, brca_down_c2, brca_down_c3)
write.table(brca_all_down, file='new_ipa_table/brca_ctla_all_downregulated.tsv', row.names=FALSE, quote=FALSE, sep='\t')


brca_both_c1 <- brca_ctla %>% filter(
                           (abs(LogFC_Group2_vs_Group1) > 0.5 & P_Group2_vs_Group1 < 0.005)
                       )


brca_both_c2 <- brca_ctla %>% filter(
                           (abs(LogFC_Group2_vs_Group3) > 0.5 & P_Group2_vs_Group3 < 0.005)
                       )


brca_both_c3 <- brca_ctla %>% filter(
                           (abs(LogFC_Group1_vs_Group3) > 0.5 & P_Group1_vs_Group3 < 0.005)
                       )


brca_both_c1 <- brca_both_c1 %>% select(gene, LogFC_Group2_vs_Group1, P_Group2_vs_Group1)
brca_both_c2 <- brca_both_c2 %>% select(gene, LogFC_Group2_vs_Group3, P_Group2_vs_Group3)
brca_both_c3 <- brca_both_c3 %>% select(gene, LogFC_Group1_vs_Group3, P_Group1_vs_Group3)
brca_all_both <- bind_rows(brca_both_c1, brca_both_c2, brca_both_c3)
write.table(brca_all_both, file='new_ipa_table/brca_ctla_all_up_and_downregulated.tsv', row.names=FALSE, quote=FALSE, sep='\t')


#############################


gitr_up_c1 <- brca_gitr %>% filter(
                           (LogFC_Group2_vs_Group1 > 0.5 & P_Group2_vs_Group1 < 0.005)
                       )

gitr_up_c2 <- brca_gitr %>% filter(
                           (LogFC_Group2_vs_Group3 > 0.5 & P_Group2_vs_Group3 < 0.005)
                       )

gitr_up_c3 <- brca_gitr %>% filter(
                           (LogFC_Group1_vs_Group3 > 0.5 & P_Group1_vs_Group3 < 0.005)
                       )

gitr_up_c1 <- gitr_up_c1 %>% select(gene, LogFC_Group2_vs_Group1, P_Group2_vs_Group1)
gitr_up_c2 <- gitr_up_c2 %>% select(gene, LogFC_Group2_vs_Group3, P_Group2_vs_Group3)
gitr_up_c3 <- gitr_up_c3 %>% select(gene, LogFC_Group1_vs_Group3, P_Group1_vs_Group3)
gitr_all_up <- bind_rows(gitr_up_c1, gitr_up_c2, gitr_up_c3)
write.table(gitr_all_up, file='new_ipa_table/brca_gitr_all_upregulated.tsv', row.names=FALSE, quote=FALSE, sep='\t')


gitr_down_c1 <- brca_gitr %>% filter(
                           (LogFC_Group2_vs_Group1 < -0.5 & P_Group2_vs_Group1 < 0.005)
                       )

gitr_down_c2 <- brca_gitr %>% filter(
                           (LogFC_Group2_vs_Group3 < -0.5 & P_Group2_vs_Group3 < 0.005)
                       )

gitr_down_c3 <- brca_gitr %>% filter(
                           (LogFC_Group1_vs_Group3 < -0.5 & P_Group1_vs_Group3 < 0.005)
                       )

gitr_down_c1 <- gitr_down_c1 %>% select(gene, LogFC_Group2_vs_Group1, P_Group2_vs_Group1)
gitr_down_c2 <- gitr_down_c2 %>% select(gene, LogFC_Group2_vs_Group3, P_Group2_vs_Group3)
gitr_down_c3 <- gitr_down_c3 %>% select(gene, LogFC_Group1_vs_Group3, P_Group1_vs_Group3)
gitr_all_down <- bind_rows(gitr_down_c1, gitr_down_c2, gitr_down_c3)
write.table(gitr_all_down, file='new_ipa_table/brca_gitr_all_downregulated.tsv', row.names=FALSE, quote=FALSE, sep='\t')



gitr_both_c1 <- brca_gitr %>% filter(
                           (abs(LogFC_Group2_vs_Group1) > 0.5 & P_Group2_vs_Group1 < 0.005)
                       )


gitr_both_c2 <- brca_gitr %>% filter(
                           (abs(LogFC_Group2_vs_Group3) > 0.5 & P_Group2_vs_Group3 < 0.005)
                       )


gitr_both_c3 <- brca_gitr %>% filter(
                           (abs(LogFC_Group1_vs_Group3) > 0.5 & P_Group1_vs_Group3 < 0.005)
                       )

gitr_both_c1 <- gitr_both_c1 %>% select(gene, LogFC_Group2_vs_Group1, P_Group2_vs_Group1)
gitr_both_c2 <- gitr_both_c2 %>% select(gene, LogFC_Group2_vs_Group3, P_Group2_vs_Group3)
gitr_both_c3 <- gitr_both_c3 %>% select(gene, LogFC_Group1_vs_Group3, P_Group1_vs_Group3)
gitr_all_both <- bind_rows(gitr_down_c1, gitr_down_c2, gitr_down_c3)
write.table(gitr_all_both, file='new_ipa_table/brca_gitr_all_up_and_downregulated.tsv', row.names=FALSE, quote=FALSE, sep='\t')


#####################
