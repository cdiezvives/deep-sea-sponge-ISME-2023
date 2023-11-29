#' Data analysis scripts for manuscript XXX
#' 
#' Metatranscriptome of deep-sea sponge species and locations
#' 
#' 
#' 

library(vegan)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(MASS)
library(ggrepel)
library(ggplot2)

# Load previous files (~30s load)
load('data-analysis-ready/metatrans-deep-species-locations/DATA_Peptides_fromIMG_withAnnotation_allTrans.RData')

# Load the row counts and the TMM transformation from RSEM ----
counts_tb <-
  readr::read_delim(file = 'data-analysis-ready/metatrans-deep-species-locations/trinity-output/AllinOne_matrix_RSEM.isoform.counts.matrix', delim = "\t") %>%
  dplyr::rename(trans_id = `...1`)

tmm_tb <-
  readr::read_delim(file = 'data-analysis-ready/metatrans-deep-species-locations/trinity-output/AllinOne_matrix_RSEM.isoform.TMM.EXPR.matrix', delim = "\t") %>%
  dplyr::rename(trans_id = `...1`)

# Save for use in subsequent script
save(counts_tb, tmm_tb,
     file = 'data-analysis-ready/metatrans-deep-species-locations/DATA_RNA_raw_data_allTrans.RData')

# Calculate proportion of reads per peptide size ----
relcounts <-
  img_pep_trans %>% 
  mutate(pep_prop = length_pep / length_trans) %>%
  left_join(
    counts_tb,
    by ='trans_id') %>%
  mutate(
    across(.cols = contains(c('Gbar', 'Povi', 'Gpac')),
           .fns = ~ pep_prop * .))

## Calculate transcripts per million (TPM) of peptides
relRPK <-   # reads per kilobase
  relcounts %>%   
  mutate(
    across(contains(c('Gbar', 'Povi', 'Gpac')),
           ~ . / pep_Kb))

relTPM <-   # transcripts per million  
  relRPK %>%
  mutate(
    across(.cols = contains(c("DR")),
           ~. / (sum(.)/1e6 )))

# Prepare short versions of files 
relcounts_short <-
  relcounts %>% 
  dplyr::select(peps, contains("DR"))

relTPM_short <-
  relTPM %>% 
  dplyr::select(peps, contains("DR"))

# Get means of raw counts ----
colnames(relcounts_short)

relcounts_short_means <-
  relcounts_short %>% 
  mutate(GbarDR15Mean = rowMeans(.[ ,c(2:4)])) %>%  
  mutate(GpacDR15Mean = rowMeans(.[ ,c(5:7)])) %>%  
  mutate(PoviDR15Mean = rowMeans(.[ ,c(11:13)])) %>% 
  mutate(PoviDR10Mean = rowMeans(.[ ,c(8:9)])) %>% # don't include "Povi_DR10_477", see manuscript
  mutate(PoviDR4Mean = rowMeans(.[ ,c(14:16)])) %>%  
  mutate(PoviDR9Mean = rowMeans(.[ ,c(17:19)])) %>% 
  dplyr::select(peps, contains("Mean"))

relTPM_short_means <-
  relTPM_short %>% 
  mutate(GbarDR15Mean = rowMeans(.[ ,c(2:4)])) %>%  
  mutate(GpacDR15Mean = rowMeans(.[ ,c(5:7)])) %>%  
  mutate(PoviDR15Mean = rowMeans(.[ ,c(11:13)])) %>% 
  mutate(PoviDR10Mean = rowMeans(.[ ,c(8:9)])) %>% 
  mutate(PoviDR4Mean = rowMeans(.[ ,c(14:16)])) %>%  
  mutate(PoviDR9Mean = rowMeans(.[ ,c(17:19)])) %>% 
  dplyr::select(peps, contains("Mean"))

# Add calculated counts and TPM to the peptide file with the annotation (img_pep_trans_ann) ----
## From all samples
img_pep_trans_ann_relcounts <- 
  img_pep_trans_ann %>% 
  left_join(relcounts_short, by = "peps")

img_pep_trans_ann_reltpm <- 
  img_pep_trans_ann %>% 
  left_join(relTPM_short, by = "peps")

## From the mean of groups
img_pep_trans_ann_relcounts_means <- 
  img_pep_trans_ann %>% 
  left_join(relcounts_short_means, by = "peps")

img_pep_trans_ann_reltpm_means <- 
  img_pep_trans_ann %>% 
  left_join(relTPM_short_means, by = "peps")

# Separate peptides annotated as microbial ----
prok<-c("Bacteria", "Archaea")

img_pep_trans_ann_relcounts_micros <-
  img_pep_trans_ann_relcounts %>% 
  filter(kingdom %in% prok)

img_pep_trans_ann_reltpm_micros <-
  img_pep_trans_ann_reltpm %>% 
  filter(kingdom %in% prok)

img_pep_trans_ann_relcounts_means_micros <-
  img_pep_trans_ann_relcounts_means %>% 
  filter(kingdom %in% prok)

img_pep_trans_ann_reltpm_means_micros <-
  img_pep_trans_ann_reltpm_means %>% 
  filter(kingdom %in% prok)

save(
  img_pep_trans_ann_relcounts_micros,
  img_pep_trans_ann_relcounts_means_micros,
  img_pep_trans_ann_reltpm_micros,
  img_pep_trans_ann_reltpm_means_micros,
  file = "data-analysis-ready/metatrans-deep-species-locations/DATA_Peptides_fromIMG_withAnnotation_realcount&tpm_means_Microbes.RData"
)

