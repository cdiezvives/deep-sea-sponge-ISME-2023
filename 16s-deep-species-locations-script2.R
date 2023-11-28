#' Data analysis scripts for manuscript XXX
#' 
#' 16S comparison among deep-sea sponge species and locations
#' 

library(vegan)
library(pheatmap)
library(tidyverse)
library(MASS)
library(ggplot2)
library(edgeR)
library(readxl)

####load("data-analysis-ready/16s-deep-species-locations/Cantabric18_DataForAnalyses.RData")

# Import data from mothur ----
counts <-
  readr::read_delim(
    file = 'final_sponge_asv.shared',
    delim = '\t',
    col_types = readr::cols('Group' = 'c', .default = 'i')
  )

taxa = read.table(
  file = 'sponge.trim.contigs.pcr.good.unique.good.filter.unique.precluster.abund.pick.pick.asv.ASV.cons.taxonomy',
  header = T,
  row.names = 1,
  stringsAsFactors = F
)

Infotable_tb <-
  read_excel('data-analysis-ready/16s-deep-species-locations/Infotable.xlsx', sheet = "Infotable")

Infotable <-
  data.frame(Infotable_tb ,row.names = "ID")

# Differentially abundant ASVs ----
## For Species
species_compare_objects <-
  list(
    gbar_gpac = list(
      col_ordering = c(1:3, 4:6), # To perform comparison Gbar vs. Gpac
      conditions = factor(c(rep("Gbar15", 3), rep("Gpac15", 3))), # To perform comparison Gbar vs. Gpac
      pair=c("Gbar15", "Gpac15")
    ),
    gbar_povi = list(
      col_ordering = c(1:3, 10:12), # To perform comparison Gbar vs. Povi
      conditions = factor(c(rep("Gbar15", 3), rep("Povi15", 3))), # To perform comparison Gbar vs. Povi
      pair=c("Gbar15", "Povi15")
    ),
    gpac_povi = list(
      col_ordering = c(4:6, 10:12), # To perform comparison Gpac vs. Povi
      conditions = factor(c(rep("Gpac15", 3), rep("Povi15", 3))), # To perform comparison Gpac vs. Povi
      pair = c("Gpac15", "Povi15")
    )
  )

## select the groups to compare
do_DE_analysis <- function(compare_object) {
  
  data <- counts[, compare_object$col_ordering]
  data <- data[rowSums(counts_RA > 0.001) >= 2, ]
  
  dge <-
    DGEList(counts = data,
            group = compare_object$conditions)
  
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge)
  
  # Run test for each comparison
  et <-
    exactTest(dge, pair = compare_object$pair) # To perform comparison Gbar vs. Gpac
  
  tTags <- topTags(et, n = NULL)
  result_table <- tTags$table
  result_table_sel <- result_table[result_table$FDR < 0.01,]
  
  if (nrow(result_table_sel) == 0){
    message('no differences found')
    return(NULL)
  } 
  
  data.frame(sampleA = compare_object$pair[1],
             sampleB = compare_object$pair[2],
             result_table_sel)
}

## Save each result in a table for the plots
result_table_sel_Gbar15Gpac15 <- do_DE_analysis(compare_object = species_compare_objects$gbar_gpac)
result_table_sel_Gbar15Povi15 <- do_DE_analysis(compare_object = species_compare_objects$gbar_povi)
result_table_sel_Gpac15Povi15 <- do_DE_analysis(compare_object = species_compare_objects$gpac_povi)

## Put them together
lrt_list_species <- list( result_table_sel_Gbar15Gpac15 ,
                          result_table_sel_Gbar15Povi15 ,
                          result_table_sel_Gpac15Povi15  )

names(lrt_list_species) <- c("GbarGpac", "GbarPovi", "GpacPovi")

save(lrt_list_species, file= "data-analysis-ready/16s-deep-species-locations/DATA_DA_lrt_list_species.RData")

## For locations
colnames(counts)

locations_compare_objects <-
  list(
    povi10_povi15 = list(
      col_ordering = c(7:9, 10:12), # To perform comparison povi10 vs povi15
      conditions = factor(c(rep("Povi10", 3), rep("Povi15", 3))), 
      pair=c("Povi10", "Povi15")
    ),
    povi10_povi4 = list(
      col_ordering = c(7:9, 13:15), # To perform comparison povi10 vs povi4
      conditions = factor(c(rep("Povi10", 3), rep("Povi4", 3))), 
      pair=c("Povi10", "Povi4")
    ),
    povi10_povi9 = list(
      col_ordering = c(7:9, 16:18), # To perform comparison povi10 vs povi9
      conditions = factor(c(rep("Povi10", 3), rep("Povi9", 3))), 
      pair = c("Povi10", "Povi9")
    ),
    povi15_povi4 = list(
      col_ordering = c(10:12, 13:15), # To perform comparison povi15 vs povi4
      conditions = factor(c(rep("Povi15", 3), rep("Povi4", 3))), 
      pair=c("Povi15", "Povi4")
    ),
    povi15_povi9 = list(
      col_ordering = c(10:12, 16:18), # To perform comparison povi15 vs povi9
      conditions = factor(c(rep("Povi15", 3), rep("Povi9", 3))), 
      pair=c("Povi15", "Povi9")
    ),   
    povi4_povi9 = list(
      col_ordering = c(13:15, 16:18), # To perform comparison povi4 vs povi9
      conditions = factor(c(rep("Povi4", 3), rep("Povi9", 3))), 
      pair=c("Povi4", "Povi9")
    )
  )

## Save each result in a table for the plots
result_table_sel_Povi10Povi15 <- do_DE_analysis(compare_object = locations_compare_objects$povi10_povi15)
result_table_sel_Povi10Povi4 <- do_DE_analysis(compare_object = locations_compare_objects$povi10_povi4)
result_table_sel_Povi10Povi9 <- do_DE_analysis(compare_object = locations_compare_objects$povi10_povi9)
result_table_sel_Povi15Povi4 <- do_DE_analysis(compare_object = locations_compare_objects$povi15_povi4) # no Differential ASV
result_table_sel_Povi15Povi9 <- do_DE_analysis(compare_object = locations_compare_objects$povi15_povi9) # no Differential ASV
result_table_sel_Povi4Povi9 <- do_DE_analysis(compare_object = locations_compare_objects$povi4_povi9) # no Differential ASV

## Put them together
lrt_list_locations <- list(result_table_sel_Povi10Povi15 ,
                          result_table_sel_Povi10Povi4 ,
                          result_table_sel_Povi10Povi9  )

names(lrt_list_locations) <- c("Povi10Povi15", "Povi10Povi4", "Povi10Povi9")

save(lrt_list_locations, file= "data-analysis-ready/16s-deep-species-locations/DATA_DA_lrt_list_locations.RData")

# Heatmap species DA-ASVs ----
anno_row <-
  seq_along(lrt_list_species) %>%   
  lapply(function(i){
    tibble(things = row.names(lrt_list_species[[i]]),
           contrast = names(lrt_list_species)[i])
  }) %>% 
  do.call(bind_rows, .) %>% 
  mutate(sig = 'diff') %>% 
  pivot_wider(names_from = contrast, values_from = sig, values_fill = NA) %>% 
  as.data.frame()

row.names(anno_row) <- anno_row$things
anno_row <- anno_row[,-1]
head(anno_row)

all_sigs <-
  lapply(lrt_list_species, row.names) %>% 
  unlist() %>% 
  unique()

colnames(counts)
counts.sig <- counts_RA[all_sigs, ]

## Plot
colnames(counts.sig)
counts.sig.spe <- counts_RA[all_sigs, c(1:6,10:12) ]
Infotable.sig.spe <- Infotable[ colnames(counts.sig.spe), ]

RdBu_ramp <-
  colorRampPalette(colors = RColorBrewer::brewer.pal(11, 'RdBu'))

pheatmap(log2(counts.sig.spe +1),
         clustering_method = 'average',
         cluster_cols = F,
         border_color =NA,
         fontsize_row = 5,
         show_rownames = F,
         show_colnames = T,
         color = c(rev(RdBu_ramp(99))),
         annotation_col = dplyr::select(Infotable.sig.spe, 'Species_short'),
         annotation_row = anno_row)

# Heatmap locations DA-ASVs ----
anno_row <-
  seq_along(lrt_list_locations) %>%   
  lapply(function(i){
    tibble(things = row.names(lrt_list_locations[[i]]),
           contrast = names(lrt_list_locations)[i])
  }) %>% 
  do.call(bind_rows, .) %>% 
  mutate(sig = 'diff') %>% 
  pivot_wider(names_from = contrast, values_from = sig, values_fill = NA) %>% 
  as.data.frame()

row.names(anno_row) <- anno_row$things
anno_row <- anno_row[,-1]
head(anno_row)

all_sigs <-
  lapply(lrt_list_locations, row.names) %>% 
  unlist() %>% 
  unique()

colnames(counts)
counts.sig <- counts_RA[all_sigs, ]

## Plot
colnames(counts.sig)
counts.sig.loc <- counts_RA[all_sigs, c(7:18) ]
Infotable.sig.loc <- Infotable[ colnames(counts.sig.loc), ]

RdBu_ramp <-
  colorRampPalette(colors = RColorBrewer::brewer.pal(11, 'RdBu'))

pheatmap(log2(counts.sig.loc +1),
         clustering_method = 'average',
         cluster_cols = F,
         border_color =NA,
         fontsize_row = 5,
         show_rownames = F,
         show_colnames = T,
         color = c(rev(RdBu_ramp(99))),
         annotation_col = dplyr::select(Infotable.sig.loc, 'Lance'),
         annotation_row = anno_row)
