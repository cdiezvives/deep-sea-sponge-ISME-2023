#' Data analysis scripts for manuscript XXX
#' 
#' Metatranscriptome of deep-sea sponge species and locations
#' 

library(edgeR)
library(readxl)
library(pheatmap)
library(tidyverse)

# Load previous files
load("data-analysis-ready/metatrans-deep-species-locations/DATA_RNA_raw_data_allTrans.RData")
load("data-analysis-ready/metatrans-deep-species-locations/DATA_Peptides_fromIMG_withAnnotation_allTrans.RData")

# select microbes
counts <- data.frame(counts_tb, row.names = "trans_id")

colnames(counts)
counts_micros_spe <- counts [ img_tax_micros_trans, c(1:6,10,11,12)]
colnames(counts_micros_spe)

counts_micros_loc <- counts [ img_tax_micros_trans, c(7:18)]
colnames(counts_micros_loc)

# EdgeR (species) ----
species_compare_objects <-
  list(
    gbar_gpac = list(
      col_ordering = c(1:3, 4:6), # To perform comparison Gbar vs. Gpac
      conditions = factor(c(rep("Gbar15", 3), rep("Gpac15", 3))), # To perform comparison Gbar vs. Gpac
      pair=c("Gbar15", "Gpac15")
    ),
    gbar_povi = list(
      col_ordering = c(1:3, 7:9), # To perform comparison Gbar vs. Povi
      conditions = factor(c(rep("Gbar15", 3), rep("Povi15", 3))), # To perform comparison Gbar vs. Povi
      pair=c("Gbar15", "Povi15")
    ),
    gpac_povi = list(
      col_ordering = c(4:6, 7:9), # To perform comparison Gpac vs. Povi
      conditions = factor(c(rep("Gpac15", 3), rep("Povi15", 3))), # To perform comparison Gpac vs. Povi
      pair = c("Gpac15", "Povi15")
    )
  )

## select the groups to compare
do_DE_analysis <- function(compare_object) {
  
  data <- counts_micros_spe[, compare_object$col_ordering]
  data <- data[rowSums(cpm(data) > 10) >= 2, ]
  
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
  data.frame(sampleA = compare_object$pair[1],
             sampleB = compare_object$pair[2],
             result_table_sel)
}


## Save each result in a table for the plots
result_table_sel_Gbar15Gpac15 <-do_DE_analysis(compare_object = species_compare_objects$gbar_gpac)
result_table_sel_Gbar15Povi15 <- do_DE_analysis(compare_object = species_compare_objects$gbar_povi)
result_table_sel_Gpac15Povi15 <- do_DE_analysis(compare_object = species_compare_objects$gpac_povi)

## Put them together
lrt_list_species <- list( result_table_sel_Gbar15Gpac15 ,
                          result_table_sel_Gbar15Povi15 ,
                          result_table_sel_Gpac15Povi15  )
names(lrt_list_species) <- c("GbarGpac", "GbarPovi", "GpacPovi")


save(lrt_list_species, file= "data-analysis-ready/metatrans-deep-species-locations/DATA_lrt_list_species.RData")

# FIGURE S2 ----
## Example results Species ----
load("data-analysis-ready/metatrans-deep-species-locations/DATA_lrt_list_species.RData")
load("data-analysis-ready/metatrans-deep-species-locations/DATA_Download_IMG_ready.RData")
load("data-analysis-ready/metatrans-deep-species-locations/DATA_Peptides_fromIMG_withAnnotation_realcount&tpm_means_Microbes.RData")

## Separate samples from the 3 species in DR15
Species_names <-
  img_pep_trans_ann_relcounts_means_micros %>%
  dplyr::select(peps, contains("DR15")) %>%
  mutate(sum = rowSums(.[,2:4])) %>%
  filter(sum>0) %>%
  pull(peps)

## Select species in the tpm file
img_pep_trans_ann_reltpm_micros_spe <-
  img_pep_trans_ann_reltpm_micros %>%
  filter(peps %in% Species_names)

img_pep_trans_ann_reltpm_means_micros_spe <-
  img_pep_trans_ann_reltpm_means_micros %>%
  filter(peps %in% Species_names)

## Example with Gbar vs. Gpac (repeat with the other comparisons)
positive <- lrt_list_species$GbarGpac %>%  filter(logFC > 0)
negative <- lrt_list_species$GbarGpac %>%  filter(logFC < 0)

listas <- c(rownames(positive), rownames(negative))

img_pep_trans_ann_tb_DE <-
  img_pep_trans_ann_reltpm_means_micros_spe %>% 
  filter( trans_id %in% listas) %>% 
  mutate(sign = ifelse(trans_id %in% rownames(positive), yes = "POS", no = "NEG"))

img_pep_trans_ann_tb_DE %>% 
  count(sign) # Often each transcripts includes more than one peptide, which increases this dataset

img_pep_trans_ann_tb_DE <-
  img_pep_trans_ann_tb_DE %>% 
  mutate(comparison = "GbarGpac")

## Count how many peptides we have in each KEGG entry for positive and negative values 
Kegg_all_POS <-  
  img_pep_trans_ann_tb_DE %>%
  filter(sign == "POS") %>% 
  count(KEGG_description) 

Kegg_all_NEG <-  
  img_pep_trans_ann_tb_DE %>%
  filter(sign == "NEG") %>% 
  count(KEGG_description)

## join them
Kegg_all_POS_NEG <-
  Kegg_all_POS %>% 
  full_join(Kegg_all_NEG, by ="KEGG_description")

## Add the information of the KEGG levels
img_pep_trans_ann_tb_DE_path <-  # this file increases because one KO can appear in more than one pathway
  img_pep_trans_ann_tb_DE %>% 
  left_join(path_link , by = c("ko_id" = "KO"), relationship = "many-to-many")

## Summarize by pathway numbers
path_all_POS <- 
  img_pep_trans_ann_tb_DE_path %>%  
  filter(sign == "POS") %>% 
  count(Path_desc)

path_all_NEG <- 
  img_pep_trans_ann_tb_DE_path %>%
  filter(sign == "NEG") %>% 
  count(Path_desc)

path_all_NEG_POS <-
  path_all_NEG %>% 
  full_join(path_all_POS, by ="Path_desc")

colnames(path_all_NEG_POS) <- c("Path_desc", "NEG", "POS")
  
img_pep_trans_ann_tb_DE_path %>% 
  count(sign)

## Add level A and B of pathways 
path_all_NEG_POS <-
  path_all_NEG_POS %>% 
  data.frame(LevelA = path_link$LevelA[match(path_all_NEG_POS$Path_desc, path_link$Path_desc)])

path_all_NEG_POS <-
  path_all_NEG_POS %>% 
  data.frame(LevelB = path_link$LevelB[match(path_all_NEG_POS$Path_desc, path_link$Path_desc)])

## Rename
path_DE <-
  path_all_NEG_POS %>% 
  filter (!is.na(Path_desc)) 

head(path_DE)
rownames(path_DE) <- path_DE$Path_desc
path_DE <-path_DE[ ,-1]
head(path_DE)

## Select pathway groups (levelB) to show
selected <-
  c("B_Carbohydrate metabolism",
    "B_Energy metabolism",
    "B_Lipid metabolism",
    "B_Amino acid metabolism",
    "B_Metabolism of other amino acids",
    "B_Glycan biosynthesis and metabolism",
    "B_Metabolism of cofactors and vitamins",
    "B_Biosynthesis of other secondary metabolites",
    "B_Xenobiotics biodegradation and metabolism",
    "B_Membrane transport",
    "B_Signal transduction",
    "B_Signaling molecules and interaction",
    "B_Cellular community - prokaryotes",
    "B_Cell motility")

## Bubble Plot
plot_long <-
path_DE %>% 
  as_tibble(rownames= "DE") %>% 
  filter(LevelB %in% selected) %>% 
  mutate(
    across(where(is.numeric), .fns = replace_na, replace = 0)) %>% 
  dplyr::select(-c("LevelA") ) %>% 
  pivot_longer(-c(DE, LevelB), names_to="species", values_to = "freq")


plot_long %>% 
  arrange(desc(LevelB)) %>% 
  mutate(DE = forcats::fct_inorder(DE)) %>% 
  ggplot( aes(x= species, y= DE, colour=LevelB,
              size = ifelse(freq== 0, NA, freq))) +  
  geom_point() +
  scale_size_continuous(range=c(0, 8))+
  labs(size = "RelAb") +
  theme_light()

