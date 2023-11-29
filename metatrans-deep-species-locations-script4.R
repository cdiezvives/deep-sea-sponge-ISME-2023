#' Data analysis scripts for manuscript XXX
#' 
#' Metatranscriptome of deep-sea sponge species and locations
#' 

library(vegan)
library(pheatmap)
library(RColorBrewer)
library (tidyverse)
library(MASS)
library(ggrepel)
library(ggplot2)
library(readxl)

# Prepare files ----
## Load tables of modules with reactions
table1 <-
  read_excel('data-analysis-ready/metatrans-deep-species-locations/Tables_reactions.xlsx', sheet = "Energy")
table2 <-
  read_excel('data-analysis-ready/metatrans-deep-species-locations/Tables_reactions.xlsx', sheet = "Cfixation")
table3 <-
  read_excel('data-analysis-ready/metatrans-deep-species-locations/Tables_reactions.xlsx', sheet = "Co&Vit")
table4 <-
  read_excel('data-analysis-ready/metatrans-deep-species-locations/Tables_reactions.xlsx', sheet = "OtherC")

## Load previous files
load("data-analysis-ready/metatrans-deep-species-locations/DATA_Download_IMG_ready.RData")
load("data-analysis-ready/metatrans-deep-species-locations/DATA_Peptides_fromIMG_withAnnotation_realcount&tpm_means_Microbes.RData")

## Rename to shorten the name
relTPM_ann <- img_pep_trans_ann_reltpm_micros
colnames(relTPM_ann)

## Keep entries of pmoABC-amoABC that are only Methyloccocales
noPmo <-  
  relTPM_ann %>% 
  dplyr::filter(ko_id %in% c("K10944", "K10945", "K10946")) %>% ## pmoA-amoA
  dplyr::filter(!grepl(x = lineage, pattern = "^Bacteria;Proteobacteria;Gammaproteobacteria;Methylococcales")) %>% 
  pull(peps)

## remove alcohol dehydrogenase (K00114) because it has double annotation with methanol dehydrogenase k14028
relTPM_ann %>% 
  filter (ko_id == "K14028") # This one is not included in our lists

many <-
  relTPM_ann %>% 
  filter (ko_id == "K00114") %>% ## This one is 1906 times
  pull(peps)

## Check when k00144 has a second annotation as K14028
misssing <-
  img_ko %>% 
  filter (KO_ID == "K14028") %>% # K14028 is double annotated in 21 instances
  pull(pep_id)

nomdh <- setdiff(many, misssing) # the ones to remove (1885)

## Remove those two sets of entries in our annotated file
relTPM_ann <- 
  relTPM_ann %>%
  dplyr::filter(!peps %in% nomdh ) %>%  # removes 1885
  dplyr::filter(!peps %in% noPmo ) # removes 89

## Remove columns of no interest now and sample DR10_477
relTPM_ann_less <- 
  as_tibble(relTPM_ann) %>% 
  dplyr::select(-Povi_DR10_477) %>% 
  dplyr::select(peps, genes, ko_id, KEGG_description,  contains(c("_DR")))

## Join info from the modules. There are many-to-many matches, keep them all
relTPM_ann_less_join <-
  relTPM_ann_less %>% 
  left_join(module_link,
            by = c('ko_id' = 'KO'), relationship = "many-to-many")  

# FIGURE 5 ----
## Level C Modules
join_char_sum_num <- function(data){
  if (class(data) == 'character') {
    paste(data, collapse = '|')
  } else if (class(data) == 'numeric') {
    sum(data,  na.rm = TRUE)
  }
}

relTPM_ann_less_join_LevelC_long <-
  relTPM_ann_less_join %>% 
  dplyr::select(c(LevelC, contains("DR"))) %>% # , KO eliminar o no
  group_by(LevelC) %>% 
  summarise(
    across(.cols = everything(),
           .fns = ~join_char_sum_num(.))) %>% 
  pivot_longer(cols= -LevelC,
               names_to = 'Sample',  	
               values_to = 'value') %>%
  mutate(species = substring(Sample, 1,9)) 

top_32levCspe <-
  relTPM_ann_less_join_LevelC_long %>% 
  group_by(LevelC) %>% 
  summarise(avg_sumTPM = mean(value)) %>%  
  arrange(desc(avg_sumTPM)) %>% 
  filter(!LevelC == "NA") %>% 
  slice(1:32) %>% 
  pull(LevelC)

## Plot as a mean of all samples
relTPM_ann_less_join_LevelC_long %>% 
  filter(LevelC %in% top_32levCspe) %>% 
  mutate(LevelC = factor(LevelC,
                         levels= rev(top_32levCspe))) %>% 
  group_by(LevelC) %>% 
  summarise(mean= mean(value), sd= sd(value)) %>% 
  arrange(desc(mean)) %>% 
  ggplot(aes(x = log(mean), y = LevelC)) + 
  geom_bar(stat="identity", linewidth=0.6) +
  geom_errorbar(aes(y=LevelC, xmin = log(mean-sd), xmax = log(mean+sd), width=0.3)) +
  theme_light()

## Plot the top Modules of each Level C
## Group by Module
relTPM_ann_less_join_summbymap_long <-   ## esto vuelve a hacer lo mismo que arrba, solo que lo deja en long
  relTPM_ann_less_join %>% 
  pivot_longer(cols = starts_with(c('Gbar', 'Gpac', 'Povi')),
               names_to = 'Sample',  	 # los headers van a una columna llamada Sample (nuevo nombre)
               values_to = 'value') %>% 
  mutate(species = substring(Sample, 1,9)) %>% 
  group_by(Module_desc, species, Sample) %>% ## lo agrupa por modulo y lo sumara
  summarise(sumTPM= sum(value)) 

## Select top 32 Modules by abundance
top_32modspe <-
  relTPM_ann_less_join_summbymap_long %>% 
  group_by(Module_desc) %>% 
  summarise(avg_sum_tpm = mean(sumTPM)) %>%  
  arrange(desc(avg_sum_tpm)) %>% 
  filter(!is.na(Module_desc)) %>% 
  slice(1:32) %>% 
  pull(Module_desc)

top_32modspe_tb <- as_tibble(top_32modspe)

relation <-
  top_32modspe_tb %>% 
  data.frame(mod_C = module_link$LevelC[match(top_32modspe_tb$value, module_link$Module_desc)])

relation2 <-
  as_tibble(relation %>% 
              mutate(value= factor(value,
                                   levels= (top_32modspe))) %>% 
              mutate(mod_C = factor(mod_C,
                                    levels= (top_32levCspe))))  

## order as the Level C
relTPM_ann_less_join_summbymap_long_top32_withC <-
  relTPM_ann_less_join_summbymap_long %>% 
  filter(Module_desc %in% top_32modspe)  %>% 
  left_join(relation, by = c("Module_desc" = "value" ))

## Plot removing two partial cycles
relTPM_ann_less_join_summbymap_long_top32_withC %>% 
  filter(!Module_desc == "M00011_Citrate cycle, second carbon oxidation, 2-oxoglutarate => oxaloacetate") %>% 
  filter(!Module_desc == "M00620_Incomplete reductive citrate cycle, acetyl-CoA => oxoglutarate") %>% 
  group_by(Module_desc, mod_C) %>% 
  summarise(mean= mean(sumTPM), sd= sd(sumTPM), .groups = "drop") %>% 
  mutate(mod_C = factor(mod_C,
                        levels= rev(top_32levCspe))) %>%   
  arrange(mod_C) %>% 
  mutate(Module_desc = forcats::fct_inorder(Module_desc)) %>% 
  ggplot(aes(x = log(mean), y = Module_desc, fill = mod_C)) + 
    geom_bar(stat="identity", size=0.6) +
    geom_errorbar(aes(y=Module_desc, xmin = log(mean-sd), xmax = log(mean+sd), width=0.3)) +
    theme_light()

# FIGURE 6 ----
## Mean expression of selected Modules (by Barplots)
relTPM_ann_mean <- as_tibble(img_pep_trans_ann_reltpm_means_micros) # rename

## Sum up peptides with same KO annotation
relTPM_ann_mean_less <-
  relTPM_ann_mean %>% 
  dplyr::filter(!peps %in% nomdh ) %>%  # remove 1885
  dplyr::filter(!peps %in% noPmo ) %>% # remove 89
  filter(!KEGG_description == "NA") %>% 
  filter(!KEGG_description == "NA_NA") %>% 
  dplyr::select(ko_id,KEGG_description, contains("DR")) 

relTPM_ann_mean_less_byKO <-
  relTPM_ann_mean_less %>% 
  group_by(ko_id) %>% 
  summarise(
    across(.cols = everything(),
           .fns = ~join_char_sum_num(.)))

## Add the values to the table
table1b <-
  left_join(table1 , relTPM_ann_mean_less_byKO, by="ko_id")

## Now we add the KOs within same module and reaction
table1b_byReac <-
  table1b %>% 
  dplyr::select(ko_id, Mod_reaction, contains("DR")) %>% 
  group_by(Mod_reaction) %>% 
  summarise(
    across(.cols = everything(),
           .fns = ~join_char_sum_num(.)))
orden <- 
  table1 %>% 
  distinct(Mod_reaction) %>% 
  pull(Mod_reaction)

## Heatmap
infotab <- 
  table1 %>% ## from the original table1
  dplyr::select(Mod_reaction, Module_desc, LevelB) %>% 
  distinct() %>% 
  tibble::column_to_rownames(var="Mod_reaction") %>%  
  as.data.frame() 

table1b_byReac %>%  
  dplyr::filter(Mod_reaction != "M0000a_R1") %>% 
  mutate_at(.vars=3:8, .fun= ~log(.+1))  %>% 
  dplyr::select(-ko_id) %>% 
  mutate(Mod_reaction = factor(Mod_reaction,
                               levels = orden)) %>% 
  arrange(Mod_reaction) %>% 
  tibble::column_to_rownames(var="Mod_reaction") %>%  
  as.data.frame() %>% 
  t() %>% 
  pheatmap( cluster_rows = F, cluster_cols = F,
            fontsize_row = 8, 
            annotation_col = infotab[, c("Module_desc","LevelB"), drop=F],
            show_rownames = T,
            show_colnames = T)

# FIGURE 7 ----
## Example with Geodia pachydermata
relTPM_ann_mean_gpac <-
  relTPM_ann_mean %>% 
  dplyr::select("peps", "ko_id" , "KEGG_description", "lineage",  "phyla", "GpacDR15Mean")

## Plot phylum level but class level for Protebacteria
relTPM_ann_mean_gpac <- 
  relTPM_ann_mean_gpac %>% 
  mutate(lineage=sub("Bacteria;Proteobacteria;","Bacteria;Proteobacteria_",lineage))

relTPM_ann_mean_gpac <- 
  relTPM_ann_mean_gpac %>% 
  mutate(phyla_class = strsplit(lineage, ';') %>% # this separates each line in a list
           sapply(function(split) paste(split[1:2], collapse = ';'))) 

## Sum up the KOs with same taxonomy
relTPM_ann_mean_gpac_byKO <-
  relTPM_ann_mean_gpac %>% 
  filter(!KEGG_description == "NA") %>% 
  filter(!KEGG_description == "NA_NA") %>% 
  dplyr::select(ko_id, phyla_class, contains("DR")) %>% 
  group_by(  ko_id, phyla_class) %>% 
  summarise(
    across(.cols = everything(),
           .fns = ~join_char_sum_num(.)))

## simplify the Module table
table1_simple <- 
  table1 %>% 
  dplyr::select("ko_id", "Desc_and_gene") 

table1_simple_tax2 <-
  left_join(table1_simple , relTPM_ann_mean_gpac_byKO , by="ko_id" ,relationship = "many-to-many") %>% 
  filter(!is.na(phyla_class))

## check that there are not more than 2
table1_simple_tax2 %>% 
  count(Desc_and_gene, phyla_class) %>% filter(n == 2) %>% print(n =Inf)

table1_simple_tax2_wider_df <- 
  table1_simple_tax2 %>% 
  dplyr::select("Desc_and_gene", "phyla_class", "GpacDR15Mean") %>%  ## elegimos name sort en luhar de KO, para el plot
  pivot_wider(values_from = GpacDR15Mean, names_from = phyla_class) %>% 
  data.frame(row.names = 'Desc_and_gene')

table1_simple_tax2_wider_df[is.na(table1_simple_tax2_wider_df)] <- 0

## Remove reactions with no expression
table1_simple_tax2_wider_df <-
  table1_simple_tax2_wider_df %>%
  filter(rowSums(table1_simple_tax2_wider_df) > 0)

## transpose
table1_simple_tax2_wider_df_inv <- t(table1_simple_tax2_wider_df)
orderO <- order(rownames(table1_simple_tax2_wider_df_inv))

RdBu_ramp <-
  colorRampPalette(colors = RColorBrewer::brewer.pal(11, 'RdBu'))

pheatmap(log10(table1_simple_tax2_wider_df_inv[orderO,]+1), 
         color = c(rev(RdBu_ramp(99))),
         cluster_rows = F, cluster_cols = F, fontsize_row = 5)   

# FIGURE S4 ----
## P. ang and H. cae are not included in the plots here because the data belongs to another research team
## Split the peptides corresponding to each species
Gbar <- 
  relTPM_ann_mean %>% 
  filter(GbarDR15Mean >0) %>% 
  filter(!is.na(ko_id)) %>% 
  dplyr::select(peps, ko_id)

Gpac <- 
  relTPM_ann_mean %>% 
  filter(GpacDR15Mean >0) %>% 
  filter(!is.na(ko_id)) %>% 
  dplyr::select(peps, ko_id)

Povi <- 
  relTPM_ann_mean %>% 
  filter(PoviDR15Mean >0) %>% 
  filter(!is.na(ko_id)) %>% 
  dplyr::select(peps, ko_id)

## Select desired Module groups  
ko_pickmod <-
  module_link %>% # loaded before
  filter(LevelB %in% c("B_Carbohydrate metabolism", "B_Energy metabolism", "B_Metabolism of cofactors and vitamins","B_Xenobiotics biodegradation")) %>% 
  pull(KO)

Gbar_pickmod <- Gbar %>%  filter(ko_id %in% ko_pickmod)
Gpac_pickmod <- Gpac %>%  filter(ko_id %in% ko_pickmod)
Povi_pickmod <- Povi %>%  filter(ko_id %in% ko_pickmod)

## Save that lists and upload them in KEGG Mapper (https://www.genome.jp/kegg/mapper/reconstruct.html)
write_tsv(Gbar_pickmod , file="data-analysis-ready/metatrans-deep-species-locations/IMG-mapper/List_Gbar_KO_mapper_pickmod.txt", col_names = F)
write_tsv(Gpac_pickmod , file="data-analysis-ready/metatrans-deep-species-locations/IMG-mapper/List_Gpac_KO_mapper_pickmod.txt", col_names = F)
write_tsv(Povi_pickmod , file="data-analysis-ready/metatrans-deep-species-locations/IMG-mapper/List_Povi_KO_mapper_pickmod.txt", col_names = F)

## Then download the tables from KEGG Mapper and prepare them in one as the following example:
tablemap <- read_delim("data-analysis-ready/metatrans-deep-species-locations/IMG-mapper/File_Table_KeegMapper.csv", delim = ",", col_names=T)
ord_spe <- c("Gbar", "Gpac", "Povi")

## Select one line (related with level C) for each plot
tablemap_wd <- 
  tablemap %>% 
  filter(LevelC %in% c("Nitrogen metabolism", "Methane metabolism", "Photosynthesis", "Sulfur metabolism")) %>% # for "Energy metabolism" we selected some of the modules
  #filter(LevelC %in% "Carbon fixation") %>% 
  #filter(LevelC == "Central carbohydrate metabolism") %>%  
  #filter(LevelC == "Other carbohydrate metabolism") %>%  
  #filter(LevelC == "Cofactor and vitamin metabolism") %>% 
  #filter(LevelB == "Xenobiotics biodegradation") %>%  
  #filter(LevelB == "Lipid metabolism") %>%  
  dplyr::select(Group, Module_desc, Percent) %>% 
  tidyr::pivot_wider(names_from = Group, values_from = Percent) %>% 
  replace(is.na(.), 0)

tablemap_wd_df <- 
  tablemap_wd %>% 
  data.frame(row.names = "Module_desc")

RdBu_ramp <- colorRampPalette(colors = RColorBrewer::brewer.pal(9, 'YlGnBu'))

## plot heatmap
pheatmap(tablemap_wd_df [ ,ord_spe] , 
         color = c((RdBu_ramp(99))),
         border_color =NA,
         cluster_rows = F, 
         cluster_cols = F,
         fontsize_row = 5)



