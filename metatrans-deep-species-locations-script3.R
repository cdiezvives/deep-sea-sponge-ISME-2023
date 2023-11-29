#' Data analysis scripts for manuscript XXX
#' 
#' Metatranscriptome of deep-sea sponge species and locations
#' 

library(vegan)
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(MASS)
library(ggrepel)
library(ggplot2)
library(eulerr)
library(ComplexUpset) #https://krassowski.github.io/complex-upset/articles/Examples_Python.html#3-adjusting-intersection-size-
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
#library(pairwiseAdonis)

# Load previous files
load('data-analysis-ready/metatrans-deep-species-locations/DATA_Peptides_fromIMG_withAnnotation_allTrans.RData')
load('data-analysis-ready/metatrans-deep-species-locations/DATA_RNA_raw_data_allTrans.RData')

# Load infotable of samples ----
Infotable_tb <- read_excel("data-analysis-ready/metatrans-deep-species-locations/Infotable.xlsx")

Infotable_df <- 
  Infotable_tb %>% 
  data.frame(row.names = "Sample")
  

# FIGURE 2: Phylogeny and Beta Diversity ----
## Select only microbes 
tmm_micros_tb <- 
  tmm_tb %>% 
  filter(trans_id %in% img_tax_micros_trans)

counts_micros_tb <- 
  counts_tb %>% 
  filter(trans_id %in% img_tax_micros_trans)

tmm_micros  <- data.frame(tmm_micros_tb, row.names = "trans_id")
counts_micros  <- data.frame(counts_micros_tb, row.names = "trans_id")
  
## Remove sample "Povi_DR10_477"
tmm_micros_17 <- tmm_micros[ , -9]

Infotable_df_17 <- Infotable_df[ !rownames(Infotable_df) =='DR10_477', ]
identical(colnames(tmm_micros_17), Infotable_df_17$ID)

## Calculate bray-Curtis dissimilarity and plot ordination of samples
bc17 <- vegdist(log2(t(tmm_micros_17)+1), method="bray")  # With log 2
mod17 <-
  betadisper(d = bc17, 
             group = Infotable_df_17$codif,
             type = c("median"),
             bias.adjust = FALSE,
             sqrt.dist = FALSE,
             add = FALSE
  )

plot(mod17) 
boxplot(mod17)
anova(mod17)  ### testing the dispersion
permutest(mod17, pairwise = TRUE, permutations = 99)

## Adonis between groups
adonis2(bc17 ~ codif , data = Infotable_df_17)
source('functions/adonis.pairwise.custom.r')
adonis.pairwise(bc17 ~ codif, data = Infotable_df_17, compare = 'codif')

## Dendrogram 
hc <- hclust(bc17, method = "complete") 
plot(hc,  hang = -1)

## Barplot at phylum level
## Get taxonomy of the longest peptide
img_pep_trans_ann_order <-
  img_pep_trans_ann %>%
  mutate(phyla2=paste(kingdom, phyla, sep=";")) %>%
  filter(!is.na(lineage)) %>%
  arrange(genes, length_pep)

img_pep_trans_ann_order_unique <-
  img_pep_trans_ann_order %>%
  distinct(trans_id, .keep_all=T)

## Add phyla to transcripts counts
counts_phyla <-
  counts_tb %>%
  data.frame(phyla=img_pep_trans_ann_order_unique$phyla2[match( counts_tb$trans_id, img_pep_trans_ann_order_unique$trans_id)]) %>%
  data.frame(kingdom=img_pep_trans_ann_order_unique$kingdom[match( counts_tb$trans_id, img_pep_trans_ann_order_unique$trans_id)])

head(counts_phyla)

counts_phyla <-
  counts_phyla %>%
  filter(phyla!="<NA>") %>%
  filter(phyla!="NA;NA")

## Function to add values of the same group
join_char_sum_num <- function(data){
  if (class(data) == 'character') {
    paste(data, collapse = '|')
  } else if (class(data) == 'numeric') {
    sum(data,  na.rm = TRUE)
  }
}

colnames(counts_phyla)

## Sum counts by phya
counts_phyla_byphyla <-
  counts_phyla %>% 
  as_tibble() %>% 
  dplyr::select(phyla, contains("DR") ) %>% 
  group_by(phyla) %>% 
  summarise(
    across(.cols = everything(),
           .fns = ~join_char_sum_num(.)))

counts_phyla_byphyla_df <-
  counts_phyla_byphyla  %>% 
  filter(!is.na(phyla)) %>%  
  data.frame( row.names ="phyla") 

## transform to percentages
sampleTotals <- colSums(counts_phyla_byphyla_df)
counts_phyla_byphyla_df_RA <- t(t(counts_phyla_byphyla_df) / sampleTotals * 100)  # lo pasamos a relatove abundance (las presencia de transc. es lo mismo que en amplicon)

## Group low abundant transcripts into "others" category
counts_phyla_byphyla_df_RA.imp  <- (counts_phyla_byphyla_df_RA)
counts_phyla_byphyla_df_RA.imp <- counts_phyla_byphyla_df_RA.imp[rowSums(counts_phyla_byphyla_df_RA.imp) > 5, ]
other <- counts_phyla_byphyla_df_RA[rowSums(counts_phyla_byphyla_df_RA) <= 5, , drop=F]
counts_phyla_byphyla_df_RA.imp.all <- rbind(counts_phyla_byphyla_df_RA.imp, "Bacteria;Others<5%" = colSums(other)  )

## Use the sample order as shown before in the dendrogram
target <- hc$labels[hc$order]

dim(counts_phyla_byphyla_df_RA.imp.all)
yep <- colorRampPalette(brewer.pal(12, "Paired"))(25) 
par(mar = c(4,2,1,10), xpd = T)  
barplot(counts_phyla_byphyla_df_RA.imp.all[ , (target)], horiz = F, space=0, col = yep,las=2, cex.names = 0.7, cex.axis = 0.7, ylim=c(100,0))
legend(x=20, y= 10, cex = 0.5, legend = rownames(counts_phyla_byphyla_df_RA.imp.all[,]), pch = 19, col=yep)

# FIGURE 3: Eurlerr plots of the transcripts ----
## Filter dataset by > 1 TPM 
otuFilter = rowSums(tmm_micros > 1) >= 1 
table(otuFilter)

tmm_micros_fil <- tmm_micros[ otuFilter ,  ]
counts_micros_fil <- counts_micros[ otuFilter ,  ]

## Transform to presence/absence
tmm_micros_PA_fil <- as.data.frame(ifelse(tmm_micros_fil > 0,1,0))
tmm_micros_PA_fil_tb <- 
  as_tibble(tmm_micros_PA_fil, rownames = "trans")

tmm_micros_PA_fil_long <-
  tmm_micros_PA_fil_tb %>%
  pivot_longer(cols=-c(trans), names_to = "sample", values_to = "Abund" ) 

## Count samples by Species and Location
tmm_micros_PA_fil_long_group <- 
  tmm_micros_PA_fil_long %>%
  mutate(group = substr(sample, 1, 9)) %>% 
  group_by(group, trans) %>%
  summarise(sum = sum (Abund))

tmm_micros_PA_fil_group <- 
  tmm_micros_PA_fil_long_group %>%
  pivot_wider(id_cols = trans, names_from = group, values_from = sum)

tmm_micros_PA_fil_group_db <-
  tmm_micros_PA_fil_group %>% 
  data.frame(row.names="trans")

colSums(tmm_micros_PA_fil_group_db[ rowSums(tmm_micros_PA_fil_group_db)>0 , ])

PA_gropus <- as.data.frame(ifelse(tmm_micros_PA_fil_group_db > 0,1,0)) 
table(rowSums(PA_gropus)>=6)

gbar15 <- as.logical(tmm_micros_PA_fil_group_db[,1])
gpac15 <- as.logical(tmm_micros_PA_fil_group_db[,2])
povi15 <- as.logical(tmm_micros_PA_fil_group_db[,4])

## plot eurlerr graph
logi_final<-cbind(gbar15, gpac15, povi15)

fit2 <- euler(logi_final)

plot(fit2,
     quantities = TRUE,
     lty = 1:3,
     labels = list(font = 4))

## Upset plots of the transcripts
upset(data = PA_gropus, 
      intersect =  c("Gbar_DR15", "Gpac_DR15", "Povi_DR15", "Povi_DR10","Povi_DR9_","Povi_DR4_"),
      mode = 'inclusive_intersection',
      set_sizes = upset_set_size() +
        theme(axis.text.x = element_text(angle=90)),
      base_annotations = list('Size'=(intersection_size(counts=TRUE, mode='inclusive_intersection',
                                                      text_colors = c(on_background='black', on_bar='black'),
                                                      text_mapping=aes(label=paste0(round(
                                                        !!get_size_mode('inclusive_intersection')/!!get_size_mode('inclusive_union') * 100), '%'))))),
      width_ratio = 0.2,
      min_degree=2,
      max_degree=3) # increase to 6 to see all combinations


