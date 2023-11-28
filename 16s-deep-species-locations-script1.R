#' Data analysis scripts for manuscript XXX
#' 
#' 16S comparison among deep-sea sponge species and locations
#' 

library(vegan)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(MASS)
library("ggrepel")
library(ggplot2)
library(eulerr)
library(edgeR)
library(readxl)
library(ComplexUpset) #https://krassowski.github.io/complex-upset/articles/Examples_Python.html#3-adjusting-intersection-size-

####load("data-analysis-ready/16s-deep-species-locations/Cantabric18_DataForAnalyses.RData")

# Import data from mothur ----
counts <-
  readr::read_delim(
    file = 'final_sponge_asv.shared',
    delim = '\t',
    col_types = readr::cols('Group' = 'c', .default = 'i')
  )


taxa <-
  read.table(
    file = 'sponge.trim.contigs.pcr.good.unique.good.filter.unique.precluster.abund.pick.pick.asv.ASV.cons.taxonomy',
    header = T,
    row.names = 1,
    stringsAsFactors = F
  )

Infotable_tb <-
  read_excel('data-analysis-ready/16s-deep-species-locations/Infotable.xlsx', sheet = "Infotable")

Infotable <-
  data.frame(Infotable_tb , row.names = "ID")

# sort names of counts
counts <- counts[ , sort(colnames(counts))]
colnames(counts)

identical(colnames(counts), rownames(Infotable))

## Relative abundance and presence/absence transformation
counts_RA <-
  t(t(counts) / colSums(counts) * 100)

colSums(counts_RA)

counts_PA <- counts
counts_PA[ counts_PA > 0 ] <- 1

# Rarefaction ----
set.seed(1010)

counts_rare <- rrarefy(t(counts), sample = min(colSums(counts)))
counts_rare <- as.data.frame(t(counts_rare))
colSums(counts_rare)

## Remove empty ASV after rarefaction
counts_rare <- counts_rare[ rowSums(counts_rare) > 0, ]
counts_rare_RA <- as.data.frame(t(t (counts_rare) / colSums(counts_rare) *100))

# Dendrogram & Barplot ----
# Dendrogram at ASV level
counts_bc <- vegdist(log2(t(counts_rare_RA)+1), method="bray")  # With log 2 ## USO ESTE

hc <- hclust(counts_bc, method = "complete") 
plot(rev(hc),  hang = -1)

## save the order for the barplot
orderD <- hc$labels[hc$order]

## barplot at level 2 to simplify colors
level2.counts.RelAb = apply(counts_RA, 2, function(x) tapply(x, taxa$level2taxa, sum)) 
level2.counts.RelAb_imp <- level2.counts.RelAb[rowSums(level2.counts.RelAb)>1,]
other <- level2.counts.RelAb[rowSums(level2.counts.RelAb)<= 1, , drop=F]
level2.counts.RelAb_imp.all <- rbind(level2.counts.RelAb_imp, "Bacteria;Others" = colSums(other)  )
dim(level2.counts.RelAb_imp.all)

yep <- colorRampPalette(brewer.pal(12, "Paired"))(22) 
par(mar = c(3,5,1,9), xpd = T)       

barplot(as.matrix(level2.counts.RelAb_imp.all[, rev(orderD)]),
        col=yep, las = 2, cex.axis = 1, cex.names=0.7, horiz=F, space=0) #, space=0)

legend(x=20, y= 100, cex = 0.5, legend = rownames(level2.counts.RelAb_imp.all), pch = 19, col=yep)
dev.off()

# Ordination ----
## PCoA
mod <- betadisper(counts_bc, Infotable$codif, 
                  type = c("median"), bias.adjust = FALSE,sqrt.dist = FALSE, add = FALSE)
boxplot(mod)
plot(mod, label = T)
mod$eig ## not the same

anova(mod)
permutest(mod, pairwise = TRUE, permutations = 99)

# Shared and specific ASVs (Venn diagram and Upset plots) ----
## Keep abundant ASVs
otuFilter = rowSums(counts_rare_RA > 0.01) >= 1
table(otuFilter)
counts_rare_RA_fil <- counts_rare_RA[ otuFilter ,  ]
counts_rare_RA_fil_PA <- as.data.frame(ifelse(counts_rare_RA_fil>0,1,0))
colSums(counts_rare_RA_fil_PA)
colSums(counts_rare_RA_fil)

## Prepare data set by group
counts_PA <- counts_rare_RA_fil_PA  ## rename
counts_PA_tb <- 
  as_tibble(counts_PA, rownames = "asv")

counts_PA_long <-
  counts_PA_tb %>%
  pivot_longer(cols=-c(asv), names_to = "sample", values_to = "Abund" ) 

counts_PA_long_group <- 
  counts_PA_long %>%
  data.frame(group = Infotable$codif[match(row.names(Infotable), counts_PA_long$sample)]) %>%
  group_by(group, asv) %>%
  summarise(sum = sum (Abund))

counts_PA_group <- 
  counts_PA_long_group %>%
  pivot_wider(id_cols = asv, names_from = group, values_from = sum)

counts_PA_group_db <-
  counts_PA_group %>% 
  data.frame(row.names="asv")

head(counts_PA_group_db)

## Eulerr diagrams 
PA_groups <- as.data.frame(ifelse(counts_PA_group_db > 0,1,0)) 
table(rowSums(PA_groups)>=6)

gbar15 <- as.logical(counts_PA_group_db[,1])
gpac15 <- as.logical(counts_PA_group_db[,2])
povi10 <- as.logical(counts_PA_group_db[,3])
povi15 <- as.logical(counts_PA_group_db[,4])
povi4 <- as.logical(counts_PA_group_db[,5])
povi9 <- as.logical(counts_PA_group_db[,6])

logi_final<-cbind(gbar15, gpac15, povi15)
fit2 <- euler(logi_final)

plot(fit2,
     quantities = TRUE,
     lty = 1:3,
     labels = list(font = 4))

## Upset plots
upset(
  PA_groups,
  c(
    "Gbar_DR15",
    "Gpac_DR15",
    "Povi_DR15",
    "Povi_DR10",
    "Povi_DR9",
    "Povi_DR4"
  ),
  mode = 'inclusive_intersection',
  set_sizes = upset_set_size() +
    theme(axis.text.x = element_text(angle = 90)),
  width_ratio = 0.2,
  min_degree = 2,
  max_degree = 6,
  sort_sets = FALSE,
  base_annotations = list('Size' = (
    intersection_size(
      counts = TRUE,
      mode = 'inclusive_intersection',
      text_colors = c(on_background =
                        'black', on_bar = 'black'),
      text_mapping = aes(label =
                           paste0(
                             round(
                               !!get_size_mode('inclusive_intersection') / !!get_size_mode('inclusive_union') * 100
                             ), '%'
                           ))
    )
  ))
)

