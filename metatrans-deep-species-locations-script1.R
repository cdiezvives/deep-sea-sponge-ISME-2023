#' Data analysis scripts for manuscript XXX
#' 
#' Metatranscriptome of deep-sea sponge species and locations
#' 
#' Wrangling IMG and trinity data for data analysis in R.
#' Creates two RData files to use with the next script
#' 

library(vegan)
library(pheatmap)
library(RColorBrewer)
library (tidyverse)
library(MASS)
library(ggrepel)
library(ggplot2)

# Load results from IMG ----
img_tax <- read_delim("data-analysis-ready/metatrans-deep-species-locations/IMG-download/Ga0496847_gene_phylogeny.tsv", delim = "\t", col_names=F)
img_ko <- read_delim('data-analysis-ready/metatrans-deep-species-locations/IMG-download/Ga0496847_ko.tsv', delim = "\t", col_names=F)
img_cog <- read_delim('data-analysis-ready/metatrans-deep-species-locations/IMG-download/Ga0496847_cog.gff', delim = "\t", col_names=F)
img_pfam <- read_delim("data-analysis-ready/metatrans-deep-species-locations/IMG-download/Ga0496847_pfam.gff", delim = "\t", col_names=F)
img_ec <- read_delim("data-analysis-ready/metatrans-deep-species-locations/IMG-download//Ga0496847_ec.tsv", delim = "\t", col_names=F)
img_map <- read_delim('data-analysis-ready/metatrans-deep-species-locations/IMG-download/Ga0496847_contig_names_mapping.tsv', delim = "\t", col_names=F)

## Add headers to files
colnames(img_tax) <- c("pep_id",	"homolog_gene_oid",	"homolog_taxon_oid",	"percent_identity",	"lineage")
colnames(img_ko) <- c("pep_id",	"mg_ko_flag",	"ko_term",	"percent_identity",	"query_start",	"query_end", "subj_start",	"subject_end",	"evalue",	"bit_score",	"align_lenght")
colnames(img_cog) <- c("pep_id",	"emptycolumn", "cog_id",	"dontknow", 	"align_length", "percent_identity",	"query_start",	"query_end", 	"someinfo")
colnames(img_pfam) <- c("pep_id",	"program","pfam_id",	"query_start",	"query_end", "percent_identity",		"subj_start",	"subject_end",	"someinfo")
colnames(img_ec) <- c("pep_id",	"img_ko_flag",	"EC(fromKEGG)",	"percent_identity",	"query_start"	,"query_end",	"subj_start",	"subject_end", "evalue",	"bit_score",	"align_lenght")

# Load extra information ----
cog_map <- read_delim("data-analysis-ready/metatrans-deep-species-locations/mapping-files/File_cog_map.csv", delim = ",", col_names=T)
pfam_map <- read_delim("data-analysis-ready/metatrans-deep-species-locations/mapping-files/File_pfam_map.csv", delim = ",", col_names=T)
ko_desc <- read_delim("data-analysis-ready/metatrans-deep-species-locations/mapping-files/File_ko_desc.csv", delim = ",", col_names=T)
module_link <- read_delim("data-analysis-ready/metatrans-deep-species-locations/mapping-files/File_module_link.csv", delim = ",", col_names=T)
path_link <- read_delim("data-analysis-ready/metatrans-deep-species-locations/mapping-files/File_path_link.csv", delim = ",", col_names=T)

# Wrangle data ----
## Split phyla in img_tax phyla
img_tax_counts <- 
  img_tax %>% 
  mutate(phyla = strsplit(lineage, ';') %>% # this separates each line in a list
           sapply(function(split) paste(split[1], collapse = ';'))) 

## Modify ko_term column
img_ko <- 
  img_ko %>% 
  mutate(KO_ID = substr(ko_term, 4, 9))   ## 833,700

## Separate information from last column
img_pfam2 <- 
  img_pfam %>% 
  separate(someinfo, c("ID", "Name", "fake_percent_id", "alignment_length", "e-value", "model_start", "model_end"), sep= ";") 

## Save data for next step
save(
  img_cog,
  img_ko,
  img_ec,
  img_tax,
  img_pfam2,
  img_map,
  img_tax_counts,
  pfam_map,
  cog_map,
  pfam_map,
  ko_desc,
  module_link,
  path_link,
  file = "data-analysis-ready/metatrans-deep-species-locations/DATA_Download_IMG_ready.RData"
)

# Obtain length of transcripts and peptides ----
## Obtain length of peptides from the amino acid sequence file from IMG output
"data-analysis-ready/metatrans-deep-species-locations/IMG-download/Ga0496847_proteins.faa"

# ```{bash}
# [get the file Ga0496847_proteins.faa, make it single line fasta]
# #!/usr/bin/perl -w
# downloaded from http://www.bioinformatics-made-simple.com
# change multiline fasta to single line fasta
# use strict;
# 
# my $input_fasta=$ARGV[0];
# open(IN,"<$input_fasta") || die ("Error opening $input_fasta $!");
# 
# my $line = <IN>; 
# print $line;
# 
# while ($line = <IN>)
# {
#   chomp $line;
#   if ($line=~m/^>/) { print "\n",$line,"\n"; }
#   else { print $line; }
# }
# 
# print "\n";
# 
# #usage: perl singleline.pl Ga0496847_proteins.faa > Ga0496847_proteins_sl.faa

# awk '/^>/ {print; next; } { seqlen = length($0); print seqlen}' Ga0496847_proteins_sl.faa > Ga0496847_proteins_sl_lenghts.txt ## Print the lenght of the sequences
# cut -d'#' -f 1 Ga0496847_proteins_sl_lenghts.txt > Ga0496847_proteins_sl_lenghts2.txt # keep first column
# awk '{printf "%s%s",$0,NR%2?"\t":RS}' Ga0496847_proteins_sl_lenghts2.txt > Ga0496847_proteins_sl_lenghts2_col.txt # move pair lines to columns
# sed -i -e 's@>@@g' Ga0496847_proteins_sl_lenghts2_col.txt  ## remove symbol ">"
# sed -i -e 's@ @@g' Ga0496847_proteins_sl_lenghts2_col.txt  #remove spaces
# ```

img_pep <- 
  read_delim('data-analysis-ready/metatrans-deep-species-locations/IMG-download/Ga0496847_proteins_sl_lenghts2_col.txt', col_names = F)

## Obtain length of transcripts from nucleotide sequences file from Trinity output
"data-analysis-ready/metatrans-deep-species-locations/Users_data_metatransc/Trinity_18samp.Trinity.fasta"

# ```{bash}
# grep ">" Trinity_18samp.Trinity.fasta > Trinity_18samp.Trinity.names.txt
# cut -d " " -f 1,2  Trinity_18samp.Trinity.names.txt > Trinity_18samp.Trinity.lenght.txt
# sed -i -e 's@>@@g' Trinity_18samp.Trinity.lenght.txt ## remove symbol ">"
# sed -i -e 's@len=@@g' Trinity_18samp.Trinity.lenght.txt ## remove symbol ">"
# ```  

trans_length <-
  read_delim('data-analysis-ready/metatrans-deep-species-locations/trinity-output/Trinity_18samp.Trinity.lenght.txt', col_names = F)

## Add values to peptide table
img_pep <-
  img_pep %>% 
  mutate(genes = substring(X1, 1,17)) 

colnames(img_pep) <- c("peps",	"length_pep",	"genes")

colnames(trans_length) <- c("trans_id",	"length_trans")

trans_length <- 
  trans_length %>% 
  mutate(trans_Kb = length_trans/1000)

## match the transcripts to the peptide file
img_pep_trans <- 
  img_pep %>% 
  data.frame(trans_id = img_map$X1[match(img_pep$genes , img_map$X2)]) %>% 
  as_tibble()

img_pep_trans <-
  img_pep_trans %>% 
  mutate(pep_Kb= length_pep/1000) %>% 
  dplyr::select(peps, genes, length_pep, pep_Kb, trans_id) 

## Add the transcript length
img_pep_trans <-  
  img_pep_trans %>% 
  left_join(trans_length, by = 'trans_id')

# Prepare combined metadata file: Add annotations to each peptide ----
## add KOs     
img_pep_trans_ann <-
  img_pep_trans %>% 
  data.frame(ko_id = img_ko$KO_ID[match(img_pep_trans$peps, img_ko$pep_id)])

head(img_pep_trans_ann)

img_pep_trans_ann <-   # add KO descriptor
  img_pep_trans_ann %>% 
  left_join(ko_desc, by = c( "ko_id"="Knumber"))

img_pep_trans_ann <- # join KO_id and descriptor
  img_pep_trans_ann %>% 
  mutate(KO_descrip =paste(ko_id, KEGG_description, sep = "_"))

## add COGs
img_pep_trans_ann <-
  img_pep_trans_ann %>% 
  data.frame(cog_id = img_cog$cog_id[match(img_pep_trans_ann$peps, img_cog$pep_id)])

img_pep_trans_ann <-
  img_pep_trans_ann %>% 
  left_join(cog_map, by = c( "cog_id" = "COG")) # add COG descriptors

img_pep_trans_ann <-
  img_pep_trans_ann %>% 
  mutate(COG_desc =paste(cog_id,annotation, sep = "_")) # join COG_id and descriptor

## add EC
img_pep_trans_ann <-
  img_pep_trans_ann %>% 
  data.frame(ec_id = img_ec$`EC(fromKEGG)`[match(img_pep_trans_ann$peps, img_ec$pep_id)])

## add Pfam
img_pep_trans_ann <-
  img_pep_trans_ann %>% 
  data.frame(pfam_id = img_pfam$pfam_id[match(img_pep_trans_ann$peps, img_pfam$pep_id)]) 

img_pep_trans_ann <-
  img_pep_trans_ann %>% 
  data.frame(pfam_names = pfam_map$PFAM_desc[match(img_pep_trans_ann$pfam_id, pfam_map$pfam_ID)])

img_pep_trans_ann <-
  img_pep_trans_ann %>% 
  mutate(PFAM_desc =paste(pfam_id,pfam_names, sep = "_")) # join PFAM_id and descriptor

## Add taxonomy from IMG
img_pep_trans_ann <-  
  img_pep_trans_ann %>% 
  data.frame(lineage = img_tax$lineage[match(img_pep_trans_ann$peps, img_tax$pep_id)])

img_pep_trans_ann <- 
  img_pep_trans_ann %>% 
  mutate(kingdom = strsplit(lineage, ';') %>%  # separate kingdom
           sapply(function(split) paste(split[1], collapse = ';'))) 

img_pep_trans_ann <- 
  img_pep_trans_ann %>% 
  mutate(phyla = strsplit(lineage, ';') %>% # separate phyla
           sapply(function(split) paste(split[2], collapse = ';'))) 

head(img_pep_trans_ann)

# Separate features by origin ----
prok <- c("Bacteria", "Archaea")

img_tax_micros <- 
  img_tax %>% 
  mutate(phyla = strsplit(lineage, ';') %>% # this separates each line in a list
           sapply(function(split) paste(split[1], collapse = ';'))) %>% 
  filter(phyla %in% prok) %>% 
  mutate(genes = substring(pep_id, 1,17)) 

img_tax_micros <-
  img_tax_micros %>% 
  left_join(img_map, by = c("genes" = "X2"))  

img_tax_micros_peps <- 
  img_tax_micros %>% 
  distinct(pep_id) %>% 
  pull(pep_id)

img_tax_micros_genes <- 
  img_tax_micros %>% 
  distinct(genes) %>% 
  pull(genes)

img_tax_micros_trans <- 
  img_tax_micros %>% 
  distinct(X1) %>% 
  pull(X1)

# Save data for next step
save(
  img_pep_trans,
  img_pep_trans_ann,
  img_tax_micros_peps,
  img_tax_micros_genes,
  img_tax_micros_trans,
  file = "data-analysis-ready/metatrans-deep-species-locations/DATA_Peptides_fromIMG_withAnnotation_AllTrans.RData"
)
