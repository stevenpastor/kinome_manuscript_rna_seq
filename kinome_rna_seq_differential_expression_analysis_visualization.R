####################

# Differential Response to ALK inhibition Through G2/M-Mediated Kinome Reprogramming in ALK-driven Neuroblastoma
# last updated: 2024-06-05
# Author: Steven Pastor
# contact: pastors@chop.edu

####################


# required libraries:
library(tximport)
library(EnsDb.Hsapiens.v86)
library(DESeq2)
library(dplyr)
library(tidyverse)
library(ggpattern)
library(msigdbr)
library(fgsea)
library(tidyr)
library(tibble)
library(ggplot2)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(reshape2)


## 1. DEA

# sample 453: raw counts from star+rsem:
dir_tmp <- 'counts/453'
sample_files <- grep("*.genes", list.files(dir_tmp), value=T)

# load (10 files) into txi object and edit colnames to remove .genes suffix:
txi.rsem <- tximport(paste(dir_tmp, sample_files, sep='/'), type = "rsem", txIn = FALSE, txOut = FALSE)
new_cols <- gsub(pattern = ".genes", replacement = "", sample_files)
colnames(txi.rsem$counts) <- new_cols

# Preemptively obtain gene symbols in addition to ensembl ids:
ensembl.genes <- rownames(txi.rsem$counts)
ensembl.genes <- gsub("\\..*","",ensembl.genes)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))

# metadata: sample info table created as a file, so refer to it for conditions that will test later (deseq2):
s = read.table('helper_files/453/sample.info', header=T)
s$name <- gsub("\\.genes", "", s$name)
rownames(s) <- s$name
s$condition = factor(s$condition)

# Build DESeq object from counts and metadata:
dds = DESeqDataSetFromMatrix(countData = round(txi.rsem$counts), colData = s, design = ~ condition)
dds = DESeq(dds)

# pairwise results: CZ=0, LR=1, V=2
# i know there is a better way and can actually obtain the condition names
# but I was provided the metadata file: re-create it better for the future
# s/t it is not 0/1/2 and has the actual drugs
cz_vs_lr_rna <- results(dds, lfcThreshold=0, alpha=0.05, contrast=c("condition", "0", "1"))
cz_vs_v_rna <- results(dds, lfcThreshold=0, alpha=0.05, contrast=c("condition", "0", "2"))
lr_vs_v_rna <- results(dds, lfcThreshold=0, alpha=0.05, contrast=c("condition", "1", "2"))

# again, ensembl ids here so get gene names:
the_counts <- data.frame(txi.rsem$counts)
colnames(the_counts) <- gsub('^X', '', gsub('\\.', '-', colnames(the_counts)))
rownames(the_counts) <- gsub("\\..*","",rownames(the_counts))
the_counts$GENEID <- rownames(the_counts)
df_joined <- left_join(the_counts, geneIDs1, by="GENEID")
the_counts_final <- left_join(the_counts, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

# remove NA anywhere:
the_counts_final_tmp <- the_counts_final %>% drop_na()

# checking a couple genes for sanity:
#the_counts_final_tmp[the_counts_final_tmp$SYMBOL == "CTDSPL2" ,]
#the_counts[the_counts$GENEID == "ENSG00000137770",]

# save the counts:
#write.table(cz_vs_lr_rna, 'degs/cz_vs_lr_rna_453.txt', sep='\t', quote=F)
#write.table(cz_vs_v_rna, 'degs/cz_vs_v_rna_453.txt', sep='\t', quote=F)
#write.table(lr_vs_v_rna, 'degs/lr_vs_v_rna_453.txt', sep='\t', quote=F)

# again, not efficient as should ingest all samples at once but
# previously only had 453 and not 1643
# for code continuity, breaking them up here

# repeat the same process for sample 1643:
dir_tmp <- 'counts/1643'
sample_files <- grep("*.genes", list.files(dir_tmp), value=T)

# genes loaded as:
txi.rsem <- tximport(paste(dir_tmp, sample_files, sep='/'), type = "rsem", txIn = FALSE, txOut = FALSE)
new_cols <- gsub(pattern = ".genes", replacement = "", sample_files)
colnames(txi.rsem$counts) <- new_cols

# get gene symbols:
ensembl.genes <- rownames(txi.rsem$counts)
ensembl.genes <- gsub("\\..*","",ensembl.genes)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))

# metadata: sample info table created as a file, so refer to it for conditions that will test later (deseq2):
s = read.table('helper_files/1643/sample.info', header=T)
s$name <- gsub("\\.genes", "", s$name)
rownames(s) <- s$name
s$condition = factor(s$condition)

# Build DESeq object
dds.1643 = DESeqDataSetFromMatrix(countData = round(txi.rsem$counts), colData = s, design = ~ condition)
dds.1643 = DESeq(dds.1643)

# pairwise results: CZ=0, LR=1, V=2
cz_vs_lr_rna <- results(dds.1643, lfcThreshold=0, alpha=0.05, contrast=c("condition", "0", "1"))
cz_vs_v_rna <- results(dds.1643, lfcThreshold=0, alpha=0.05, contrast=c("condition", "0", "2"))
lr_vs_v_rna <- results(dds.1643, lfcThreshold=0, alpha=0.05, contrast=c("condition", "1", "2"))

# get gene names:
the_counts <- data.frame(txi.rsem$counts)
colnames(the_counts) <- gsub('^X', '', gsub('\\.', '-', colnames(the_counts)))
rownames(the_counts) <- gsub("\\..*","",rownames(the_counts))
the_counts$GENEID <- rownames(the_counts)
df_joined <- left_join(the_counts, geneIDs1, by="GENEID")
the_counts_final <- left_join(the_counts, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

the_counts_final_tmp <- the_counts_final %>% drop_na()

# save the counts:
#write.table(cz_vs_lr_rna, 'degs/cz_vs_lr_rna_1643.txt', sep='\t', quote=F)
#write.table(cz_vs_v_rna, 'degs/cz_vs_v_rna_1643.txt', sep='\t', quote=F)
#write.table(lr_vs_v_rna, 'degs/lr_vs_v_rna_1643.txt', sep='\t', quote=F)


## 2. next, obtain degs from above (if necessary) and perform various operations for figures

# for speed, resume here by reading in the 6 files of the degs per sample-condition:
# you will need to edit the results_dir to point to wherever you cloned this repo:
results_dir <- 'degs/'
cz_vs_lr_rna_453 <- read.table(paste0(results_dir, 'cz_vs_lr_rna_453.txt'))
cz_vs_lr_rna_1643 <- read.table(paste0(results_dir, 'cz_vs_lr_rna_1643.txt'))
cz_vs_v_rna_453 <- read.table(paste0(results_dir, 'cz_vs_v_rna_453.txt'))
cz_vs_v_rna_1643 <- read.table(paste0(results_dir, 'cz_vs_v_rna_1643.txt'))
lr_vs_v_rna_453 <- read.table(paste0(results_dir, 'lr_vs_v_rna_453.txt'))
lr_vs_v_rna_1643 <- read.table(paste0(results_dir, 'lr_vs_v_rna_1643.txt'))

# add in gene ids to side and write out:
rownames(cz_vs_lr_rna_453) <- gsub("\\..*","",rownames(cz_vs_lr_rna_453))
cz_vs_lr_rna_453$GENEID <- rownames(cz_vs_lr_rna_453)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = cz_vs_lr_rna_453$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
cz_vs_lr_rna_453_final <- left_join(cz_vs_lr_rna_453, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(cz_vs_v_rna_453) <- gsub("\\..*","",rownames(cz_vs_v_rna_453))
cz_vs_v_rna_453$GENEID <- rownames(cz_vs_v_rna_453)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = cz_vs_v_rna_453$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
cz_vs_v_rna_453_final <- left_join(cz_vs_v_rna_453, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(lr_vs_v_rna_453) <- gsub("\\..*","",rownames(lr_vs_v_rna_453))
lr_vs_v_rna_453$GENEID <- rownames(lr_vs_v_rna_453)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = lr_vs_v_rna_453$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
lr_vs_v_rna_453_final <- left_join(lr_vs_v_rna_453, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(cz_vs_lr_rna_1643) <- gsub("\\..*","",rownames(cz_vs_lr_rna_1643))
cz_vs_lr_rna_1643$GENEID <- rownames(cz_vs_lr_rna_1643)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = cz_vs_lr_rna_1643$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
cz_vs_lr_rna_1643_final <- left_join(cz_vs_lr_rna_1643, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(cz_vs_v_rna_1643) <- gsub("\\..*","",rownames(cz_vs_v_rna_1643))
cz_vs_v_rna_1643$GENEID <- rownames(cz_vs_v_rna_1643)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = cz_vs_v_rna_1643$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
cz_vs_v_rna_1643_final <- left_join(cz_vs_v_rna_1643, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(lr_vs_v_rna_1643) <- gsub("\\..*","",rownames(lr_vs_v_rna_1643))
lr_vs_v_rna_1643$GENEID <- rownames(lr_vs_v_rna_1643)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = lr_vs_v_rna_1643$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
lr_vs_v_rna_1643_final <- left_join(lr_vs_v_rna_1643, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

# filter each by padj < 0.05:
cz_vs_lr_rna_453_filt <- as.data.frame(cz_vs_lr_rna_453) %>% dplyr::filter(padj < 0.05)
cz_vs_v_rna_453_filt <- as.data.frame(cz_vs_v_rna_453) %>% dplyr::filter(padj < 0.05)
lr_vs_v_rna_453_filt <- as.data.frame(lr_vs_v_rna_453) %>% dplyr::filter(padj < 0.05)
cz_vs_lr_rna_1643_filt <- as.data.frame(cz_vs_lr_rna_1643) %>% dplyr::filter(padj < 0.05)
cz_vs_v_rna_1643_filt <- as.data.frame(cz_vs_v_rna_1643) %>% dplyr::filter(padj < 0.05)
lr_vs_v_rna_1643_filt <- as.data.frame(lr_vs_v_rna_1643) %>% dplyr::filter(padj < 0.05)

# again, add in gene ids:
rownames(cz_vs_lr_rna_453_filt) <- gsub("\\..*","",rownames(cz_vs_lr_rna_453_filt))
cz_vs_lr_rna_453_filt$GENEID <- rownames(cz_vs_lr_rna_453_filt)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = cz_vs_lr_rna_453_filt$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
cz_vs_lr_rna_453_filt_final <- left_join(cz_vs_lr_rna_453_filt, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(cz_vs_v_rna_453_filt) <- gsub("\\..*","",rownames(cz_vs_v_rna_453_filt))
cz_vs_v_rna_453_filt$GENEID <- rownames(cz_vs_v_rna_453_filt)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = cz_vs_v_rna_453_filt$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
cz_vs_v_rna_453_filt_final <- left_join(cz_vs_v_rna_453_filt, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(lr_vs_v_rna_453_filt) <- gsub("\\..*","",rownames(lr_vs_v_rna_453_filt))
lr_vs_v_rna_453_filt$GENEID <- rownames(lr_vs_v_rna_453_filt)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = lr_vs_v_rna_453_filt$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
lr_vs_v_rna_453_filt_final <- left_join(lr_vs_v_rna_453_filt, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(cz_vs_lr_rna_1643_filt) <- gsub("\\..*","",rownames(cz_vs_lr_rna_1643_filt))
cz_vs_lr_rna_1643_filt$GENEID <- rownames(cz_vs_lr_rna_1643_filt)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = cz_vs_lr_rna_1643_filt$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
cz_vs_lr_rna_1643_filt_final <- left_join(cz_vs_lr_rna_1643_filt, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(cz_vs_v_rna_1643_filt) <- gsub("\\..*","",rownames(cz_vs_v_rna_1643_filt))
cz_vs_v_rna_1643_filt$GENEID <- rownames(cz_vs_v_rna_1643_filt)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = cz_vs_v_rna_1643_filt$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
cz_vs_v_rna_1643_filt_final <- left_join(cz_vs_v_rna_1643_filt, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(lr_vs_v_rna_1643_filt) <- gsub("\\..*","",rownames(lr_vs_v_rna_1643_filt))
lr_vs_v_rna_1643_filt$GENEID <- rownames(lr_vs_v_rna_1643_filt)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = lr_vs_v_rna_1643_filt$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
lr_vs_v_rna_1643_filt_final <- left_join(lr_vs_v_rna_1643_filt, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

# from padj<0.05 filtered degs, remove NA
cz_vs_lr_rna_453_filt_final <- cz_vs_lr_rna_453_filt_final %>% drop_na(SYMBOL)
cz_vs_v_rna_453_filt_final <- cz_vs_v_rna_453_filt_final %>% drop_na(SYMBOL)
lr_vs_v_rna_453_filt_final <- lr_vs_v_rna_453_filt_final %>% drop_na(SYMBOL)
cz_vs_lr_rna_1643_filt_final <- cz_vs_lr_rna_1643_filt_final %>% drop_na(SYMBOL)
cz_vs_v_rna_1643_filt_final <- cz_vs_v_rna_1643_filt_final %>% drop_na(SYMBOL)
lr_vs_v_rna_1643_filt_final <- lr_vs_v_rna_1643_filt_final %>% drop_na(SYMBOL)

# # number of de genes at padj < 0.05:
# dim(cz_vs_lr_rna_453_filt_final)[1]
# dim(cz_vs_v_rna_453_filt_final)[1]
# dim(lr_vs_v_rna_453_filt_final)[1]
# dim(cz_vs_lr_rna_1643_filt_final)[1]
# dim(cz_vs_v_rna_1643_filt_final)[1]
# dim(lr_vs_v_rna_1643_filt_final)[1]

# want to see the shared degs between the conditions and same sample:
# join them by external_gene_name; first, add contrast to each:
lr_vs_v_rna_453_filt_final$contrast <- "lr_vs_v"
cz_vs_v_rna_453_filt_final$contrast <- "cz_vs_v"
lrcz_453_joined <- inner_join(lr_vs_v_rna_453_filt_final, cz_vs_v_rna_453_filt_final, "SYMBOL")
dim(lrcz_453_joined)[1] # 133

lr_vs_v_rna_1643_filt_final$contrast <- "lr_vs_v"
cz_vs_v_rna_1643_filt_final$contrast <- "cz_vs_v"
# there will be one instance of multi-matching:
lrcz_1643_joined <- inner_join(lr_vs_v_rna_1643_filt_final, cz_vs_v_rna_1643_filt_final, "SYMBOL")
dim(lrcz_1643_joined) # 1217

# separate the columns into the single drugs now for downstream processes:
cz_tmp_453 <- lrcz_453_joined[, c("SYMBOL", "contrast.y", "log2FoldChange.y")]
lr_tmp_453 <- lrcz_453_joined[, c("SYMBOL", "contrast.x", "log2FoldChange.x")]
colnames(cz_tmp_453) <- c("SYMBOL", "contrast", "log2FoldChange")
colnames(lr_tmp_453) <- c("SYMBOL", "contrast", "log2FoldChange")

cz_tmp_1643 <- lrcz_1643_joined[, c("SYMBOL", "contrast.y", "log2FoldChange.y")]
lr_tmp_1643 <- lrcz_1643_joined[, c("SYMBOL", "contrast.x", "log2FoldChange.x")]
colnames(cz_tmp_1643) <- c("SYMBOL", "contrast", "log2FoldChange")
colnames(lr_tmp_1643) <- c("SYMBOL", "contrast", "log2FoldChange")

# row bind them:
lrcz_453_joined_melted_tmp <- rbind(cz_tmp_453, lr_tmp_453)
lrcz_1643_joined_melted_tmp <- rbind(cz_tmp_1643, lr_tmp_1643)

# sort by L2FC CZ:
lrcz_453_joined_melted_tmp <- lrcz_453_joined_melted_tmp[order(lrcz_453_joined_melted_tmp$log2FoldChange),]
lrcz_1643_joined_melted_tmp <- lrcz_1643_joined_melted_tmp[order(lrcz_1643_joined_melted_tmp$log2FoldChange),]

# as an aside/curiosity: shared genes between 453 and 1643:
dim(inner_join(lrcz_1643_joined, lrcz_453_joined, by="SYMBOL"))

################################################################

## Was originally figure 4 part c but changed since the code was created
# thus, this is currently mis-named but keeping naming convention for now
# should be figure 5 now in the actual manuscript
# kinases provided by Smita:
figure_4_c_kinases <- read.table('helper_files/list_of_kinases.txt', header=T)

# check number of kinases:
length(figure_4_c_kinases$Gene_names)

# # check in the 6 sample-conditions:
# dim(cz_vs_v_rna_453_filt_final[cz_vs_v_rna_453_filt_final$SYMBOL %in% figure_4_c_kinases$Gene_names,])
# dim(cz_vs_v_rna_1643_filt_final[cz_vs_v_rna_1643_filt_final$SYMBOL %in% figure_4_c_kinases$Gene_names,])
# dim(lr_vs_v_rna_453_filt_final[lr_vs_v_rna_453_filt_final$SYMBOL %in% figure_4_c_kinases$Gene_names,])
# dim(lr_vs_v_rna_1643_filt_final[lr_vs_v_rna_1643_filt_final$SYMBOL %in% figure_4_c_kinases$Gene_names,])
# dim(cz_vs_lr_rna_453_filt_final[cz_vs_lr_rna_453_filt_final$SYMBOL %in% figure_4_c_kinases$Gene_names,])
# dim(cz_vs_lr_rna_1643_filt_final[cz_vs_lr_rna_1643_filt_final$SYMBOL %in% figure_4_c_kinases$Gene_names,])

# filter rna-seq for these (gene symbols) so can plot them
# manually add contrasts in prettier manner
cz_vs_v_rna_453_plotme <- cz_vs_v_rna_453_filt_final[cz_vs_v_rna_453_filt_final$SYMBOL %in% figure_4_c_kinases$Gene_names,][,c("SYMBOL", "log2FoldChange")]
lr_vs_v_rna_453_plotme <- lr_vs_v_rna_453_filt_final[lr_vs_v_rna_453_filt_final$SYMBOL %in% figure_4_c_kinases$Gene_names,][,c("SYMBOL", "log2FoldChange")]
cz_vs_v_rna_453_plotme$contrast <- "Crizotinib vs Vehicle"
lr_vs_v_rna_453_plotme$contrast <- "Lorlatinib vs Vehicle"
plot_me_453 <- rbind(cz_vs_v_rna_453_plotme, lr_vs_v_rna_453_plotme)

cz_vs_v_rna_1643_plotme <- cz_vs_v_rna_1643_filt_final[cz_vs_v_rna_1643_filt_final$SYMBOL %in% figure_4_c_kinases$Gene_names,][,c("SYMBOL", "log2FoldChange")]
lr_vs_v_rna_1643_plotme <- lr_vs_v_rna_1643_filt_final[lr_vs_v_rna_1643_filt_final$SYMBOL %in% figure_4_c_kinases$Gene_names,][,c("SYMBOL", "log2FoldChange")]
cz_vs_v_rna_1643_plotme$contrast <- "Crizotinib vs Vehicle"
lr_vs_v_rna_1643_plotme$contrast <- "Lorlatinib vs Vehicle"
plot_me_1643 <- rbind(cz_vs_v_rna_1643_plotme, lr_vs_v_rna_1643_plotme)

# what if I put them all into one figure:
cz_vs_v_rna_453_plotme$sample_contrast <- as.factor("COG-N-453 Crizotinib vs Vehicle")
lr_vs_v_rna_453_plotme$sample_contrast <- as.factor("COG-N-453 Lorlatinib vs Vehicle")
cz_vs_v_rna_1643_plotme$sample_contrast <- as.factor("NB1643 Crizotinib vs Vehicle")
lr_vs_v_rna_1643_plotme$sample_contrast <- as.factor("NB1643 Lorlatinib vs Vehicle")

plot_me_both <- rbind(cz_vs_v_rna_453_plotme, lr_vs_v_rna_453_plotme, cz_vs_v_rna_1643_plotme, lr_vs_v_rna_1643_plotme)

# separate figure by sample and force 0 l2fc for missing ones?
tmp_one <- merge(cz_vs_v_rna_453_plotme, lr_vs_v_rna_453_plotme, by="SYMBOL", all=T)
tmp_two <- merge(cz_vs_v_rna_1643_plotme, lr_vs_v_rna_1643_plotme, by="SYMBOL", all=T)
all_outer_joined <- merge(tmp_one, tmp_two, by="SYMBOL", all=T)
all_outer_joined <- all_outer_joined[, c("SYMBOL", 
                                         "log2FoldChange.x.x", "sample_contrast.x.x",
                                         "log2FoldChange.y.x", "sample_contrast.y.x", 
                                         "log2FoldChange.x.y", "sample_contrast.x.y",
                                         "log2FoldChange.y.y", "sample_contrast.y.y")]

# input symbols as separate col:
all_symbols <- all_outer_joined$SYMBOL

# separate out:
tmp_453_cz <- all_outer_joined[all_outer_joined$sample_contrast.x.x == "COG-N-453 Crizotinib vs Vehicle", ][, c("log2FoldChange.x.x", "sample_contrast.x.x")]
colnames(tmp_453_cz) <- c("log2FoldChange", "sample_contrast")
tmp_453_cz["log2FoldChange"][is.na(tmp_453_cz["log2FoldChange"])] <- 0
tmp_453_cz$sample_contrast <- as.factor('COG-N-453 Crizotinib vs Vehicle')

tmp_453_lr <- all_outer_joined[all_outer_joined$sample_contrast.y.x == "COG-N-453 Lorlatinib vs Vehicle", ][, c("log2FoldChange.y.x", "sample_contrast.y.x")]
colnames(tmp_453_lr) <- c("log2FoldChange", "sample_contrast")
tmp_453_lr["log2FoldChange"][is.na(tmp_453_lr["log2FoldChange"])] <- 0
tmp_453_lr$sample_contrast <- as.factor('COG-N-453 Lorlatinib vs Vehicle')

tmp_1643_cz <- all_outer_joined[all_outer_joined$sample_contrast.x.y == "NB1643 Crizotinib vs Vehicle", ][, c("log2FoldChange.x.y", "sample_contrast.x.y")]
colnames(tmp_1643_cz) <- c("log2FoldChange", "sample_contrast")
tmp_1643_cz["log2FoldChange"][is.na(tmp_1643_cz["log2FoldChange"])] <- 0
tmp_1643_cz$sample_contrast <- as.factor('NB1643 Crizotinib vs Vehicle')

tmp_1643_lr <- all_outer_joined[all_outer_joined$sample_contrast.y.y == "NB1643 Lorlatinib vs Vehicle", ][, c("log2FoldChange.y.y", "sample_contrast.y.y")]
colnames(tmp_1643_lr) <- c("log2FoldChange", "sample_contrast")
tmp_1643_lr["log2FoldChange"][is.na(tmp_1643_lr["log2FoldChange"])] <- 0
tmp_1643_lr$sample_contrast <- as.factor('NB1643 Lorlatinib vs Vehicle')

a <- cbind(all_symbols, tmp_453_cz)
b <- cbind(all_symbols, tmp_453_lr)
c <- cbind(all_symbols, tmp_1643_cz)
d <- cbind(all_symbols, tmp_1643_lr)

samples_to_plot_453 <- rbind(a,b)
samples_to_plot_1643 <- rbind(c,d)

## only genes in 453 or 1643, also sorted by expression
samples_to_plot_453_not_0 <- samples_to_plot_453[samples_to_plot_453$log2FoldChange != 0, ]
samples_to_plot_1643_not_0 <- samples_to_plot_1643[samples_to_plot_1643$log2FoldChange != 0, ]

# better colors and horizontal log2 fold changes:
samples_to_plot_453_not_0$all_symbols <- factor(samples_to_plot_453_not_0$all_symbols, levels = unique(samples_to_plot_453_not_0$all_symbols))
samples_to_plot_1643_not_0$all_symbols <- factor(samples_to_plot_1643_not_0$all_symbols, levels = unique(samples_to_plot_1643_not_0$all_symbols))

mycolors <- c("blue", "red")

samples_to_plot_453_not_0 <- samples_to_plot_453_not_0[order(-samples_to_plot_453_not_0$log2FoldChange, samples_to_plot_453_not_0$all_symbols),]
samples_to_plot_1643_not_0 <- samples_to_plot_1643_not_0[order(-samples_to_plot_1643_not_0$log2FoldChange, samples_to_plot_1643_not_0$all_symbols),]

## ignore for now: plotting these:
# ggplot(samples_to_plot_453_not_0, aes(y = reorder(all_symbols, log2FoldChange), x = log2FoldChange, fill = factor(sample_contrast))) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.y = element_text(size = 5.5), legend.position="top") +
#   ylab("Kinase Gene Names") +
#   geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
#   xlim(-2, 3.5) +
#   scale_fill_manual(values=mycolors) +
#   guides(fill=guide_legend(title=""))
# 
# ggplot(samples_to_plot_1643_not_0, aes(y = reorder(all_symbols, log2FoldChange), x = log2FoldChange, fill = sample_contrast)) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.y = element_text(size = 5.5), legend.position="top") +
#   ylab("Kinase Gene Names") +
#   geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
#   xlim(-2, 3.5) +
#   scale_fill_manual(values=mycolors) +
#   guides(fill=guide_legend(title=""))


## fig c: no adj cutoffs
# again, this is incorrectly named - this is figure 5a

# here is a list of genes of interest from Smita
# manually vectors
# check to see if Smita's genes are in each of the 6 sample-conditions:
smitas_1643_genes <- c("AURKA","PLK1","CDK2","CDK1","PLK4","NEK2","BUB1B",
                       "PKMYT1","PLK3","BUB1","PBK","RET","MELK","AURKB","PLK2")
smitas_453_genes <- c("WEE1","CHEK2","BUB1B","AURKA","BUB1","CHEK1","PLK4","NEK2",
                      "PLK3","PLK1","CDK1","MELK","CDK2","PBK","PKMYT1","AURKB")

# resume for convenience:
results_dir <- 'degs/'
cz_vs_lr_rna_453 <- read.table(paste0(results_dir, 'cz_vs_lr_rna_453.txt'))
cz_vs_lr_rna_1643 <- read.table(paste0(results_dir, 'cz_vs_lr_rna_1643.txt'))
cz_vs_v_rna_453 <- read.table(paste0(results_dir, 'cz_vs_v_rna_453.txt'))
cz_vs_v_rna_1643 <- read.table(paste0(results_dir, 'cz_vs_v_rna_1643.txt'))
lr_vs_v_rna_453 <- read.table(paste0(results_dir, 'lr_vs_v_rna_453.txt'))
lr_vs_v_rna_1643 <- read.table(paste0(results_dir, 'lr_vs_v_rna_1643.txt'))

# add SYMBOLs to each of the 6:
rownames(cz_vs_lr_rna_453) <- gsub("\\..*","",rownames(cz_vs_lr_rna_453))
cz_vs_lr_rna_453$GENEID <- rownames(cz_vs_lr_rna_453)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = cz_vs_lr_rna_453$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
cz_vs_lr_rna_453_final <- left_join(cz_vs_lr_rna_453, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(cz_vs_v_rna_453) <- gsub("\\..*","",rownames(cz_vs_v_rna_453))
cz_vs_v_rna_453$GENEID <- rownames(cz_vs_v_rna_453)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = cz_vs_v_rna_453$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
cz_vs_v_rna_453_final <- left_join(cz_vs_v_rna_453, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(lr_vs_v_rna_453) <- gsub("\\..*","",rownames(lr_vs_v_rna_453))
lr_vs_v_rna_453$GENEID <- rownames(lr_vs_v_rna_453)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = lr_vs_v_rna_453$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
lr_vs_v_rna_453_final <- left_join(lr_vs_v_rna_453, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(cz_vs_lr_rna_1643) <- gsub("\\..*","",rownames(cz_vs_lr_rna_1643))
cz_vs_lr_rna_1643$GENEID <- rownames(cz_vs_lr_rna_1643)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = cz_vs_lr_rna_1643$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
cz_vs_lr_rna_1643_final <- left_join(cz_vs_lr_rna_1643, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(cz_vs_v_rna_1643) <- gsub("\\..*","",rownames(cz_vs_v_rna_1643))
cz_vs_v_rna_1643$GENEID <- rownames(cz_vs_v_rna_1643)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = cz_vs_v_rna_1643$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
cz_vs_v_rna_1643_final <- left_join(cz_vs_v_rna_1643, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(lr_vs_v_rna_1643) <- gsub("\\..*","",rownames(lr_vs_v_rna_1643))
lr_vs_v_rna_1643$GENEID <- rownames(lr_vs_v_rna_1643)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = lr_vs_v_rna_1643$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
lr_vs_v_rna_1643_final <- left_join(lr_vs_v_rna_1643, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

# now, filter for them and plot them:
cz_vs_v_rna_453_final <- cz_vs_v_rna_453_final[cz_vs_v_rna_453_final$SYMBOL %in% smitas_453_genes,]
lr_vs_v_rna_453_final <- lr_vs_v_rna_453_final[lr_vs_v_rna_453_final$SYMBOL %in% smitas_453_genes,]
cz_vs_v_rna_1643_final <- cz_vs_v_rna_1643_final[cz_vs_v_rna_1643_final$SYMBOL %in% smitas_1643_genes,]
lr_vs_v_rna_1643_final <- lr_vs_v_rna_1643_final[lr_vs_v_rna_1643_final$SYMBOL %in% smitas_1643_genes,]

cz_vs_v_rna_453_final$contrast <- "Crizotinib vs Vehicle"
lr_vs_v_rna_453_final$contrast <- "Lorlatinib vs Vehicle"
cz_vs_v_rna_1643_final$contrast <- "Crizotinib vs Vehicle" 
lr_vs_v_rna_1643_final$contrast <- "Lorlatinib vs Vehicle"

cz_vs_v_rna_453_final <- cz_vs_v_rna_453_final[match(smitas_453_genes, cz_vs_v_rna_453_final$SYMBOL), ]
lr_vs_v_rna_453_final <- lr_vs_v_rna_453_final[match(smitas_453_genes, lr_vs_v_rna_453_final$SYMBOL), ]
cz_vs_v_rna_1643_final <- cz_vs_v_rna_1643_final[match(smitas_1643_genes, cz_vs_v_rna_1643_final$SYMBOL), ]
lr_vs_v_rna_1643_final <- lr_vs_v_rna_1643_final[match(smitas_1643_genes, lr_vs_v_rna_1643_final$SYMBOL), ]

fig_c_453_plotme <- rbind(cz_vs_v_rna_453_final, lr_vs_v_rna_453_final)
fig_c_1643_plotme <- rbind(cz_vs_v_rna_1643_final, lr_vs_v_rna_1643_final)

# 453 Keep the order:
fig_c_453_plotme <- fig_c_453_plotme[, c("log2FoldChange", "padj", "SYMBOL", "contrast")]
fig_c_453_plotme$SYMBOL <- as.factor(fig_c_453_plotme$SYMBOL)
idx <- order(c(seq_along(smitas_453_genes), seq_along(smitas_453_genes)))
tmp <- unlist(c(smitas_453_genes,smitas_453_genes))[idx]
tmp2 <- c("Crizotinib vs Vehicle", "Lorlatinib vs Vehicle","Crizotinib vs Vehicle", "Lorlatinib vs Vehicle","Crizotinib vs Vehicle", "Lorlatinib vs Vehicle",
          "Crizotinib vs Vehicle", "Lorlatinib vs Vehicle","Crizotinib vs Vehicle", "Lorlatinib vs Vehicle","Crizotinib vs Vehicle", "Lorlatinib vs Vehicle",
          "Crizotinib vs Vehicle", "Lorlatinib vs Vehicle","Crizotinib vs Vehicle", "Lorlatinib vs Vehicle","Crizotinib vs Vehicle", "Lorlatinib vs Vehicle",
          "Crizotinib vs Vehicle", "Lorlatinib vs Vehicle","Crizotinib vs Vehicle", "Lorlatinib vs Vehicle","Crizotinib vs Vehicle", "Lorlatinib vs Vehicle",
          "Crizotinib vs Vehicle", "Lorlatinib vs Vehicle","Crizotinib vs Vehicle", "Lorlatinib vs Vehicle","Crizotinib vs Vehicle", "Lorlatinib vs Vehicle",
          "Crizotinib vs Vehicle", "Lorlatinib vs Vehicle")
tmp3 <- as.data.frame(cbind(tmp, tmp2))
colnames(tmp3) <- c("SYMBOL", "contrast")
tmp3$SYMBOL <- as.factor(tmp3$SYMBOL)
fig_c_453_plotme_final <- left_join(tmp3, fig_c_453_plotme, by = c('SYMBOL', 'contrast'))
fig_c_453_plotme_final$SYMBOL <- as.factor(fig_c_453_plotme_final$SYMBOL)
fig_c_453_plotme_final$contrast <- as.factor(fig_c_453_plotme_final$contrast)
fig_c_453_plotme_final$SYMBOL <- factor(fig_c_453_plotme_final$SYMBOL,levels=unique(fig_c_453_plotme_final$SYMBOL))

# 1643 Keep the order:
fig_c_1643_plotme <- fig_c_1643_plotme[, c("log2FoldChange", "padj", "SYMBOL", "contrast")]
fig_c_1643_plotme$SYMBOL <- as.factor(fig_c_1643_plotme$SYMBOL)
idx <- order(c(seq_along(smitas_1643_genes), seq_along(smitas_1643_genes)))
tmp <- unlist(c(smitas_1643_genes,smitas_1643_genes))[idx]
tmp2 <- c("Crizotinib vs Vehicle", "Lorlatinib vs Vehicle","Crizotinib vs Vehicle", "Lorlatinib vs Vehicle","Crizotinib vs Vehicle", "Lorlatinib vs Vehicle",
          "Crizotinib vs Vehicle", "Lorlatinib vs Vehicle","Crizotinib vs Vehicle", "Lorlatinib vs Vehicle","Crizotinib vs Vehicle", "Lorlatinib vs Vehicle",
          "Crizotinib vs Vehicle", "Lorlatinib vs Vehicle","Crizotinib vs Vehicle", "Lorlatinib vs Vehicle","Crizotinib vs Vehicle", "Lorlatinib vs Vehicle",
          "Crizotinib vs Vehicle", "Lorlatinib vs Vehicle","Crizotinib vs Vehicle", "Lorlatinib vs Vehicle","Crizotinib vs Vehicle", "Lorlatinib vs Vehicle",
          "Crizotinib vs Vehicle", "Lorlatinib vs Vehicle","Crizotinib vs Vehicle", "Lorlatinib vs Vehicle","Crizotinib vs Vehicle", "Lorlatinib vs Vehicle")
tmp3 <- as.data.frame(cbind(tmp, tmp2))
colnames(tmp3) <- c("SYMBOL", "contrast")
tmp3$SYMBOL <- as.factor(tmp3$SYMBOL)
fig_c_1643_plotme_final <- left_join(tmp3, fig_c_1643_plotme, by = c('SYMBOL', 'contrast'))
fig_c_1643_plotme_final$SYMBOL <- as.factor(fig_c_1643_plotme_final$SYMBOL)
fig_c_1643_plotme_final$contrast <- as.factor(fig_c_1643_plotme_final$contrast)
fig_c_1643_plotme_final$SYMBOL <- factor(fig_c_1643_plotme_final$SYMBOL,levels=unique(fig_c_1643_plotme_final$SYMBOL))

# intersect 453 and 1643
fig_c_453_plotme_final$sample <- "COG-N-453x"
fig_c_1643_plotme_final$sample <- "NB-1643"

findme <- intersect(fig_c_453_plotme_final$SYMBOL, fig_c_1643_plotme_final$SYMBOL)

fig_c_1643_plotme_final_inter <- fig_c_1643_plotme_final[fig_c_1643_plotme_final$SYMBOL %in% findme,]
fig_c_453_plotme_final_inter <- fig_c_453_plotme_final[fig_c_453_plotme_final$SYMBOL %in% findme,]

both_plot_me <- rbind(fig_c_1643_plotme_final_inter, fig_c_453_plotme_final_inter)
both_plot_me$contrast_and_sample <- paste0(both_plot_me$sample, " ", both_plot_me$contrast)

# then, Y is L2FC and X is each gene while 1643 is solid red and blue and 453 is dashed red and blue

tiff("figures/figure_5a.tiff", units="in", width=16, height=8, res=300)
ggplot(both_plot_me, aes(fill=contrast_and_sample, y=SYMBOL, x=log2FoldChange, pattern=contrast_and_sample)) + 
  geom_bar_pattern(position="dodge", stat="identity", pattern_spacing = 0.01, pattern_angle = 45,pattern_frequency = 5) +
  theme_bw() +
  scale_fill_manual(values=c("blue","red","blue", "red")) +
  scale_pattern_manual(values=c('none', 'none', 'stripe', 'stripe')) +
  theme(text = element_text(size=16))
dev.off()

## was Smita's figure d but is now 5b

# again, Smita's genes of interest:
d_genes <- c("FER", "TNK2", "PTK2B", "PTK2", "ALK")

# resume for convenience:
results_dir <- 'degs/'
cz_vs_lr_rna_453 <- read.table(paste0(results_dir, 'cz_vs_lr_rna_453.txt'))
cz_vs_lr_rna_1643 <- read.table(paste0(results_dir, 'cz_vs_lr_rna_1643.txt'))
cz_vs_v_rna_453 <- read.table(paste0(results_dir, 'cz_vs_v_rna_453.txt'))
cz_vs_v_rna_1643 <- read.table(paste0(results_dir, 'cz_vs_v_rna_1643.txt'))
lr_vs_v_rna_453 <- read.table(paste0(results_dir, 'lr_vs_v_rna_453.txt'))
lr_vs_v_rna_1643 <- read.table(paste0(results_dir, 'lr_vs_v_rna_1643.txt'))

# add SYMBOLs to each of the 6:
rownames(cz_vs_v_rna_453) <- gsub("\\..*","",rownames(cz_vs_v_rna_453))
cz_vs_v_rna_453$GENEID <- rownames(cz_vs_v_rna_453)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = cz_vs_v_rna_453$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
cz_vs_v_rna_453_final <- left_join(cz_vs_v_rna_453, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(lr_vs_v_rna_453) <- gsub("\\..*","",rownames(lr_vs_v_rna_453))
lr_vs_v_rna_453$GENEID <- rownames(lr_vs_v_rna_453)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = lr_vs_v_rna_453$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
lr_vs_v_rna_453_final <- left_join(lr_vs_v_rna_453, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(cz_vs_v_rna_1643) <- gsub("\\..*","",rownames(cz_vs_v_rna_1643))
cz_vs_v_rna_1643$GENEID <- rownames(cz_vs_v_rna_1643)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = cz_vs_v_rna_1643$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
cz_vs_v_rna_1643_final <- left_join(cz_vs_v_rna_1643, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(lr_vs_v_rna_1643) <- gsub("\\..*","",rownames(lr_vs_v_rna_1643))
lr_vs_v_rna_1643$GENEID <- rownames(lr_vs_v_rna_1643)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = lr_vs_v_rna_1643$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
lr_vs_v_rna_1643_final <- left_join(lr_vs_v_rna_1643, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

# find the genes:
cz_vs_v_rna_453_final <- cz_vs_v_rna_453_final[cz_vs_v_rna_453_final$SYMBOL %in% d_genes,]
lr_vs_v_rna_453_final <- lr_vs_v_rna_453_final[lr_vs_v_rna_453_final$SYMBOL %in% d_genes,]
cz_vs_v_rna_1643_final <- cz_vs_v_rna_1643_final[cz_vs_v_rna_1643_final$SYMBOL %in% d_genes,]
lr_vs_v_rna_1643_final <- lr_vs_v_rna_1643_final[lr_vs_v_rna_1643_final$SYMBOL %in% d_genes,]

cz_vs_v_rna_453_final$contrast <- "Crizotinib vs Vehicle"
lr_vs_v_rna_453_final$contrast <- "Lorlatinib vs Vehicle"
cz_vs_v_rna_1643_final$contrast <- "Crizotinib vs Vehicle"
lr_vs_v_rna_1643_final$contrast <- "Lorlatinib vs Vehicle"

plot_me_453_d <- rbind(cz_vs_v_rna_453_final, lr_vs_v_rna_453_final)
plot_me_1643_d <- rbind(cz_vs_v_rna_1643_final, lr_vs_v_rna_1643_final)

# force order: 
plot_me_453_d$SYMBOL <- factor(plot_me_453_d$SYMBOL,levels=unique(plot_me_453_d$SYMBOL))
plot_me_1643_d$SYMBOL <- factor(plot_me_1643_d$SYMBOL,levels=unique(plot_me_1643_d$SYMBOL))

plot_me_453_d <- rbind(plot_me_453_d[3,], plot_me_453_d[8,], plot_me_453_d[1,], plot_me_453_d[6,],
                       plot_me_453_d[2,], plot_me_453_d[7,], plot_me_453_d[4,], plot_me_453_d[9,],
                       plot_me_453_d[5,], plot_me_453_d[10,])

plot_me_1643_d <- rbind(plot_me_1643_d[3,], plot_me_1643_d[8,], plot_me_1643_d[1,], plot_me_1643_d[6,],
                        plot_me_1643_d[2,], plot_me_1643_d[7,], plot_me_1643_d[4,], plot_me_1643_d[9,],
                        plot_me_1643_d[5,], plot_me_1643_d[10,])

###########################
## stop here if need padj and use the below to get them together with sample names:
#plot_me_453_d$contrast_and_sample <- paste0("COG-N-453x", " ", plot_me_453_d$contrast)
#plot_me_1643_d$contrast_and_sample <- paste0("NB-1643", " ", plot_me_1643_d$contrast)
#both_plot_me_d <- rbind(plot_me_453_d, plot_me_1643_d)

#write.table(both_plot_me, 'helper_files/figure_5a_pasted.txt', sep='\t', quote=F, row.names = F)
#write.table(both_plot_me_d, 'helper_files/figure_5b_pasted.txt', sep='\t', quote=F, row.names = F)
###########################

plot_me_453_d <- plot_me_453_d[, c("log2FoldChange", "SYMBOL", "contrast")]
plot_me_1643_d <- plot_me_1643_d[, c("log2FoldChange", "SYMBOL", "contrast")]

plot_me_453_d$contrast <- factor(plot_me_453_d$contrast,levels=unique(plot_me_453_d$contrast))
plot_me_1643_d$contrast <- factor(plot_me_1643_d$contrast,levels=unique(plot_me_1643_d$contrast))

plot_me_453_d$SYMBOL <- factor(plot_me_453_d$SYMBOL,levels=d_genes)
plot_me_1643_d$SYMBOL <- factor(plot_me_1643_d$SYMBOL,levels=d_genes)

## both on same plot:
plot_me_453_d$contrast_and_sample <- paste0("COG-N-453x", " ", plot_me_453_d$contrast)
plot_me_1643_d$contrast_and_sample <- paste0("NB-1643", " ", plot_me_1643_d$contrast)
both_plot_me_d <- rbind(plot_me_453_d, plot_me_1643_d)

tiff("figures/figure_5b.tiff", units="in", width=16, height=8, res=300)
ggplot(both_plot_me_d, aes(fill=contrast_and_sample, y=SYMBOL, x=log2FoldChange, pattern=contrast_and_sample)) + 
  geom_bar_pattern(position="dodge", stat="identity", pattern_spacing = 0.01, pattern_angle = 45,pattern_frequency = 5) +
  theme_bw() +
  scale_fill_manual(values=c("blue","red","blue", "red")) +
  scale_pattern_manual(values=c('none', 'none', 'stripe', 'stripe')) +
  theme(text = element_text(size=16))
dev.off()

# supplemental data
# all genes down-regulated across all sample-conditions

# resume from base results:
results_dir <- 'degs/'
cz_vs_lr_rna_453 <- read.table(paste0(results_dir, 'cz_vs_lr_rna_453.txt'))
cz_vs_lr_rna_1643 <- read.table(paste0(results_dir, 'cz_vs_lr_rna_1643.txt'))
cz_vs_v_rna_453 <- read.table(paste0(results_dir, 'cz_vs_v_rna_453.txt'))
cz_vs_v_rna_1643 <- read.table(paste0(results_dir, 'cz_vs_v_rna_1643.txt'))
lr_vs_v_rna_453 <- read.table(paste0(results_dir, 'lr_vs_v_rna_453.txt'))
lr_vs_v_rna_1643 <- read.table(paste0(results_dir, 'lr_vs_v_rna_1643.txt'))

# add SYMBOLs to each of the 6:
rownames(cz_vs_lr_rna_453) <- gsub("\\..*","",rownames(cz_vs_lr_rna_453))
cz_vs_lr_rna_453$GENEID <- rownames(cz_vs_lr_rna_453)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = cz_vs_lr_rna_453$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
cz_vs_lr_rna_453_final <- left_join(cz_vs_lr_rna_453, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(cz_vs_v_rna_453) <- gsub("\\..*","",rownames(cz_vs_v_rna_453))
cz_vs_v_rna_453$GENEID <- rownames(cz_vs_v_rna_453)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = cz_vs_v_rna_453$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
cz_vs_v_rna_453_final <- left_join(cz_vs_v_rna_453, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(lr_vs_v_rna_453) <- gsub("\\..*","",rownames(lr_vs_v_rna_453))
lr_vs_v_rna_453$GENEID <- rownames(lr_vs_v_rna_453)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = lr_vs_v_rna_453$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
lr_vs_v_rna_453_final <- left_join(lr_vs_v_rna_453, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(cz_vs_lr_rna_1643) <- gsub("\\..*","",rownames(cz_vs_lr_rna_1643))
cz_vs_lr_rna_1643$GENEID <- rownames(cz_vs_lr_rna_1643)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = cz_vs_lr_rna_1643$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
cz_vs_lr_rna_1643_final <- left_join(cz_vs_lr_rna_1643, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(cz_vs_v_rna_1643) <- gsub("\\..*","",rownames(cz_vs_v_rna_1643))
cz_vs_v_rna_1643$GENEID <- rownames(cz_vs_v_rna_1643)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = cz_vs_v_rna_1643$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
cz_vs_v_rna_1643_final <- left_join(cz_vs_v_rna_1643, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

rownames(lr_vs_v_rna_1643) <- gsub("\\..*","",rownames(lr_vs_v_rna_1643))
lr_vs_v_rna_1643$GENEID <- rownames(lr_vs_v_rna_1643)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys = lr_vs_v_rna_1643$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
lr_vs_v_rna_1643_final <- left_join(lr_vs_v_rna_1643, geneIDs1, by="GENEID") %>%
  mutate(GENEID=ifelse(is.na(SYMBOL), GENEID, SYMBOL))

# only consider all vs vehicle:
cz_vs_v_rna_453_final$sample_and_contrast <- paste0("COG-N-453x", " ", "Crizotinib vs Vehicle")
lr_vs_v_rna_453_final$sample_and_contrast <- paste0("COG-N-453x", " ", "Lorlatinib vs Vehicle")
cz_vs_v_rna_1643_final$sample_and_contrast <- paste0("NB-1643", " ", "Crizotinib vs Vehicle")
lr_vs_v_rna_1643_final$sample_and_contrast <- paste0("NB-1643", " ", "Lorlatinib vs Vehicle")

# filter for down-reg only, intersect, and combine:
cz_vs_v_rna_453_final_down <- na.omit(cz_vs_v_rna_453_final[cz_vs_v_rna_453_final$log2FoldChange < 0,])
lr_vs_v_rna_453_final_down <- na.omit(lr_vs_v_rna_453_final[lr_vs_v_rna_453_final$log2FoldChange < 0,])
cz_vs_v_rna_1643_final_down <- na.omit(cz_vs_v_rna_1643_final[cz_vs_v_rna_1643_final$log2FoldChange < 0,])
lr_vs_v_rna_1643_final_down <- na.omit(lr_vs_v_rna_1643_final[lr_vs_v_rna_1643_final$log2FoldChange < 0,])

cz_vs_v_rna_453_final_down <- cz_vs_v_rna_453_final_down %>% dplyr::filter(padj < 0.05)
lr_vs_v_rna_453_final_down <- lr_vs_v_rna_453_final_down %>% dplyr::filter(padj < 0.05)
cz_vs_v_rna_1643_final_down <- cz_vs_v_rna_1643_final_down %>% dplyr::filter(padj < 0.05)
lr_vs_v_rna_1643_final_down <- lr_vs_v_rna_1643_final_down %>% dplyr::filter(padj < 0.05)

a <- intersect(cz_vs_v_rna_453_final_down$SYMBOL, lr_vs_v_rna_453_final_down$SYMBOL)
b <- intersect(cz_vs_v_rna_1643_final_down$SYMBOL, lr_vs_v_rna_1643_final_down$SYMBOL)

combo_down <- intersect(a,b)

combo_down_bind <- rbind(cz_vs_v_rna_453_final_down, lr_vs_v_rna_453_final_down,
                         cz_vs_v_rna_1643_final_down, lr_vs_v_rna_1643_final_down)

combo_down_bind_plot <- combo_down_bind[combo_down_bind$SYMBOL %in% combo_down,]

# stacked barchart:
tiff("figures/supplemental_rna_seq_all_shared_downreg.tiff", units="in", width=16, height=8, res=300)
ggplot(combo_down_bind_plot, aes(fill=sample_and_contrast, y=SYMBOL, x=log2FoldChange, pattern=sample_and_contrast)) + 
  geom_bar_pattern(position="dodge", stat="identity", pattern_spacing = 0.01, pattern_angle = 45,pattern_frequency = 5) +
  theme_bw() +
  scale_fill_manual(values=c("blue","red","blue", "red")) +
  scale_pattern_manual(values=c('none', 'none', 'stripe', 'stripe')) +
  theme(text = element_text(size=16))
dev.off()


## figure 5: MYCN amplification in the 6 samples
# use vst or rlog function from DESeq2 to compensate for the effect of different library sizes and place into log2 scale

# obtain log2 counts
vsd <- vst(dds,blind=TRUE)
vsd.1643 <- vst(dds.1643, blind=TRUE)

vsd_counts <- vsd@assays@data@listData[[1]]
vsd.1643_counts <- vsd.1643@assays@data@listData[[1]]

dds_counts <- dds@assays@data@listData[[1]]
dds.1643_counts <- dds.1643@assays@data@listData[[1]]

# MYCN ENSG00000134323, ENSG00000134323.10:
vsd_counts_mycn <- vsd_counts[rownames(vsd_counts) == "ENSG00000134323.10",]
vsd.1643_counts_mycn <- vsd.1643_counts[rownames(vsd.1643_counts) == "ENSG00000134323.10",]

dds_counts_mycn <- dds_counts[rownames(dds_counts) == "ENSG00000134323.10",]
dds.1643_counts_mycn <- dds.1643_counts[rownames(dds.1643_counts) == "ENSG00000134323.10",]

# want collapse the replicates:
vsd_counts_mycn_df <- data.frame("values" = vsd_counts_mycn, "group" = names(vsd_counts_mycn), row.names = NULL)
a <- unlist(strsplit(vsd_counts_mycn_df$group, "_"))
vsd_counts_mycn_df$group <- a[seq(2, length(a), 3)]

# re-order like Smita has it in other figures:
vsd_counts_mycn_df <- rbind(vsd_counts_mycn_df[vsd_counts_mycn_df$group=="Vehicle",], 
                            vsd_counts_mycn_df[vsd_counts_mycn_df$group=="Crizotinib",], 
                            vsd_counts_mycn_df[vsd_counts_mycn_df$group=="Lorlatinib",])
vsd_counts_mycn_df$group <- factor(vsd_counts_mycn_df$group, levels=c("Vehicle", "Crizotinib", "Lorlatinib"))

vsd.1643_counts_mycn_df <- data.frame("values" = vsd.1643_counts_mycn, "group" = names(vsd.1643_counts_mycn), row.names = NULL)
a <- unlist(strsplit(vsd.1643_counts_mycn_df$group, "_"))
vsd.1643_counts_mycn_df$group <- a[seq(2, length(a), 3)]

# re-order like Smita has it:
vsd.1643_counts_mycn_df <- rbind(vsd.1643_counts_mycn_df[vsd.1643_counts_mycn_df$group=="Vehicle",], 
                                 vsd.1643_counts_mycn_df[vsd.1643_counts_mycn_df$group=="Crizotinib",], 
                                 vsd.1643_counts_mycn_df[vsd.1643_counts_mycn_df$group=="Lorlatinib",])
vsd.1643_counts_mycn_df$group <- factor(vsd.1643_counts_mycn_df$group, levels=c("Vehicle", "Crizotinib", "Lorlatinib"))

# dotplot:
vsd_counts_mycn_df_2 <- vsd_counts_mycn_df %>% group_by(group) %>% 
  mutate(min=min(values,na.rm=T), max=max(values,na.rm=T), avg=mean(values,na.rm=T))

vsd.1643_counts_mycn_df_2 <- vsd.1643_counts_mycn_df %>% group_by(group) %>% 
  mutate(min=min(values,na.rm=T), max=max(values,na.rm=T), avg=mean(values,na.rm=T))

# normalize Vehicle to 1:
# VEHICLE: find mean, subtract FROM mean for each value, add one to each of these
# CZ/LR: subtract FROM VEHICLE mean (V mean - value), add one to each, 1 divided by these values (1/value)
mean_subtract <- vsd_counts_mycn_df_2[vsd_counts_mycn_df_2$group == "Vehicle",]$avg[1]
vsd_counts_mycn_df_2$subt <- mean_subtract - vsd_counts_mycn_df_2$values
vsd_counts_mycn_df_2$add_one <- vsd_counts_mycn_df_2$subt + 1

tmp_add <- c()
for (i in 1:length(vsd_counts_mycn_df_2$group)) {
  if (vsd_counts_mycn_df_2$group[i] == "Vehicle") {
    tmp_add <- c(tmp_add, vsd_counts_mycn_df_2$add_one[i])
  }
  else {
    tmp_add <- c(tmp_add, 1/vsd_counts_mycn_df_2$add_one[i])
  }
}
vsd_counts_mycn_df_2$final_plotme <- tmp_add

# 1643:
mean_subtract <- vsd.1643_counts_mycn_df_2[vsd.1643_counts_mycn_df_2$group == "Vehicle",]$avg[1]
vsd.1643_counts_mycn_df_2$subt <- mean_subtract - vsd.1643_counts_mycn_df_2$values
vsd.1643_counts_mycn_df_2$add_one <- vsd.1643_counts_mycn_df_2$subt + 1

tmp_add <- c()
for (i in 1:length(vsd.1643_counts_mycn_df_2$group)) {
  if (vsd.1643_counts_mycn_df_2$group[i] == "Vehicle") {
    tmp_add <- c(tmp_add, vsd.1643_counts_mycn_df_2$add_one[i])
  }
  else {
    tmp_add <- c(tmp_add, 1/vsd.1643_counts_mycn_df_2$add_one[i])
  }
}
vsd.1643_counts_mycn_df_2$final_plotme <- tmp_add

# add back in min, max, avgs:
vsd_counts_mycn_df_2_new <- vsd_counts_mycn_df_2 %>% 
  group_by(group) %>% mutate(min=min(final_plotme,na.rm=T), max=max(final_plotme,na.rm=T), avg=mean(final_plotme,na.rm=T))

vsd.1643_counts_mycn_df_2_new <- vsd.1643_counts_mycn_df_2 %>% 
  group_by(group) %>% mutate(min=min(final_plotme,na.rm=T), max=max(final_plotme,na.rm=T), avg=mean(final_plotme,na.rm=T))

## plot these
tiff("figures/figure_5c.tiff", units="in", width=16, height=8, res=300)
ggplot(vsd_counts_mycn_df_2_new, aes(x=reorder(group, final_plotme), y=final_plotme, fill=group)) + 
  geom_dotplot(binaxis='y', stackdir='center', stroke=NA, stackratio=0.1, dotsize=0.2) +
  geom_errorbar(aes(ymin = avg, ymax = max, color=factor(group)), width = 0.5,size=0.8) +
  geom_errorbar(aes(ymin = min, ymax = avg, color=factor(group)), width = 0.5,size=0.8) +
  geom_errorbar(aes(ymin = avg, ymax = avg), width = 0.9,size=0.8) +
  scale_x_discrete(limits = levels(vsd_counts_mycn_df_2_new$group)) +
  geom_point(aes(shape = group, color=group), size = 4) +
  scale_shape_manual(values = c(16, 15, 17)) +
  scale_fill_manual(values=c("black", "blue","red")) +
  scale_color_manual(values=c("black", "blue","red")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylim(0.0, 1.2)
dev.off()

tiff("figures/figure_5d.tiff", units="in", width=16, height=8, res=300)
ggplot(vsd.1643_counts_mycn_df_2_new, aes(x=reorder(group, final_plotme), y=final_plotme, fill=group)) + 
  geom_dotplot(binaxis='y', stackdir='center', stroke=NA, stackratio=0.1, dotsize=0.2) +
  geom_errorbar(aes(ymin = avg, ymax = max, color=factor(group)), width = 0.5,size=0.8) +
  geom_errorbar(aes(ymin = min, ymax = avg, color=factor(group)), width = 0.5,size=0.8) +
  geom_errorbar(aes(ymin = avg, ymax = avg), width = 0.9,size=0.8) +
  scale_x_discrete(limits = levels(vsd.1643_counts_mycn_df_2_new$group)) +
  geom_point(aes(shape = group, color=group), size = 4) +
  scale_shape_manual(values = c(16, 15, 17)) +
  scale_fill_manual(values=c("black", "blue","red")) +
  scale_color_manual(values=c("black", "blue","red")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylim(0.0, 1.2)
dev.off()

