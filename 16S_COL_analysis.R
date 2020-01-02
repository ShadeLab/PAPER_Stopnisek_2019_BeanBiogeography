#*************************************************
#IMPORTANT: run first the 16S_US_analysis.R script
#*************************************************

library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(RColorBrewer)
library(tidyverse)
library(gridExtra)
library(limma)
library(maps)
library(reshape2)
library(broom)
library(ggpubr)
library(phyloseq)
library(DESeq2)
library(stats4)
library(minpack.lm)
library(Hmisc)
library(tidytext)
library(stringr)
library(ggrepel)
library(ggExtra)

setwd('~/Documents/git/PAPER_Stopnisek_2019_BeanBiogeography/')

###**********************************************************************************
#Colombia comparison

#reference based and de-novo OTU picking
otu_col <- read.table('Data/OTU_table_col.txt', sep='\t', row.names = 1, header=T)

map_col <- read.table('Data/map_col.txt', row.names = 1, sep='\t', header=T)

otu_col <- otu_col[,order(colnames(otu_col))]
# Order the samples of the map the same way
map_col=map_col[order(rownames(map_col)),]
# Check to make sure they all match with each other
rownames(map_col)==colnames(otu_col)

otu_col_bulk <- otu_col[,map_col$Group == 'Bulk']
otu_col_bulk <- otu_col_bulk[rowSums(otu_col_bulk)>0,]
otu_col_rhizo <- otu_col[,map_col$Group != 'Bulk']
otu_col_rhizo <- otu_col_rhizo[rowSums(otu_col_rhizo)>0,]
set.seed(51)
otu_col_rare <- t(rrarefy(t(otu_col_rhizo), min(colSums(otu_col_rhizo)))) 
set.seed(83)
otu_col_rare_bulk<- t(rrarefy(t(otu_col_bulk), min(colSums(otu_col)))) 

#col occ_abun plot
#remove the bulk samples
otu_col_rare_rhizo <- otu_col_rare[rowSums(otu_col_rare)>0,]
col_otu_PA_rhizo <- 1*((otu_col_rare_rhizo>0)==1)
col_otu_PA_rhizo <- col_otu_PA_rhizo[rowSums(col_otu_PA_rhizo)>0,]
Occ_col <- rowSums(col_otu_PA_rhizo)/ncol(col_otu_PA_rhizo)

col_rare_rhizo_rel <- decostand(otu_col_rare_rhizo, method = 'total', MARGIN = 2)
com_abund_col_rhizo <- rowSums(col_rare_rhizo_rel)/ncol(col_rare_rhizo_rel)
col_rhizo_df_occ <- data.frame(otu=names(Occ_col), occ=Occ_col) 
col_rhizo_df_abun <- data.frame(otu=names(com_abund_col_rhizo), abun=log10(com_abund_col_rhizo))
occ_abun_col_rhizo <- left_join(col_rhizo_df_occ, col_rhizo_df_abun)

length(Occ_col[Occ_col==1]) #848 OTUs makes the Colombian core

#making a list of names from the US study - only the reference picked OTUs
US_rhizo_test <- as.data.frame(otu_rare_rhizo)
US_rhizo_test$names <- rownames(US_rhizo_test)
US_ref_names <- US_rhizo_test %>% 
  filter(!str_detect(names, 'OTU_dn'))
US_ref_otus <- US_ref_names$names

Col_rhizo_test <- as.data.frame(otu_col_rare_rhizo)
Col_rhizo_test$names <- rownames(Col_rhizo_test)
Col_ref_names <- Col_rhizo_test %>% 
  filter(!str_detect(names, 'OTU_dn'))
Col_ref_otus <- Col_ref_names$names

#making a list of col de-novo OTUs that have 100% similarity to US de-novo OTUs (results from blastn)
blast_out <- read.table('Data/blast_results_denovo.txt', header=T)
col_us_match <- as.character(blast_out$Br_subject)

#combining the name lists
names_matching <- append(US_ref_otus,col_us_match)
col_names_matching <- append(Col_ref_otus, col_us_match)

occ_abun_col_rhizo$found <- 'Colombia unique (n=8493)'
occ_abun_col_rhizo$found[occ_abun_col_rhizo$otu %in% names_matching] <- 'shared with US (n=3359)'

########
# Fig 3A

ggplot(data=occ_abun_col_rhizo, aes(x=abun, y=occ, bg=found)) +
  geom_point(size=3, pch=21, alpha=.5) +
  theme_bw()+
  scale_fill_manual(values=c('black', 'white')) +
  labs(x=paste("log10(mean abundace)\n (n=", nrow(occ_abun_col_rhizo)," OTUs)", sep=''), y= paste("Occupancy (n=", length(colnames(col_otu_PA_rhizo))," samples)",sep=''), 
       bg= NULL) +
  theme(plot.title = element_text(hjust = 0.5, size=12),
        legend.position = 'top',
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  guides(fill = guide_legend(override.aes = list(alpha=1)))

US_occ_1 <- as.vector(combined_occ_data$otu[combined_occ_data$occ== 1])
COL_occ_1 <- as.vector(occ_abun_col_rhizo$otu[occ_abun_col_rhizo$occ==1])

#number of OTUs sahred between the cores?
#first removig denovo OTUs that are not shared
core_US_otus <- US_occ_1[US_occ_1 %in% names_matching] 
core_COL_otus <- COL_occ_1[COL_occ_1 %in% names_matching] 

core_otus <- core_COL_otus[core_COL_otus %in% core_US_otus] #48 OTUs shared between the Occ=1

##Investigating 48 OTUs found in both US and Colombian cores
global_core <- core_otus
core_subset_OTU_table <- otu_us_rare_rhizo[rownames(otu_us_rare_rhizo) %in% global_core,]

combined_occ_data_v2=combined_occ_data
combined_occ_data_v2$bean_found[combined_occ_data_v2$otu %in% global_core] <- 'core'

otu_rare_rhizo <- otu_us_rare[,map_combined$soil=='rhizosphere']
otu.rhizo.rel.abun <- decostand(otu_rare_rhizo, method="total", MARGIN=2) #calculating relative abundance

core_tmp <- data.frame(otu=as.factor(rownames(otu.rhizo.rel.abun)),otu.rhizo.rel.abun) %>%
  gather(sample_ID, abun, -otu) %>%
  left_join(map_combined[,c('sample_ID','pH', 'state', 'bean', 'plot', 'site')], by='sample_ID') %>%
  left_join(tax_v1, by='otu') 

core_tmp_v2 <- core_tmp[core_tmp$otu %in% global_core, ] %>%
  group_by(Phylum,Order,site) %>%
  dplyr::summarise(
    n=sum(abun)/length(unique(sample_ID))) %>%
  arrange(desc(n)) 

col_rare_rhizo_rel <- decostand(otu_col_rare_rhizo, method = 'total', MARGIN = 2)
core_tmp_col <- data.frame(otu=as.factor(rownames(col_rare_rhizo_rel)),col_rare_rhizo_rel) %>%
  gather(sample_ID, abun, -otu)
core_tmp_col$country <- 'Colombia'

core_tmp_col_v2 <- core_tmp_col[as.character(core_tmp_col$otu) %in% global_core,]
tax_core <- tax_v1[tax_v1$otu %in% global_core, ]
core_tmp_col_v2$otu <- as.character(core_tmp_col_v2$otu)
core_tmp_col_v2$otu

core_col_tmp_v3 <- core_tmp_col_v2 %>%
  left_join(tax_core, by='otu') 

core_col_tmp_v4 <- core_col_tmp_v3 %>%
  group_by(country,Phylum,Order) %>%
  summarise(
    n=sum(abun)/length(country)) %>%
  arrange(desc(n))

core_col_tmp_v4$site <- 'Colombia'
core_tmp_v2$country <- 'US'
global_core_df <- rbind(core_tmp_v2, core_col_tmp_v4)

core_col <- core_col_tmp_v3
core_col$site <- 'Colombia'
core_tmp_short <- core_tmp[core_tmp$otu %in% global_core,][,-c(4,5,6,7,9,10,11,19)]
core_tmp_short$country <- 'US'

core_tmp_v2 <- core_tmp[core_tmp$otu %in% global_core, ] %>%
  group_by(Phylum,Order) %>%
  dplyr::summarise(
    n=sum(abun)/length(unique(sample_ID)),
    std=sqrt(sum((abun-mean(abun))^2/(length(unique(sample_ID))-1)))) %>%
  arrange(desc(n))

core_tmp_v3 <- core_tmp[core_tmp$otu %in% global_core, ] %>%
  mutate(country='US') %>%
  group_by(Phylum,Class,Order, Genus) %>%
  dplyr::summarise(
    n=sum(abun)/length(country),
    std=sqrt(sum((abun-mean(abun))^2/(length(country)-1)))) %>%
  arrange(desc(n)) 

core_col_tmp_v4 <- core_col_tmp_v3 %>%
  group_by(Phylum,Order) %>%
  dplyr::summarise(
    n=sum(abun)/length(sample_ID),
    std=sqrt(sum((abun-mean(abun))^2/(length(sample_ID)-1)))) %>%
  arrange(desc(n))

core_tmp_v2$site <- 'US'
core_tmp_v2$country <- 'US'
core_col_tmp_v4$site <- 'Colombia'
core_col_tmp_v4$country <- 'Colombia'
joined_core_taxa <- rbind(core_col_tmp_v4, core_tmp_v2)

core_v1 <- core_tmp[core_tmp$otu %in% global_core,] %>%
  mutate(country='US') %>%
  mutate(new_names = if_else((Order == "uncultured bacterium" | Order == "Ambiguous_taxa" | is.na(Order)), Phylum, Order),
         new_names = if_else((new_names == 'uncultured Acidobacteria bacterium'), 'Acidobacteria bacterium', new_names))

core_col_v1 <- core_col_tmp_v3 %>%
  mutate(new_names = if_else((Order == "uncultured bacterium" | Order == "Ambiguous_taxa" | is.na(Order)), Phylum, Order),
         new_names = if_else((new_names == 'uncultured Acidobacteria bacterium'), 'Acidobacteria bacterium', new_names))

core_v1 <- core_v1[,-c(4,5,6,7,8,9,10,11,18)]
core_col_v1 <- core_col_v1[,-c(5,6,7,14)]

core_all_v1 <- rbind(core_col_v1,core_v1)

core_all_v1 <- core_all_v1 %>%
  mutate(new_names=as.factor(new_names),
         new_names=fct_reorder(new_names, abun))

x = tapply(core_all_v1$abun, core_all_v1$new_names, function(x) mean(x))
x = sort(x, TRUE)
core_all_v1$new_names = factor(as.character(core_all_v1$new_names), levels=names(x))

########
# Fig 3C

auto_plot_fig <- ggplot(core_all_v1,aes(x=new_names, y=abun, fill=Phylum, color=country)) +
  theme_bw()+
  scale_color_manual(values=c('black', 'grey')) +
  labs(y='Relative abundance', x='Order') +
  scale_fill_manual(values=tax_colors)+
  geom_boxplot() +
  theme(legend.position = c(.7, .72),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.background=element_blank()) +
  coord_flip()

coreTaxaNo <- core_all_v1 %>%
  group_by(Phylum, new_names) %>%
  summarise(n_taxa=length(unique(otu))) %>%
  ggplot(aes(x=new_names, n_taxa, fill=Phylum)) +
  geom_bar(stat='identity') +
  ylab('Number of taxa') +
  xlab('') +
  scale_fill_manual(values=tax_colors)+
  xlab('') +
  theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = 'none',
        legend.text=element_text(size=6),
        legend.title = element_text(size=8))+
  scale_y_continuous(breaks=seq(0,9,3))+
  coord_flip()

grid.arrange(auto_plot_fig, coreTaxaNo, widths=c(2.5,.7))


######
#Venn
######
matchBLAST <- col_us_match

allUS <- as.data.frame(otu_rare_rhizo) %>%
  mutate(otu=rownames(otu_rare_rhizo),
         otu2=if_else((otu %in% matchBLAST), otu, paste(otu,'us', sep='.')),
         otu3=if_else(grepl('^OTU.*us$', otu2), otu2, otu))
rownames(allUS) <- NULL
allUS <- (allUS)[,-c(31, 32)]

coreUS <- as.data.frame(otu_rare_rhizo) %>%
  mutate(otu=rownames(otu_rare_rhizo)) %>%
  filter(otu %in% US_occ_1) %>%
  mutate(otu2=if_else((otu %in% matchBLAST), otu, paste(otu,'us', sep='.')),
         otu3=if_else(grepl('^OTU.*us$', otu2), otu2, otu)) 
rownames(coreUS) <- coreUS$otu3
coreUS <- (coreUS)[,-c(31, 32, 33)] 
coreUS_column<- data.frame((1*(rowSums(coreUS)>0)))
coreUS_column$otu3 <- rownames(coreUS_column)
names(coreUS_column)[1] <- 'UScore'
rownames(coreUS_column) <- NULL

allCOL <- as.data.frame(otu_col_rare_rhizo) %>%
  mutate(otu=rownames(otu_col_rare_rhizo),
         otu2=if_else((otu %in% matchBLAST), otu, paste(otu,'col', sep='.')),
         otu3=if_else(grepl('^OTU.*col$', otu2), otu2, otu))
rownames(allCOL) <- NULL
allCOL <- (allCOL)[,-c(33, 34)]

coreCOL <- as.data.frame(otu_col_rare_rhizo) %>%
  mutate(otu=rownames(otu_col_rare_rhizo)) %>%
  filter(otu %in% COL_occ_1) %>%
  mutate(otu2=if_else((otu %in% matchBLAST), otu, paste(otu,'col', sep='.')),
         otu3=if_else(grepl('^OTU.*col$', otu2), otu2, otu)) 
rownames(coreCOL) <- coreCOL$otu3
coreCOL <- (coreCOL)[,-c(33, 34, 35)]

coreCOL_column<- data.frame((1*(rowSums(coreCOL)>0)))
coreCOL_column$otu3 <- rownames(coreCOL_column)
names(coreCOL_column)[1] <- 'COLcore'
rownames(coreCOL_column) <- NULL

allOTUs<- full_join(allUS, allCOL)
allOTUs<- full_join(allOTUs, coreCOL_column)
allOTUs<- full_join(allOTUs, coreUS_column)

allOTUs[is.na(allOTUs)] <- 0
rownames(allOTUs) <- allOTUs$otu3
allOTUs$otu3 <- NULL
US_taxa_noCore <- allOTUs[,c(1:30)]
COL_taxa_noCore <- allOTUs[,c(31:62)]
US_taxa_onlyCore <- data.frame(otu=rownames(allOTUs), coreUS=allOTUs[,63])
COL_taxa_onlyCore <- data.frame(otu=rownames(allOTUs), coreCOL=allOTUs[,64])
cores <- cbind(US_taxa_onlyCore[,2],COL_taxa_onlyCore[,2])

US_taxa_noCore_venn <- 1*(rowSums(US_taxa_noCore)>0)
COL_taxa_noCore_venn <- 1*(rowSums(COL_taxa_noCore)>0)

venn_df <- cbind(US_taxa_noCore_venn,cores)
venn_df <- cbind(venn_df, COL_taxa_noCore_venn)
venn_df[rowSums(venn_df)==4,]

colnames(venn_df) <- c("US", "COL core", "US core", "COL")
venn_df <- venn_df[rowSums(venn_df)>0,]
v_count=vennCounts(venn_df)

########
# Fig 3B

vennDiagram(v_count) 
