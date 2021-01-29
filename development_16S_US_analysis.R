#Core taxa in the development study
library(RColorBrewer)
library(viridis)
library(metagMisc)
library(microbiome)
library(reshape2)
library(ggpubr)
library(broom)
library(dplyr)
library(phyloseq)
library(decontam)
library(vegan)
library(tidyverse)
library(lemon)
library(scales)
library(ape)
library(heatmap)

setwd('~/Documents/git/PAPER_Stopnisek_2019_BeanBiogeography/')

taxa_global <- read.table('Data/global_core.txt', sep='\t')
global_core<- as.character(taxa_global$x)
taxa_US <- read.table('Data/US_core.txt', sep='\t')
core_US_otus<- as.character(taxa_US$x)

otu=read.table('Data/OTU_table_denovo.txt', header=T, row.names = 1)
map=readRDS('Data/map.rds')
rownames(map) <- map$ID
tax=readRDS('Data/taxOTUS.rds')

# Removing OTU contaminants using decontam package
#creating a phyloseq object
OTU=otu_table(as.matrix(otu), taxa_are_rows = T)
TAX=tax_table(as.matrix(tax))
MAP=sample_data(map)

otuPhyloseq=phyloseq(OTU,MAP,TAX)

#Filtering mitochondria, Chloroplast and Unclassified taxa
otuPhyloseq_filt <- otuPhyloseq %>%
  subset_taxa(Kingdom != "Unassigned")

otuPhyloseq_filt <- otuPhyloseq_filt %>%
  subset_taxa((Family != "Mitochondria") | is.na(Family))

otuPhyloseq_filt <- otuPhyloseq_filt %>%
  subset_taxa((Class != "Chloroplast") | is.na(Class))

filtered_otus <- phyloseq_to_df(otuPhyloseq_filt, addtax=T, addtot=T)

#-------------------------------------------------------------------------
#Using decontam package to remove contamination by prevalence
sample_data(otuPhyloseq_filt)$is.neg <- sample_data(otuPhyloseq_filt)$Type == "control"

contamdf.prev <- isContaminant(otuPhyloseq_filt, method="prevalence", neg="is.neg", threshold=0.5)
keepOTU <- rownames(contamdf.prev[contamdf.prev$contaminant=='FALSE',])
omit.OTUs <- rownames(contamdf.prev[contamdf.prev$contaminant=='TRUE',])

otuPhyloseq_filt_decont <- prune_taxa(keepOTU,otuPhyloseq_filt)
otuPhyloseq_filt_decont_samp <- subset_samples(otuPhyloseq_filt_decont, Type != 'control')

map_final <- read.csv('Data/metadata_table.csv', header=T, row.names = 1)
otu_final <- phyloseq_to_df(otuPhyloseq_filt_decont_samp, addtax=F, addtot=T)
otu_final <- otu_final[complete.cases(otu_final),]

tax_final <- phyloseq_to_df(otuPhyloseq_filt_decont_samp, addtax=T, addtot=F)
tax_final <- tax_final[1:8]

rownames(otu_final) <- otu_final$OTU
otu_final$Total <- NULL
otu_final$OTU <- NULL

#ordering samples
otu_final <- otu_final[,order(colnames(otu_final))]
map_final=map_final[order(rownames(map_final)),]
tax_filt <- tax[rownames(tax) %in% rownames(otu_final),]

PowersoilMRC <- as.character(map$ID[map$Site=='MRC' & map$IsolationMethod=='PowerSoil'])

#-------------------------------------------------------------------------
#Rarefiying data to 15000 reads per sample 
set.seed(007)

OTU.rare <- t(rrarefy(t(otu_final), 15000)) 
OTU.rare <- OTU.rare[,colSums(OTU.rare)>14999]
map=map_final[map_final$ID %in% colnames(OTU.rare),]

OTU.rare <- OTU.rare[,order(colnames(OTU.rare))]
map=map[order(map$ID),]

colnames(OTU.rare)==map$ID

# Core taxa abundance
map <- map[!(map$ID %in% PowersoilMRC),]
OTU.rare <- OTU.rare[,!(colnames(OTU.rare) %in% PowersoilMRC)]
rel.abun <- decostand(OTU.rare, method = 'total', MARGIN = 2)
#sum(rownames(OTU.rare) %in% global_core) #48 out of 48 global core taxa are included in the rarefied data

# Creating rarefaction curves
curve_colors_rare <- rep("darkgreen", ncol(OTU.rare))
curve_colors_rare[map$Compartment=="rhizosphere"] <- "burlywood"
rarecurve(t(OTU.rare), step=1000, label = FALSE, col = curve_colors_rare)


global_abun_df<- data.frame(otu=rownames(rel.abun), rel.abun) %>%
  gather(ID, abun, -otu) %>%
  left_join(map) %>%
  filter(otu %in% global_core,
         !is.na(Timepoint)) %>%
  group_by(ID,Timepoint, Site, Compartment) %>%
  summarise(coreNo=sum(abun),
            timeNo=length(unique(ID)),
            rel=coreNo/timeNo) 

global_abun_plot <- ggplot(global_abun_df, aes(x=as.factor(Timepoint), y=rel, color=Compartment, fill=Compartment)) +
  geom_boxplot(position=position_dodge(1), outlier.shape = NA, alpha=.5)+  
  geom_point(position = position_jitterdodge())+
  scale_color_manual(values = c('#009E73',"#E69F00")) +
  scale_fill_manual(values = c('#009E73',"#E69F00")) +
  scale_x_discrete(labels=c("V2", "V5",'flowering','pod filling','senescence', 'dry')) +
  theme_classic()+
  ylim(0,.68) +
  facet_grid(~Site) +
  labs(y='Relative abundance (%)', x=NULL) +
  theme(strip.text.x = element_blank(), 
        legend.position = 'none',
        axis.text.x = element_text(angle=60, hjust = 1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  stat_compare_means(label = "p.signif")

t_test = global_abun_df %>%
  group_by(Timepoint, Site) %>%
  do(tidy(t.test(.$rel~.$Compartment)))

alpha_aov <- data.frame(global_abun_df) %>%
  group_by(Timepoint, Site) %>%
  do(tidy(aov(.$rel~.$Compartment)))

global_core_occupancy <- data.frame(otu=rownames(rel.abun), rel.abun) %>%
  gather(ID, abun, -otu) %>%
  left_join(map) %>%
  filter(otu %in% global_core,
         !is.na(Timepoint)) %>%
  group_by(otu,Timepoint, Site, Compartment) %>%
  mutate(occ=if_else(abun>0,1, 0)) %>%
  summarise(occ_sum=sum(occ),
            timeNo=length(unique(ID)),
            occ_end=occ_sum/timeNo,
            occupancy=if_else(occ_end==1, '1', '<1'))

tax_v2 <- tax_filt
tax_v2$otu <- rownames(tax_v2)

z_sc2=data.frame(otu=rownames(rel.abun), rel.abun) %>%
  gather(ID, abun, -otu) %>%
  left_join(map) %>%
  left_join(tax_v2) %>%
  filter(otu %in% global_core,
         !is.na(Timepoint)) %>%
  group_by(Site, Timepoint, Compartment, otu, Phylum) %>%
  summarise(relabun_mean=mean(abun))%>%
  group_by(Site, Compartment, otu, Phylum) %>%
  mutate(z_score=scale(relabun_mean)) %>%
  arrange(Phylum)

z_score_df <- left_join(z_sc2, global_core_occupancy)

z_score_df$otu<- factor(z_score_df$otu, levels=unique(z_score_df$otu), ordered=T)

taxa_clas <- z_score_df%>%
  left_join(tax_v2)%>%
  group_by(otu,Phylum, Family,Genus)%>%
  summarise(n_count=length(otu))

global_occ_plot=ggplot(z_score_df,aes(x=as.factor(Timepoint), y=otu, color=z_score)) +
  geom_point(size=3) +
  theme_classic()+
  scale_size_continuous(range=c(0,4))+
  labs(size='Occupancy', color="z-score", y='Global core OTUs', x='Timepoint')+
  #scale_color_viridis(option = "D", direction = -1) +
  scale_colour_gradient2(low = muted("red"), mid='grey80', high = muted("blue"))+
  facet_grid(~Site*Compartment) +
  theme(strip.text.x = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.position = 'bottom',
        legend.box = 'vertical')

us_core_df <- data.frame(otu=rownames(rel.abun), rel.abun) %>%
  gather(ID, abun, -otu) %>%
  left_join(map) %>%
  filter(otu %in% core_US_otus,
         !is.na(Timepoint),
         !otu %in% global_core) %>%
  group_by(ID,Timepoint, Site, Compartment) %>%
  summarise(coreNo=sum(abun),
            timeNo=length(unique(ID)),
            rel=coreNo/timeNo)

us_core_plot <- ggplot(us_core_df, aes(x=as.factor(Timepoint), y=rel, color=Compartment, fill=Compartment)) +
  geom_boxplot(position=position_dodge(1), outlier.shape = NA, alpha=.5)+  
  geom_point(position = position_jitterdodge())+
  scale_color_manual(values = c('#009E73',"#E69F00")) +
  scale_fill_manual(values = c('#009E73',"#E69F00")) +  
  scale_x_discrete(labels=c("V2", "V5",'flowering','pod filling','senescence', 'dry')) +
  theme_classic()+
  ylim(0,.68) +
  facet_grid(~Site) +
  labs(y='Relative abundance (%)', x=NULL) +
  theme(strip.text.x = element_blank(),
        axis.text.x = element_text(angle=60, hjust = 1),
        legend.title=element_text(size=8), 
        legend.text=element_text(size=8),
        legend.position = c(.2,.8),
        legend.background = element_rect('transparent'),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  stat_compare_means(label = "p.signif")

#######
# Fig 3

ggarrange(global_abun_plot, us_core_plot,labels = c("A", "B"), ncol = 2, nrow = 1)

########
# Fig S6

ggplot(global_core_occupancy,aes(x=as.factor(Timepoint), y=otu, color=occ_end, size=occ_end)) +
  geom_point() +
  theme_classic()+
  scale_size_continuous(range=c(0,4))+
  labs(size=NULL, color="Occupancy", y='Global core OTUs', x='Timepoint')+
  scale_color_viridis(option = "D", direction = -1) +
  facet_grid(~Site*Compartment) +
  theme(strip.text.x = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))

global_core_occupancy %>%
  filter(occ_end == 0)

global_core_occupancy %>%
  group_by(Site, Compartment, Timepoint, timeNo) %>%
  summarise(occ1=sum(occ_end==1),
            occ.6=sum(occ_end<.5)) %>%
  pivot_longer(values_to = "occ",
               cols=c(occ1, occ.6)) %>%
  ggplot(aes(x=Timepoint, y = occ, group=name, col= name)) +
  geom_line(col='grey70') +
  geom_point() +
  theme_pubr() +
  labs(title='Number of OTUs with occupancy 1',subtitle = "The total number of OTUs is 48", y="# OTUs") +
  geom_hline(yintercept = 48, color='grey60', size=.5, linetype='dashed')+
  facet_grid(~Site*Compartment) +
  scale_color_manual(name="Occupancy",labels=c("<60%", "100%"), values = c('#FF6C90', 'black')) +
  theme(legend.position = c(0.1, .5),
        legend.justification = c(0.1, .5))
  
#'  DESeq2 for detecting rhizoplane enriched OTUs
#'  For that we need first to remove all unpaired samples and samples that have low counts
otu_dev <- otu_final[,colnames(otu_final) %in% colnames(OTU.rare)]
otu_dev <- otu_dev[,order(colnames(otu_dev))]

rownames(map) <- map$ID
map_dev=map[rownames(map) %in% colnames(OTU.rare),]
map_dev=map_dev[order(rownames(map)),]

map_tmp <- as.data.frame(map_dev %>%
                           mutate(group=paste(Site, Genotype, Stage, Plot, Plant, sep='.')) %>%
                           filter(Type == 'sample')%>%
                           group_by(group) %>%
                           mutate(n=n_distinct(ID)) %>%
                           filter(n==2))

map_tmp %>%
  filter(Site == 'MRC') %>%
  group_by(Timepoint) %>%
  summarise(n=n_distinct(ID))

#' Producing 214 paired samples in the map file
#' Removing also from the OTU table
otu_dev_v2 <- otu_final[,colnames(otu_final) %in% map_tmp$ID]
otu_dev=otu_dev_v2
otu_dev <- otu_dev[rowSums(otu_dev)>0,]

tax_filt <- tax[rownames(tax) %in% rownames(otu_dev),]
rownames(map_tmp) <- map_tmp$ID

OTU=otu_table(as.matrix(otu_dev), taxa_are_rows = T)
TAX=tax_table(as.matrix(tax_filt))
MAP=sample_data(map_tmp)
otuPhyloseq=phyloseq(OTU,MAP,TAX)

#' Plants in flowering stage (V2) growig at MRC site 
physeq_mrc_V2 = subset_samples(otuPhyloseq, Site == 'MRC' & Stage == 'V2')
diagdds_bean = phyloseq_to_deseq2(physeq_mrc_V2, ~ Compartment)
diagdds_bean = DESeq2::DESeq(diagdds_bean, test="Wald", fitType="parametric")
res_bean = DESeq2::results(diagdds_bean, cooksCutoff = FALSE)

alpha = 0.05
sigtab_bean = res_bean[which(res_bean$padj < alpha),]
sigtab_bean = cbind(as(sigtab_bean, "data.frame"), as(tax_table(physeq_mrc_V2)[rownames(sigtab_bean),], "matrix"))

rhizoplaneV2M <- sigtab_bean[rownames(sigtab_bean) %in% global_core,]
rhizoplaneV2M$site <- 'MRC'
rhizoplaneV2M$time <- 1
rhizoplaneV2M$OTU <- rownames(rhizoplaneV2M)
rownames(rhizoplaneV2M) <- NULL

#' Plants in flowering stage (V2) growig at SVERC site 
physeq_sverc_V2 = subset_samples(otuPhyloseq, Site == 'SVERC' & Stage == 'V2')
diagdds_bean = phyloseq_to_deseq2(physeq_sverc_V2, ~ Compartment)
diagdds_bean = DESeq2::DESeq(diagdds_bean, test="Wald", fitType="parametric")
res_bean = DESeq2::results(diagdds_bean, cooksCutoff = FALSE)

alpha = 0.05
sigtab_bean = res_bean[which(res_bean$padj < alpha), ]
sigtab_bean = cbind(as(sigtab_bean, "data.frame"), as(tax_table(physeq_sverc_V2)[rownames(sigtab_bean), ], "matrix"))

rhizoplaneV2S <- sigtab_bean[rownames(sigtab_bean) %in% global_core & sigtab_bean$log2FoldChange<0,]
rhizoplaneV2S$site <- 'SVERC'
rhizoplaneV2S$time <- 1
rhizoplaneV2S$OTU <- rownames(rhizoplaneV2S)
rownames(rhizoplaneV2S) <- NULL

#' Plants in flowering stage (V5) growig at MRC site 
physeq_mrc_V5 = subset_samples(otuPhyloseq, Site == 'MRC' & Stage == 'V5')
diagdds_bean = phyloseq_to_deseq2(physeq_mrc_V5, ~ Compartment)
diagdds_bean = DESeq2::DESeq(diagdds_bean, test="Wald", fitType="parametric")
res_bean = DESeq2::results(diagdds_bean, cooksCutoff = FALSE)

alpha = 0.05
sigtab_bean = res_bean[which(res_bean$padj < alpha) , ]
sigtab_bean = cbind(as(sigtab_bean, "data.frame"), as(tax_table(physeq_mrc_V5)[rownames(sigtab_bean), ], "matrix"))

rhizoplaneV5M <- sigtab_bean[rownames(sigtab_bean) %in% global_core,]
rhizoplaneV5M$site <- 'MRC'
rhizoplaneV5M$time <- 2
rhizoplaneV5M$OTU <- rownames(rhizoplaneV5M)
rownames(rhizoplaneV5M) <- NULL

#' Plants in flowering stage (V5) growig at SVERC site 
physeq_sverc_V5 = subset_samples(otuPhyloseq, Site == 'SVERC' & Stage == 'V5')
diagdds_bean = phyloseq_to_deseq2(physeq_sverc_V5, ~ Compartment)
diagdds_bean = DESeq2::DESeq(diagdds_bean, test="Wald", fitType="parametric")
res_bean = DESeq2::results(diagdds_bean, cooksCutoff = FALSE)

alpha = 0.05
sigtab_bean = res_bean[which(res_bean$padj < alpha), ]
sigtab_bean = cbind(as(sigtab_bean, "data.frame"), as(tax_table(physeq_sverc_V5)[rownames(sigtab_bean), ], "matrix"))

rhizoplaneV5S <- sigtab_bean[rownames(sigtab_bean) %in% global_core & sigtab_bean$log2FoldChange<0,]
rhizoplaneV5S$site <- 'SVERC'
rhizoplaneV5S$time <- 2
rhizoplaneV5S$OTU <- rownames(rhizoplaneV5S)
rownames(rhizoplaneV5S) <- NULL

#' Plants in flowering stage (Timepoint 3) growig at MRC site 
physeq_mrc_flow = subset_samples(otuPhyloseq, Site == 'MRC' & Timepoint == 3)
diagdds_bean = phyloseq_to_deseq2(physeq_mrc_flow, ~ Compartment)
diagdds_bean = DESeq2::DESeq(diagdds_bean, test="Wald", fitType="parametric")
res_bean = DESeq2::results(diagdds_bean, cooksCutoff = FALSE)

alpha = 0.05
sigtab_bean = res_bean[which(res_bean$padj < alpha) , ]
sigtab_bean = cbind(as(sigtab_bean, "data.frame"), as(tax_table(physeq_mrc_flow)[rownames(sigtab_bean), ], "matrix"))

rhizoplaneFlowM <- sigtab_bean[rownames(sigtab_bean) %in% global_core,]
rhizoplaneFlowM$site <- 'MRC'
rhizoplaneFlowM$time <- 3
rhizoplaneFlowM$OTU <- rownames(rhizoplaneFlowM)
rownames(rhizoplaneFlowM) <- NULL

#' Plants in flowering stage (Timepoint 3) growig at SVERC site 
physeq_sverc_flow = subset_samples(otuPhyloseq, Site == 'SVERC' & Timepoint == 3)
diagdds_bean = phyloseq_to_deseq2(physeq_sverc_flow, ~ Compartment)
diagdds_bean = DESeq2::DESeq(diagdds_bean, test="Wald", fitType="parametric")
res_bean = DESeq2::results(diagdds_bean, cooksCutoff = FALSE)

alpha = 0.05
sigtab_bean = res_bean[which(res_bean$padj < alpha), ]
sigtab_bean = cbind(as(sigtab_bean, "data.frame"), as(tax_table(physeq_sverc_flow)[rownames(sigtab_bean), ], "matrix"))

rhizoplaneFlowS <- sigtab_bean[rownames(sigtab_bean) %in% global_core & sigtab_bean$log2FoldChange<0,]
rhizoplaneFlowS$site <- 'SVERC'
rhizoplaneFlowS$time <- 3
rhizoplaneFlowS$OTU <- rownames(rhizoplaneFlowS)
rownames(rhizoplaneFlowS) <- NULL

#' Plants in pod pod filling (Timepoint 4) growig at MRC site 
physeq_mrc_pod = subset_samples(otuPhyloseq, Site== 'MRC' & Timepoint == 4)
diagdds_bean = phyloseq_to_deseq2(physeq_mrc_pod, ~ Compartment)
diagdds_bean = DESeq2::DESeq(diagdds_bean, test="Wald", fitType="parametric")
res_bean = DESeq2::results(diagdds_bean, cooksCutoff = FALSE)

alpha = 0.05
sigtab_bean = res_bean[which(res_bean$padj < alpha) , ]
sigtab_bean = cbind(as(sigtab_bean, "data.frame"), as(tax_table(physeq_mrc_pod)[rownames(sigtab_bean), ], "matrix"))

rhizoplanePodM <- sigtab_bean[rownames(sigtab_bean) %in% global_core,]
rhizoplanePodM$site <- 'MRC'
rhizoplanePodM$time <- 4
rhizoplanePodM$OTU <- rownames(rhizoplanePodM)
rownames(rhizoplanePodM) <- NULL

#' Plants in pod filling stage (Timepoint 4) growig at SVERC site 
physeq_sverc_pod = subset_samples(otuPhyloseq, Site== 'SVERC' & Timepoint == 4)
diagdds_bean = phyloseq_to_deseq2(physeq_sverc_pod, ~ Compartment)
diagdds_bean = DESeq2::DESeq(diagdds_bean, test="Wald", fitType="parametric")
res_bean = DESeq2::results(diagdds_bean, cooksCutoff = FALSE)

alpha = 0.05
sigtab_bean = res_bean[which(res_bean$padj < alpha), ]
sigtab_bean = cbind(as(sigtab_bean, "data.frame"), as(tax_table(physeq_sverc_pod)[rownames(sigtab_bean), ], "matrix"))

rhizoplanePodS <- sigtab_bean[rownames(sigtab_bean) %in% global_core,]
rhizoplanePodS$site <- 'SVERC'
rhizoplanePodS$time <- 4
rhizoplanePodS$OTU <- rownames(rhizoplanePodS)
rownames(rhizoplanePodS) <- NULL

#' Plants in senescenc (Timepoint 5 & 6) growig at MRC site 
physeq_mrc_senesc = subset_samples(otuPhyloseq, Site== 'MRC' & Timepoint == 5)
diagdds_bean = phyloseq_to_deseq2(physeq_mrc_senesc, ~ Compartment)
diagdds_bean = DESeq2::DESeq(diagdds_bean, test="Wald", fitType="parametric")
res_bean = DESeq2::results(diagdds_bean, cooksCutoff = FALSE)

alpha = 0.05
sigtab_bean = res_bean[which(res_bean$padj < alpha) , ]
sigtab_bean = cbind(as(sigtab_bean, "data.frame"), as(tax_table(physeq_mrc_senesc)[rownames(sigtab_bean), ], "matrix"))

rhizoplaneSenescM <- sigtab_bean[rownames(sigtab_bean) %in% global_core,]
rhizoplaneSenescM$time <- 5
rhizoplaneSenescM$site <- 'MRC'
rhizoplaneSenescM$OTU <- rownames(rhizoplaneSenescM)
rownames(rhizoplaneSenescM) <- NULL

physeq_mrc_dry = subset_samples(otuPhyloseq, Site== 'MRC' & Timepoint == 6)
diagdds_bean = phyloseq_to_deseq2(physeq_mrc_dry, ~ Compartment)
diagdds_bean = DESeq2::DESeq(diagdds_bean, test="Wald", fitType="parametric")
res_bean = DESeq2::results(diagdds_bean, cooksCutoff = FALSE)

alpha = 0.05
sigtab_bean = res_bean[which(res_bean$padj < alpha) , ]
sigtab_bean = cbind(as(sigtab_bean, "data.frame"), as(tax_table(physeq_mrc_dry)[rownames(sigtab_bean), ], "matrix"))

rhizoplaneDryM <- sigtab_bean[rownames(sigtab_bean) %in% global_core,]
rhizoplaneDryM$time <- 6
rhizoplaneDryM$site <- 'MRC'
rhizoplaneDryM$OTU <- rownames(rhizoplaneDryM)
rownames(rhizoplaneDryM) <- NULL

#' Plants in senescenc (Timepoint 5 & 6) growig at SVERC site 
physeq_sverc_senesc = subset_samples(otuPhyloseq, Site== 'SVERC' & Timepoint %in% c(5,6))
diagdds_bean = phyloseq_to_deseq2(physeq_sverc_senesc, ~ Compartment)
diagdds_bean = DESeq2::DESeq(diagdds_bean, test="Wald", fitType="parametric")
res_bean = DESeq2::results(diagdds_bean, cooksCutoff = FALSE)

alpha = 0.05
sigtab_bean = res_bean[which(res_bean$padj < alpha), ]
sigtab_bean = cbind(as(sigtab_bean, "data.frame"), as(tax_table(physeq_sverc_senesc)[rownames(sigtab_bean), ], "matrix"))

rhizoplaneSenescS <- sigtab_bean[rownames(sigtab_bean) %in% global_core & sigtab_bean$log2FoldChange<0,]
rhizoplaneSenescS$time <- 6
rhizoplaneSenescS$site <- "SVERC"
rhizoplaneSenescS$OTU <- rownames(rhizoplaneSenescS)
rownames(rhizoplaneSenescS) <- NULL

DESeq_results <- rbind(rhizoplaneV2M, rhizoplaneV2S, rhizoplaneV5M, rhizoplaneV5S, rhizoplaneFlowM, rhizoplaneFlowS, rhizoplanePodM, rhizoplanePodS, rhizoplaneSenescM, rhizoplaneDryM, rhizoplaneSenescS)

DESeq_results %>%
  group_by(site, time) %>%
  summarise(n_distinct(OTU)) %>%
  arrange(-desc(time))

DESeq_results_MRC <- rbind(rhizoplaneV2M, rhizoplaneV5M, rhizoplaneFlowM, rhizoplanePodM, rhizoplaneSenescM, rhizoplaneDryM)

DESeq2_MRC_tmp = DESeq_results_MRC %>%
  mutate(DESeq2_predict=if_else(log2FoldChange>0, 1, -1) ) %>%
  mutate(time=paste("timepoint_", time, sep=''))%>%
  select(OTU, time, log2FoldChange) %>%
  pivot_wider(names_from = time, values_from=log2FoldChange, values_fill = 0)

data.frame(OTU=rownames(rel.abun), rel.abun) %>%
  gather(ID, abun, -OTU) %>%
  left_join(map) %>%
  filter(OTU %in% rhizoplaneV2M$OTU,
         Site == 'MRC',
         Type == "sample",
         Timepoint==1) %>%
  left_join(rhizoplaneV2M)%>%
  mutate(predict=if_else(log2FoldChange>0, 'up', 'down')) %>%
  ggplot(aes(x=Compartment,y=abun))+ 
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter() +
  facet_wrap(predict~OTU, scales = 'free_y')

#' Sloan neutral model
#' Creating seperate datsets for rhizopshere and rhizoplane
map_rhizo = map[map$Compartment == 'rhizosphere',]
map_plane= map[map$Compartment == 'rhizoplane',]

OTU.rare.rhizo=OTU.rare[,colnames(OTU.rare) %in% as.character(map_rhizo$ID)]
OTU.rare.rhizo=OTU.rare.rhizo[rowSums(OTU.rare.rhizo)>0,]
OTU.rare.plane=OTU.rare[,colnames(OTU.rare) %in% as.character(map_plane$ID)]
OTU.rare.plane=OTU.rare.plane[rowSums(OTU.rare.plane)>0,]

#' Rhizosplane dataset:
spp=t(OTU.rare.plane)

#' Models for the whole community
obs.np=sncm.fit(spp, taxon=F, stats=FALSE, pool=NULL)
sta.np=sncm.fit(spp, taxon=F, stats=TRUE, pool=NULL)

above.pred.plane=sum(obs.np$freq > (obs.np$pred.upr), na.rm=TRUE)/sta.np$Richness
below.pred.plane=sum(obs.np$freq < (obs.np$pred.lwr), na.rm=TRUE)/sta.np$Richness

ap = obs.np$freq > (obs.np$pred.upr)
bp = obs.np$freq < (obs.np$pred.lwr)

Predict <- obs.np
Predict$otu <- rownames(Predict)

Predic_Plane=Predict %>%
  filter(otu %in% global_core) %>%
  left_join(tax_v2) %>%
  mutate(prediction_plane=if_else(freq<pred.lwr, "below", 'neutral'),
         prediction_plane=if_else(freq>pred.upr, "above", prediction_plane)) %>%
  select(-c(p, bino.pred, bino.lwr, bino.upr, y)) %>%
  rename(freq.plane=freq)%>%
  rename(freq.pred.plane=freq.pred) %>%
  rename(pred.lwr.plane=pred.lwr) %>%
  rename(pred.upr.plane=pred.upr)

#' Rhizosphere dataset:
spp=t(OTU.rare.rhizo)

#Models for the whole community
obs.np=sncm.fit(spp, taxon=F, stats=FALSE, pool=NULL)
sta.np=sncm.fit(spp, taxon=F, stats=TRUE, pool=NULL)

above.pred.plane=sum(obs.np$freq > (obs.np$pred.upr), na.rm=TRUE)/sta.np$Richness
below.pred.plane=sum(obs.np$freq < (obs.np$pred.lwr), na.rm=TRUE)/sta.np$Richness

ap = obs.np$freq > (obs.np$pred.upr)
bp = obs.np$freq < (obs.np$pred.lwr)

Predict <- obs.np
Predict$otu <- rownames(Predict)

Predic_Rhizo=Predict %>%
  filter(otu %in% global_core) %>%
  left_join(tax_v2) %>%
  mutate(prediction_rhizo=if_else(freq>pred.upr, "above", 'neutral'),
         prediction_rhizo=if_else(freq<pred.lwr, "below", prediction_rhizo)) %>%
  select(-c(p, bino.pred, bino.lwr, bino.upr, y)) %>%
  rename(freq.rhizo=freq)%>%
  rename(freq.pred.rhizo=freq.pred) %>%
  rename(pred.lwr.rhizo=pred.lwr) %>%
  rename(pred.upr.rhizo=pred.upr)

Development_Model <- left_join(Predic_Rhizo, Predic_Plane)
Development_Model <- left_join(Development_Model, Predic_Biogeo)
write_delim(Development_Model, path = '~/Desktop/CoreTaxaModel.txt')

#' **********************
#' Model considering only MRC sites due to usage of same DNA isolation methods
#' **********************

map_plane_mrc= map[map$Compartment == 'rhizoplane' &  map$Site == 'MRC',]
OTU.rare.plane.mrc=OTU.rare[,colnames(OTU.rare) %in% as.character(map_plane_mrc$ID)]
OTU.rare.plane.mrc=OTU.rare.plane.mrc[rowSums(OTU.rare.plane.mrc)>0,]

spp=t(OTU.rare.plane.mrc)
obs.np=sncm.fit(spp, taxon=F, stats=FALSE, pool=NULL)
sta.np=sncm.fit(spp, taxon=F, stats=TRUE, pool=NULL)

above.pred.plane.mrc=sum(obs.np$freq > (obs.np$pred.upr), na.rm=TRUE)/sta.np$Richness
below.pred.plane.mrc=sum(obs.np$freq < (obs.np$pred.lwr), na.rm=TRUE)/sta.np$Richness

ap = obs.np$freq > (obs.np$pred.upr)
bp = obs.np$freq < (obs.np$pred.lwr)

Predict <- obs.np
Predict$otu <- rownames(Predict)

Predic_Plane_MRC =Predict %>%
  filter(otu %in% global_core) %>%
  left_join(tax_v2) %>%
  mutate(prediction_plane=if_else(freq>pred.upr, "above", 'neutral'),
         prediction_plane=if_else(freq<pred.lwr, "below", prediction_plane)) %>%
  select(-c(p, bino.pred, bino.lwr, bino.upr, y)) %>%
  rename(freq.plane=freq)%>%
  rename(freq.pred.plane=freq.pred) %>%
  rename(pred.lwr.plane=pred.lwr) %>%
  rename(pred.upr.plane=pred.upr)


map_rhizo_mrc= map[map$Compartment == 'rhizosphere' &  map$Site == 'MRC',]
OTU.rare.rhizo.mrc=OTU.rare[,colnames(OTU.rare) %in% as.character(map_rhizo_mrc$ID)]
OTU.rare.rhizo.mrc=OTU.rare.rhizo.mrc[rowSums(OTU.rare.rhizo.mrc)>0,]

spp=t(OTU.rare.rhizo.mrc)
obs.np=sncm.fit(spp, taxon=F, stats=FALSE, pool=NULL)
sta.np=sncm.fit(spp, taxon=F, stats=TRUE, pool=NULL)

above.pred.rhizo.mrc=sum(obs.np$freq > (obs.np$pred.upr), na.rm=TRUE)/sta.np$Richness
below.pred.rhizo.mrc=sum(obs.np$freq < (obs.np$pred.lwr), na.rm=TRUE)/sta.np$Richness

ap = obs.np$freq > (obs.np$pred.upr)
bp = obs.np$freq < (obs.np$pred.lwr)

Predict <- obs.np
Predict$otu <- rownames(Predict)

Predic_Rhizo_MRC =Predict %>%
  filter(otu %in% global_core) %>%
  left_join(tax_v2) %>%
  mutate(prediction_rhizo=if_else(freq>pred.upr, "above", 'neutral'),
         prediction_rhizo=if_else(freq<pred.lwr, "below", prediction_rhizo)) %>%
  select(-c(p, bino.pred, bino.lwr, bino.upr, y)) %>%
  rename(freq.rhizo=freq)%>%
  rename(freq.pred.rhizo=freq.pred) %>%
  rename(pred.lwr.rhizo=pred.lwr) %>%
  rename(pred.upr.rhizo=pred.upr)

CoreTaxaModel_MRC = left_join(Predic_Rhizo_MRC, Predic_Plane_MRC) %>%
  left_join(Predic_Biogeo) %>%
  select(otu, Kingdom, Phylum, Class, Order, Family, Genus, Species, prediction_biogeo,prediction_rhizo,prediction_plane)

TableModelDESeq_MRC_tmp=full_join(CoreTaxaModel_MRC,data.frame(DESeq2_MRC_tmp), by=c('otu'='OTU'))

write_delim(CoreTaxaModel_MRC, path='~/Desktop/CoreTaxaModel_MRC.txt')
write_delim(TableModelDESeq_MRC_tmp, path='~/Desktop/CoreTaxaModelandDESeq2_MRC_v3.txt')

long_df <- TableModelDESeq_MRC_tmp %>%
  pivot_longer(timepoint_1:timepoint_6, names_to = "time", values_to = "DESeq") %>%
  mutate(DESeq=if_else(is.na(DESeq), 0, DESeq)) %>%
  arrange(desc(Phylum))

long_df$otu<- factor(long_df$otu, levels=unique(long_df$otu), ordered=T)

deseqPlot<- ggplot(long_df, aes(x=time, y=otu)) +
  geom_tile(aes(fill = DESeq), colour = "white") + 
  scale_fill_gradient2(low = "#E69F00", high = "#009E73", mid='white')+ 
  scale_x_discrete(expand = c(0, 0), labels=c("V2", "V5", "flowering", "pod filling", "senescence", 'dry')) +
  theme(axis.text.x = element_text(angle=60, hjust = 1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.title.x = element_blank())

tax_colors <-  c('Acidobacteria'='#3cb44b','Actinobacteria'='#a9a9a9' ,'Bacteroidetes'='#aaffc3',
                 'BRC1'='#4363d8','Chloroflexi'='#f58231','Euryarchaeota'='#911eb4','FBP'='#42d4f4',
                 'Firmicutes'='#f032e6','Gemmatimonadetes'='#bfef45','Nitrospirae'='#469990', 
                 'Planctomycetes'='#e6beff','Proteobacteria'= '#e6194B','Thaumarchaeota'='#808000','Verrucomicrobia'='#000075')

phyloPlot<- ggplot(long_df, aes(x=1,y=otu, fill=Phylum)) +
  geom_tile(colour = "white", ) + 
  scale_fill_manual(values = tax_colors)+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))

long_df_v2 <- TableModelDESeq_MRC_tmp %>%
  pivot_longer(prediction_biogeo:prediction_plane, names_to = "prediction", values_to = "model_fit") %>%
  arrange(desc(Phylum))

long_df_v2$otu<- factor(long_df_v2$otu, levels=unique(long_df_v2$otu), ordered=T)

predictionPlot<- ggplot(long_df_v2, aes(x=prediction, y=otu)) +
  geom_tile(aes(fill = as.factor(model_fit)), colour = "white") + 
  scale_fill_manual(values = c("dodgerblue3", "firebrick4",'white')) + 
  scale_x_discrete(expand = c(0, 0), labels=c("Biogeography\n(n=30)", "Rhizoplane\nMRC (n=62)",'Rhizosphere\nMRC (n=63)')) +
  theme(axis.text.x = element_text(angle=60, hjust = 1),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.ticks.y = element_blank())+
  labs(fill='Model fit')

ggarrange(deseqPlot, predictionPlot,phyloPlot, ncol = 3, widths = c(1,.6,.6)) 

install.packages('patchwork')
library(patchwork)
deseqPlot + predictionPlot + phyloPlot



CoreTaxaModel_MRC %>%
  filter(otu %in% c('JF429082.1.1514', 'FM209319.1.1474', 'GU550579.1.1244', 'JF429082.1.1514', 'FM209322.1.1496','FR753119.1.1437',
                    'GU940705.1.1318', 'FR853451.1.1447', 'FN546860.1.1411', 'HG423553.1.1283', 'FQ658643.1.1351','GQ389063.1.1478',
                    'FR675947.1.1459', 'GQ115602.1.1307'))

rhizoCoreMatrix <- OTU.rare.rhizo.mrc[rownames(OTU.rare.rhizo.mrc) %in% global_core,]
planeCoreMatrix <- OTU.rare.plane.mrc[rownames(OTU.rare.plane.mrc) %in% global_core,]


rhizoCoreOcc <- 1*(rowSums(rhizoCoreMatrix)>0)

