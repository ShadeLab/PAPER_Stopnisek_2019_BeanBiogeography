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

global_abun_df<- data.frame(otu=rownames(rel.abun), rel.abun) %>%
  gather(ID, abun, -otu) %>%
  left_join(map) %>%
  filter(otu %in% global_core,
         !is.na(Timepoint)) %>%
  group_by(ID,Timepoint, Site, Compartment) %>%
  summarise(coreNo=sum(abun),
            timeNo=length(unique(ID)),
            rel=coreNo/timeNo) 

global_abun_plot <- ggplot(global_abun_df, aes(x=as.factor(Timepoint), y=rel*100, color=Compartment, fill=Compartment)) +
  geom_boxplot(position=position_dodge(1), outlier.shape = NA, alpha=.5)+  
  geom_point(position = position_jitterdodge())+
  scale_color_manual(values = c("#E69F00", '#009E73')) +
  scale_fill_manual(values = c("#E69F00", '#009E73')) +
  theme_classic()+
  facet_grid(~Site) +
  labs(y='Relative abundance (%)', x=NULL) +
  theme(strip.text.x = element_blank(), 
        legend.position = 'none',
        #axis.text.x = element_blank(),
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

us_core_plot <- ggplot(us_core_df, aes(x=as.factor(Timepoint), y=rel*100, color=Compartment, fill=Compartment)) +
  geom_boxplot(position=position_dodge(1), outlier.shape = NA, alpha=.5)+  
  geom_point(position = position_jitterdodge())+
  scale_color_manual(values = c("#E69F00", '#009E73')) +
  scale_fill_manual(values = c("#E69F00", '#009E73')) +
  theme_classic()+
  facet_grid(~Site) +
  labs(y='Relative abundance (%)', x='Timepoint') +
  theme(strip.text.x = element_blank(),
        legend.title=element_text(size=8), 
        legend.text=element_text(size=8),
        legend.position = c(.2,.8),
        legend.background = element_rect('transparent'),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  stat_compare_means(label = "p.signif")

#######
# Fig 4

ggarrange(ggarrange(global_abun_plot, us_core_plot, 
                    labels = c("A", "B"), ncol = 1, nrow = 2),
          global_occ_plot, labels=c("","C"),ncol=2, widths = c(.6,1)) 

########
# Fig S5

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
