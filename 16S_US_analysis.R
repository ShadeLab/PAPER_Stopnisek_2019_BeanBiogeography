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

tax_colors <-  c('Acidobacteria'='#3cb44b','Actinobacteria'='#a9a9a9' ,'Bacteroidetes'='#aaffc3',
                 'BRC1'='#4363d8','Chloroflexi'='#f58231','Euryarchaeota'='#911eb4','FBP'='#42d4f4',
                 'Firmicutes'='#f032e6','Gemmatimonadetes'='#bfef45','Nitrospirae'='#469990', 
                 'Planctomycetes'='#e6beff','Proteobacteria'= '#e6194B','Thaumarchaeota'='#808000','Verrucomicrobia'='#000075')

setwd('~/Documents/git/PAPER_Stopnisek_2019_BeanBiogeography/')
######################################################
# Dataset preparation

otu_us <- read.table('Data/OTU_table_US.txt', sep='\t', row.names = 1, header=T)
map_combined <- read.table('Data/map.txt', row.names = 1, sep='\t', header=T)

#remove lines where taxonomy includes mitochondria
otu_us <- otu_us[!grepl("Mitochondria", otu_us$taxonomy),]
otu_us <- otu_us[!grepl("Chloroplast", otu_us$taxonomy),]

tax_otu <- otu_us['taxonomy']
tax_otu <- tax_otu %>%
  separate(taxonomy, into=c("Kingdom", "Phylum", "Class", 
                            "Order", "Family", "Genus", "Species"), sep="; ", remove=F)
tax_otu[2:8] <- lapply(tax_otu[2:8], function(x) gsub(".*__", "", x))

otu_us['taxonomy'] <- NULL
otu_us['SVERCc1'] <- NULL
map_combined <- map_combined[-1,] 
# Order the samples
otu_us <- otu_us[,order(colnames(otu_us))]
# Order the samples of the map the same way
map_combined=map_combined[order(rownames(map_combined)),]
# Check to make sure they all match with each other
rownames(map_combined)==colnames(otu_us)

map_combined$sample_ID <- rownames(map_combined)


#########
# Fig S1A

data.frame(otu=as.factor(rownames(otu_us)),otu_us) %>%
  gather(sample_ID, abun, -otu) %>%
  left_join(map_combined, by='sample_ID') %>%
  group_by(sample_ID) %>%
  summarise(readCount=sum(abun)) %>%
  arrange(desc(readCount)) %>%
  ggplot(aes(x=sample_ID,y=readCount)) +
  geom_bar(stat='identity')+
  geom_hline(yintercept=min(colSums(otu_us)), linetype="dashed", color = "red") +
  theme_bw()+
  xlab('Sample name') +
  ylab('read number') +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))

#########
# Fig S1B
#Rarefying data to the sample with lowest reads (31255)
rarecurve(t(otu_us), step=1000, label = FALSE, xlim=c(0,31255))

######################################################
# Rarefying data
# whole dataset
set.seed(3)
otu_us_rare <- t(rrarefy(t(otu_us), min(colSums(otu_us)))) 
otu_us_rare <- otu_us_rare[rowSums(otu_us_rare)>0,]
otu_rare <- otu_us_rare

# rhizosphere only
otu_us_rhizo <- otu_us[,map_combined$soil=='rhizosphere']
set.seed(33)
otu_us_rhizo_rare <- t(rrarefy(t(otu_us_rhizo), min(colSums(otu_us_rhizo)))) 

######################################################
# Alpha Diversity

#calculating indices:
#richness
s <- specnumber(otu_us_rare,MARGIN=2)
# shannon
h <- vegan::diversity(t(otu_us_rare), "shannon")
#pielou
pielou=h/log(s)

map.div <- map_combined
map.div$Richness <- s
map.div$Shannon <- h
map.div$Pielou <- pielou 

map.alpha <- melt(map.div, id.vars=c("sample_name","bean", "plot", "state", 'soil', 'site', 'irrigation', 'fertilization', 'pH', 'OM', 'Nitrogen', 'P'), 
                  measure.vars=c("Richness", "Shannon", "Pielou"))

map.alpha %>%
  filter(variable=='Richness') %>%
  group_by(variable, site) %>%
  summarise(mean_var=mean(value))

ggplot(map.alpha[map.alpha$soil=='rhizosphere',], aes(y=value,x=bean, color=bean))+
  scale_color_manual(values = c('darkorange', 'black')) + 
  facet_wrap( ~ variable, scales = "free_y") + 
  theme_bw()+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))+
  geom_boxplot() 

ggplot(map.alpha[map.alpha$soil=='rhizosphere',], aes(y=value, x=site))+
  facet_wrap( ~ variable, scales = "free_y") + 
  theme_bw()+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))+
  geom_boxplot() 

t.test(value~bean, data=map.alpha[map.alpha$variable=="Richness" & map.alpha$soil=='rhizosphere',])
tidy(aov(value~site, data=map.alpha[map.alpha$variable=="Richness",]))

t_test_rich = map.alpha %>%
  filter(soil=='rhizosphere') %>%
  group_by(variable) %>%
  do(tidy(t.test(.$value~.$bean)))

alpha_aov <- data.frame(map.alpha) %>%
  group_by(variable) %>%
  do(tidy(aov(.$value~.$pH + .$soil + .$bean + .$fertilization + .$irrigation + .$site + .$bean*.$site + .$pH*.$site)))

bartlett.test(value~site, data = map.alpha[map.alpha$variable=='Richness' & map.alpha$soil=='rhizosphere',])
bartlett.test(value~fertilization, data = map.alpha[map.alpha$variable=='Richness'  & map.alpha$soil=='rhizosphere',])
bartlett.test(value~bean, data = map.alpha[map.alpha$variable=='Richness'  & map.alpha$soil=='rhizosphere',])
bartlett.test(value~soil, data = map.alpha[map.alpha$variable=='Richness'  & map.alpha$soil=='rhizosphere',])

######################################################
#Beta diversity
#PCoA
bean.dist <- vegdist(t(otu_us_rare), method="bray")
bean.pcoa <- cmdscale(bean.dist, eig=TRUE)

map_short <- map_combined[,-c(1:23,26)]
bean_envfit <- envfit(bean.pcoa, map_short)

ax1.bean <- bean.pcoa$eig[1]/sum(bean.pcoa$eig)
ax2.bean <- bean.pcoa$eig[2]/sum(bean.pcoa$eig)

ax.bean=ax1.bean+ax2.bean

bean_plot_colors <- rep("#e41a1c", nrow(map_combined))
bean_plot_colors[map_combined$site=="WA"] <- "#377eb8"
bean_plot_colors[map_combined$site=="CO"] <- "#4daf4a"
bean_plot_colors[map_combined$site=="NE"] <- "#984ea3"
bean_plot_colors[map_combined$site=="SVREC"] <- "darkorange"

bean_plot_shapes <- rep(21, nrow(map_combined))
bean_plot_shapes[map_combined$bean=="CELRK"] <- 23
bean_plot_shapes[map_combined$soil=="bulk"] <- 22

#siteVSsoil PCoA
plot(bean.pcoa$points[,1], bean.pcoa$points[,2], cex=1.5, bg=bean_plot_colors, pch=bean_plot_shapes,
     xlab=paste("PCoA1: ",100*round(ax1.bean,3),"% var. explained",sep=""), 
     ylab=paste("PCoA2: ",100*round(ax2.bean,3),"% var. explained",sep=""))
plot(bean_envfit, p=0.001, col='black')
legend(.17, -.08, legend=c('CO','MRF','SVREC','NE', 'WA'), col=c("#4daf4a","#e41a1c", "darkorange","#984ea3","#377eb8"), 
       box.lty=0, cex= .8, pch = 20, title='Location')
legend(-.03, -.06, legend = c('CELRK', 'Eclipse', 'Root\nassociated'), col='black',pch=c(23,21, 22), cex=.8, box.col = 0, title='Genotype')

adonis(bean.dist~map_combined$bean, strata=map_combined$site)
adonis(bean.dist~map_combined$site, strata=map_combined$bean)
adonis(bean.dist~map_combined$site/map_combined$bean)
adonis(bean.dist~map_combined$irrigation) 
adonis(bean.dist~map_combined$fertilization)
adonis(bean.dist~map_combined$pH) 
adonis(bean.dist~map_combined$site)
adonis(bean.dist~map_combined$soil)
adonis(bean.dist~map_combined$bean)
adonis(bean.dist~map_combined$P)
adonis(bean.dist~map_combined$Nitrogen)
adonis(bean.dist~map_combined$OM)
adonis(bean.dist~map_combined$NO3)
adonis(bean.dist~map_combined$NH4)


######################################################
#Taxonomic distribution 
tax_v1 <- tax_otu
# creating column where last taxonomic identifier is attached to the OTU id 
features <- c(sprintf("OTU%05d", seq(1,28976)),"label")
tax_v1$otu_id <- features
tax_v1 <- tax_v1 %>% select(otu_id, everything())

lastValue <- function(x) tail(x[!is.na(x)], 1)
last_taxons<- apply(tax_v1, 1, lastValue)
tax_v1$last_taxon <- last_taxons
tax_v1$final_names <- paste(tax_v1$otu_id, tax_v1$last_taxon, sep='-')

tax_v1$otu <- rownames(tax_v1)

otu.rel.abun <- decostand(otu_us_rare, method="total", MARGIN=2) #calculating relative abundance
top.otus = as.data.frame(head(sort(rowSums(otu.rel.abun), decreasing = T),10))
abund.top.otus <- otu.rel.abun[rownames(otu.rel.abun) %in% rownames(top.otus),]

tmp_v3 <- data.frame(otu=as.factor(rownames(otu.rel.abun)),otu.rel.abun) %>%
  gather(sample_ID, abun, -otu) %>%
  left_join(map_combined[,c('sample_ID','pH', 'state', 'bean', 'plot', 'site', 'soil')], by='sample_ID') %>%
  left_join(tax_v1, by='otu') %>%
  group_by(Phylum,site) %>% #include bean if needed to split the data into genotypes
  summarise(
    n=sum(abun)/length(unique(sample_ID))) %>%
  arrange(desc(n))

temp_test <- data.frame(otu=as.factor(rownames(otu.rel.abun)),otu.rel.abun) %>%
  gather(sample_ID, abun, -otu) %>%
  left_join(map_combined[,c('sample_ID','pH', 'state', 'bean', 'plot', 'site', 'soil')], by='sample_ID') %>%
  left_join(tax_v1, by='otu') %>%
  filter(soil != 'bulk') %>%
  group_by(Phylum, sample_ID) %>%
  summarise(n=sum(abun),
            sample_n=length(unique(sample_ID)),
            rel_abun_phylum=n/sample_n) %>%
  mutate(
    tax_names=ifelse(rel_abun_phylum>=.01, Phylum, 'other'))
#filter(Phylum == 'Proteobacteria')
arrange(desc(rel_abun_phylum)) %>%
  group_by(Phylum)%>%
  summarise(min=min(n),
            max=max(n)) %>%
  arrange(desc(max))

###########
#Figure S3A

ggplot(temp_test,aes(sample_ID, n, fill=tax_names)) +
  geom_bar(stat='identity') +
  theme_bw()+
  labs(y="Relative abundace", x= NULL, fill='Phylum') +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=6),
        legend.position = 'bottom') +
  guides(fill=guide_legend(ncol=3)) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))

######################################################
# Identifying taxa differently abundant between the plant genotypes (DESeq2)
taxonomy_phyloseq<- tax_v1[rownames(tax_v1) %in% rownames(otu_us),]
dim(otu_us)
OTU = otu_table(otu_us, taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxonomy_phyloseq))
SAM = sample_data(map_combined)
taxa_names(OTU)
physeq <- merge_phyloseq(phyloseq(OTU, TAX), SAM)

map_rhizo_us <- map_combined[map_combined$soil!='bulk',]

physeq_rhizo = subset_samples(physeq, soil!= 'bulk')
diagdds_bean = phyloseq_to_deseq2(physeq_rhizo, ~ bean)
diagdds_bean = DESeq2::DESeq(diagdds_bean, test="Wald", fitType="parametric")
res_bean = DESeq2::results(diagdds_bean, cooksCutoff = FALSE)
alpha = 0.05

sigtab_bean = res_bean[which(res_bean$padj < alpha), ]
sigtab_bean = cbind(as(sigtab_bean, "data.frame"), as(tax_table(physeq_rhizo)[rownames(sigtab_bean), ], "matrix"))

sigtab_bean$Genus <- as.character(sigtab_bean$Genus)

sigtabgen_bean <- sigtab_bean %>%
  mutate(new_names = if_else((Genus == "uncultured bacterium" | Genus == "Ambiguous_taxa" | Genus == "uncultured"| is.na(Genus)), 'other', Genus))

######################################################
#Venn diagram
#make subset OTU tables
map_combined$state <- as.character(map_combined$state)
CO_otu <- otu_us_rare[,map_combined$location=="Fort_Collins_CO"& map_combined$soil=="rhizosphere"] 
MI_otu <- otu_us_rare[,map_combined$location=="Montcalm_MI"& map_combined$soil=="rhizosphere"]
NE_otu <- otu_us_rare[,map_combined$location=="Scotts_Bluff_County_NE"& map_combined$soil=="rhizosphere"]
WA_otu <- otu_us_rare[,map_combined$location=="Othello_WA"& map_combined$soil=="rhizosphere"]
SA_otu <- otu_us_rare[,map_combined$location=="Saginaw_MI" & map_combined$soil=="rhizosphere" ]

CELRK_otu <- otu_us_rare[,map_combined$bean=="CELRK"]
Eclipse_otu <- otu_us_rare[,map_combined$bean=="Eclipse"]

rhizo_otu <- otu_us_rare[,map_combined$soil =="rhizosphere"]
otu_us_rare_rhizo=rhizo_otu
CELRK_rhizo_otu <- otu_us_rare_rhizo[,map_rhizo_us$bean=="CELRK" & map_rhizo_us$soil=="rhizosphere"]
Eclipse_rhizo_otu <- otu_us_rare[,map_rhizo_us$bean=="Eclipse" & map_rhizo_us$soil=="rhizosphere"]
# make presence absence list from soil and plant into 1 & 0
CO_venn <- 1*(rowSums(CO_otu)>0)
MI_venn <- 1*(rowSums(MI_otu)>0)
NE_venn <- 1*(rowSums(NE_otu)>0)
WA_venn <- 1*(rowSums(WA_otu)>0)
SA_venn <- 1*(rowSums(SA_otu)>0)

# combine vectors into a matrix cbind = column bind (r bind = row bind)
venndata <- cbind(CO_venn,MI_venn)
venndata <- cbind(venndata,NE_venn)
venndata <- cbind(venndata,WA_venn)
venndata <- cbind(venndata,SA_venn)

#location venn 
colnames(venndata) <- c ("CO", "MRF", "NE", "WA", "SVREC")
venndata <- venndata[rowSums(venndata)>0,]
v=vennCounts(venndata)
v2=round(v[,"Counts"]/sum(v[,"Counts"]),2)

###########
#Figure S3A
Fig3B <- vennDiagram(v) 

#Percent shared by site
CO_present <- CO_otu[rowSums(CO_otu)>1,]
MI_present <- MI_otu[rowSums(MI_otu)>1,]
NE_present <- MI_otu[rowSums(NE_otu)>1,]
WA_present <- MI_otu[rowSums(WA_otu)>1,]
SA_present <- MI_otu[rowSums(SA_otu)>1,]

CO_share <- 2173/length(rownames(CO_present))
MI_share <- 2173/length(rownames(MI_present))
NE_share <- 2173/length(rownames(NE_present))
WA_share <- 2173/length(rownames(WA_present))
SA_share <- 2173/length(rownames(SA_present))

mean(c(CO_share,MI_share,NE_share, WA_share, SA_share))


######################################################
#Occupancy abundance relationship

#selecting variety unique OTUs
CELRK_rhizo_otu <- CELRK_rhizo_otu[rowSums(CELRK_rhizo_otu)>0,]
Eclipse_rhizo_otu <- Eclipse_rhizo_otu[rowSums(Eclipse_rhizo_otu)>0,]

CELRK_uniq <- CELRK_rhizo_otu[!(rownames(CELRK_rhizo_otu) %in% rownames(Eclipse_rhizo_otu)),]
Eclipse_uniq <- Eclipse_rhizo_otu[!(rownames(Eclipse_rhizo_otu) %in% rownames(CELRK_rhizo_otu)),]

#remove the bulk samples
otu_rare_rhizo <- otu_us_rare[,map_combined$soil=='rhizosphere']
otu_rare_rhizo <- otu_rare_rhizo[rowSums(otu_rare_rhizo)>0,]
bean_otu_PA_rhizo <- 1*((otu_rare_rhizo>0)==1)
bean_otu_PA_rhizo <- bean_otu_PA_rhizo[rowSums(bean_otu_PA_rhizo)>0,]
Occ <- rowSums(bean_otu_PA_rhizo)/ncol(bean_otu_PA_rhizo)

com_abund_bean_rhizo <- rowSums(otu_rare_rhizo)/ncol(otu_rare_rhizo)
rhizo_abund <- decostand(otu_rare_rhizo, method='total', MARGIN=2)
com_abund_bean_rhizo <- rowSums(rhizo_abund)/ncol(rhizo_abund)

bean_rhizo_df_occ <- data.frame(otu=names(Occ), occ=Occ) 
bean_rhizo_df_abun <- data.frame(otu=names(com_abund_bean_rhizo), abun=log10(com_abund_bean_rhizo))
occ_abun_bean_rhizo <- left_join(bean_rhizo_df_occ, bean_rhizo_df_abun)

occ_abun_bean_rhizo$bean_found <- 'shared'
occ_abun_bean_rhizo$bean_found[occ_abun_bean_rhizo$otu %in% rownames(CELRK_uniq)] <- 'CELRK'
occ_abun_bean_rhizo$bean_found[occ_abun_bean_rhizo$otu %in% rownames(Eclipse_uniq)] <- 'Eclipse'

size_occ<- data.frame(otu=as.factor(rownames(bean_otu_PA_rhizo)),bean_otu_PA_rhizo) %>%
  gather(sample_ID, abun_rhizo, -otu) %>%
  left_join(map_combined[,c('sample_ID','pH', 'site', 'bean', 'plot', 'soil')], by='sample_ID') %>%
  left_join(tax_v1, by='otu') %>%
  #filter(soil!='bulk') %>%
  group_by(otu, site, Family, Genus, final_names) %>%
  mutate(sum_abun=sum(abun_rhizo),
         n_rep=length(abun_rhizo),
         rel_occ=sum_abun/n_rep,
         is_present=1*((sum_abun>0)==1))

#calculating presence across sites
size_occ %>%
  group_by(otu, site) %>%
  dplyr::summarize(n=sum(is_present),
                   presence= 1*((n>0)==1)) %>%
  group_by(otu) %>%
  dplyr::summarize(
    total_presence=sum(presence)
  ) -> tmp_occ

combined_occ_data <- left_join(tmp_occ, occ_abun_bean_rhizo)

combined_occ_data%>%
  #filter(total_presence==1) %>%
  group_by(total_presence,bean_found) %>%
  summarise(counts=length(otu)) %>%
  group_by(total_presence) %>%
  mutate(tot_counts=sum(counts),
         rel_contribution=counts/tot_counts)

##*********************************
#Sloan neutral model 
spp=t(otu_rare_rhizo)
taxon=as.vector(size_occ$taxonomy)

#Models for the whole community
obs.np=sncm.fit(spp, taxon=F, stats=FALSE, pool=NULL)
sta.np=sncm.fit(spp, taxon, stats=TRUE, pool=NULL)
sta.np.16S <- sta.np

above.pred=sum(obs.np$freq > (obs.np$pred.upr), na.rm=TRUE)/sta.np$Richness
below.pred=sum(obs.np$freq < (obs.np$pred.lwr), na.rm=TRUE)/sta.np$Richness

ap = obs.np$freq > (obs.np$pred.upr)
bp = obs.np$freq < (obs.np$pred.lwr)

Predict <- obs.np
Predict$otu <- rownames(Predict)
UpPredictOTU <- unique(Predict$otu[ap==TRUE])
DownPedictedOTU <- unique(Predict$otu[bp==TRUE])

Predic_Biogeo=Predict %>%
  filter(otu %in% global_core) %>%
  left_join(tax_v2) %>%
  mutate(prediction_biogeo=if_else(freq>pred.upr, "above", 'neutral'),
         prediction_biogeo=if_else(freq<pred.lwr, "below", prediction_biogeo)) %>%
  select(-c(p, bino.pred, bino.lwr, bino.upr, y)) 

########
# Fig 2A

ggplot(data=combined_occ_data, aes(x=abun, y=occ)) +
    theme_bw()+
    geom_point(pch=21,  size=3, aes(fill=bean_found)) +
    geom_line(color='black', data=obs.np, size=2, aes(y=obs.np$freq.pred, x=log10(obs.np$p))) +
    geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.upr, x=log10(obs.np$p)))+
    geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.lwr, x=log10(obs.np$p)))+
    scale_fill_manual(aes(breaks = bean_found), values=c('darkorange','black', 'white')) +
    xlim(-6,-1.5)+
    labs(x=paste("log10(mean abundance)\n (n=",sta.np.16S$Richness," OTUs)", sep=''), y=paste("Occupancy (n=",sta.np.16S$Samples," samples)",sep=''), fill='Genotype') +
    theme(legend.position="top",
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    guides(fill = guide_legend(override.aes = list(alpha=1)),
           fill = guide_legend(title=NULL))


######################################################
#Core community 
length(Occ[Occ==1]) #258 with occupancy = 1
Occ_bact <- Occ
together_occ16S <- otu.rel.abun[rownames(otu.rel.abun) %in% occ_abun_bean_rhizo$otu[occ_abun_bean_rhizo$occ==1],] 

occ_16S <- data.frame(otu=as.factor(rownames(together_occ16S)),together_occ16S) %>%
  gather(sample_ID, abun, -otu) %>%
  left_join(map_combined[,c('sample_name','sample_ID','site','bean', 'soil')], by='sample_ID') %>%
  left_join(tax_v1, by='otu') %>%
  filter(soil != 'bulk') %>%
  group_by(sample_name,site, Phylum, Class, otu) %>%
  mutate(
    n=abun/length(abun))

core_species <- data.frame(data.frame(otu=as.factor(rownames(together_occ16S)),together_occ16S)) %>%
  gather(sample_ID, abun, -otu) %>%
  left_join(map_combined[c('sample_ID','pH', 'state', 'bean', 'sample_name','soil')], by='sample_ID') %>%
  left_join(tax_v1, by='otu') %>%
  filter(soil != 'bulk') %>%
  mutate(new_names = if_else((Order == "uncultured bacterium" | Order == "Ambiguous_taxa" | is.na(Order)), Phylum, Order),
         new_names = if_else((new_names == 'uncultured Acidobacteria bacterium'), 'Acidobacteria bacterium', new_names)) %>%
  arrange(desc(abun))

x = tapply(core_species$abun, core_species$new_names, function(x) mean(x))
x = sort(x, TRUE)
core_species$new_names = factor(as.character(core_species$new_names), levels=names(x))

########
# Fig 2C

US_core_taxa <- ggplot(core_species,aes(new_names, abun, fill=Phylum)) +
  theme_bw()+
  geom_boxplot() +
  scale_fill_manual(values=tax_colors)+
  labs(y='Relative abundance', x='Order') +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = c(.7,.7),
        legend.title = element_text(size=10),
        axis.text.y = element_text(size=8),
        legend.text = element_text(size=8),
        legend.key = element_rect(colour = 'transparent', fill = 'transparent', size = 0.2, linetype='dashed'))+
  coord_flip()

coreUStaxaCounts <- core_species %>%
  group_by(Phylum,new_names) %>%
  summarise(n_taxa=length(unique(otu))) %>%
  ggplot(aes(x=new_names, n_taxa, fill=Phylum)) +
  geom_bar(stat='identity') +
  ylab('Number of taxa') +
  xlab('') +
  scale_fill_manual(values=tax_colors)+
  theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = 'none',
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  coord_flip()

grid.arrange(US_core_taxa, coreUStaxaCounts, widths=c(2.5,.7))


#' Distance decay analysis 

df.sites=data.frame(name=c("CO", "NE", "MRF", "WA", "SVREC"),
                    lat=c(40.5,41.9, 43.3, 46.8, 43.4),
                    lon=c(-104.8,-103.8, -85.1, -121, -83.7))
site_distance_matrix <- round(GeoDistanceInMetresMatrix(df.sites))/1000

DISTdf <- data.frame(site1=as.factor(rownames(site_distance_matrix)),site_distance_matrix) %>%
  gather(site2, dist, -site1) %>%
  mutate(combined=paste(site1, site2, sep='-'))

BCdf <- data.frame(sample1=as.factor(rownames(as.matrix(bean.dist))),as.matrix(bean.dist)) %>%
  gather(sample2, BC, -sample1)

BCdf_its <- data.frame(sample1=as.factor(rownames(as.matrix(bean_its.dist))),as.matrix(bean_its.dist)) %>%
  gather(sample2, BC, -sample1)

BCdf=BCdf_its

BCdf$sample1=as.character(BCdf$sample1)
BCdf$sample2=as.character(BCdf$sample2)

BCdf$site1=ifelse(grepl("NE", BCdf$sample1), "NE", BCdf$sample1)
BCdf$site1=ifelse(grepl("SVERC", BCdf$sample1), "SVREC", BCdf$site1)
BCdf$site1=ifelse(grepl("MRF", BCdf$sample1), "MRF", BCdf$site1)
BCdf$site1=ifelse(grepl("WA", BCdf$sample1), "WA", BCdf$site1)
BCdf$site1=ifelse(grepl("CO", BCdf$sample1), "CO", BCdf$site1)

BCdf$site2=ifelse(grepl("NE", BCdf$sample2), "NE", BCdf$sample2)
BCdf$site2=ifelse(grepl("SVERC", BCdf$sample2), "SVREC", BCdf$site2)
BCdf$site2=ifelse(grepl("MRF", BCdf$sample2), "MRF", BCdf$site2)
BCdf$site2=ifelse(grepl("WA", BCdf$sample2), "WA", BCdf$site2)
BCdf$site2=ifelse(grepl("CO", BCdf$sample2), "CO", BCdf$site2)

BCdf$combined=paste(BCdf$site1, BCdf$site2, sep="-")

BCdistDF <- BCdf %>% 
  left_join(DISTdf, by='combined') %>%
  filter(dist!=0.000)
m <- lm(log(1-BCdistDF$BC) ~ log(BCdistDF$dist), BCdistDF)
summary(m)

BCdistPlot <- BCdf %>% 
  left_join(DISTdf, by='combined') %>%
  filter(dist!=0.000) %>%
  ggplot(aes(x=log(dist), y=log(1-BC))) +
  geom_point() +
  labs(x='Geographic distance (ln(km))', y='Community similarity (ln(BC))') +
  stat_smooth(method = "lm", size = .8,level = .95) +
  theme_bw()
