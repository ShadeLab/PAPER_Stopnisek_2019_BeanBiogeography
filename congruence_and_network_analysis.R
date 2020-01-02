#******************************************************
# For running the code below 16S_US_analysis.R and 
# ITS_analysis.R need to be executed!
#******************************************************

#analysis of congruence between ITS and 16S datsets

otu_its.noS1 <- otu_its[,-25]
otu_its.noS1.rare <- t(rrarefy(t(otu_its.noS1), min(colSums(otu_its.noS1))))
otu_its.noS1.rare <- otu_its.noS1.rare[rowSums(otu_its.noS1.rare)>0,]
bean_its_noS1.dist <- vegdist(t(otu_its.noS1.rare), method='bray')

proITSvs16S <- protest(bean_its_noS1.dist, bean.dist)
plot(proITSvs16S, kind=2)
# 
procruTEMP <- procrustes(bean_its_noS1.dist, bean.dist)
plot(procruTEMP, scores='lc', symmetric=T, choices=c(1:3))

#######################################################################################
#MENA network analysis 
MENA_ZP <- read.csv('Data/ZiPi_bean_0.88_new.csv', header = T)
core_US_all <- append(taxa_US, core_ITS)
sum(!MENA_ZP$ID %in% core_US_all)

zp_plot <-MENA_ZP %>%
  mutate(otu=ID) %>%
  mutate(core=if_else(otu %in% core_US_all, 'US core (n=75 OTUs)', 'no (n=497 OTUs)')) %>%
  left_join(tax_v1[,c('Phylum', "Genus", 'otu')], by='otu') %>%
  ggplot(aes(y=Zi, x=Pi, color=core)) +
  scale_color_manual(values = c('black', '#009900')) +
  geom_hline(yintercept=2.5, linetype="dashed")+
  geom_vline(xintercept=c(.62), linetype="dashed")+
  geom_point() +
  labs(y='Zi (within-module connectivity)', x='Pi (among-module connectivity)', color='Membership') +
  theme_classic()+
  theme(legend.position = c(.4,.1),
        legend.background=element_blank())

########
# Fig 5B
ggMarginal(zp_plot, groupColour = TRUE, groupFill = TRUE)


modules_No <- MENA_ZP
names(modules_No)[2] <- 'otu'

########
# Fig 5C

data.frame(otu=as.factor(rownames(otu.rhizo.rel.abun)),otu.rhizo.rel.abun) %>%
    gather(sample_ID, abun, -otu) %>%
    left_join(modules_No, by='otu') %>%
    filter(!(is.na(No..module))) %>%
    left_join(map_combined[,c('sample_ID','pH', 'state', 'bean', 'plot', 'site')], by='sample_ID') %>%
    left_join(tax_v1[,c('Kingdom','Phylum','Family' ,"Genus", 'otu', 'final_names')], by='otu') %>%
    filter(No..module %in% c(0,1,2,3))%>%
    group_by(site,sample_ID, No..module) %>%
    summarise(sum_abun=sum(abun),
              n_samples=length(sample_ID),
              rel_abun=sum_abun/n_samples) %>%
    ggplot(aes(x=site, y=rel_abun, color=site)) +
    geom_boxplot(alpha=.4, outlier.shape = NA) +
    geom_jitter() +
    scale_color_manual(values=c("#4daf4a","#e41a1c", "#984ea3","darkorange","#377eb8"))+
    labs(y='Relative abundance', x='Sampling location') +
    facet_grid(~No..module, scales = "free_y") +
    theme_classic()+
    theme(legend.position = 'none',
          strip.text.x = element_text(size=10),
          strip.background = element_rect(colour="transparent", fill="transparent")))


bact_fung_inter<- modules_No %>%
  left_join(tax_v1, by='otu') %>%
  filter(otu %in% c(unique(its_core_taxa$otu_id), core_US_otus)) %>%
  mutate(Kingdom=if_else(is.na(Kingdom), 'Fungi', Kingdom)) %>%
  group_by(No..module, Kingdom) %>%
  summarise(n=length(otu))

mena_int<- read.csv('Data/interaction_mena.csv', sep=',', header=F)
archaeaOTU <- tax_v1$otu[tax_v1$Kingdom=='Archaea']

coreOTUs <- c(unique(its_core_taxa$otu_id), core_US_otus)
interTAXA <- mena_int %>% separate(V1, c('start','connect', 'end','value'), sep='([\\(\\)\\=])', remove=T) %>%
  mutate(start=str_trim(start, side = "both"),
         end=str_trim(end, side = "both"),
         start.memb=if_else(start %in% coreOTUs, 1,0),
         end.memb=if_else(end %in% coreOTUs, 1, 0),
         core.int=start.memb+end.memb,
         domain.start= if_else(start %in% tax_v1$otu, '16S', 'ITS'),
         domain.end= if_else(end %in% tax_v1$otu, '16S', 'ITS'),
         fungus=if_else(start %in% tax_its$otu_id, 1,0),
         fungus=if_else(end %in% tax_its$otu_id,1,fungus),
         archaea=if_else(start %in% archaeaOTU, 1,0),
         archaea=if_else(end %in% archaeaOTU, 1, archaea),
         bacteria= if_else(start %in% tax_v1$otu, 1, 0),
         bacteria= if_else(end %in% tax_v1$otu, 1, bacteria),
         bact_fung=if_else(bacteria+fungus==2, 1, 0),
         bact_arch=if_else(bacteria+archaea==2, 1, 0),
         arch_fung=if_else(archaea+fungus==2, 1, 0),
         core_core=if_else(core.int>1,1,0))

tax_16S <- tax_v1[,-c(1,2)]
names(tax_16S)[10] <- 'otu_id'
tax_combined <- rbind(tax_16S, tax_its)
alluv.df <- interTAXA %>%
  filter(start %in% coreOTUs) %>%
  left_join(tax_combined[,c('Kingdom','Phylum','Order','Family','Genus' ,'otu_id')], by=c('start'='otu_id')) %>%
  left_join(tax_combined[,c('Kingdom','Phylum','Order','Genus', 'otu_id')], by=c('end'='otu_id')) %>%
  group_by(Kingdom.x,Phylum.x,Order.x,Family,Genus.x, Kingdom.y, Phylum.y, Genus.y) %>%
  summarise(n_connect=length(start)) %>%
  group_by(Genus.x)%>%
  filter(n_connect>1) %>%
  spread(Genus.y, n_connect) %>%
  replace(., is.na(.), 0)
