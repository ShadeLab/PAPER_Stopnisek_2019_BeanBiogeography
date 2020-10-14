#************************
###   ITS analysis
#########################

otu_its <- read.table('Data/OTU_table_ITS.txt', sep='\t', row.names = 1, header=T)
map_its <- read.table('Data/map_ITS.txt', row.names = 1, sep='\t', header=T)
tax_its <- read.csv('Data/consensus_taxonomy_latest.csv', header=T, na.strings=c("","NA"), row.names=1)

# Order the samples
otu_its <- otu_its[,order(colnames(otu_its))]
# Order the samples of the map the same way
map_its=map_its[order(rownames(map_its)),]
# Check to make sure they all match with each other
rownames(map_its)==colnames(otu_its)

#Rarefying data to the sample with lowest reads (22716 )
#rarecurve(t(otu_its), step=100, label = FALSE, xlim=c(0,22716))

#No of reads per sample
readNo_ITS <- data.frame(otu=as.factor(rownames(otu_its)),otu_its) %>%
  gather(sample_name, abun, -otu) %>%
  left_join(map_its, by='sample_name') %>%
  group_by(sample_name) %>%
  summarise(readCount=sum(abun)) %>%
  arrange(desc(readCount)) %>%
  ggplot(aes(x=sample_name,y=readCount)) +
  geom_bar(stat='identity')+
  geom_hline(yintercept=min(colSums(otu_its)), linetype="dashed", color = "red") +
  theme_bw()+
  xlab('Sample name') +
  ylab('read number') +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))

#Rarefying the data
set.seed(107)
otu_its_rare <- t(rrarefy(t(otu_its), min(colSums(otu_its))))
otu_its_rare <- otu_its_rare[rowSums(otu_its_rare)>0,]

# ######################################################
# # Alpha Diversity
s_its <- specnumber(otu_its_rare,MARGIN=2)
h_its <- vegan::diversity(t(otu_its_rare), "shannon")
pielou_its=h_its/log(s_its)
# 
map.div.its <- map_its
map.div.its$Richness <- s_its
map.div.its$Shannon <- h_its
map.div.its$Pielou <- pielou_its

map.alpha.its <- melt(map.div.its, id.vars=c("sample_name","bean", "plot", "state", 'soil', 'site', 'fertilization', 'irrigation', 'pH'),
                  measure.vars=c("Richness", "Shannon", "Pielou"))

map.alpha.its %>%
  filter(variable=='Richness') %>%
  group_by(variable) %>%
  summarise(mean_var=mean(value))

map.alpha.its %>%
  filter(soil=='rhizosphere') %>%
  group_by(variable) %>%
  summarise(mean_var=mean(value))

bartlett.test(value~site, data = map.alpha.its)
bartlett.test(value~bean, data = map.alpha.its)
bartlett.test(value~pH, data = map.alpha.its)
bartlett.test(value~fertilization, data = map.alpha.its)
bartlett.test(value~bean, data = map.alpha.its)

t.test(value~bean, data=map.alpha.its[map.alpha.its$variable=="Shannon" & map.alpha.its$soil=='rhizosphere',])
tidy(aov(value~site, data=map.alpha.its[map.alpha.its$variable=="Richness",]))
summary(aov(value~site, data=map.alpha.its[map.alpha.its$variable=="Richness",]))

ggplot(map.alpha.its[map.alpha.its$soil=='rhizosphere',], aes(y=value,x=bean, color=bean))+
  scale_color_manual(values = c('darkorange', 'black')) + 
  facet_wrap( ~ variable, scales = "free_y") + 
  theme_bw()+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))+
  geom_boxplot()

ggplot(map.alpha.its[map.alpha.its$soil=='rhizosphere' & map.alpha.its$variable=="Richness",], aes(y=value, x=bean, color=bean))+
  scale_color_manual(values = c('darkorange', 'black')) + 
  theme_bw()+
  geom_boxplot() +
  xlab('Genotype')+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = 'none')

# ######################################################
#Beta diversity
#PCoA
bean_its.dist <- vegdist(t(otu_its_rare), method="bray")
bean_its.pcoa <- cmdscale(bean_its.dist, eig=TRUE)

map_short_its <- map_its[,-c(1:23,26)]
bean_its_envfit <- envfit(bean_its.pcoa, map_short_its)

ax1.bean_its <- bean_its.pcoa$eig[1]/sum(bean_its.pcoa$eig)
ax2.bean_its <- bean_its.pcoa$eig[2]/sum(bean_its.pcoa$eig)
ax3.bean_its <- bean_its.pcoa$eig[3]/sum(bean_its.pcoa$eig)

bean_plot_colors <- rep("#e41a1c", nrow(map_its))
bean_plot_colors[map_its$site=="WA"] <- "#377eb8"
bean_plot_colors[map_its$site=="CO"] <- "#4daf4a"
bean_plot_colors[map_its$site=="NE"] <- "#984ea3"
bean_plot_colors[map_its$site=="SVREC"] <- "darkorange"

bean_plot_shapes <- rep(21, nrow(map_its))
bean_plot_shapes[map_its$bean=="CELRK"] <- 23
bean_plot_shapes[map_its$soil=="bulk"] <- 22

plot(bean_its.pcoa$points[,1], bean_its.pcoa$points[,2], cex=1.5, bg=bean_plot_colors, pch=bean_plot_shapes,
     xlab=paste("PCoA1: ",100*round(ax1.bean_its,3),"% var. explained",sep=""), 
     ylab=paste("PCoA2: ",100*round(ax2.bean_its,3),"% var. explained",sep=""))
plot(bean_its_envfit, p=0.001, col='black')
legend(-.15, -.1, legend=c('CO','MRF','SVREC','NE', 'WA'), col=c("#4daf4a","#e41a1c", "darkorange","#984ea3","#377eb8"), 
       box.lty=0, cex= .8, pch = 20, title='Location')
legend(-.4, -.2, legend = c('CELRK', 'Eclipse'), col='black',pch=c(23,21), cex=.8, box.col = 0, title='Soil')

adonis(bean_its.dist~map_its$bean) #no diff between the genotypes
adonis(bean_its.dist~map_its$soil) #no differences between bulk and rhizo
adonis(bean_its.dist~map_its$site) #p=0.001
adonis(bean_its.dist~map_its$pH) #p=0.001
adonis(bean_its.dist~map_its$irrigation) 
adonis(bean_its.dist~map_its$fertilization) 
adonis(bean_its.dist~map_its$P)
adonis(bean_its.dist~map_its$Nitrogen)
adonis(bean_its.dist~map_its$NH4)
adonis(bean_its.dist~map_its$NO3)
adonis(bean_its.dist~map_its$OM)
adonis(bean_its.dist~map_its$site/map_its$bean) #no diff between the genotypes
adonis(bean_its.dist~map_its$bean*map_its$site) #no diff between the genotypes

#are co-occurring species different when you only compare species at the same sites?
adonis(bean_its.dist~map_its$bean, strata=map_its$site)
#Are sites different when you only comparing the same species?
adonis(bean_its.dist~map_its$site, strata=map_its$bean)


# ######################################################
# #Taxonomic distribution 
map_its$sample_ID <- rownames(map_its)
tax_its$otu_id <- rownames(tax_its)
tax_its <- tax_its %>% select(otu_id, everything())

lastValue <- function(x) tail(x[!is.na(x)], 1)
last_taxons<- apply(tax_its, 1, lastValue)
tax_its$last_taxon <- last_taxons
tax_its$final_names <- paste(tax_its$otu_id, tax_its$last_taxon, sep='-')

otu_its_rare_rhizo <- otu_its_rare[,map_its$soil=='rhizosphere']
otu.rel.abun.ITS <- decostand(otu_its_rare_rhizo, method="total", MARGIN=2) #calculating relative abundance

tmp_its_v3 <- data.frame(otu_id=as.factor(rownames(otu.rel.abun.ITS)),otu.rel.abun.ITS) %>%
  gather(sample_ID, abun, -otu_id) %>%
  left_join(map_its[,c('sample_ID','pH', 'state', 'bean', 'plot', 'site', 'fertilization', 'irrigation')], by='sample_ID') %>%
  left_join(tax_its, by='otu_id') %>%
  group_by(Phylum) %>%
  summarise(
    n=sum(abun)/length(unique(sample_ID))) %>%
  arrange(desc(n)) %>%
  group_by(Phylum)%>%
  summarise(min=min(n),
            max=max(n)) %>%
  arrange(desc(max))
  
tax_its$Phylum <- as.character(tax_its$Phylum)

site_ITS <- data.frame(otu_id=as.factor(rownames(otu.rel.abun.ITS)),otu.rel.abun.ITS) %>%
  gather(sample_ID, abun, -otu_id) %>%
  left_join(map_its[,c('sample_ID','pH', 'state', 'bean', 'plot', 'site', 'fertilization', 'irrigation')], by='sample_ID') %>%
  left_join(tax_its, by='otu_id') %>%
  group_by(Phylum, sample_ID) %>%
  summarise(
    n=sum(abun)/length(unique(sample_ID))) %>%
  mutate(tax_names=ifelse(n>=.01, Phylum, 'other'),
         tax_names=ifelse(is.na(tax_names), 'other',tax_names )) #%>%
  # ggplot(aes(y=n, x=sample_ID, fill= tax_names)) +
  # geom_bar(stat='identity') +
  # theme_bw()+
  # labs(y="Relative abundace", x= NULL, fill='Phylum') +
  # theme(plot.title = element_text(hjust = 0.5),
  #       legend.text=element_text(size=6),
  #       legend.position = 'bottom',
  #       axis.text.x=element_text(angle = 45, hjust = 1)) +
  # guides(fill=guide_legend(ncol=3)) 

All_ITS <- data.frame(otu_id=as.factor(rownames(otu.rel.abun.ITS)),otu.rel.abun.ITS) %>%
  gather(sample_ID, abun, -otu_id) %>%
  left_join(map_its[,c('sample_ID','pH', 'state', 'bean', 'plot', 'site', 'fertilization', 'irrigation')], by='sample_ID') %>%
  left_join(tax_its, by='otu_id') %>%
  group_by(Phylum) %>%
  summarise(
    n=sum(abun)/length(unique(sample_ID))) %>%
  mutate(tax_names=ifelse(n>=.01, Phylum, 'other'),
         tax_names=ifelse(is.na(tax_names), 'other',tax_names ),
         sample_ID='all') #%>%
  # ggplot(aes(y=n, x=as.character(sample_ID), fill= tax_names)) +
  # geom_bar(stat='identity') +
  # theme_bw()+
  # labs(y="Relative abundace", x= NULL, fill='Phylum') +
  # theme(plot.title = element_text(hjust = 0.5),
  #       legend.text=element_text(size=6),
  #       legend.position = 'bottom') +
  # guides(fill=guide_legend(ncol=3)) 

#########
# FigS 3C

ITS_dt <- rbind(as.data.frame(site_ITS), as.data.frame(All_ITS)) %>%
  mutate(sep=if_else(sample_ID == 'all', 'all', 'site')) %>%
  ggplot(aes(y=n, x=as.character(sample_ID), fill=tax_names)) +
  geom_bar(stat='identity') +
  facet_grid(~sep, scales='free_x', space='free_x')+
  theme_bw()+
  labs(y="Relative abundance", x= NULL, fill='Phylum') +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=6),
        legend.position = 'bottom',
        axis.text.x=element_text(angle = 45, hjust = 1),
        strip.text.x = element_blank()) +
  guides(fill=guide_legend(ncol=3)) 

######################################################
#Venn diagram
#make subset OTU tables
map_its$state <- as.character(map_its$state)
CO_otu <- otu_its_rare[,map_its$site=="CO"]
MRF_otu <- otu_its_rare[,map_its$site=="MRF"]
NE_otu <- otu_its_rare[,map_its$site=="NE"]
WA_otu <- otu_its_rare[,map_its$site=="WA"]
SVREC_otu <- otu_its_rare[,map_its$site=="SVREC"]

CELRK_otu <- otu_its_rare[,map_its$bean=="CELRK"]
Eclipse_otu <- otu_its_rare[,map_its$bean=="Eclipse"]

# make presence absence list from soil and plant into 1 & 0
CO_venn <- 1*(rowSums(CO_otu)>0)
MRF_venn <- 1*(rowSums(MRF_otu)>0)
NE_venn <- 1*(rowSums(NE_otu)>0)
WA_venn <- 1*(rowSums(WA_otu)>0)
SVREC_venn <- 1*(rowSums(SVREC_otu)>0)

CELRK_venn <- 1*(rowSums(CELRK_otu)>0)
Eclipse_venn <- 1*(rowSums(Eclipse_otu)>0)

# combine vectors into a matrix cbind = column bind (r bind = row bind)
venndata <- cbind(CO_venn,MRF_venn)
venndata <- cbind(venndata,NE_venn)
venndata <- cbind(venndata,WA_venn)
venndata <- cbind(venndata,SVREC_venn)

# venndata_bean <- cbind(CELRK_venn,Eclipse_venn)
colnames(venndata) <- c("CO", "MRF", "NE", "CO", 'SVREC')
venndata <- venndata[rowSums(venndata)>0,]
v=vennCounts(venndata)
v2=round(v[,"Counts"]/sum(v[,"Counts"]),2)

#########
# FigS 3D

vennDiagram(v) 

#Shared by site:
CO_present <- CO_otu[rowSums(CO_otu)>1,]
MI_present <- MI_otu[rowSums(MI_otu)>1,]
NE_present <- MI_otu[rowSums(NE_otu)>1,]
WA_present <- MI_otu[rowSums(WA_otu)>1,]
SA_present <- MI_otu[rowSums(SA_otu)>1,]

CO_share <- 70/length(rownames(CO_present))
MI_share <- 70/length(rownames(MI_present))
NE_share <- 70/length(rownames(NE_present))
WA_share <- 70/length(rownames(WA_present))
SA_share <- 70/length(rownames(SA_present))

mean_shared=mean(c(CO_share,MI_share,NE_share, WA_share, SA_share))

######################################################
#Occupancy abundance figure
#selecting variety unique OTUs

CELRK_otu <- CELRK_otu[rowSums(CELRK_otu)>0,]
Eclipse_otu <- Eclipse_otu[rowSums(Eclipse_otu)>0,]

CELRK_uniq <- CELRK_otu[!(rownames(CELRK_otu) %in% rownames(Eclipse_otu)),]
Eclipse_uniq <- Eclipse_otu[!(rownames(Eclipse_otu) %in% rownames(CELRK_otu)),]

#Occupancy abundance figure - combined data
otu_ITS_rare_bean <- otu_its_rare[rowSums(otu_its_rare)>0,]
bean_otu_PA <- 1*((otu_ITS_rare_bean>0)==1)
bean_otu_PA <- bean_otu_PA[rowSums(bean_otu_PA)>0,]
Occ <- rowSums(bean_otu_PA)/ncol(bean_otu_PA)
com_abund_bean <- rowSums(otu_ITS_rare_bean)/ncol(otu_ITS_rare_bean)

bean_df_occ <- data.frame(otu=names(Occ), occ=Occ)
bean_df_abun <- data.frame(otu=names(com_abund_bean), abun=log10(com_abund_bean))
bean_df_abun$found <- 'shared'
bean_df_abun$found[bean_df_abun$otu %in% rownames(CELRK_uniq)] <- 'CELRK'
bean_df_abun$found[bean_df_abun$otu %in% rownames(Eclipse_uniq)] <- 'Eclipse'

bean_df_abun$col[bean_df_abun$found == 'CELRK'] <- 'darkorange'
bean_df_abun$col[bean_df_abun$found == 'Eclipse'] <- 'black'
bean_df_abun$col[bean_df_abun$found == 'shared'] <- 'lightgray'

occ_abun_bean <- left_join(bean_df_abun, bean_df_occ)

size_occ_its<- data.frame(otu_id=as.factor(rownames(bean_otu_PA)),bean_otu_PA) %>%
  gather(sample_ID, abun_rhizo, -otu_id) %>%
  left_join(map_its[,c('sample_ID', 'site', 'bean', 'soil')], by='sample_ID') %>%
  left_join(tax_its, by='otu_id') %>%
  #filter(soil!='bulk') %>%
  group_by(otu_id, site, Family, Genus, final_names) %>%
  mutate(sum_abun=sum(abun_rhizo),
         n_rep=length(abun_rhizo),
         rel_occ=sum_abun/n_rep,
         is_present=1*((sum_abun>0)==1))

#calculating presence across sites
size_occ_its %>%
  group_by(otu_id, site) %>%
  dplyr::summarize(n=sum(is_present),
                   presence= 1*((n>0)==1)) %>%
  group_by(otu_id) %>%
  dplyr::summarize(
    total_presence=sum(presence)
  ) -> tmp_occ_its

names(occ_abun_bean)[1] <- 'otu_id'
combined_occ_data_its <- left_join(tmp_occ_its, occ_abun_bean)

#Sloan neutral model 
spp=t(otu_its_rare_rhizo)
taxon=as.vector(rownames(otu_its_rare_rhizo))

#Models for the whole community
obs.np=sncm.fit(spp, taxon, stats=FALSE, pool=NULL)
sta.np=sncm.fit(spp, taxon, stats=TRUE, pool=NULL)
sta.np.ITS <- sta.np

above.pred=sum(obs.np$freq > (obs.np$pred.upr), na.rm=TRUE)/sta.np$Richness
below.pred=sum(obs.np$freq < (obs.np$pred.lwr), na.rm=TRUE)/sta.np$Richness

ap = obs.np$freq > (obs.np$pred.upr)
bp = obs.np$freq < (obs.np$pred.lwr)

########
# Fig 2B

ggplot(data=occ_abun_bean, aes(x=abun, y=occ)) +
    theme_bw()+
    geom_point(pch=21,  size=3, aes(fill=found)) +
    geom_line(color='black', data=obs.np, size=2, aes(y=obs.np$freq.pred, x=log10(obs.np$p)+4.4)) +
    geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.upr, x=log10(obs.np$p)+4.4))+
    geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.lwr, x=log10(obs.np$p)+4.4))+
    scale_fill_manual(aes(breaks = found), values=c('darkorange','black', 'white')) +
    labs(x=paste("log10(mean abundace)\n (n=", sta.np.ITS$Richness," OTUs)", sep=''), y=paste("Occupancy (n=",sta.np.ITS$Samples," samples)", sep='')) +
    theme(legend.position="top",
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    guides(fill = guide_legend(override.aes = list(alpha=1)),
           fill = guide_legend(title=NULL))

######################################################
#Core community
#taking into consideration the 1505 otu shared accross the locations
length(occ_abun_bean$otu[occ_abun_bean$occ==1]) #13
core_ITS <- as.character(occ_abun_bean$otu[occ_abun_bean$occ==1])

together_occ1 <- otu.rel.abun.ITS[rownames(otu.rel.abun.ITS) %in% occ_abun_bean$otu[occ_abun_bean$occ==1],]

occ_ITS <- data.frame(otu_id=as.factor(rownames(together_occ1)),together_occ1) %>%
  gather(sample_ID, abun, -otu_id) %>%
  left_join(map_its[,c('sample_name','sample_ID','site','bean')], by='sample_ID') %>%
  left_join(tax_its, by='otu_id') %>%
  group_by(site, Phylum, Class, otu_id, final_names) %>%
  summarise(n=sum(abun)/length(site))

plot_core_ITS <- data.frame(data.frame(otu_id=as.factor(rownames(together_occ1)),together_occ1)) %>%
  gather(sample_ID, abun, -otu_id) %>%
  left_join(map_its[,c('sample_ID','pH', 'state', 'bean', 'plot', 'sample_name', 'soil')], by='sample_ID') %>%
  left_join(tax_its, by='otu_id') %>%
  filter(soil != 'bulk') %>%
  group_by(Order) %>%
  dplyr::summarise(
    n=sum(abun)/length(unique(sample_name))) %>%
  arrange(desc(n))

plot_core_ITS$Order2 <- factor(plot_core_ITS$Order, levels = plot_core_ITS$Order[order(plot_core_ITS$n)])

core_all_v1 <- core_all_v1 %>%
  mutate(new_names=as.factor(new_names),
         new_names=fct_reorder(new_names, abun))

its_core_taxa <- data.frame(data.frame(otu_id=as.factor(rownames(together_occ1)),together_occ1)) %>%
  gather(sample_ID, abun, -otu_id) %>%
  left_join(map_its[,c('sample_ID','pH', 'state', 'bean', 'plot', 'sample_name', 'soil')], by='sample_ID') %>%
  left_join(tax_its, by='otu_id') %>%
  filter(soil != 'bulk') %>%
  mutate(Order=as.character(Order),
         new_names = if_else(is.na(Order), 'other', Order),
         Phylum=if_else(is.na(Phylum), 'Fungi', Phylum)) %>%
  arrange(desc(abun))

x = tapply(its_core_taxa$abun, its_core_taxa$new_names, function(x) mean(x))
x = sort(x, TRUE)
its_core_taxa$new_names = factor(as.character(its_core_taxa$new_names), levels=names(x))

core_ITS_US <- ggplot(its_core_taxa,aes(x=new_names, y=abun, fill=Phylum)) +
  theme_bw()+
  labs(y='Relative abundance', x='Order', fill='Phylum') +
  geom_boxplot() +
  scale_fill_manual(values=c('#e6194B', '#3cb44b', '#a9a9a9'))+
  theme(legend.position = c(.7,.8),
        legend.text = element_text(size=8),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()

coreTaxaNo_ITS <-its_core_taxa %>%
  group_by(Phylum,new_names) %>%
  summarise(n_taxa=length(unique(otu_id))) %>%
  ggplot(aes(x=new_names, n_taxa, fill=Phylum)) +
  scale_fill_manual(values=c('#e6194B', '#3cb44b', '#a9a9a9'))+
  geom_bar(stat='identity') +
  ylab('Number of taxa') +
  xlab('') +
  labs(fill='Phylum')+
  theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = 'none',
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  coord_flip()

########
# Fig 2D

grid.arrange(core_ITS_US, coreTaxaNo_ITS,  widths = c(2.5,.7))



its_core_taxa %>%
  group_by(Phylum,new_names,Family,Genus,otu_id, Species) %>%
  summarise(meanTax=mean(abun)) %>%
  arrange(desc(meanTax))

 
