#****************************************************************************
#ZOTUs
######
zotu_us <- read.table('Data/ZOTU_filtered_table.txt', sep='\t', row.names = 1, header=T)
map_combined <- read.table('Data/map.txt', row.names = 1, sep='\t', header=T)

#remove lines where taxonomy includes mitochondria
zotu_us <- zotu_us[!grepl("Mitochondria", zotu_us$taxonomy),]
zotu_us <- zotu_us[!grepl("Chloroplast", zotu_us$taxonomy),]

tax_zotu <- zotu_us['taxonomy']
tax_zotu <- tax_zotu %>%
  separate(taxonomy, into=c("Kingdom", "Phylum", "Class", 
                            "Order", "Family", "Genus", "Species"), sep="; ", remove=F)
tax_zotu[2:8] <- lapply(tax_zotu[2:8], function(x) gsub(".*__", "", x))

zotu_us['taxonomy'] <- NULL
zotu_us['SVERCc1'] <- NULL

# Order the samples
zotu_us <- zotu_us[,order(colnames(zotu_us))]
# Check to make sure they all match with each other
map_combined <- map_combined[order(rownames(map_combined)),]
map_combined <- map_combined[-25,]
rownames(map_combined)==colnames(zotu_us)

set.seed(4)
zotu_us_rare <- t(rrarefy(t(zotu_us), min(colSums(zotu_us)))) 
zotu_us_rare <- zotu_us_rare[rowSums(zotu_us_rare)>0,]

#PCoA
bean.dist.zotu <- vegdist(t(zotu_us_rare), method="bray")
bean.zotu.pcoa <- cmdscale(bean.dist.zotu, eig=TRUE)

map_short <- map_combined[,-c(1:23,26)]
bean_envfit <- envfit(bean.pcoa, map_short)

ax1.bean <- bean.pcoa$eig[1]/sum(bean.pcoa$eig)
ax2.bean <- bean.pcoa$eig[2]/sum(bean.pcoa$eig)

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
legend(-.03, -.08, legend = c('CELRK', 'Eclipse', 'Bulk soil'), col='black',pch=c(23,21, 22), cex=.8, box.col = 0, title='Genotype')

#Occ_abun
zotu_us_rare_rhizo <- zotu_us_rare[,map_combined$soil == 'rhizosphere']
zotu_rare <-  zotu_us_rare_rhizo[rowSums(zotu_us_rare_rhizo)>0,]
zotu_rare_PA <- 1*((zotu_rare>0)==1)
zotu_rare_PA <- zotu_rare_PA[rowSums(zotu_rare_PA)>0,]
Occ_zotu <- rowSums(zotu_rare_PA)/ncol(zotu_rare_PA)

zotu_rare_rel <- decostand(zotu_rare, method = 'total', MARGIN = 2)
com_abund_zotu <- rowSums(zotu_rare_rel)/ncol(zotu_rare_rel)
zotu_df_occ <- data.frame(otu=names(Occ_zotu), occ=Occ_zotu) 
zotu_df_abun <- data.frame(otu=names(com_abund_zotu), abun=log10(com_abund_zotu))
occ_abun_zotu <- left_join(zotu_df_occ, zotu_df_abun)

length(Occ_zotu[Occ_zotu==1])

ggplot(data=occ_abun_zotu, aes(x=abun, y=occ)) +
  geom_point() +
  theme_bw()+
  #scale_fill_manual(values=c('black', 'white')) +
  labs(x=paste("log10(mean abundace)\n (n=", nrow(occ_abun_zotu)," OTUs)", sep=''), y= paste("Occupancy (n=", length(colnames(zotu_rare_PA))," samples)",sep=''), 
       bg= NULL) +
  theme(plot.title = element_text(hjust = 0.5, size=12),
        legend.position = 'top',
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Mantel ZOTU vs OTU US study
mantel(bean.dist.zotu, bean.dist)
protest(bean.pcoa, bean.zotu.pcoa)


# Investigating ecotypes in the core 
#**********************************
# To do in command line
# # Using seqtk tool (https://github.com/lh3/seqtk) for extract fastq reads that are found in the global core OTUs
# module load seqtk
# seqtk subseq merged_fastq_us/merged.fq fastqCOREids.lst > seqtk_COREfastq.fq
# Using FASTX tool for finding unique sequences 
# ml GCC/5.4.0-2.26  OpenMPI/1.10.3
# ml FASTX-Toolkit/0.0.14
# fastx_collapser -v -i seqtk_COREfastq.fq -o FASTX_output.fq
# Input: 421369 sequences (representing 421369 reads)
# Output: 113074 sequences (representing 421369 reads)
# Search for matches using blastn
# makeblastdb -in FASTX_output.fn -dbtype nucl -out ref
# blastn -db ref -query seqtk_COREfastq.fn -out blastCORE.txt -outfmt 6 -perc_identity 100 -max_target_seqs 1
#**********************************
#Extract compressed file OTU_fastq_map.zip in Data folder!

US_fastq <- read.table("Data/OTU_fastq_map.txt", sep='\t')
BLAST.coreReads <- read.table('Data/blastCORE.txt')
BLAST.coreReads <- BLAST.coreReads[,c(1:4)]
names(BLAST.coreReads) <- c('readID', 'uniqueRead', 'matchPercent', 'length')
US_fastq <- US_fastq[,c(4,9,10)]
names(US_fastq) <- c('percentMatch','readID','refID')

map_v2=map_combined
map_v2$sample_ID <- rownames(map_v2)

subOTUsiteSpecific <- left_join(BLAST.coreReads, US_fastq) %>%
  tidyr::separate(readID, c('sample_ID', 'read.No'), '[.]', remove = F) %>%
  left_join(map_v2[,c('state','soil', 'sample_ID','site')], by='sample_ID') %>%
  filter(soil=='rhizosphere') %>%
  group_by(refID, uniqueRead) %>%
  summarise(presenceSite=length(unique(site)))

subOTUoccupancy<- subOTUsiteSpecific %>%
  group_by(refID) %>%
  summarise(nUniques=length(unique(uniqueRead)),
            occ5=sum(presenceSite==5),
            occ_less_5=sum(presenceSite<5)) %>%
  mutate(name=paste(refID,' (',nUniques,',',occ5,')', sep='')) %>%
  melt(id.vars=c('refID','name','nUniques'), measure.vars =c('occ5', 'occ_less_5')) %>%
  mutate(relCount=value/nUniques) 

COREsubOTUtable <- read.table("Data/ZOTU_subOTU_UNOISE_table.txt", header = T, row.names = 1)
COREsubOTUmap <- read.table("Data/ZOTUS_subOTU_UNOISE_map.txt", header=F)
names(COREsubOTUmap) <- c('readID', 'ZOTUid')

subOTU_match <- left_join(BLAST.coreReads, US_fastq) %>%
  tidyr::separate(readID, c('sample_ID', 'read.No'), '[.]', remove = F) %>%
  left_join(map_v2[,c('state','soil', 'sample_ID', 'site')], by='sample_ID') %>%
  left_join(COREsubOTUmap) %>%
  na.omit()

temp <- subOTU_match %>%
  group_by(refID, ZOTUid) %>%
  summarise(core=if_else(length(unique(site))>4, 1,0),
            countZOTU=length(readID))%>% 
  group_by(refID) %>%
  mutate(allReads=sum(countZOTU)) %>%
  mutate(relabun=countZOTU/allReads) %>%
  filter(core>0)%>%
  group_by(refID) %>%
  summarise(maxRel=max(relabun))

temp_v2 <- subOTU_match %>%
  group_by(refID, ZOTUid) %>%
  summarise(core=if_else(length(unique(site))>4, 1,0),
            countZOTU=length(readID))%>% 
  group_by(refID) %>%
  mutate(allReads=sum(countZOTU)) %>%
  mutate(relabun=countZOTU/allReads,
         maxRel=relabun)

NoCOREzotus<- temp_v2%>%
  group_by(refID)%>%
  summarise(NoCoreZOTU=sum(core))

########
# Fig S4

temp_v2 %>%
  left_join(plotData[,c('refID', 'name')]) %>%
  filter(!is.na(name)) %>%
  ggplot(aes(name, maxRel, col=as.factor(core)))+
  geom_point() +
  labs(x=NULL, y='Relative abundance', col='Occupancy') +
  theme_bw()+
  theme(axis.text = element_text(color='grey40'),
        axis.title = element_text(color='grey40'),
        axis.text.x = element_text(angle=90, size=10, hjust=1, vjust = 0.5, color='grey40'))

