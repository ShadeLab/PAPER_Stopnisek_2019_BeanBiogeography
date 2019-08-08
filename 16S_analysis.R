#Using Peter Rosenmai's code for calculating the Distance matrix for geographic points (https://eurekastatistics.com/calculating-a-distance-matrix-for-geographic-points-using-r/)
ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.
  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
  
  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }
  
  n.geopoints <- nrow(df.geopoints)
  
  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints
  
  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})
  
  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  
  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  
  return(mat.distances)
}


tax_colors <-  c('Acidobacteria'='#3cb44b','Actinobacteria'='#a9a9a9' ,'Bacteroidetes'='#aaffc3',
                 'BRC1'='#4363d8','Chloroflexi'='#f58231','Euryarchaeota'='#911eb4','FBP'='#42d4f4',
                 'Firmicutes'='#f032e6','Gemmatimonadetes'='#bfef45','Nitrospirae'='#469990', 
                 'Planctomycetes'='#e6beff','Proteobacteria'= '#e6194B','Thaumarchaeota'='#808000','Verrucomicrobia'='#000075')

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
######################################################
# Dataset preparation

otu_us <- read.table('Data/OTU_table_US.txt', sep='\t', row.names = 1, header=T)
map_combined <- read.table('Data/map.txt', row.names = 1, sep='\t', header=T)

#remove lines where taxonomy includes mitochondria
otu_us <- otu_us[!grepl("Mitochondria", otu_us$taxonomy),]
otu_us <- otu_us[!grepl("Chloroplast", otu_us$taxonomy),]

tax_us <- otu_us['taxonomy']
tax_us <- tax_us %>%
  separate(taxonomy, into=c("Kingdom", "Phylum", "Class", 
                            "Order", "Family", "Genus", "Species"), sep="; ", remove=F)
tax_us[2:8] <- lapply(tax_us[2:8], function(x) gsub(".*__", "", x))

otu_us['taxonomy'] <- NULL         #removing taxonomy column from the OTU table
map_combined <- map_combined[-1,]  #removing sample SVERCc1 because of low read coverage 
otu_us['SVERCc1'] <- NULL          #removing sample SVERCc1 because of low read coverage

# Order the samples
otu_us <- otu_us[,order(colnames(otu_us))]
# Order the samples of the map the same way
map_combined=map_combined[order(rownames(map_combined)),]
# Check to make sure they all match with each other
rownames(map_combined)==colnames(otu_us)

map_combined$sample_ID <- rownames(map_combined)

#Rarefying data to the sample with lowest reads (31255)
# setEPS()
# postscript('rarefaction_curves_16S.eps', width=4, height = 3.8)
# rarecurve(t(otu_us), step=1000, label = FALSE, xlim=c(0,31255))
# dev.off()

#How many samples per fertilization treatment
ggplot(map_combined, aes(x=fertilization, fill=site)) + 
  theme_bw()+
  geom_histogram(stat='count', position = 'dodge')

# How many reads per sample?
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

# How many reads per growing location?
data.frame(otu=as.factor(rownames(otu_us)),otu_us) %>%
  gather(sample_ID, abun, -otu) %>%
  left_join(map_combined, by='sample_ID') %>%
  group_by(site,bean) %>%
  summarise(sampleNo=length(unique(sample_ID))) %>%
  arrange(desc(sampleNo)) %>%
  ggplot(aes(x=site,y=sampleNo, fill=bean)) +
  geom_bar(stat='identity', position='dodge')+
  theme_bw()+
  xlab('Growing location') +
  ylab('Counts') +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))

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

#########
#Figure 1
#########
fig1 <- ggplot(map.alpha[map.alpha$soil=='rhizosphere' & map.alpha$variable=="Richness",], aes(y=value, x=bean, color=bean))+
  scale_color_manual(values = c('darkorange', 'black')) + 
  theme_bw()+
  geom_boxplot() +
  xlab('Genotype')+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = 'none')

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
ax3.bean <- bean.pcoa$eig[3]/sum(bean.pcoa$eig)
ax4.bean <- bean.pcoa$eig[4]/sum(bean.pcoa$eig)

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
legend(-.03, -.08, legend = c('CELRK', 'Eclipse', 'Bulk soil'), col='black',pch=c(23,21, 22), cex=.8, box.col = 0, title='Genotype')


#are co-occurring species different when you only compare species at the same sites?
adonis(bean.dist~map_combined$bean, strata=map_combined$site)
#Are sites different when you only comparing the same species?
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


#Dissimilarity between the samples
df.sites=data.frame(name=c("CO", "NE", "MRF", "WA", "SVREC"),
                    lat=c(40.5,41.9, 43.3, 46.8, 43.4),
                    lon=c(-104.8,-103.8, -85.1, -121, -83.7))
site_distance_matrix <- round(GeoDistanceInMetresMatrix(df.sites))/1000

DISTdf <- data.frame(site1=as.factor(rownames(site_distance_matrix)),site_distance_matrix) %>%
  gather(site2, dist, -site1) %>%
  mutate(combined=paste(site1, site2, sep='-'))

BCdf <- data.frame(sample1=as.factor(rownames(as.matrix(bean.dist))),as.matrix(bean.dist)) %>%
  gather(sample2, BC, -sample1)

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

BCdistPlot <- BCdf %>% 
  left_join(DISTdf, by='combined') %>%
  filter(dist!=0.000) %>%
  ggplot(aes(x=log(dist), y=log(1-BC))) +
  geom_point() +
  labs(x='Geographic distance (ln(km))', y='Community similarity (ln(BC))') +
  stat_smooth(method = "lm", size = .8,level = .95)+
  theme_bw()

BCdistDF <- BCdf %>% 
  left_join(DISTdf, by='combined') %>%
  filter(dist!=0.000)
m <- lm(log(BCdistDF$BC) ~ log(BCdistDF$dist), BCdistDF)
summary(m)

######################################################
#Taxonomic distribution 
tax_v1 <- tax_us  

# creating column where last taxonomic identifier is attached to the OTU id 
features <- c(sprintf("OTU%05d", seq(1,28976)),"label")
tax_v1$otu_id <- features
tax_v1 <- tax_v1 %>% select(otu_id, everything())

lastValue <- function(x) tail(x[!is.na(x)], 1)
last_taxons<- apply(tax_v1, 1, lastValue)
tax_v1$last_taxon <- last_taxons
tax_v1$final_names <- paste(tax_v1$otu_id, tax_v1$last_taxon, sep='-')

tax_v1$otu <- rownames(tax_v1)

otu_us_rare_rhizo <- otu_us_rare
otu.rel.abun <- decostand(otu_us_rare_rhizo, method="total", MARGIN=2) #calculating relative abundance
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

##########
#Figure S2
##########
site_16S <- ggplot(temp_test,aes(sample_ID, n, fill=tax_names)) +
  geom_bar(stat='identity') +
  theme_bw()+
  labs(y="Relative abundace", x= NULL, fill='Phylum') +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=6),
        legend.position = 'bottom') +
  guides(fill=guide_legend(ncol=3)) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))


#######
#Using DESeq to determine taxa responing to the different fertilization types

taxonomy_phyloseq<- tax_v1[rownames(tax_v1) %in% rownames(otu_us),]
dim(otu_us)
OTU = otu_table(otu_us, taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxonomy_phyloseq))
SAM = sample_data(map_combined)
taxa_names(OTU)
physeq <- merge_phyloseq(phyloseq(OTU, TAX), SAM)

#No Vs Synthetic
physeq_synthetic = subset_samples(physeq, fertilization != "manure")
diagdds = phyloseq_to_deseq2(physeq_synthetic, ~ fertilization)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.001

sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_synthetic)[rownames(sigtab), ], "matrix"))
dim(sigtab)

sigtab$Genus <- as.character(sigtab$Genus)

sigtabgen <- sigtab %>%
  mutate(new_names = if_else((Genus == "uncultured bacterium" | Genus == "Ambiguous_taxa" | Genus == "uncultured"| is.na(Genus)), 'other', Genus))

# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$new_names, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$new_names = factor(as.character(sigtabgen$new_names), levels=names(x))
sigtabgen %>%
  mutate(regulation=if_else(sigtabgen$log2FoldChange>0, 'increased', 'decreased')) %>%
  filter(Phylum == 'Proteobacteria') %>%
  ggplot(aes(x=new_names, y=log2FoldChange, color=regulation)) + 
  geom_point(size=2) + 
  scale_color_manual(values=c('decreased'='#EFC000FF', 'increased'='#868686FF'))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  labs(color='Synthetic Vs\nNo fertilizer')+
  xlab('Genus')+
  coord_flip()
dim(sigtabgen)

syntVSno <- sigtabgen[sigtabgen$Phylum=='Proteobacteria',]

#No Vs Synthetic
physeq_manure = subset_samples(physeq, fertilization != "synthetic")
diagdds = phyloseq_to_deseq2(physeq_manure, ~ fertilization)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.001

sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_manure)[rownames(sigtab), ], "matrix"))
dim(sigtab)

sigtabgen <- sigtab %>%
  mutate(new_names = if_else((Genus == "uncultured bacterium" | Genus == "Ambiguous_taxa" | Genus == "uncultured"| is.na(Genus)), 'other', Genus))

sigtab$Genus <- as.character(sigtab$Genus)
sigtabgen <- sigtab %>%
  mutate(new_names = if_else((Genus == "uncultured bacterium" | Genus == "Ambiguous_taxa" | Genus == "uncultured"| is.na(Genus)), 'other', Genus))
dim(sigtabgen)

# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$new_names, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$new_names = factor(as.character(sigtabgen$new_names), levels=names(x))
sigtabgen %>%
  mutate(regulation=if_else(sigtabgen$log2FoldChange>0, 'increased', 'decreased')) %>%
  #filter(Phylum == 'Proteobacteria') %>%
  ggplot(aes(x=new_names, y=log2FoldChange, color=regulation)) + 
  geom_point(size=2) + 
  scale_color_manual(values=c('increased'='#EFC000FF', 'decreased'='#0073C2FF'))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  labs(color='No fertilizer\nVs manure')+
  xlab('Genus')+
  coord_flip()

manureVSno <- sigtabgen[sigtabgen$Phylum=='Proteobacteria',]

#Synthetic Vs Manure
physeq_fert = subset_samples(physeq, fertilization != "no fertilization")
diagdds = phyloseq_to_deseq2(physeq_fert, ~ fertilization)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.001

sigtab = res[which(res$padj < alpha), ]

sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_fert)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)

sigtab$Genus <- as.character(sigtab$Genus)
sigtabgen <- sigtab %>%
  mutate(new_names = if_else((Genus == "uncultured bacterium" | Genus == "Ambiguous_taxa" | Genus == "uncultured"| is.na(Genus)), 'other', Genus))

# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$new_names, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$new_names = factor(as.character(sigtabgen$new_names), levels=names(x))
sigtabgen %>%
  mutate(regulation=if_else(sigtabgen$log2FoldChange>0, 'increased', 'decreased')) %>%
  filter(Phylum == 'Proteobacteria') %>%
  ggplot(aes(x=new_names, y=log2FoldChange, color=regulation)) + 
  geom_point(size=2) + 
  scale_color_manual(values=c('increased'='#868686FF', 'decreased'='#0073C2FF'))+
  #theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=6))+
  labs(color='Synthetic\nVs manure')+
  xlab('Genus')+
  coord_flip()

manureVSsynth <- sigtabgen[sigtabgen$Phylum=='Proteobacteria',]

#***********************************************************************
rhizo_increased<- sigtabgen%>%
  filter(Genus == 'Rhizobium')

otu_us_rare_rhizo <- otu_us_rare[,map_combined$soil=='rhizosphere']
oturhizo.rel.abun <- decostand(otu_us_rare_rhizo, method="total", MARGIN=2) #calculating relative abundance

rhizo_tmp <- data.frame(otu=as.factor(rownames(oturhizo.rel.abun)),oturhizo.rel.abun) %>%
  gather(sample_ID, abun, -otu) %>%
  left_join(map_combined[,c('sample_ID','pH', 'state', 'bean', 'plot', 'site', 'fertilization')], by='sample_ID') %>%
  left_join(tax_v1, by='otu') %>%
  filter(final_names %in% rhizo_increased$final_names & abun!= 0)

ggplot(rhizo_tmp,aes(x=factor(fertilization,levels=c('no fertilization', 'synthetic', 'manure')), y=abun, group=otu)) +
  geom_point() +
  theme_bw()+
  labs(y="Relative abundance of\n Rhizobium sp.", x= "") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=6),
        legend.position = 'bottom')

manureVSsynth$otu <- as.character(manureVSsynth$otu)
manureVSno$otu <- as.character(manureVSno$otu)
syntVSno$otu <- as.character(syntVSno$otu)

deseqSelected <- c(manureVSsynth$otu, manureVSno$otu, syntVSno$otu)
deseqSelected <- unique(deseqSelected)

# Conduct indicator species analysis - need processing power! 
# (done on high performance computational cluster at MSU)
#indval_input <- otu_us_rare[rownames(otu_us_rare) %in% deseqSelected,]
#ind.fert <- indval(indval_input, map_combined$fertilization)

indval_results <- read.csv('Data/indvalspec_result.csv', header=T)

#Proportion of taxa responding to fertilization brouped by proteobacterial classes
indval_results %>%
  left_join(tax_v1, by='otu') %>%
  group_by(Phylum, fert) %>%
  mutate(n_treat=length(fert)) %>%
  group_by(Phylum,Class, fert) %>%
  summarise(count=length(otu)/unique(n_treat)) %>%
  ggplot(aes(x=Class,y=count)) +
  geom_bar(stat = 'identity') +
  labs(x=NULL, y='Proportion of taxa') +
  facet_grid(rows = vars(fert))+
  theme_bw() +
  theme(strip.text.x = element_text(size=12),
        strip.background = element_rect(colour="transparent", fill="transparent"),
        axis.text.x = element_text(angle=60, hjust = 1))

# Realtive abundance of taxa responding to fertilization grouped by proteobacterial classes
indvalOTUS <- data.frame(otu=rownames(otu.rel.abun), otu.rel.abun) %>%
  gather(sample_ID, abun, -otu) %>%
  filter(otu %in% indval_results$otu) %>%
  left_join(tax_v1, by='otu') %>%
  left_join(map_combined, by='sample_ID') %>%
  filter(soil != 'bulk') %>%
  mutate(new_names2=if_else(is.na(Genus)| Genus == 'uncultured' | Genus == 'Ambiguous_taxa' | Genus == 'uncultured bacterium', Class, Genus),
         class_v2=if_else(Class == 'Alphaproteobacteria', 'Alpha', Class),
         class_v2=if_else(class_v2 == 'Betaproteobacteria', 'Beta', class_v2),
         class_v2=if_else(class_v2 == 'Gammaproteobacteria', 'Gamma', class_v2),
         class_v2=if_else(class_v2 == 'Deltaproteobacteria', 'Delta', class_v2)) %>%
  #filter(new_names2 == 'Dyella') %>%
  group_by(new_names2) %>%
  mutate(count_taxa=length(unique(otu)),
         new_names3=paste(new_names2,' (',count_taxa,')',sep='')) %>%
  ggplot(aes(x=new_names3, y=abun, fill=new_names3)) +
  geom_boxplot() +
  facet_grid(fertilization~class_v2, scales = 'free_x', space='free')+
  labs(x=NULL, y=NULL) +
  theme_bw()+
  theme(legend.position = 'none',
        strip.text = element_blank(),
        strip.background = element_rect(colour="transparent", fill="transparent"),
        axis.text.x = element_text(angle=90, hjust = 1, vjust = .5),
        axis.text.y = element_blank())

# Presence-absence of taxa responding to fertilization grouped by proteobacterial classes
indvalPresence <- data.frame(otu=rownames(otu.rel.abun), otu.rel.abun) %>%
  gather(sample_ID, abun, -otu) %>%
  filter(otu %in% indval_results$otu) %>%
  left_join(tax_v1, by='otu') %>%
  left_join(map_combined, by='sample_ID') %>%
  left_join(indval_results, by='otu') %>%
  filter(soil != 'bulk') %>%
  mutate(new_names2=if_else(is.na(Genus)| Genus == 'uncultured' | Genus == 'Ambiguous_taxa' | Genus == 'uncultured bacterium', Class, Genus),
         class_v2=if_else(Class == 'Alphaproteobacteria', 'Alpha', Class),
         class_v2=if_else(class_v2 == 'Betaproteobacteria', 'Beta', class_v2),
         class_v2=if_else(class_v2 == 'Gammaproteobacteria', 'Gamma', class_v2),
         class_v2=if_else(class_v2 == 'Deltaproteobacteria', 'Delta', class_v2)) %>%
  group_by(class_v2,new_names2,fert) %>%
  summarise(count_taxa=length(unique(otu))) %>%
  mutate(fert = fct_relevel(fert, "synthetic", "no fertilization", "manure")) %>%
  ggplot(aes(size=count_taxa, x=new_names2, y=fert, fill=new_names2)) +
  geom_point(pch=21)+
  facet_grid(~class_v2, scales = 'free_x', space='free')+
  labs(x=NULL, y=NULL)+
  scale_y_discrete()+
  theme_bw()+
  theme(legend.position = 'none',
        strip.text.x = element_text(size=12),
        strip.background = element_rect(colour="transparent", fill="transparent"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

#joining figures
#grid.arrange(arrangeGrob(indvalPresence,indvalOTUS, ncol = 1, heights = c(.2,1.2)))

#otu_rare_rhizo_20 <- otu_us_rare_rhizo[rowSums(otu_us_rare_rhizo)>20,]

# Identifying taxa differently abundant between the plant genotypes (DESeq2)
map_rhizo_us <- map_combined[map_combined$soil!='bulk',]

physeq_rhizo = subset_samples(physeq, soil!= 'bulk')
diagdds_bean = phyloseq_to_deseq2(physeq_rhizo, ~ bean)
diagdds_bean = DESeq(diagdds_bean, test="Wald", fitType="parametric")
res_bean = results(diagdds_bean, cooksCutoff = FALSE)
alpha = 0.05

sigtab_bean = res_bean[which(res_bean$padj < alpha), ]
sigtab_bean = cbind(as(sigtab_bean, "data.frame"), as(tax_table(physeq_rhizo)[rownames(sigtab_bean), ], "matrix"))

sigtab_bean$Genus <- as.character(sigtab_bean$Genus)

sigtabgen_bean <- sigtab_bean %>%
  mutate(new_names = if_else((Genus == "uncultured bacterium" | Genus == "Ambiguous_taxa" | Genus == "uncultured"| is.na(Genus)), 'other', Genus))
#write.csv(sigtabgen_bean, 'TableS3.csv')

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
Fig1S <- vennDiagram(v) #2173 shared between sites - CORE?

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
#Occupancy abundance figure
#selecting variety unique OTUs
CELRK_rhizo_otu <- CELRK_rhizo_otu[rowSums(CELRK_rhizo_otu)>0,]
Eclipse_rhizo_otu <- Eclipse_rhizo_otu[rowSums(Eclipse_rhizo_otu)>0,]

CELRK_uniq <- CELRK_rhizo_otu[!(rownames(CELRK_rhizo_otu) %in% rownames(Eclipse_rhizo_otu)),]
Eclipse_uniq <- Eclipse_rhizo_otu[!(rownames(Eclipse_rhizo_otu) %in% rownames(CELRK_rhizo_otu)),]

#Occupancy abundance figure - combined data
otu_rare_bean <- rhizo_otu[rowSums(rhizo_otu)>0,]
bean_otu_PA <- 1*((otu_rare_bean>0)==1)
bean_otu_PA <- bean_otu_PA[rowSums(bean_otu_PA)>0,]
Occ <- rowSums(bean_otu_PA)/ncol(bean_otu_PA)
com_abund_bean <- rowSums(rhizo_otu)/ncol(rhizo_otu)

bean_df_occ <- data.frame(otu=names(Occ), occ=Occ) 
bean_df_abun <- data.frame(otu=names(com_abund_bean), abun=log10(com_abund_bean))
bean_df_abun$found <- 'shared'
bean_df_abun$found[bean_df_abun$otu %in% rownames(CELRK_uniq)] <- 'CELRK'
bean_df_abun$found[bean_df_abun$otu %in% rownames(Eclipse_uniq)] <- 'Eclipse'

bean_df_abun$col[bean_df_abun$found == 'CELRK'] <- 'darkorange'
bean_df_abun$col[bean_df_abun$found == 'Eclipse'] <- 'black'
bean_df_abun$col[bean_df_abun$found == 'shared'] <- 'lightgray'

bean_df_abun_new <- left_join(bean_df_abun, tax_v1)

occ_abun_bean <- left_join(bean_df_abun_new, bean_df_occ)

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

# nrow(combined_occ_data[combined_occ_data$bean_found=='shared',]) 
# nrow(combined_occ_data[combined_occ_data$bean_found=='CELRK',]) 
# nrow(combined_occ_data[combined_occ_data$bean_found=='Eclipse',]) 
# nrow(combined_occ_data)

##*********************************
#Sloan neutral model 
spp=t(otu_rare_rhizo)
taxon=as.vector(size_occ$taxonomy)

#Models for the whole community
obs.np=sncm.fit(spp, taxon, stats=FALSE, pool=NULL)
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

# sum(UpPredictOTU %in% US_occ_1)
# sum(DownPedictedOTU %in% US_occ_1)

plot(x=log10(obs.np$p), y=obs.np$freq, xlab="Log Abundance", ylab="Occurrence Frequency")
points(x=log10(obs.np$p[ap==TRUE]), y=obs.np$freq[ap==TRUE], col="#2D7DD2", pch=19)
points(x=log10(obs.np$p[bp==TRUE]), y=obs.np$freq[bp==TRUE], col="#F45D01", pch=19)
lines(obs.np$freq.pred~log10(obs.np$p), col="grey80", lty=1, lwd=6)
lines(obs.np$pred.upr~log10(obs.np$p), col="grey80", lty=1, lwd=3)
lines(obs.np$pred.lwr~log10(obs.np$p), col="grey80", lty=1, lwd=3)

#colour coding of the uniqueness and overplotting the Sloan model

(occ_abun_bean_fig <- ggplot(data=combined_occ_data, aes(x=abun, y=occ)) +
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
           fill = guide_legend(title=NULL)))

######################################################
#Core community 
#taking into consideration the 2734 otu shared accross the locations
length(Occ[Occ==1]) #258 with occupancy = 1
Occ_bact <- Occ
length(occ_abun_bean_rhizo$otu[occ_abun_bean_rhizo$occ==1]) #258
together_occ16S <- otu.rel.abun[rownames(otu.rel.abun) %in% occ_abun_bean_rhizo$otu[occ_abun_bean_rhizo$occ==1],] 

occ_16S <- data.frame(otu=as.factor(rownames(together_occ16S)),together_occ16S) %>%
  gather(sample_ID, abun, -otu) %>%
  left_join(map_combined[,c('sample_name','sample_ID','site','bean', 'soil')], by='sample_ID') %>%
  left_join(tax_v1, by='otu') %>%
  filter(soil != 'bulk') %>%
  group_by(sample_name,site, Phylum, Class, otu) %>%
  mutate(
    n=abun/length(abun))
#Many taxa in the core are novel and don't match to any ref in SILVA db (at Genus or even Class level)

species_level_core <- data.frame(otu=as.factor(rownames(together_occ16S)),together_occ16S) %>%
  gather(sample_ID, abun, -otu) %>%
  left_join(map_combined[,c('sample_name','sample_ID','site','bean', 'soil')], by='sample_ID') %>%
  left_join(tax_v1, by='otu') %>%
  filter(soil != 'bulk') %>%
  group_by(site,Phylum) %>%
  summarise(
    n=sum(abun)/length(unique(sample_ID))) %>%
  mutate(rel_abun=n/sum(n)) %>%
  arrange(desc(rel_abun)) %>%
  ggplot(aes(site, rel_abun, fill= Phylum)) +
  theme_bw()+
  geom_bar(stat='identity', colour='black') +
  labs(y="Relative abundace", x= "Sample locations") +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=8))

plot_core_16S <- data.frame(data.frame(otu=as.factor(rownames(together_occ16S)),together_occ16S)) %>%
  gather(sample_ID, abun, -otu) %>%
  left_join(map_combined[c('sample_ID','pH', 'state', 'bean', 'sample_name','soil')], by='sample_ID') %>%
  left_join(tax_v1, by='otu') %>%
  mutate(new_names = if_else((Order == "uncultured bacterium" | Order == "Ambiguous_taxa" | is.na(Order)), Phylum, Order),
         new_names = if_else((new_names == 'uncultured Acidobacteria bacterium'), 'Acidobacteria bacterium', new_names)) %>%
  filter(soil != 'bulk') %>%
  group_by(Phylum,new_names) %>%
  dplyr::summarise(
    n=sum(abun)/length(unique(sample_name))) %>%
  arrange(desc(n))

plot_core_16S$Order2 <- factor(plot_core_16S$new_names, levels = plot_core_16S$new_names[order(plot_core_16S$n)])

ggplot(plot_core_16S,aes(Order2, n, fill=Phylum)) +
  theme_bw()+
  geom_bar(stat = 'identity', color='black')+
  labs(y='Combined relative abundance', x='Order') +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=8),
        legend.position = 'none')+
  coord_flip()

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

#Colombia Sloan model
spp=t(otu_col_rare)
taxon=as.vector(rownames(otu_col_rare))

#Models for the whole community
obs.np=sncm.fit(spp, taxon, stats=FALSE, pool=NULL)
sta.np=sncm.fit(spp, taxon, stats=TRUE, pool=NULL)
sta.np.16S.COL <- sta.np

above.pred=sum(obs.np$freq > (obs.np$pred.upr), na.rm=TRUE)/sta.np$Richness
below.pred=sum(obs.np$freq < (obs.np$pred.lwr), na.rm=TRUE)/sta.np$Richness

ap = obs.np$freq > (obs.np$pred.upr)
bp = obs.np$freq < (obs.np$pred.lwr)

plot(x=log10(obs.np$p), y=obs.np$freq, xlab="Log Abundance", ylab="Occurrence Frequency")
points(x=log10(obs.np$p[ap==TRUE]), y=obs.np$freq[ap==TRUE], col="red", pch=19)
points(x=log10(obs.np$p[bp==TRUE]), y=obs.np$freq[bp==TRUE], col="blue", pch=19)
lines(obs.np$freq.pred~log10(obs.np$p), col="yellow", lty=1, lwd=6)
lines(obs.np$pred.upr~log10(obs.np$p), col="yellow", lty=1, lwd=3)
lines(obs.np$pred.lwr~log10(obs.np$p), col="yellow", lty=1, lwd=3)

#How many taxa in the Colombian core (Occ = 1) is found in the US dataset?
sum(occ_abun_col_rhizo$otu[occ_abun_col_rhizo$occ==1] %in% names_matching) #434
sum(combined_occ_data$otu[combined_occ_data$occ==1] %in% col_names_matching) #147

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

core_plot <- core_tmp_v2 %>%
  group_by(site) %>%
  mutate(n=n/sum(n)) %>%
  arrange(desc(n)) %>%
  ggplot(aes(x=site, y=n, fill= Order)) +
  geom_bar(stat='identity', colour='black') +
  theme_bw()+
  labs(y="Relative abundace", x= "Site") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=6),
        legend.position = 'bottom')

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

###Combined both datasets and plotting the rel abund of the global core taxa
core_col_tmp_v4$site <- 'Colombia'
core_tmp_v2$country <- 'US'
global_core_df <- rbind(core_tmp_v2, core_col_tmp_v4)

####
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

view_core <- core_all_v1 %>%
  group_by(Order,final_names) %>%
  summarise(otu_sum=length(otu))

core_all_v1 <- core_all_v1 %>%
  mutate(new_names=as.factor(new_names),
         new_names=fct_reorder(new_names, abun))

x = tapply(core_all_v1$abun, core_all_v1$new_names, function(x) mean(x))
x = sort(x, TRUE)
core_all_v1$new_names = factor(as.character(core_all_v1$new_names), levels=names(x))

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


core_all_v1 %>%
  #filter(Order == 'Micrococcales') %>%
  group_by(Order,Genus)%>%
  summarise(n=mean(abun)) %>%
  arrange(desc(-n))

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

colnames(venn_df) <- c("US", "US core", "COL core", "COL")
venn_df <- venn_df[rowSums(venn_df)>0,]
v_count=vennCounts(venn_df)
v_count_2=round(v_count[,"Counts"]/sum(v_count[,"Counts"]),2)
Fig1S <- vennDiagram(v_count) 


#How many global core taxa are uncultured?
tax_v1[tax_v1$otu %in% global_core,] %>%
  mutate(newNames=if_else(is.na(Genus), 'uncultured', Genus),
         newNames=if_else(newNames %in% c('uncultured bacterium','uncultured Acidobacteria bacterium', 'Ambiguous_taxa'), 'uncultured', newNames)) %>%
  ggplot(aes(newNames, fill=Phylum)) +
  geom_histogram(stat='count') +
  coord_flip() +
  labs(x='Genus') +
  theme_classic()+
  theme(legend.position = c(.5,.5))


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

# subOTUids <- US_fastq$V9[US_fastq$V10 %in% global_core]
# write.table(subOTUids,'~/Desktop/fatqCOREids.txt', sep='\t')

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

length(unique(subOTU_match$ZOTUid))

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

subOTU_maxAbund <- left_join(temp, temp_v2[c('refID','core', 'maxRel')])

plotData <- subOTU_match %>%
  group_by(refID) %>%
  summarise(nZOTUs=length(unique(ZOTUid))) %>%
  left_join(subOTU_maxAbund) %>%
  left_join(NoCOREzotus) %>%
  filter(core>0)%>%
  mutate(name=paste(refID," (",NoCoreZOTU,"/",nZOTUs,")", sep='')) 

plotData$name_new <- factor(plotData$name, levels = plotData$name[order(plotData$maxRel)])

plot_temp_data=temp_v2

plot_temp_data %>%
  left_join(plotData[,c('refID', 'name')]) %>%
  filter(!is.na(name)) %>%
  ggplot(aes(name, maxRel, col=as.factor(core)))+
  geom_point() +
  labs(x=NULL, y='Relative abundance', col='Occupancy') +
  theme_bw()+
  theme(axis.text = element_text(color='grey40'),
        axis.title = element_text(color='grey40'),
        axis.text.x = element_text(angle=90, size=10, hjust=1, vjust = 0.5, color='grey40'))

