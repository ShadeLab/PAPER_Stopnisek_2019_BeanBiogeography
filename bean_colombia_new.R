col_otu_new <- read.table('../bean_microbiota/new_Col/OTU_table.txt', header=T, sep='\t', row.names = 1)
map_col_new <- read.table('../bean_microbiota/new_Col/map.txt', head=T, sep='\t')
tax_col_new <- read.table('../bean_microbiota/new_Col/full_rep_set.nr_v128.wang.taxonomy', header=F, sep='\t')
dim(tax_col_new)

names(tax_col_new) <- c('otu','taxonomy')
tax_col_new <- tax_col_new %>%
  separate(taxonomy, into=c("Kingdom", "Phylum", "Class", 
                            "Order", "Family", "Genus", "Species"), sep=";", remove=T)

tax_col_new <- tax_col_new[!grepl("Mitochondria", tax_col_new$Family),]
tax_col_new <- tax_col_new[!grepl("Chloroplast", tax_col_new$Class),]
tax_col_new <- tax_col_new[!grepl("unknown", tax_col_new$Kingdom),]
tax_col_new <- tax_col_new[!grepl("Eukaryota", tax_col_new$Kingdom),]

tax_col_new$otu <- as.character(tax_col_new$otu)
col_otu_new_filt <- col_otu_new[rownames(col_otu_new) %in% tax_col_new$otu,]
tax_col_new <- tax_col_new[tax_col_new$otu %in% rownames(col_otu_new),]

dim(tax_col_new)
dim(col_otu_new_filt)

colSums(col_otu_new_filt)
set.seed(085)
Col_rare <- t(rrarefy(t(col_otu_new_filt), min(colSums(col_otu_new_filt)))) 
Col_rare <- Col_rare[rowSums(Col_rare)>0,]
dim(Col_rare)
tax_col_filt <- tax_col_new[tax_col_new$otu %in% rownames(Col_rare),]

map_agri <- map_col_new[map_col_new$study_Project=='old',]

count(map_col_new$Plant_status == 'Bulk_soil')
unique(map_col_new$ID)
core_taxa_Cel_new <- data.frame(otu=rownames(Col_rare),Col_rare)%>%
  gather(ID,abun, -otu) %>%
  left_join(map_col_new) %>%
  filter(Plant_status != 'Bulk_soil',
         Soil_type == 'Forest_soil') %>%
  filter(otu %in% global_core) %>%
  mutate(abun=if_else(abun>0, 1, 0)) %>%
  group_by(otu) %>%
  summarise(frequency=sum(abun)/length(unique(ID)))
  
summary(core_taxa_Cel_new)
write.csv(core_taxa_Cel_new,'../../../Desktop/Col_new_core_occupancy_rhizosphre_forest.csv')  
nrow(tax_col_new[tax_col_new$otu %in% global_core,])
