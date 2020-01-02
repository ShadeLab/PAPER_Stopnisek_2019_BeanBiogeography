taxa_global <- read.table('Data/global_core.txt', sep='\t')
global_core<- as.character(taxa_global$x)

col_otu_new <- read.table('Data/OTU_table_COL.v2.txt', header=T, sep='\t', row.names = 1)
map_col_new <- read.table('Data/map_COL.v2.txt', head=T, sep='\t')
tax_col_new <- read.table('Data/tax_COL.v2.taxonomy', header=F, sep='\t')

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

set.seed(085)
Col_rare <- t(rrarefy(t(col_otu_new_filt), min(colSums(col_otu_new_filt)))) 
Col_rare <- Col_rare[rowSums(Col_rare)>0,]
tax_col_filt <- tax_col_new[tax_col_new$otu %in% rownames(Col_rare),]

map_agri <- map_col_new[map_col_new$study_Project=='old',]

unique(map_col_new$ID)

core_taxa_Col_v2_forest <- data.frame(otu=rownames(Col_rare),Col_rare)%>%
  gather(ID,abun, -otu) %>%
  left_join(map_col_new) %>%
  filter(Plant_status != 'Bulk_soil',
         Soil_type == 'Forest_soil') %>%
  filter(otu %in% global_core) %>%
  mutate(abun=if_else(abun>0, 1, 0)) %>%
  group_by(otu) %>%
  summarise(frequency=sum(abun)/length(unique(ID)))

core_taxa_Col_v2_agricultural <- data.frame(otu=rownames(Col_rare),Col_rare)%>%
  gather(ID,abun, -otu) %>%
  left_join(map_col_new) %>%
  filter(Plant_status != 'Bulk_soil',
         Soil_type != 'Forest_soil') %>%
  filter(otu %in% global_core) %>%
  mutate(abun=if_else(abun>0, 1, 0)) %>%
  group_by(otu) %>%
  summarise(frequency=sum(abun)/length(unique(ID)))
  
