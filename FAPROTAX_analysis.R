UpPredictOTU 
DownPedictedOTU 
global_core
OTUus<- otu_rare_rhizo[rowSums(otu_rare_rhizo)>0,]
neutral <- rownames(OTUus)[!(rownames(OTUus) %in% c(UpPredictOTU, DownPedictedOTU, global_core))]


OTUtable <- read.table('US/OTU_table_final_US.txt', sep='\t', row.names = 1, header=T)
OTUtable$SVERCc1 <- NULL

OTU_up <- OTUtable[rownames(OTUtable)%in%UpPredictOTU,]
OTU_upTAX <- OTU_up[39]
OTU_up$taxonomy <- NULL

OTU_down <- OTUtable[rownames(OTUtable) %in% DownPedictedOTU,]
OTU_downTAX <- OTU_down[39]
OTU_down$taxonomy <- NULL

OTU_core <- OTUtable[rownames(OTUtable) %in% global_core,]
OTU_coreTAX <- OTU_core[39]
OTU_core$taxonomy <- NULL

OTU_neutral <- OTUtable[rownames(OTUtable) %in% neutral,]
OTU_neutralTAX <- OTU_neutral[39]
OTU_neutral$taxonomy <- NULL

write.table(OTU_up, '../../../Desktop/OTU_up.txt', sep='\t')
write.table(OTU_up, '../../../Desktop/OTU_upTAX_v2.txt', sep='\t')

write.table(OTU_down, '../../../Desktop/OTU_down.txt', sep='\t')
write.table(OTU_down, '../../../Desktop/OTU_downTAX_v2.txt', sep='\t')

write.table(OTU_neutral, '../../../Desktop/OTU_neutral.txt', sep='\t')
write.table(OTU_neutral, '../../../Desktop/OTU_neutralTAX_v2.txt', sep='\t')

write.table(OTU_core, '../../../Desktop/OTU_core.txt', sep='\t')
write.table(OTU_core, '../../../Desktop/OTU_coreTAX_v2.txt', sep='\t')

CORE_norm<- read.table('../../../Desktop/functional_table_CORE.txt', sep='\t', header =T, row.names = 1)
UP_norm<- read.table('../../../Desktop/functional_table_UP.txt', sep='\t', header =T, row.names = 1)
DOWN_norm<- read.table('../../../Desktop/functional_table_DOWN.txt', sep='\t', header =T, row.names = 1)
NEUTRAL_norm<- read.table('../../../Desktop/functional_table_NEUTRAL.txt', sep='\t', header =T, row.names = 1)

library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)

CORE_norm <- CORE_norm[-91,]
core.fun.rel.abun <- decostand(CORE_norm, method="total", MARGIN=2) 
CORE_df <- data.frame(fun=as.factor(rownames(core.fun.rel.abun)), core.fun.rel.abun) %>%
  gather(sample_ID, abun, -fun) %>%
  mutate(dataset='core') %>%
  left_join(map_v2[,c('sample_ID', 'site', 'soil')], by='sample_ID') %>%
  filter(soil != 'bulk') %>%
  group_by(fun, dataset) %>%
  summarise(rel_abun=sum(abun)/length(sample_ID)) %>%
  arrange(desc(rel_abun))

UP_norm <- UP_norm[-91,]
up.fun.rel.abun <- decostand(UP_norm, method="total", MARGIN=2) 
UP_df <- data.frame(fun=as.factor(rownames(up.fun.rel.abun)), up.fun.rel.abun) %>%
  gather(sample_ID, abun, -fun) %>%
  mutate(dataset='up') %>%
  left_join(map_v2[,c('sample_ID', 'site', 'soil')], by='sample_ID') %>%
  filter(soil != 'bulk') %>%
  group_by(fun, dataset) %>%
  summarise(rel_abun=sum(abun)/length(sample_ID)) %>%
  arrange(desc(rel_abun))

DOWN_norm <- DOWN_norm[-91,]
down.fun.rel.abun <- decostand(DOWN_norm, method="total", MARGIN=2) 
DOWN_df <- data.frame(fun=as.factor(rownames(down.fun.rel.abun)), down.fun.rel.abun) %>%
  gather(sample_ID, abun, -fun) %>%
  mutate(dataset='down') %>%
  left_join(map_v2[,c('sample_ID', 'site', 'soil')], by='sample_ID') %>%
  filter(soil != 'bulk') %>%
  group_by(fun, dataset) %>%
  summarise(rel_abun=sum(abun)/length(sample_ID)) %>%
  arrange(desc(rel_abun))

NEUTRAL_norm <- NEUTRAL_norm[-91,]
neutral.fun.rel.abun <- decostand(NEUTRAL_norm, method="total", MARGIN=2) 
NEUTRAL_df <- data.frame(fun=as.factor(rownames(neutral.fun.rel.abun)), neutral.fun.rel.abun) %>%
  gather(sample_ID, abun, -fun) %>%
  mutate(dataset='neutral') %>%
  left_join(map_v2[,c('sample_ID', 'site', 'soil')], by='sample_ID') %>%
  filter(soil != 'bulk') %>%
  group_by(fun, dataset) %>%
  summarise(rel_abun=sum(abun)/length(sample_ID)) %>%
  arrange(desc(rel_abun))

faprotax_df <- rbind(CORE_df, UP_df, DOWN_df, NEUTRAL_df)

filter_fun <- rbind(CORE_df, UP_df, DOWN_df, NEUTRAL_df) %>%
  group_by(fun) %>%
  summarise(sum=sum(rel_abun)) %>%
  filter(sum>0)

faprotax_df <- faprotax_df[as.character(faprotax_df$fun) %in% as.character(filter_fun$fun),]

faprotax_df$Percentage <- round(faprotax_df$rel_abun * 100, 2)
faprotax_df[faprotax_df$Percentage == 0, 'Percentage'] <- NA
faprotax_df$dataset <- factor(faprotax_df$dataset)
levels(faprotax_df$dataset) <- c("core", "up", "down", "neutral")

ggplot(faprotax_df, aes(x=dataset, y=fun, fill=dataset, size=Percentage)) +
  geom_point(shape=21, alpha=.5) +
  #facet_grid(~dataset, scales='free_x') +
  scale_fill_brewer(guide='none', palette='Dark2') +
  scale_size_area(max_size = 10, name=element_text("Functional Abundance (%)")) +
  ylab("Functional Group") +
  xlab(NULL)+
  theme_bw()+
  theme(legend.position='bottom', legend.title = element_text(size=10))

faprotax_test=faprotax_df
faprotax_test[faprotax_test$rel_abun == 0, 'rel_abun'] <- NA
data_wide <- spread(faprotax_test, dataset, rel_abun)

core_diff=log2(data_wide[,2]/data_wide[,5])
down_diff=log2(data_wide[,4] / data_wide[,5])
up_diff=log2(data_wide[,3] / data_wide[,5])

diff_df=cbind(core_diff, up_diff,down_diff)
diff_df$fun <- data_wide$fun

diff_data <- gather(diff_df, dataset, diff_abundance, core:down)

diff_data_filt=diff_data %>%
  mutate(diff_abundance=if_else(is.na(diff_abundance), 0, diff_abundance)) %>%
  group_by(fun) %>%
  filter(sum(diff_abundance)!=0,
         !(fun %in% c('human_gut', 'mammal_gut', 'human_pathogens_all', 'human_pathogens_pneumonia')))

ggplot(diff_data_filt[diff_data_filt$diff_abundance!=0,],aes(x=dataset, y=diff_abundance, fill=dataset, group=dataset)) +
  #geom_point(position = position_dodge(width=0.9))+
  geom_bar(position = 'dodge', stat = 'identity')+
  geom_hline(yintercept=0)+
  facet_wrap(~fun, ncol = 5)+
  theme_classic() +
  labs(x=NULL, y='Differential abundance to neutral taxa (log2)')

dot_plot <- diff_df %>%
  filter(fun %in% unique(as.character(diff_data_filt$fun))) %>%
  mutate(core=if_else(is.na(core), 0, core),
         up=if_else(is.na(up), 0, up),
         down=if_else(is.na(down), 0, down)) %>%
  arrange(core) 

df_dotplot<- as.data.frame(dot_plot) %>%
  mutate(fun=as.character(fun),
         fun = fct_inorder(fun))
df_dotplot[df_dotplot$core == 0, 'core'] <- NA
df_dotplot[df_dotplot$up == 0, 'up'] <- NA
df_dotplot[df_dotplot$down == 0, 'down'] <- NA

onlyCORE <- df_dotplot[c(1,4)]
onlyCORE <- onlyCORE[!is.na(onlyCORE$core),]

ggplot(onlyCORE) +
  
  theme_classic() +  
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.title.x=element_text(size=8, color = 'grey30')) +
  
  geom_point(aes(x = -3, y = df_dotplot$fun),size = 0, col = "white") + 
  
  # add the horizontal discipline lines
  geom_hline(yintercept = 1:paste(length(dot_plot$fun)), col = "grey80", linetype = 2, alpha=.3) +
  geom_vline(xintercept = 0, col="grey80") +
  geom_vline(xintercept = c(-4,-3,-2,-1,1,2,3,4,5), col="grey80", linetype = 3) +
  
  geom_point(aes(x = core, y = fun), size = 7, col = "#97CC04", alpha=.8) +
  geom_point(aes(x = up, y = fun), size = 7, col = "#F45D01", alpha=.8) +
  geom_point(aes(x = down, y = fun), size = 7, col = "#2D7DD2", alpha=.8) +
  
  geom_point(aes(x = 4, y = 4.5), size = 4, col = "#97CC04") +
  geom_point(aes(x = 4, y = 3), size = 4, col = "#F45D01") +
  geom_point(aes(x = 4, y = 1.5), size = 4, col = "#2D7DD2") +

  scale_x_continuous(breaks=seq(-4,5,1)) +
  #adding text
  geom_text(aes(x = up, y = fun, 
                label = paste0(round(up, 2))),
            col = "white", size=2.5) +
  geom_text(aes(x = core, y = fun, 
                label = paste0(round(core, 2))),
            col = "white",size=2.5) +
  geom_text(aes(x = down, y = fun, 
                label = paste0(round(down, 2))),
            col = "white", size=2.5) +
  # add a label above the first two points
  #geom_text(aes(x = x, y = y, label = label, col = label),
  #          data.frame(x = c(2.47+.4, 1.27-.5,1.68), y = 27, 
  #                     label = c("core", "up", 'down')), size = 4) +
  #scale_color_manual(values = c("#97CC04", "#F45D01", "#2D7DD2"), guide = "none")+
  xlab("Predicted functions abundance\nexpressed as log2-fold difference\nto taxa fitting to the neutral model")


ggplot(onlyCORE) +
  
  theme_classic() +  
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.title.x=element_text(size=10, color = 'grey30')) +
  
  geom_point(aes(x = -1, y = onlyCORE$fun),size = 0, col = "white") + 
  
  # add the horizontal discipline lines
  geom_hline(yintercept = 1:paste(length(onlyCORE$fun)), col = "grey80", linetype = 2) +
  geom_vline(xintercept = 0, col="grey30", size=1.2) +
  geom_vline(xintercept = c(-1,1,2), col="grey30", linetype = 3) +
  
  geom_point(aes(x = core, y = fun), size = 12, col = "#97CC04") +
  
  scale_x_continuous(breaks=seq(-4,5,1)) +
  #adding text
  geom_text(aes(x = core, y = fun, 
                label = paste0(round(core, 2))),
            col = "white",size=4) +
xlab("Predicted functions abundance expressed as log2-fold difference\nto taxa fitting to the neutral model")


