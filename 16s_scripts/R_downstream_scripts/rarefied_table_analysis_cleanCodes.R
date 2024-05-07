
###########Preparation################
#load required packages
packages_to_load = c('stringr','phyloseq','vegan','tidyr','ggplot2','qiime2R',
                     'dplyr','ANCOMBC','usedist','rioja','ggpubr','spaa','ape',
                     'cowplot','ggforestplot','GGally','scales','metagMisc','ggsci','rstatix','wesanderson','RColorBrewer','jakR')
lapply(packages_to_load,require, character.only = T)

#generating color palettes
col1 = wes_palette("Rushmore1")
col2 = wes_palette("GrandBudapest2")
col3 = wes_palette('Darjeeling2')
col21 <- rev(c("tomato1","darkblue","turquoise1","lightblue","darkred","mediumblue","purple","bisque",
               "greenyellow","yellow","violetred2","darkgreen","darkgoldenrod1","deeppink3","cadetblue4",
               "orchid2","seagreen3","purple4","dodgerblue2","red","gray27"))

###########load variables##########
PFIG = FAS = 'R_analysis_output/rarefied_table'
#dir.create(PFIG,recursive = T)

#set random seed
set.seed(666)

############metadata pre-processing###################
metadata = read.csv('analysis_from_bcdPrm/metadata_parsed_validated.tsv',sep = '\t')
summary(metadata)
#keep the column type row
metadata_typerow = metadata[1,]
metadata = metadata[-1,]
metadata = select(metadata,c(`sample.id`,`Condition`,`Day.numeric`,`Day.categorical`,`Condition_day`,`Replicate`,`Experiment`))
metadata = metadata[order(metadata$Condition_day),]
metadata$Replicate_index = sequence(table(metadata$Condition_day))
metadata$Day.numeric = as.numeric(metadata$Day.numeric)
# function to collapse days into phases
collapse_phase = function(df) {
  df$phase = ifelse(df$Day.numeric == 0 , 'pre-infection',
                    ifelse(df$Day.numeric <= 5, 'early infection', 'late infection'))
  return(df)
}
metadata = collapse_phase(metadata)
write.table(metadata, paste(FAS,'metadata_selected_rows.csv',sep = '/'),sep = '\t', quote = F,row.names = F)

metadata %>%
  group_by(Day.categorical,Condition) %>%
  summarize(Sample_count = n()) -> metadata.summary

############input data using phyloseq package###############
#Create a phyloseq object
pr<-qza_to_phyloseq(
  features="analysis_from_bcdPrm/PE_bcdPrm_core-metrics-results/rarefied_table.qza",
  tree="analysis_from_bcdPrm/PE_bcdPrm_rooted-tree.qza",
  taxonomy="analysis_from_bcdPrm/PE_bcdPrm_taxonomy.qza",
  metadata=paste(FAS,'metadata_selected_rows.csv',sep = '/'))
pr
asv = data.frame(otu_table(pr)) #note that the asv is the raw table which has not been normalized by either rows or columns
tax = data.frame(tax_table(pr))
metadata= data.frame(sample_data(pr))

# parse metadata
summary(metadata)
metadata$Day.categorical = as.numeric(metadata$Day.categorical)
metadata$Day.categorical = ifelse(nchar(metadata$Day.categorical) == 1, paste0('0',metadata$Day.categorical),metadata$Day.categorical)
metadata$Replicate_index = as.character(metadata$Replicate_index)
#####################ASV table pre-processing##############################
#order the asv table based on the total counts across all samples
asv$sum = apply(asv,1,sum)
asv.o = asv[order(asv$sum,decreasing = T),]
head(asv.o)

#remove the ghost ASVs
asv.o %>% dplyr::filter(sum >0) %>% as.data.frame -> asv.of # nothing is removed
nrow(asv.of) #the total number of ASVs after rarefaction = 3716
sum(asv.of) #the total number of reads after rarefaction = 12591096

#simplify the ASV ID
nums <- sprintf('%0.4d', 1:nrow(asv.of)) ##creating string "00001", "00002"。。。
ASVnames <- paste("ASV",nums, sep="")
ASVconversion <- data.frame(ASV_ID = rownames(asv.of),ASVnames, stringsAsFactors = F) #make df with two columns
rownames(asv.of) = ASVconversion$ASVnames
write.csv(x = ASVconversion, file = paste(FAS,"/ASVconversion.csv", sep = ''),
          row.names = T)

# parse metadata
asv.of = select(asv.of, -sum)
asv.of = asv.of[,order(names(asv.of))]
metadata.o = metadata[order(rownames(metadata)),]
summary(names(asv.of) == rownames(metadata.o)) # sanity check
#calculate frequency -- normalize by sample
asv.off = decostand(asv.of,2,method = 'total')
colSums(asv.off) #col sum = 1

#############taxonomic table pre-processing################
#remove the ghost ASVs in the taxonomy table
tax = subset(tax, rownames(tax) %in% ASVconversion$ASV_ID) # nothing is removed
rownames(tax) = ASVconversion$ASVnames[match(rownames(tax),ASVconversion$ASV_ID)]
tax.o = tax[order(rownames(tax)),]

######################Taxonomic bar plot#####################
# function to collapse tax/asv table to a desired taxanomy rank
collapse_tax = function(rank='Phylum',ft=9){
  # rank='Phylum'
  # ft=9
  #aggregate the ASV table to the phylum level
  tax.o %>%
    dplyr::select(rank) %>%
    merge(x = ., y = asv.off, all.y = T, by = 0) -> asv.offp
  rownames(asv.offp) <- asv.offp$Row.names
  asv.offp = asv.offp[,-c(1)]
  #change NA in Phylum column to unknown
  asv.offp[[rank]][is.na(asv.offp[[rank]])] = 'Unknown'
  #select the desired taxa rank, and collapse onto the rank
  if (rank == 'Kingdom'){
    asv.offp = aggregate(. ~ Kingdom,asv.offp, sum)
  }  else if (rank == 'Phylum') {
    asv.offp = aggregate(. ~ Phylum,asv.offp, sum)
  } else if (rank == 'Class') {
    asv.offp = aggregate(. ~ Class,asv.offp, sum)
  } else if (rank == 'Order') {
    asv.offp = aggregate(. ~ Order,asv.offp, sum)
  } else if (rank == 'Family') {
    asv.offp = aggregate(. ~ Family,asv.offp, sum)
  } else if (rank == 'Genus') {
    asv.offp = aggregate(. ~ Genus,asv.offp, sum)
  } else if (rank == 'Species') {
    asv.offp = aggregate(. ~ Species,asv.offp, sum)
  }
  # sanity check
  if (all(round(colSums(asv.offp[,-c(1)])*100) == 100)) {
    rownames(asv.offp) = asv.offp[[rank]]
    asv.offp = select(asv.offp, select = -rank)
    #retain the ft most abundant taxa; calculate the sum of all the others
    asv.offp = asv.offp[order(rowSums(asv.offp), decreasing = T),]  
    asv.offps = bind_rows(asv.offp[c(1:ft),],colSums(asv.offp[c(c(ft+1):nrow(asv.offp)),]))
    rownames(asv.offps)[c(ft+1)] = 'All others'
    if (all(round(colSums(asv.offps)*100) ==100)) {
      rownames(asv.offps) -> asv.offps[[rank]]
      asv.offps.s = gather(asv.offps,"SD","frequency",-{{rank}})
      #merge the stacked ASV table with a subset of metadata to include some categorization columns
      metadata.o %>%
        dplyr::select(`Condition`,`Day.categorical`,`Replicate_index`) %>%
        merge(x = ., y = asv.offps.s, all.x = T, by.x = 0, by.y = "SD") -> asv.offps.s
      
      #coerce the phylum column to ordered factors for plotting
      asv.offps.s[[rank]] = factor(asv.offps.s[[rank]], levels = c(rownames(asv.offps)))
      asv.offps.s$Row.names = as.character(asv.offps.s$Row.names)
    } else {asv.offps.s = NULL}
  } else {asv.offps.s = NULL}
  return(list(asv.offp,asv.offps,asv.offps.s))
}

#make bar plot
asv.taxa.genus = collapse_tax(rank='Genus')
ggplot(data = asv.taxa.genus[[3]], mapping = aes(x = `Replicate_index`, y = `frequency`,fill = Genus)) +
  geom_bar(stat="identity",color = 'black') +
  scale_fill_manual(values = col21) +
  facet_grid(`Condition`~`Day.categorical`,scales = "free_x") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.3,'cm'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10, face = 'bold', color = 'black'),
        axis.text.y = element_text(size=15, face = 'bold', color = 'black'),
        axis.title=element_text(size=18,face="bold"),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"),
        strip.text.y = element_text(size = 15, color = "black", face = "bold"),
        strip.background = element_rect(color="black", fill="grey95", size=1.5, linetype="solid"))+
  xlab("Replicates") + ylab('Relative Abundance') +
  guides(fill=guide_legend(nrow=2))
ggsave(paste(PFIG,"/",'Taxonomy_top_9_genera_barplot.png',sep = ""),width = 10, height = 5, device = 'png')

#######################Alpha diversity analysis############################
# function to read in the diversity metrix
read_alpha_diversity = function(path_to_qza) {
  alpha = read_qza(path_to_qza)
  alpha <- alpha$data %>% tibble::rownames_to_column("SampleID") 
  alpha = merge(x = alpha,
                y = metadata.o,
                by.x = 'SampleID',
                by.y = 0,
                all.Y = T)
  return(alpha)
}
#shannon diversity matrix

shannon = read_alpha_diversity("analysis_from_bcdPrm/PE_bcdPrm_core-metrics-results/shannon_vector.qza")
evenness = read_alpha_diversity("analysis_from_bcdPrm/PE_bcdPrm_core-metrics-results/evenness_vector.qza")
faithpd = read_alpha_diversity("analysis_from_bcdPrm/PE_bcdPrm_core-metrics-results/faith_pd_vector.qza")
obsasv = read_alpha_diversity("analysis_from_bcdPrm/PE_bcdPrm_core-metrics-results/observed_features_vector.qza")

(plot.shannon=ggboxplot(shannon, x = 'Condition', y = 'shannon_entropy', fill = 'Condition', 
                        size = 1, add = "jitter")+
    xlab("Phase") + ylab('Shannon index') + ggtitle('Shannon Index') +
    coord_cartesian(ylim = c(0,max(shannon[,2]+0.5)))+
    theme(text = element_text(size = 22), legend.position = 'none')+
    scale_fill_lancet(alpha=0.8) +
    stat_compare_means(method = "wilcox.test",label.x = 0.7,label.y = 0.05,size = 6)) 


(plot.evenness=ggboxplot(evenness, x = 'Condition', y = 'pielou_evenness', fill = 'Condition', 
                         size = 1, add = "jitter")+
    xlab("Phase") + ylab('Pielou evenness') + ggtitle('Pielou Evenness') +
    coord_cartesian(ylim = c(0,max(evenness[,2]+0.5)))+
    theme(text = element_text(size = 22), legend.position = 'none')+
    scale_fill_lancet(alpha=0.8) +
    stat_compare_means(method = "wilcox.test",label.x = 0.7,label.y = 0.05,size = 6)) 


(plot.obsasv=ggboxplot(obsasv, x = 'Condition', y = 'observed_features', fill = 'Condition', 
                       size = 1, add = "jitter")+
    xlab("Phase") + ylab('Observed features') + ggtitle('Observed Features') +
    coord_cartesian(ylim = c(0,max(obsasv[,2]+10)))+
    theme(text = element_text(size = 22), legend.position = 'none')+
    scale_fill_lancet(alpha=0.8) +
    stat_compare_means(method = "wilcox.test",label.x = 0.7,label.y = 0.5,size = 6)) 


(plot.faithpd=ggboxplot(faithpd, x = 'Condition', y = 'faith_pd', fill = 'Condition', 
                        size = 1, add = "jitter")+
    xlab("Phase") + ylab('Faith PD') + ggtitle('Faith PD') +
    coord_cartesian(ylim = c(0,max(faithpd[,2]+0.5)))+
    theme(text = element_text(size = 22), legend.position = 'none')+
    scale_fill_lancet(alpha=0.8) +
    stat_compare_means(method = "wilcox.test",label.x = 0.7,label.y = 0.5,size = 6)) 

plot_grid(plot.shannon,plot.evenness,plot.obsasv,plot.faithpd,
          nrow = 2) 
ggsave(paste(PFIG,"/",'alpha_diversity.png',sep = ""),width = 25, height = 20, device = 'png',limitsize = FALSE)

###########################PCoA plot using the Weighted UniFrac distance###############
#read the distance matrix output by QIIME2
dm = read_qza("analysis_from_bcdPrm/PE_bcdPrm_core-metrics-results/weighted_unifrac_distance_matrix.qza")
#select the rows we need
dm = dist_subset(dm$data,rownames(metadata.o))
#some distance measures may result in negative eigenvalues. In that case, add a correction:
pcoa.wunifrac = pcoa(dm,correction = 'cailliez')
#transform the pcoa.wunifrac object into a df for making it plottable with ggplot2
coordinates.pcoa.wunifrac = data.frame(pcoa.wunifrac$vectors,check.names = FALSE)

#plot for all samples
#Adonis test
metadata.o$Condition = as.factor(metadata.o$Condition)
all(rownames(metadata.o) == names(dm))
adonis.r = adonis2(dm~Condition,metadata.o,permut = 1000)
adonis.p = round(adonis.r$`Pr(>F)`[1],digit = 4)

betadisp.r = permutest(betadisper(dm,metadata.o$Condition),permutations = 1000)
betadisp.p = round(betadisp.r$tab[[6]][1],digit = 4)
#plot
color = metadata.o$Condition
(pcoa=ggplot()+
    geom_point(data=subset(coordinates.pcoa.wunifrac[,c(1,2)]),
               mapping = aes(x = `Axis.1`, y = `Axis.2`, colour = color),
               alpha = 0.8, size = 3) + 
    labs(title = 'Weighted UniFrac Distance',
         subtitle = paste0('p(adonis) = ', adonis.p, ', p(betadisp) = ', betadisp.p))+
    stat_ellipse(data=subset(coordinates.pcoa.wunifrac[,c(1,2)]),
                 mapping = aes(x = `Axis.1`, y = `Axis.2`, colour = color),
                 type = 'norm') +
    xlab(paste('PC1 (',label_percent()(pcoa.wunifrac$values$Rel_corr_eig[1]),')',sep = '')) +
    ylab(paste('PC2 (',label_percent()(pcoa.wunifrac$values$Rel_corr_eig[2]),')',sep = '')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = 'black'),
          legend.title=element_blank(),
          legend.position = 'bottom',
          legend.text = element_text(size=15))+
    scale_color_jama())
ggsave(paste(PFIG,"/",'pcoa_weightedUniFrac.coloredByCondition.png',sep = ""),width = 4, height = 4, device = 'png')

###############Session information############
sessionInfo()
devtools::session_info()
