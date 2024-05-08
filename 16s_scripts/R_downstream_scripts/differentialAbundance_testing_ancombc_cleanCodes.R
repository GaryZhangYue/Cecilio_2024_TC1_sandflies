
###########Preparation################
#load required packages
packages_to_load = c('stringr','phyloseq','vegan','tidyr','ggplot2','qiime2R',
                     'dplyr','ANCOMBC','usedist','rioja','ggpubr','spaa','ape',
                     'cowplot','ggforestplot','GGally','scales','metagMisc','ggsci','rstatix','wesanderson','RColorBrewer','jakR','iCAMP')

lapply(packages_to_load,require, character.only = T)

#generating color palettes
col1 = wes_palette("Rushmore1")
col2 = wes_palette("GrandBudapest2")
col3 = wes_palette('Darjeeling2')
col21 <- rev(c("tomato1","darkblue","turquoise1","lightblue","darkred","mediumblue","purple","bisque",
               "greenyellow","yellow","violetred2","darkgreen","darkgoldenrod1","deeppink3","cadetblue4",
               "orchid2","seagreen3","purple4","dodgerblue2","red","gray27"))

###########load variables##########
PFIG = FAS = 'R_analysis_output/unrarefied_table'
#dir.create(PFIG,recursive = T)

#set random seed
set.seed(666)

###############################functions###############################
# manually calculate log bias-corrected absolute abundance;
lca = function(pu.d.ac,date) {
  
  samp_frac = pu.d.ac$samp_frac
  # Replace NA with 0
  samp_frac[is.na(samp_frac)] = 0 
  # Add pesudo-count (1) to avoid taking the log of 0
  log_obs_abn = log(pu.d.ac$feature_table + 1)
  # Adjust the log observed abundances
  log_corr_abn = t(t(log_obs_abn) - samp_frac)
  # calculate mean and std
  #TC1
  serratia.avg.TC1 = mean(log_corr_abn['Genus:Serratia',which(startsWith(colnames(log_corr_abn),'TC1'))])
  serratia.sd.TC1 = sd(log_corr_abn['Genus:Serratia',which(startsWith(colnames(log_corr_abn),'TC1'))])
  delftia.avg.TC1 = mean(log_corr_abn['Genus:Delftia',which(startsWith(colnames(log_corr_abn),'TC1'))])
  delftia.sd.TC1 = sd(log_corr_abn['Genus:Delftia',which(startsWith(colnames(log_corr_abn),'TC1'))])
  #CTRL
  serratia.avg.CTRL = mean(log_corr_abn['Genus:Serratia',which(startsWith(colnames(log_corr_abn),'CTRL'))])
  serratia.sd.CTRL = sd(log_corr_abn['Genus:Serratia',which(startsWith(colnames(log_corr_abn),'CTRL'))])
  delftia.avg.CTRL = mean(log_corr_abn['Genus:Delftia',which(startsWith(colnames(log_corr_abn),'CTRL'))])
  delftia.sd.CTRL = sd(log_corr_abn['Genus:Delftia',which(startsWith(colnames(log_corr_abn),'CTRL'))])
  return(data.frame(Genus = c('Genus:Serratia','Genus:Serratia','Genus:Delftia','Genus:Delftia'),
                    Day = c(date,date,date,date),
                    Condition = c('CTRL','TC1','CTRL','TC1'),
                    mean = c(serratia.avg.CTRL,serratia.avg.TC1,delftia.avg.CTRL,delftia.avg.TC1),
                    std = c(serratia.sd.CTRL,serratia.sd.TC1,delftia.sd.CTRL,delftia.sd.TC1)))
}

####################################Data preprocessing###############################
#use the un-rarefied table for ANCOMBC analysis as suggested by the authors of ANCOM-BC
pu<-qza_to_phyloseq(
  features="analysis_from_bcdPrm/PE_bcdPrm_dada2.qza",
  tree="analysis_from_bcdPrm/PE_bcdPrm_rooted-tree.qza",
  taxonomy="analysis_from_bcdPrm/PE_bcdPrm_taxonomy.qza",
  metadata='R_analysis_output/rarefied_table/metadata_selected_rows.csv')
pu.asv = data.frame(otu_table(pu),check.names = F)
pu.tax = data.frame(tax_table(pu),check.names = F)
metadata.o = data.frame(sample_data(pu),check.names = F)
pu.tre = phy_tree(pu)

#remove the ghost ASVs
dim(subset(pu.asv, rowSums(pu.asv)==0)) # no ghost ASVs

# sanity check using the a function from iCAMP
#check if sample id matches between metadata and the ASV table
pu.sampleid.check = match.name(cn.list=list(comm=pu.asv), rn.list = list(env=metadata.o))
#check if the ASV id matches across the ASV table, phylogenetic tree tips and taxonomic table
pu.sampleid.check = match.name(rn.list = list(asv=pu.asv,tax=pu.tax),tree.list=list(tree=pu.tre)) 

res.all = data.frame(NULL)
log_corr_abn.serratia_delftia <- data.frame(genus = NULL, mean = NULL, std = NULL, date = NULL)

############################### TC1 vs Control on different days ###############################
for (i in unique(metadata.o$Day.numeric)) {
  print(i)
  pu.d = subset_samples(pu, Day.numeric ==  i)
  pu.d.ac = ancombc(phyloseq = pu.d, #the input data (phyloseq object)
                    formula = 'Condition', #the character string expresses how microbial absolute abundances for each taxon depend on the variables in metadata.
                    tax_level = "Genus", #character. The taxonomic level of interest. The input data can be agglomerated at different taxonomic levels
                    p_adj_method = "holm", #character. method to adjust p-values. Default is "holm".
                    prv_cut = 0, #a numerical fraction between 0 and 1. Taxa with prevalences less than prv_cut will be excluded in the analysis. Default is 0.10.
                    lib_cut = 0, #a numerical threshold for filtering samples based on library sizes. Samples with library sizes less than lib_cut will be excluded in the analysis. Default is 0, i.e. do not discard any sample.
                    group = NULL, #character. the name of the group variable in metadata. group should be discrete. Specifying group is required for detecting structural zeros and performing global test. Default is NULL. If the group of interest contains only two categories, leave it as NULL.
                    struc_zero = F, #logical. whether to detect structural zeros based on group. Default is FALSE.
                    neg_lb = F, #logical. whether to classify a taxon as a structural zero using its asymptotic lower bound. Default is FALSE.
                    tol = 1e-5, #numeric. the iteration convergence tolerance for the E-M algorithm. Default is 1e-05.
                    max_iter = 100, #numeric. the maximum number of iterations for the E-M algorithm. Default is 100.
                    conserve = TRUE, #logical. whether to use a conservative variance estimator for the test statistic. It is recommended if the sample size is small and/or the number of differentially abundant taxa is believed to be large. Default is FALSE.
                    alpha = 0.05, #numeric. level of significance. Default is 0.05.
                    global = F, #logical. whether to perform the global test. Default is FALSE.
                    n_cl = 1, #numeric. The number of nodes to be forked. For details, see ?parallel::makeCluster. Default is 1 (no parallel computing).
                    verbose = TRUE) #logical. Whether to generate verbose output during the ANCOM-BC fitting process. Default is FALSE.
  
  res = as.data.frame(pu.d.ac$res)
  res$Day = i
  res.all = bind_rows(res.all,res)
  log_corr_abn.serratia_delftia = bind_rows(log_corr_abn.serratia_delftia,lca(pu.d.ac,i))
}

res.all.selectedCol <- res.all %>%
  select("Day","lfc.taxon","lfc.ConditionTC1","se.ConditionTC1","p_val.ConditionTC1","q_val.ConditionTC1")
names(res.all.selectedCol)[2] = 'Genus'
res.all.selectedCol.Delfia_Serratia <- res.all.selectedCol %>%
  subset(Genus %in% c("Genus:Delftia","Genus:Serratia"))
write.csv(res.all.selectedCol,paste0(PFIG,'/ancombc2_allGenera.csv'))         
write.csv(res.all.selectedCol.Delfia_Serratia,paste0(PFIG,'/ancombc2_Delfia-Serratia.csv'))

############################### Plot the log change with p value annotated ###############################
res.all.selectedCol.Delfia_Serratia$Day = factor(res.all.selectedCol.Delfia_Serratia$Day,levels = c('0','2','5','7','9','12'))

res.all.selectedCol.Delfia_Serratia$p_asterisks = 
  ifelse(res.all.selectedCol.Delfia_Serratia$p_val.ConditionTC1 > 0.05, 'ns',
         ifelse(res.all.selectedCol.Delfia_Serratia$p_val.ConditionTC1 > 0.01, "*",
                ifelse(res.all.selectedCol.Delfia_Serratia$p_val.ConditionTC1 > 0.001, '**', '***')))

res.all.selectedCol.Delfia_Serratia$q_asterisks = 
  ifelse(res.all.selectedCol.Delfia_Serratia$q_val.ConditionTC1 > 0.05, 'ns',
         ifelse(res.all.selectedCol.Delfia_Serratia$q_val.ConditionTC1 > 0.01, "*",
                ifelse(res.all.selectedCol.Delfia_Serratia$q_val.ConditionTC1 > 0.001, '**', '***')))


ggplot(data = res.all.selectedCol.Delfia_Serratia, 
       mapping = aes(x = Day,y = lfc.ConditionTC1,fill = Genus))+
  geom_bar(stat = 'identity',
           position=position_dodge())+
  geom_errorbar(mapping = aes(ymin = lfc.ConditionTC1 - se.ConditionTC1,
                              ymax = lfc.ConditionTC1 + se.ConditionTC1),
                width = 0.3,
                size = 0.4,
                position = position_dodge(width = 1.5)) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5)+
  facet_grid(.~Genus, scales = "free_x")+
  labs(x='Day', y="log(TC1/CTRL)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        text = element_text(size = 16, colour = 'black'),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"),
        strip.text.y = element_text(size = 15, color = "black", face = "bold")) +
  geom_text(mapping = aes(x = Day, y = lfc.ConditionTC1 + 1.5*sign(lfc.ConditionTC1), label = p_asterisks),
            position = position_dodge(width = 1),
            size = 3, color = 'red')+
  scale_fill_uchicago(alpha=0.8)+
  ylim(-0.5,6)
ggsave(paste(PFIG,'ancombc_delfia-serratia_pval_annot.png',sep = '/'),width = 8,height = 8, device = 'png')

############################### Plot the bias-corrected absolute abundance ###############################
log_corr_abn.serratia_delftia$Day = as.factor(log_corr_abn.serratia_delftia$Day)
test <- res.all.selectedCol.Delfia_Serratia
test$`.y.` = 'mean'
test$group1 = 'CTRL'
test$group2 = 'TC1'
test$xmin = rep(1:6,each=2)
test$xmax = rep(1:6,each=2)
names(test)[5] = 'p'
names(test)[6] = 'p.adj'
names(test)[7] = "p.adj.signif"
test$Condition = 'CTRL'

log_corr_abn.serratia_delftia$merge_anchor = paste0(log_corr_abn.serratia_delftia$Genus,'_',log_corr_abn.serratia_delftia$Day)
log_corr_abn.serratia_delftia$y.position = log_corr_abn.serratia_delftia$mean+log_corr_abn.serratia_delftia$std+0.5
test$merge_anchor  = paste0(test$Genus,'_',test$Day)
test = merge(test,subset(log_corr_abn.serratia_delftia,Condition == 'TC1',select = c(merge_anchor,y.position)),
             by = 'merge_anchor',all.x = T)
ggplot(data = log_corr_abn.serratia_delftia, 
       mapping = aes(x = Day,y = mean,fill = Condition))+
  geom_bar(stat = 'identity',
           position=position_dodge())+
  facet_grid(.~Genus) + 
  geom_errorbar(mapping = aes(ymin = mean - std,
                              ymax = mean + std),
                width = 0.3,
                size = 0.4,
                position = position_dodge(width = 0.9)) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5)+
  facet_grid(.~Genus, scales = "free_x")+
  labs(x='Day', y="log(Bias corrected absolute abundance)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'top',
        text = element_text(size = 16, colour = 'black'),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"),
        strip.text.y = element_text(size = 15, color = "black", face = "bold")) + scale_fill_uchicago(alpha=0.8) +
  stat_pvalue_manual(test)

ggsave(paste(PFIG,'ancombc_delfia-serratia_bcAbundance_pval_annot.png',sep = '/'),width = 8,height = 8, device = 'png')


###############Session information############
sessionInfo()
devtools::session_info()
