#######
#######


######## clean workspace #######

rm(list=ls())

######## load packages ########

library(plyr)
library(dplyr)
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(rstatix)
#library(devtools)
#library(CoDaSeq)
#library(zCompositions)
library(tibble)
library(ggplot2)
library(vegan)
library(agricolae)
library(picante)
library(gridExtra)
library(lme4)
library(lmerTest)
library(microbiome)
library(DT)
library(eulerr)
library(microbiomeutilities)
library(ggrepel)
#library(ranacapa)
library(sjPlot)
library(sjmisc)
library(MuMIn)
library(jtools)
library(car)
library(scales)
library(DHARMa)

######## Finally, set a random seed for reproducibility #######

set.seed(777)

######## read RDS file into R ###### 

ps = readRDS("~/Dropbox (Personal)/Ghana/03_Microbiome/R/working_data/filtered_219samples.rds")

##### check the phyloseq object #####
#sample_names(ps)
#rank_names(ps)
#sample_data(ps)
#tax_table(ps)
#(str(data.frame(sample_data(ps))))

sample_data(ps)$Site<-factor(sample_data(ps)$Site)

sample_data(ps)$sample_period<-as.factor(sample_data(ps)$sample_period)
levels(sample_data(ps)$sample_period)
levels(sample_data(ps)$sample_period) = c("Jan-Feb2012",  "Jan-Feb2011",   "Jan-Feb2012",  "Jul-Aug2011",  "Jul-Aug2012",  "Jul-Aug2011",  
                                          "Jul-Aug2012",   "Mar-Apr2011",  "Mar-Apr2012",  "Mar-Apr2011",   "Mar-Apr2012", "May-Jun2011", 
                                          "May-Jun2012",  "May-Jun2011",   "May-Jun2012",   "Nov-Dec2010",  "Nov-Dec2011",  "Nov-Dec2010",  
                                          "Nov-Dec2011",   "Sep-Oct2010",  "Sep-Oct2011",  "Sep-Oct2010",   "Sep-Oct2011")


##### subsetting for indiviuals lacking sex or age information ###### 
ps<-subset_samples(ps, AGE !="NA" & #### 1 samples = 218 samples left
                     SEX !="NA")  #### 0 samples = 218 samples left

##### double check we only are dealing with one species and rename ps  #####
ps.LinD<-subset_samples(ps, Lineage =="D") #219 samples

###### summarise
microbiome::summarize_phyloseq(ps.LinD)

mean(sample_sums(ps.LinD))

ps.LinD <- subset_taxa(ps.LinD, taxa_sums(ps.LinD)>0) 

######

meta.LinD <- as(sample_data(ps.LinD), "data.frame")

#clean up sample names
sample_names(ps) <- meta.LinD$FaecalID

ddply(meta.LinD, c("infection_status"), summarise, n=length(infection_status)) #sample sizes per group
#1       uninfected 46
#2    2Blog_monoinf 41
#3  2Bbasal_monoinf 70
#4     229E_monoinf 61


ddply(meta.LinD, c("sample_period"), summarise, n=length(sample_period)) #sample sizes per group



##### lets make a pretty compositional bar plot #####
######## now for the microbiome ########

ps.comp<-ps.LinD

##### find most common phyla
ps_phyla<-tax_glom(ps.comp, taxrank = "Phylum")

sample_data(ps.comp)<-sample_data(ps.comp)[,c("sample.ID")]
ps_phyla <- microbiome::transform(ps_phyla, "compositional") #### ps was not previously transformed - i.e. this is the first time - just to keep in mind 
phyla_per_sample<-psmelt(ps_phyla)

get_group_abundances(ps_phyla, level="Phylum", group="Phylum", transform = "compositional") %>%arrange (-mean_abundance)

#1 Firmicutes       Firmicutes              0.670        0.221  
#2 Proteobacteria   Proteobacteria          0.272        0.206  
#3 Synergistota     Synergistota            0.0215       0.0631 
#4 Actinobacteriota Actinobacteriota        0.0175       0.0274 
#5 Bacteroidota     Bacteroidota            0.0118       0.0251 
#6 Campilobacterota Campilobacterota        0.00444      0.0357 
#7 Fusobacteriota   Fusobacteriota          0.00278      0.00627

#### now common classes
ps_class<-tax_glom(ps.comp, taxrank = "Class")

sample_data(ps.comp)<-sample_data(ps.comp)[,c("sample.ID")]
ps_class <- microbiome::transform(ps_class, "compositional") #### ps was not previously transformed - i.e. this is the first time - just to keep in mind 
class_per_sample<-psmelt(ps_class)

get_group_abundances(ps_class, level="Class", group="Class", transform = "compositional") %>%arrange (-mean_abundance)

#1 Bacilli             Bacilli                 0.540         0.243    
#2 Gammaproteobacteria Gammaproteobacteria     0.249         0.199    
#3 Clostridia          Clostridia              0.128         0.158    
#4 Alphaproteobacteria Alphaproteobacteria     0.0221        0.0438   
#5 Synergistia         Synergistia             0.0215        0.0631   
#6 Actinobacteria      Actinobacteria          0.0132        0.0240   
#7 Bacteroidia         Bacteroidia             0.0118        0.0251
#...

##### Visualisation ####

#1. hierarchical clustering

unifrac_dist <- phyloseq::distance(ps.LinD, method = "Unifrac")

#Save as dendrogram
ward <- as.dendrogram(hclust(unifrac_dist, method = "ward.D2"))
#Ward D2 considers the distance between the centroids of the clusters being merged, 
#while Ward D considers the distance between the individual data points 
#and the mean of the merged cluster

library(ggdendro)

dendro.plot<-ggdendrogram(ward, rotate = T)+
  theme_classic2()+ylab("Unweighted Unifrac Distance")+
  theme(axis.line.x = element_blank())

# 2. Compositional plot 
library(microViz)

##### reorder the samples according to the dendrogram order
dendroorder<-gsub("\\.","-",labels(ward))
ps.rearranged<-ps_reorder(ps.LinD, rev(dendroorder))
row.names(sample_data(ps.rearranged)) #double-check
sample_data(ps.rearranged)[218,]
sample_data(ps.rearranged)$infection_status

comp.plot<-ps.rearranged %>%
  comp_barplot(
    tax_level = "Class", n_taxa = 7,
    bar_outline_colour = NA,
    sample_order = "asis",
    bar_width = 0.7,
    taxon_renamer = toupper
  ) + #coord_flip() + 
  theme(#axis.ticks.y = element_blank(), 
        #axis.text.y = element_blank(), 
    legend.position = "top")


#####  alpha diversity metrices and some stats ######
#####################################################

######## first we rarefy ########
minimum = min(sample_sums(ps.LinD)) ###rarefying threshold

standf = function(x, t=minimum) round(t * (x / sum(x))) #standardise by this sampling depth 

ps.LinD.rare = transform_sample_counts(ps.LinD, standf) #normalise data accordingly 

otu.table <- as.data.frame(otu_table(ps.LinD.rare))
df.pd <- pd(t(otu.table), phy_tree(ps.LinD.rare), include.root=T)
meta.LinD$Phyogenetic_diversity <- df.pd$PD
alpha = estimate_richness(ps.LinD.rare, measures = c("Observed", "Chao1")) 
meta.LinD<-cbind(meta.LinD, alpha)

### Visualisations ###

alpha.InfectionStatus<-ggplot(meta.LinD, aes(infection_status, Phyogenetic_diversity, fill=infection_status))+geom_boxplot()+
  scale_fill_manual(values = c("grey85", "#33638DFF", "#30BB75FF", "#FDE725FF"))+
  geom_point(shape=21, size=2.5)+theme_bw()+ylab("Phylogenetic Diversity")+
  theme(axis.title.x = element_text(colour = 'black', face="bold", size=12), 
        axis.title.y = element_text(colour = 'black', face="bold", size=12), 
        legend.position = "none")+ scale_x_discrete(name ="Infection Status", 
                                                    labels=c("uninfected"="uninfected",
                                                             "2Blog_monoinf"="2B",
                                                             "2Bbasal_monoinf"="2Bbasal", 
                                                             "229E_monoinf"="229E"))+
  geom_signif(data=meta.LinD,stat="signif",position="identity",
              comparisons=list(c("uninfected","2Blog_monoinf")),map_signif_level = TRUE, annotations="*", y_position = 47)+
  geom_signif(data=meta.LinD,stat="signif",position="identity",
              comparisons=list(c("2Blog_monoinf","229E_monoinf")),map_signif_level = TRUE, annotations="***", y_position = 44)+
  geom_signif(data=meta.LinD,stat="signif",position="identity",
              comparisons=list(c("2Blog_monoinf","2Bbasal_monoinf")),map_signif_level = TRUE, annotations="*", y_position = 50)


Chao1.InfectionStatus<-ggplot(meta.LinD, aes(infection_status, Chao1, fill=infection_status))+geom_boxplot()+
  scale_fill_manual(values = c("grey85", "#33638DFF", "#30BB75FF", "#FDE725FF"))+
  geom_point(shape=21, size=2.5)+theme_bw()+ylab("Chao1")+
  theme(axis.title.x = element_text(colour = 'black', face="bold", size=12), 
        axis.title.y = element_text(colour = 'black', face="bold", size=12), 
        legend.position = "none")+ scale_x_discrete(name ="Infection Status", 
                                                    labels=c("uninfected"="uninfected",
                                                             "2Blog_monoinf"="2B",
                                                             "2Bbasal_monoinf"="2Bbasal", 
                                                             "229E_monoinf"="229E")) +
  geom_signif(data=meta.LinD,stat="signif",position="identity",
              comparisons=list(c("uninfected","2Blog_monoinf")),map_signif_level = TRUE, annotations="**", y_position = 500)+
  geom_signif(data=meta.LinD,stat="signif",position="identity",
              comparisons=list(c("2Blog_monoinf","229E_monoinf")),map_signif_level = TRUE, annotations="***", y_position = 450)+
  geom_signif(data=meta.LinD,stat="signif",position="identity",
              comparisons=list(c("2Blog_monoinf","2Bbasal_monoinf")),map_signif_level = TRUE, annotations="*", y_position = 400) +
  geom_signif(data=meta.LinD,stat="signif",position="identity",
              comparisons=list(c("229E_monoinf","2Bbasal_monoinf")),map_signif_level = TRUE, annotations="*", y_position = 350)


Observed.InfectionStatus<-ggplot(meta.LinD, aes(infection_status, Observed, fill=infection_status))+geom_boxplot()+
  scale_fill_manual(values = c("grey85", "#33638DFF", "#30BB75FF", "#FDE725FF"))+
  geom_point(shape=21, size=2.5)+theme_bw()+ylab("Observed ASVs")+
  theme(axis.title.x = element_text(colour = 'black', face="bold", size=12), 
        axis.title.y = element_text(colour = 'black', face="bold", size=12), 
        legend.position = "none")+ scale_x_discrete(name ="Infection Status", 
                                                    labels=c("uninfected"="uninfected",
                                                             "2Blog_monoinf"="2B",
                                                             "2Bbasal_monoinf"="2Bbasal", 
                                                             "229E_monoinf"="229E")) +
  geom_signif(data=meta.LinD,stat="signif",position="identity",
                 comparisons=list(c("uninfected","2Blog_monoinf")),map_signif_level = TRUE, annotations="*", y_position = 475)+
  geom_signif(data=meta.LinD,stat="signif",position="identity",
              comparisons=list(c("2Blog_monoinf","229E_monoinf")),map_signif_level = TRUE, annotations="***", y_position = 425)+
  geom_signif(data=meta.LinD,stat="signif",position="identity",
              comparisons=list(c("2Blog_monoinf","2Bbasal_monoinf")),map_signif_level = TRUE, annotations="*", y_position = 375) 
  
ggarrange(Observed.InfectionStatus, Chao1.InfectionStatus, nrow=1, align="hv", 
          labels=c("A", "B"))


###### lmer models including sample period as random effect #######

complex_mod<-lmer(sqrt(Phyogenetic_diversity) # Chao1 and Observed log , PD sqrt,
                  ~ infection_status+#AGE+SEX+ Site+
                    (1|sample_period), 
                 data=meta.LinD, 
                 REML=F, #If your random effects are nested, or you have only one random effect, and if your data are balanced (i.e., similar sample sizes in each factor group) set REML to FALSE, because you can use maximum likelihood
                 na.action=na.fail) #more complex model with random effect sex or site or year ?! changes nothing 
plot(simulateResiduals(fittedModel = complex_mod, plot = F))
dredge(complex_mod) #does not retain location, age or sex
anova(complex_mod)
difflsmeans(complex_mod, test.effs = "infection_status", ddf="Satterthwaite")

##### for future visualisations #####
##### define annotations for all future annotations ####

annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  hjustvar = c(1.1),
  vjustvar = c(2))

###### CT values and alpha diversity #######
###### here we investigate whether a correlation exists between infection intensity
###### and gut microbial alpha diversity (results only for PD given)

cor.test(meta.LinD$CT_value229E, meta.LinD$Phyogenetic_diversity)

Alpha.229E<-ggplot(meta.LinD, aes(CT_value229E, Phyogenetic_diversity))+geom_point(fill="#FDE725FF",shape=21, size=2.5)+
  geom_smooth(method="lm", colour="#FDE725FF",fill="#FDE725FF")+theme_bw()+xlab("CT-value of 229E infections")+
  ylab("Faith's Phylogenetic Diversity")+
  theme(axis.title.x = element_text(colour = 'black', face="bold", size=12), 
        axis.title.y = element_blank())+
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label="R2=0.01, p=0.969"))

cor.test(meta.LinD$CT_value2B, meta.LinD$Phyogenetic_diversity)

Alpha.2B<-ggplot(meta.LinD, aes(CT_value2B, Phyogenetic_diversity))+geom_point(fill="#33638DFF",shape=21, size=2.5)+
  geom_smooth(method="lm", colour="#33638DFF",fill="#33638DFF")+theme_bw()+xlab("CT-value of 2B infections")+
  ylab("Faith's Phylogenetic Diversity")+
  theme(axis.title.x = element_text(colour = 'black', face="bold", size=12), 
        axis.title.y = element_blank())+
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label="R2=0.50, p=0.001"))

cor.test(meta.LinD$CT_value2Bbasal, meta.LinD$Phyogenetic_diversity)

Alpha.2Bbasal<-ggplot(meta.LinD, aes(CT_value2Bbasal, Phyogenetic_diversity))+geom_point(fill="#30BB75FF",shape=21, size=2.5)+
  geom_smooth(method="lm", colour="#30BB75FF",fill="#30BB75FF")+theme_bw()+xlab("CT-value of 2Bbasal infections")+
  ylab("Faith's Phylogenetic Diversity")+
  theme(axis.title.x = element_text(colour = 'black', face="bold", size=12), 
        axis.title.y = element_blank())+
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label="R2=-0.04, p=0.756"))


ggarrange(alpha.InfectionStatus, Alpha.2B, Alpha.2Bbasal, Alpha.229E, nrow=1, align="hv", 
          labels=c("C", "D", "E", "F"))


########## BETA DIVERSITY ########
##################################
###### we tax_glom to genus and remove unknown genera ######
ps.LinD.genera<-tax_glom(ps.LinD.rare, taxrank = "Genus") ### not rarefied 
ps.LinD.genera<-subset_taxa(ps.LinD.genera, !is.na(Genus) & !Genus %in% c("uncultured"))

wunifrac_dist <- phyloseq::distance(ps.LinD.genera, method = "wUnifrac")
unifrac_dist<-phyloseq::distance(ps.LinD.genera, method = "Unifrac")

adonis2(unifrac_dist ~ infection_status+Site+AGE+SEX,  
        data = meta.LinD, strata=meta.LinD$sample_period, method="unifrac", by="margin")
adonis2(wunifrac_dist ~ infection_status+Site+AGE+SEX, 
        data = meta.LinD, strata=meta.LinD$sample_period, method="wunifrac", by="margin")

ordination_Uni<-ordinate(ps.LinD.genera, method="PCoA", distance=unifrac_dist)

unweighted <- plot_ordination(ps.LinD.genera, ordination_Uni) + 
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = infection_status)) + # put ellipses around group centroids
  geom_point(aes(
    fill = infection_status), # fill color by infection status
    size = 4, # make points size 4
    pch = 21, # Make points circular with a border 
    col = "black") + # the border colour should be black
  #labs(subtitle = "a)") + # insert subtitle
  #ggtitle("Unweighted Unifrac") + # insert title
  theme_bw() + scale_fill_manual(values=c("grey85", "#33638DFF", "#30BB75FF", "#FDE725FF")) +
  labs(fill = "infection status") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', face="bold")) +
  theme(legend.position = "none", axis.title = element_text(size=12, face="bold")) 
unweighted


ordination_wUni<-ordinate(ps.LinD.genera, method="PCoA", distance=wunifrac_dist)

weighted <- plot_ordination(ps.LinD.genera, ordination_wUni) + 
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = infection_status)) + # put ellipses around group centroids
  geom_point(aes(
    fill = infection_status), # fill color by infection status
    size = 4, # make points size 4
    pch = 21, # Make points circular with a border 
    col = "black") + # the border colour should be black
  #labs(subtitle = "a)") + # insert subtitle
  #ggtitle("Unweighted Unifrac") + # insert title
  theme_bw() + scale_fill_manual(values=c("grey85", "#33638DFF", "#30BB75FF", "#FDE725FF")) +
  labs(fill = "infection status") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', face="bold"))+ # center align the plot title
  theme(legend.position = "none", axis.title = element_text(size=12, face="bold")) 
weighted


disp.wunifrac = betadisper(wunifrac_dist, meta.LinD$infection_status)
permutest(disp.wunifrac, pairwise=TRUE, permutations=1000)
disp.wunifrac.df <- data.frame(group = disp.wunifrac$group, distances = disp.wunifrac$distances)
plot.disp.core.wUni<-ggplot(disp.wunifrac.df, aes(x=group, y=distances, fill=group))+geom_boxplot()+
  geom_point(color="black", size=2, pch=21)+
  ylab("Distance from Centroid")+xlab("Infection Status")+
  scale_fill_manual(values=c("grey85", "#33638DFF", "#30BB75FF", "#FDE725FF"))+theme_bw()+
  scale_x_discrete(labels=c("uninfected" = "uninfected", "2Blog_monoinf" = "CoV-2B", "2Bbasal_monoinf" = "CoV2B-basal", "229E_monoinf" = "CoV-229E"))+
  theme(legend.position = "none", axis.title = element_text(size=12, face="bold"))

disp.unifrac = betadisper(unifrac_dist, meta.LinD$infection_status)
permutest(disp.unifrac, pairwise=TRUE, permutations=1000)
disp.unifrac.df <- data.frame(group = disp.unifrac$group, distances = disp.unifrac$distances)
plot.disp.core.Uni<-ggplot(disp.unifrac.df, aes(x=group, y=distances, fill=group))+geom_boxplot()+
  geom_point(color="black", size=2, pch=21)+
  ylab("Distance from Centroid")+xlab("Infection Status")+
  scale_fill_manual(values=c("grey85", "#33638DFF", "#30BB75FF", "#FDE725FF"))+theme_bw()+
  scale_x_discrete(labels=c("uninfected" = "uninfected", "2Blog_monoinf" = "CoV-2B", "2Bbasal_monoinf" = "CoV2B-basal", "229E_monoinf" = "CoV-229E"))+
  theme(legend.position = "none", axis.title = element_text(size=12, face="bold"))+
  geom_signif(data=disp.unifrac.df,stat="signif",position="identity",
              comparisons=list(c("uninfected","2Bbasal_monoinf")),map_signif_level = TRUE, annotations="**", y_position = 0.64)+
  geom_signif(data=disp.unifrac.df,stat="signif",position="identity",
              comparisons=list(c("uninfected","229E_monoinf")),map_signif_level = TRUE, annotations="**", y_position = 0.66)+
  geom_signif(data=disp.unifrac.df,stat="signif",position="identity",
              comparisons=list(c("2Blog_monoinf","229E_monoinf")),map_signif_level = TRUE, annotations="*", y_position = 0.68)
  

beta.all.InfectionStatus<-ggarrange(unweighted, plot.disp.core.Uni,
                                           weighted, plot.disp.core.wUni,
                                ncol=2, nrow=2, widths=c(1,1), align="hv", 
                                labels=c("A", "B", "C", "D"))


##### Location effect #####

unweighted.Loc <- plot_ordination(ps.LinD.genera, ordination_Uni) + 
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Site)) + # put ellipses around group centroids
  geom_point(aes(
    fill = Site), # fill color by infection status
    size = 4, # make points size 4
    pch = 21, # Make points circular with a border 
    col = "black") + # the border colour should be black
  #labs(subtitle = "a)") + # insert subtitle
  #ggtitle("Unweighted Unifrac") + # insert title
  theme_bw() + scale_fill_manual(values=c("#008a89", "#8affd8", "#b5004e", "#be9d1f", "#fb9b53"),
                                 labels=c("BUO1"="Buoyem1", "BUO2"="Buoyem2", "FO"="Forikrom", "KW1"="Kwamang1", "KW2"="Kwamang2")) +
  labs(fill = "Site") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', face="bold")) +
  theme(legend.position = "top", axis.title = element_text(size=12, face="bold")) 
unweighted.Loc

weighted.Loc <- plot_ordination(ps.LinD.genera, ordination_wUni) + 
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Site)) + # put ellipses around group centroids
  geom_point(aes(
    fill = Site), # fill color by infection status
    size = 4, # make points size 4
    pch = 21, # Make points circular with a border 
    col = "black") + # the border colour should be black
  #labs(subtitle = "a)") + # insert subtitle
  #ggtitle("Unweighted Unifrac") + # insert title
  theme_bw() + scale_fill_manual(values=c("#008a89", "#8affd8", "#b5004e", "#be9d1f", "#fb9b53")) +
  labs(fill = "Site") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', face="bold"))+ # center align the plot title
  theme(legend.position = "none", axis.title = element_text(size=12, face="bold")) 
weighted.Loc

beta.Location<-ggarrange(unweighted.Loc, weighted.Loc,
                                         ncol=2, nrow=1, widths=c(1,1), align="h", common.legend = T, 
                         labels=c("A", "B"))


##### lets look at CT values in relation to unifrac distances #####
#### first 2B ####
ps.LinD.genera.CoV2B<-subset_samples(ps.LinD.genera, infection_status=="2Blog_monoinf")

wunifrac_dist <- phyloseq::distance(ps.LinD.genera.CoV2B, method = "wUnifrac")
unifrac_dist<-phyloseq::distance(ps.LinD.genera.CoV2B, method = "Unifrac")

meta.LinD.genera.CoV2B <- as(sample_data(ps.LinD.genera.CoV2B), "data.frame")
meta.LinD.genera.CoV2B$unweightedUnifrac<-colMeans(as.matrix(unifrac_dist))
meta.LinD.genera.CoV2B$weightedUnifrac<-colMeans(as.matrix(wunifrac_dist))

cor.test(meta.LinD.genera.CoV2B$unweightedUnifrac, meta.LinD.genera.CoV2B$CT_value2B)
cor.test(meta.LinD.genera.CoV2B$weightedUnifrac, meta.LinD.genera.CoV2B$CT_value2B)

Beta.2B<-ggplot(meta.LinD.genera.CoV2B, aes(CT_value2B, unweightedUnifrac))+geom_point(fill="#33638DFF",shape=21, size=2.5)+
  geom_smooth(method="lm", colour="#33638DFF",fill="#33638DFF")+theme_bw()+xlab("CT-value of 2B infections")+
  ylab("Unweighted Unifrac")+
  theme(axis.title = element_text(colour = 'black', face="bold", size=12))+
  scale_y_continuous(breaks = seq(0.50, 0.80, by = .10))+
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label="R2=-0.43, p=0.005"))

##### for CoV-229E CT values #####
ps.LinD.genera.CoV229E<-subset_samples(ps.LinD.genera, infection_status=="229E_monoinf")

wunifrac_dist <- phyloseq::distance(ps.LinD.genera.CoV229E, method = "wUnifrac")
unifrac_dist<-phyloseq::distance(ps.LinD.genera.CoV229E, method = "Unifrac")

meta.LinD.genera.CoV229E <- as(sample_data(ps.LinD.genera.CoV229E), "data.frame")
meta.LinD.genera.CoV229E$unweightedUnifrac<-colMeans(as.matrix(unifrac_dist))
meta.LinD.genera.CoV229E$weightedUnifrac<-colMeans(as.matrix(wunifrac_dist))

cor.test(meta.LinD.genera.CoV229E$unweightedUnifrac, meta.LinD.genera.CoV229E$CT_value229E)
cor.test(meta.LinD.genera.CoV229E$weightedUnifrac, meta.LinD.genera.CoV229E$CT_value229E)

Beta.229E<-ggplot(meta.LinD.genera.CoV229E, aes(CT_value229E, unweightedUnifrac))+geom_point(fill="#FDE725FF",shape=21, size=2.5)+
  geom_smooth(method="lm", colour="#FDE725FF",fill="#FDE725FF")+theme_bw()+xlab("CT-value of 229E infections")+
  ylab("Averaged Unweighted Unifrac")+
  theme(axis.title = element_text(colour = 'black', face="bold", size=12), 
        axis.title.y = element_blank())+
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label="R2=-0.16, p=0.208"))

###### 2B basal  ######

ps.LinD.genera.CoV2Bbasal<-subset_samples(ps.LinD.genera, infection_status=="2Bbasal_monoinf")

wunifrac_dist <- phyloseq::distance(ps.LinD.genera.CoV2Bbasal, method = "wUnifrac")
unifrac_dist<-phyloseq::distance(ps.LinD.genera.CoV2Bbasal, method = "Unifrac")

meta.LinD.genera.CoV2Bbasal <- as(sample_data(ps.LinD.genera.CoV2Bbasal), "data.frame")
meta.LinD.genera.CoV2Bbasal$unweightedUnifrac<-colMeans(as.matrix(unifrac_dist))
meta.LinD.genera.CoV2Bbasal$weightedUnifrac<-colMeans(as.matrix(wunifrac_dist))

cor.test(meta.LinD.genera.CoV2Bbasal$unweightedUnifrac, meta.LinD.genera.CoV2Bbasal$CT_value2Bbasal)
cor.test(meta.LinD.genera.CoV2Bbasal$weightedUnifrac, meta.LinD.genera.CoV2Bbasal$CT_value2Bbasal)

Beta.2Bbasal<-ggplot(meta.LinD.genera.CoV2Bbasal, aes(CT_value2Bbasal, unweightedUnifrac))+geom_point(fill="#30BB75FF",shape=21, size=2.5)+
  geom_smooth(method="lm", colour="#30BB75FF",fill="#30BB75FF")+theme_bw()+xlab("CT-value of 2Bbasal infections")+
  ylab("Averaged Unweighted Unifrac")+
  theme(axis.title.x = element_text(colour = 'black', face="bold", size=12), 
        axis.title.y = element_blank())+
  scale_y_continuous(breaks = seq(0.50, 0.80, by = .10))+
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label="R2=0.18, p=0.133"))


ggarrange(unweighted, Beta.2B, Beta.2Bbasal, Beta.229E, 
          nrow=1, align="hv", labels=c("G", "H", "I", "J"))


###### JSDM ######
##### now we subset for only core genera #####
ps.LinD.genera<-tax_glom(ps.LinD, taxrank = "Genus") ### not rarefied 
ps.LinD.genera<-subset_taxa(ps.LinD.genera, !is.na(Genus) & !Genus %in% c("uncultured"))
ps.LinD.core<-microbiome::core(ps.LinD.genera, detection = 0, prevalence = 0.50)
genus<-data.frame(tax_table(ps.LinD.core)[,"Genus"])


my_subset <- subset(otu_table(ps.LinD.genera), 
                    rownames(otu_table(ps.LinD.genera)) 
                    %in% c("ASV-368", #Weissella
                           "ASV-386", #Lactobacillus
                           "ASV-527", #Lactococcus
                           "ASV-559", #Streptococcus
                           "ASV-784", #Ureaplasma
                           "ASV-786", #Mycoplasma
                           "ASV-1185", #Pectobacterium
                           "ASV-1327", #Morganella
                           "ASV-1348", #Acinetobacter
                           "ASV-1794", #Asaia
                           "ASV-2355", #Bartonella
                           "ASV-2758", #Fusobacterium
                           "ASV-4023", #Alistipes
                           "ASV-4369", #Empedobacter
                           "ASV-4781", #Dysgonomonas
                           "ASV-5106", #Bacteroides
                           "ASV-5637", #Mycobacterium
                           "ASV-5886", #Cutibacterium
                           "ASV-6855", #Raoultibacter
                           "ASV-7449", #Clostridium_sensu_stricto_1
                           "ASV-7571", #Staphylococcus
                           "ASV-7690", #Enterococcus
                           "ASV-7789", #Gemella
                           "ASV-7803", #Paeniclostridium
                           "ASV-8129", #Christensenellaceae_R-7_group
                           "ASV-8948", #Lachnoclostridium
                           "ASV-9968", #Incertae_Sedis
                           "ASV-10318" #Candidatus_Soleaferrea
                    ))


ps.gllvm <- merge_phyloseq(my_subset, tax_table(ps.LinD.genera), sample_data(ps.LinD.genera)
                           #, phy_tree(ps.LinD.genera)
                           )

##### linear latent variable models as joint speecies distribution model ####
library(gllvm)

##### transform data to account for its compositionality
ps.gllvm.clr <-microbiome::transform((otu_table(ps.gllvm)), transform = "clr")

genus<-data.frame(tax_table(ps.gllvm)[,"Genus"])
taxa_names(ps.gllvm)<-genus$Genus

ASV<-data.frame(t(otu_table(ps.gllvm)))

meta.LinD.core.clr <- data.frame(sample_data(ps.gllvm))
meta.LinD.core.clr$SeqDepth<-as.numeric(rescale(sample_sums(ps.gllvm), to =c(-1,1)))
meta.LinD.core.clr<-meta.LinD.core.clr[,c(217,218,34,33,27, 28, 38)]

gllvm1<-gllvm(ASV, 
              meta.LinD.core.clr, #otu table and meta data as data frame 
          formula = ~ 
            infection_status + Site + AGE + SEX + SeqDepth, 
          family = "negative.binomial", 
          row.eff = ~(1|sample_period), #random effect
          num.lv = 3)

par(mfrow = c(1, 2))
plot(gllvm1, which = 1:2, var.colors = 1, n.plot = 20) #looks good
summary(gllvm1)

##### building a graph ##### 

df <- coef(gllvm1)
est_df <- data.frame(df$Intercept)
est_df2 <- data.frame(df$Xcoef)
est_df3 <- merge(est_df, est_df2, by=0)

#reorder genera

row.names(est_df3)<-est_df3$Row.names
est_df3<-est_df3[colnames(ASV),]

#turn df into long format 

names(est_df3)[1]<-"Genus"
names(est_df3)[2]<-"Intercept"
estimates_df <- gather(est_df3, Treatment, Estimate, names(est_df3)[2]:names(est_df3)[ncol(est_df3)], factor_key=TRUE)

# extract CIs

confint_df <- data.frame(confint(gllvm1))
confint_df <- rbind(confint_df[grep("^Xcoef", rownames(confint_df)), ], confint_df[grep("^Intercept", rownames(confint_df)), ]) 

# adding the variable levels as a column

variables <- colnames(est_df3)[3:ncol(est_df3)]
variables <- c(variables, "Intercept")
variables1 <- rep(variables, nrow(est_df3))
variables2 <- variables1[order(match(variables1, variables))]

confint_df$Treatment<-variables2

# column with taxa names

confint_df$Genus <- rep(colnames(ASV), length(unique(confint_df$Treatment)))

# now have estimates and confidence intervals as seperate data frames, but they are in different formats. Need to massage them into one dataframe for plotting.

merged <- merge(estimates_df, confint_df, by =c("Treatment", "Genus"))

names(merged)[4]<-"CI_lower"
names(merged)[5]<-"CI_upper"

unique(merged$Treatment)


merged$Sig <- !data.table::between(0, merged$CI_lower, merged$CI_upper) 
merged$Genus<-as.factor(merged$Genus)

write.csv(merged, file="~/Dropbox (Personal)/Ghana/03_Microbiome/R/gllvm_results_08.4.csv")

subset_merged<-subset(merged, Treatment != "Intercept" & 
                        Treatment != "SeqDepth")

variable_names <- c('infection_status2Blog_monoinf' = "CoV-2B", 
                    'infection_status229E_monoinf' = "CoV229E", 
                    'infection_status2Bbasal_monoinf' = "CoV2Bbasal",
                    'SiteBUO2' = "BUO2", 
                    'SiteFO' = "FO",
                    "SiteKW1"= "KW1",
                    "SiteKW2" = "KW2",
                    'AGESubadult' = "subadult", 
                    'SEXM'='male')


ggplot(subset_merged, aes(y = Genus, x=Estimate, fill= Sig))+
  facet_wrap(~Treatment, scales="free_x", nrow=1, labeller = as_labeller(variable_names))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey") + 
  geom_errorbar(aes(xmin=CI_lower, xmax = CI_upper),width=0, size=1, colour="grey")+
  geom_point(size=2.5, pch=21)+ 
  theme_classic(base_size = 12) +
  scale_fill_manual(values=c("grey", "red"))+
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(face="bold"), axis.title.y = element_blank())+
  theme(axis.text.y = element_text(colour=subset_merged$Colour)) +
  scale_y_discrete(label = function(x) abbreviate(x, minlength=7, FALSE))


subset_merged2<-subset(merged, Treatment == "infection_status2Blog_monoinf" )

gllvm.plot<-ggplot(subset_merged2, aes(y = Genus, x=Estimate, fill= Sig))+
  facet_wrap(~Treatment, scales="free_x", nrow=1, labeller = as_labeller(variable_names))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey") + 
  geom_errorbar(aes(xmin=CI_lower, xmax = CI_upper),width=0, size=1, colour="grey")+
  geom_point(size=2.5, pch=21)+ 
  theme_classic(base_size = 12) +
  scale_fill_manual(values=c("grey", "red"))+
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(face="bold"), axis.title.y = element_blank())+
  theme(axis.text.y = element_text(colour=subset_merged$Colour)) +
  scale_y_discrete(label = function(x) abbreviate(x, minlength=7, FALSE))



######## prep for GAMMs ########
top_genera<-microbiomeutilities::aggregate_top_taxa2(ps, "Genus", top = 40)
print(otu_table(top_genera))[39:40,]
ps.LinD.core.clr <-microbiome::transform(otu_table(ps.LinD.core), transform = "clr")
top_genera_melt<-psmelt(ps.LinD.core.clr)
head(top_genera_melt)
meta.LinD.core<-data.frame(sample_data(ps.LinD.core))
meta.LinD.core$Seq_depth<-as.numeric(sample_sums(ps.LinD.core))
meta.LinD.core$Seq_depth_scaled <- as.numeric(scales::rescale(meta.LinD.core$Seq_depth, to = c(-1,1)))
top_hippo_genera<- merge(top_genera_melt, meta.LinD.core, by.x = "Sample", by.y = "sample.ID")
head(top_hippo_genera)

Mycoplasma<-subset(top_hippo_genera, OTU == "ASV-786") #Mycoplasma
Staphylococcus<-subset(top_hippo_genera, OTU == "ASV-7571") #Staphylococcus
Alistipes<-subset(top_hippo_genera, OTU == "ASV-4023") #Alistipes
Candidatus_Soleaferrea<-subset(top_hippo_genera, OTU == "ASV-10318") #Candidatus_Soleaferrea
Gemella<-subset(top_hippo_genera, OTU == "ASV-7789") #Gemella
Paeniclostridium<-subset(top_hippo_genera, OTU == "ASV-7803") #Paenicolstridium
Christensenellaceae_R.7_group<-subset(top_hippo_genera, OTU == "ASV-8129") #Christensenellaceae_R.7_group
Raoultibacter<-subset(top_hippo_genera, OTU == "ASV-6855") #Raoultibacter
Ureaplasma<-subset(top_hippo_genera, OTU == "ASV-784") #Ureaplasma
Pectobacterium<-subset(top_hippo_genera, OTU == "ASV-1185") #Pectobacterium


##### Mycoplasma #####
Mycoplasma.2B<-subset(Mycoplasma, CT_value2B!="NA")
Mycoplasma.uninfected<-subset(Mycoplasma, infection_status=="uninfected")
Mycoplasma.uninfected$CT_value2B<-44.0
library(sjmisc)
Mycoplasma.2B<-merge_df(Mycoplasma.2B,Mycoplasma.uninfected)

library(mgcv)

Mycoplasma_gam <- mgcv::gam(Abundance~
                              s(CT_value2B, bs="cr") +
                             Site+SEX+AGE+
                             s(Seq_depth_scaled, k=5), 
                           data=Mycoplasma.2B,
                           family = gaussian)

print(summary(Mycoplasma_gam)) 
gam.check(Mycoplasma_gam)
plot(Mycoplasma_gam)

plot.gam(Mycoplasma_gam)

Mycoplasma_plot<-tidymv::plot_smooths(
  model = Mycoplasma_gam,
  series = CT_value2B, 
  exclude_random = TRUE,
  exclude_terms = c("Site", "SEX", "AGE", "Seq_depth_scaled")
) + theme_bw() + 
  ylab("Mycoplasma") + xlab("CT-value of 2B infection") +
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.title = element_text(colour="black",size=12,face="bold"))+
geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label="F=4.81, p=0.031"))

#####
Staphylococcus.2B<-subset(Staphylococcus, CT_value2B!="NA")
Staphylococcus.uninfected<-subset(Staphylococcus, infection_status=="uninfected")
Staphylococcus.uninfected$CT_value2B<-44.0
Staphylococcus.2B<-merge_df(Staphylococcus.2B,Staphylococcus.uninfected)

##### Staphylococcus #####
Staphylococcus_gam <- mgcv::gam(Abundance~
                              Site+SEX+AGE+
                              s(CT_value2B, bs="cr") +
                              s(Seq_depth_scaled, k=10), 
                            data=Staphylococcus.2B,
                            family = gaussian)

print(summary(Staphylococcus_gam)) 
gam.check(Staphylococcus_gam)
plot(Staphylococcus_gam)


Staphylococcus_plot<-tidymv::plot_smooths(
  model = Staphylococcus_gam,
  series = CT_value2B,
  exclude_random = TRUE,
  exclude_terms = c("Site", "SEX", "AGE", "Seq_depth_scaled")
) + theme_bw() + scale_linetype_manual(values=c("solid")) + 
  ylab("Staphyloc...") + xlab("CT-value of 2B infection") +
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.title = element_text(colour="black",size=12,face="bold"))+
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label="F=5.35, p=0.024"))

##### Alistipes ######
Alistipes.2B<-subset(Alistipes, CT_value2B!="NA")
Alistipes.uninfected<-subset(Alistipes, infection_status=="uninfected")
Alistipes.uninfected$CT_value2B<-44.0
Alistipes.2B<-merge_df(Alistipes.2B,Alistipes.uninfected)

##### Alistipes #####
Alistipes_gam <- mgcv::gam(Abundance~
                               Site+SEX+
                               AGE+
                               s(CT_value2B, bs="cr") +
                             s(Seq_depth_scaled),
                             data=Alistipes.2B,
                             family = gaussian)

print(summary(Alistipes_gam)) 
gam.check(Alistipes_gam)
plot(Alistipes_gam)


Alistipes_plot<-tidymv::plot_smooths(
  model = Alistipes_gam,
  series = CT_value2B,   
  exclude_random = TRUE,
  exclude_terms = c("Site", "SEX", "AGE", "Seq_depth_scaled")
) + theme_bw() +
  ylab("Alistipes") + xlab("CT-value of 2B infection") +
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.title = element_text(colour="black",size=12,face="bold"))+
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label="F=3.32, p=0.007"))



##### Gemella ######
Gemella.2B<-subset(Gemella, CT_value2B!="NA")
Gemella.uninfected<-subset(Gemella, infection_status=="uninfected")
Gemella.uninfected$CT_value2B<-44.0
Gemella.2B<-merge_df(Gemella.2B,Gemella.uninfected)

##### Gemella.2B #####
Gemella_gam <- mgcv::gam(Abundance~
                                 Site+SEX+AGE+
                                 s(CT_value2B, bs="cr") +
                               s(Seq_depth_scaled), #+
                               #sample_period,
                               data=Gemella.2B,
                               family = gaussian)

print(summary(Gemella_gam)) 
gam.check(Gemella_gam)
plot(Gemella_gam)

Gemella_plot<-tidymv::plot_smooths(
  model = Gemella_gam,
  series = CT_value2B,
  exclude_random = TRUE,
  exclude_terms = c("Site", "SEX", "AGE", "Seq_depth_scaled"),
) + theme_bw() +
  ylab("Gemella") + xlab("CT-value of 2B infection") +
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.title = element_text(colour="black",size=12,face="bold"))+
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label="F=1.175, p=0.282"))


##### Candidatus_Soleaferrea ######
Candidatus_Soleaferrea.2B<-subset(Candidatus_Soleaferrea, CT_value2B!="NA")
Candidatus_Soleaferrea.uninfected<-subset(Candidatus_Soleaferrea, infection_status=="uninfected")
Candidatus_Soleaferrea.uninfected$CT_value2B<-44.0
Candidatus_Soleaferrea.2B<-merge_df(Candidatus_Soleaferrea.2B,Candidatus_Soleaferrea.uninfected)

##### Candidatus_Soleaferrea #####
Candidatus_Soleaferrea_gam <- mgcv::gam(Abundance~
                                               Site+SEX+AGE+
                                               s(CT_value2B, bs="cr") +
                                             s(Seq_depth_scaled), #+
                                             #sample_period,
                                             data=Candidatus_Soleaferrea.2B,
                                             family = gaussian)

print(summary(Candidatus_Soleaferrea_gam)) 
gam.check(Candidatus_Soleaferrea_gam)
plot(Candidatus_Soleaferrea_gam)


Candidatus_Soleaferrea_plot<-tidymv::plot_smooths(
  model = Candidatus_Soleaferrea_gam,
  series = CT_value2B,
  exclude_random = TRUE,
  exclude_terms = c("Site", "SEX", "AGE", "Seq_depth_scaled"),
) + theme_bw() +
  ylab("C. Soleaferrea") + xlab("CT-value of 2B infection") +
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.title = element_text(colour="black",size=12,face="bold"))+
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label="F=10.27, p=0.002"))


##### Christensenellaceae_R.7 ######
Christensenellaceae_R.7_group.2B<-subset(Christensenellaceae_R.7_group, CT_value2B!="NA")
Christensenellaceae_R.7_group.uninfected<-subset(Christensenellaceae_R.7_group, infection_status=="uninfected")
Christensenellaceae_R.7_group.uninfected$CT_value2B<-44.0
Christensenellaceae_R.7_group.2B<-merge_df(Christensenellaceae_R.7_group.2B,Christensenellaceae_R.7_group.uninfected)

##### Christensenellaceae_R.7_group #####
Christensenellaceae_R.7_group_gam <- mgcv::gam(Abundance~
                                                              Site+SEX+AGE+
                                                              s(CT_value2B, bs="cr") +
                                                            s(Seq_depth_scaled), #+
                                                            #sample_period,
                                                            data=Christensenellaceae_R.7_group.2B,
                                                            family = gaussian)

print(summary(Christensenellaceae_R.7_group_gam)) 
gam.check(Christensenellaceae_R.7_group_gam)
plot(Christensenellaceae_R.7_group_gam)


Christensenellaceae_R.7_group_plot<-tidymv::plot_smooths(
  model = Christensenellaceae_R.7_group_gam,
  series = CT_value2B,
  exclude_random = TRUE,
  exclude_terms = c("Site", "SEX", "AGE", "Seq_depth_scaled"),
) + theme_bw() +
  ylab("Christensen...") + xlab("CT-value of 2B infection") +
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.title = element_text(colour="black",size=12,face="bold"))+
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label="F=15.49, p<0.001"))


##### Raoultibacter ######
Raoultibacter.2B<-subset(Raoultibacter, CT_value2B!="NA")
Raoultibacter.uninfected<-subset(Raoultibacter, infection_status=="uninfected")
Raoultibacter.uninfected$CT_value2B<-44.0
Raoultibacter.2B<-merge_df(Raoultibacter.2B,Raoultibacter.uninfected)



##### Raoultibacter #####
Raoultibacter_gam <- mgcv::gam(Abundance~
                                    Site+SEX+AGE+
                                    s(CT_value2B, bs="cr") +
                                    s(Seq_depth_scaled), #+
                                  #sample_period,
                                  data=Raoultibacter.2B,
                                  family = gaussian)

print(summary(Raoultibacter_gam)) 
gam.check(Raoultibacter_gam)
plot(Raoultibacter_gam)


Raoultibacter_plot<-tidymv::plot_smooths(
  model = Raoultibacter_gam,
  series = CT_value2B,
  exclude_random = TRUE,
  exclude_terms = c("Site", "SEX", "AGE", "Seq_depth_scaled"),
) + theme_bw() +
  ylab("Raoultibacter") + xlab("CT-value of 2B infection") +
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.title = element_text(colour="black",size=12,face="bold"))+
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label="F=2.412, p=0.071"))



##### Paeniclostridium ######
Paeniclostridium.2B<-subset(Paeniclostridium, CT_value2B!="NA")
Paeniclostridium.uninfected<-subset(Paeniclostridium, infection_status=="uninfected")
Paeniclostridium.uninfected$CT_value2B<-44.0
Paeniclostridium.2B<-merge_df(Paeniclostridium.2B,Paeniclostridium.uninfected)


##### Paeniclostridium #####
Paeniclostridium_gam <- mgcv::gam(Abundance~
                                 Site+SEX+AGE+
                                 s(CT_value2B, bs="cr") +
                                 s(Seq_depth_scaled), #+
                               #sample_period,
                               data=Paeniclostridium.2B,
                               family = gaussian)

print(summary(Paeniclostridium_gam)) 
gam.check(Paeniclostridium_gam)
plot(Paeniclostridium_gam)


Paeniclostridium_plot<-tidymv::plot_smooths(
  model = Paeniclostridium_gam,
  series = CT_value2B,
  exclude_random = TRUE,
  exclude_terms = c("Site", "SEX", "AGE", "Seq_depth_scaled"),
) + theme_bw() +
  ylab("Paeniclost...") + xlab("CT-value of 2B infection") +
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.title = element_text(colour="black",size=12,face="bold"))+
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label="F=0.71, p=0.695"))



##### Ureaplasma ######
Ureaplasma.2B<-subset(Ureaplasma, CT_value2B!="NA")
Ureaplasma.uninfected<-subset(Ureaplasma, infection_status=="uninfected")
Ureaplasma.uninfected$CT_value2B<-44.0
Ureaplasma.2B<-merge_df(Ureaplasma.2B,Ureaplasma.uninfected)


##### Ureaplasma #####
Ureaplasma_gam <- mgcv::gam(Abundance~
                                 Site+SEX+AGE+
                                 s(CT_value2B, bs="cr") +
                                 s(Seq_depth_scaled), #+
                               #sample_period,
                               data=Ureaplasma.2B,
                               family = gaussian)

print(summary(Ureaplasma_gam)) 
gam.check(Ureaplasma_gam)
plot(Ureaplasma_gam)


Ureaplasma_plot<-tidymv::plot_smooths(
  model = Ureaplasma_gam,
  series = CT_value2B,
  exclude_random = TRUE,
  exclude_terms = c("Site", "SEX", "AGE", "Seq_depth_scaled"),
) + theme_bw() +
  ylab("Ureaplasma") + xlab("CT-value of 2B infection") +
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.title = element_text(colour="black",size=12,face="bold"))+
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label="F=0.65, p=0.534"))



##### Pectobacterium ######
Pectobacterium.2B<-subset(Pectobacterium, CT_value2B!="NA")
Pectobacterium.uninfected<-subset(Pectobacterium, infection_status=="uninfected")
Pectobacterium.uninfected$CT_value2B<-44.0
Pectobacterium.2B<-merge_df(Pectobacterium.2B,Pectobacterium.uninfected)


##### Pectobacterium #####
Pectobacterium_gam <- mgcv::gam(Abundance~
                                 Site+SEX+AGE+
                                 s(CT_value2B, bs="cr") +
                                 s(Seq_depth_scaled), #+
                               #sample_period,
                               data=Pectobacterium.2B,
                               family = gaussian)

print(summary(Pectobacterium_gam)) 
gam.check(Pectobacterium_gam)
plot(Pectobacterium_gam)


Pectobacterium_plot<-tidymv::plot_smooths(
  model = Pectobacterium_gam,
  series = CT_value2B,
  exclude_random = TRUE,
  exclude_terms = c("Site", "SEX", "AGE", "Seq_depth_scaled"),
) + theme_bw() +
  ylab("Pectobact...") + xlab("CT-value of 2B infection") +
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.title = element_text(colour="black",size=12,face="bold"))+
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label="F=9.67, p<0.001"))


gam.plots<-ggarrange(Gemella_plot, 
                     Raoultibacter_plot, 
                     Fusobacterium_plot, 
                     Candidatus_Soleaferrea_plot,
                     Ureaplasma_plot,
                     Pectobacterium_plot,
                     Paeniclostridium_plot,
                     nrow=2, ncol=4, align="hv", labels= c("A","B", "C", "D", "E", "F", "G"))

gam.plots<-ggarrange(Mycoplasma_plot, 
                     Staphylococcus_plot, 
                     Pectobacterium_plot,
                     Ureaplasma_plot,
                     Alistipes_plot, 
                     Christensenellaceae_R.7_group_plot, 
                     #Gemella_plot, 
                     #Raoultibacter_plot, 
                     Candidatus_Soleaferrea_plot,
                     Paeniclostridium_plot,
                     nrow=2, ncol=4, align="hv", labels= c("B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"))
ggarrange(gllvm.plot, gam.plots, nrow=1, ncol=2, widths=c(1,3), labels=c("A", ""))



