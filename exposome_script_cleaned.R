##########ANALYSIS OF THE EFFECTS OF PESTICIDE EXPOSURE ON THE GUT MICROBIOME OF DAG3 PARTICIPANTS
####### Last updated: 26 Jan 2022
####### Created by: Milla Gois 
###expected output: (3) plots of age/sex/employment distribution, (1) pdf with plot of alpha diversity 

#Load packages and input files: 
set.seed(123)
library(dplyr)
library(ggplot2)

phenos <- read.table("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-rgacesa/DAG3_data_ready/phenotypes/DAG3_metadata_merged_ready_v27.csv", sep=",", fill =T, header =T)
#work questionnaires - LL
work2 <- read.table("/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v1/tab_separated_labels/Questionnaire_Work.dat", sep="\t", fill =T, header =T)
#Pesticide data upload:
final_pesticide_dag3 <- read.table("/groups/umcg-lifelines/tmp01/users/umcg-mgois/exposome_project/pesticide_exposure_DAG3.txt", fill = T, header = T, sep = "\t")
#ids - DAG3
ids <- read.table("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-rgacesa/DAG3_data_ready/dag3_lld_linkage_files/pairing_files/PSEUDOIDEXT.dat", sep="\t", fill =T, header =T)
#bacterial dag3 metaphlan3
DAG3_metaphlan <- read.table("/groups/umcg-lifelines/tmp01/users/umcg-mgois/exposome_project/metaphlan3_pseudoids.txt", sep="\t", header =T)
#paticipants linkage file
participants <- read.table("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-rgacesa/DAG3_data_ready/microbiome/QC/key_DAG3_pID.tsv", fill = T, header =T, sep ="\t")

#prepare IDs df:
DAG3.ids <- ids

##CLR transformation function

library(data.table)
do_clr_externalWeighting = function(interest_matrix, core_matrix){
        if(any(interest_matrix==0)) interest_matrix = interest_matrix + min(interest_matrix[interest_matrix>0])/2
        if(any(core_matrix==0)) core_matrix = core_matrix + min(core_matrix[core_matrix>0])/2
        
        #estimate weighting parameter
        gm_mean = function(x, na.rm=TRUE){
                exp(sum(log(x), na.rm=na.rm) / length(x))
        }
        Gmean_core = apply(core_matrix, 1, gm_mean)
        
        #do transformation
        data_prepared = cbind(Gmean_core,interest_matrix)
        data_transformed = t(apply(data_prepared,1,function(x){
                log(x / x[1])[-1]
        }))
        colnames(data_transformed) = colnames(data_transformed)
        rownames(data_transformed) = rownames(data_transformed)
        data_transformed
}

##removed participant PSEUDOIDEXT == 30421653 because had two entries in "participants" 
DAG3_metaphlan <- DAG3_metaphlan[!DAG3_metaphlan$PSEUDOIDEXT == 30421653, ]

#Prepare data for CLR transformation and run the function:
DAG3_metaphlan_final <- DAG3_metaphlan
DAG3_metaphlan <- DAG3_metaphlan[,-1]
rownames(DAG3_metaphlan) <- DAG3_metaphlan$PSEUDOIDEXT
DAG3_metaphlan$PSEUDOIDEXT <- NULL

taxa_transformed = do_clr_externalWeighting(DAG3_metaphlan,DAG3_metaphlan[,grep("[.]s__",colnames(DAG3_metaphlan))])
taxa_transformed = taxa_transformed[,colSums(DAG3_metaphlan>0)>nrow(DAG3_metaphlan) * 0.05]

dag3_taxa_transf <- taxa_transformed

#subset the microbiome data for which we also have pesticide data: 
common_DAG3 <- intersect(final_pesticide_dag3$PSEUDOIDEXT, DAG3_metaphlan_final$PSEUDOIDEXT)
dag3_metaphlan_exposure <- final_pesticide_dag3[final_pesticide_dag3$PSEUDOIDEXT %in% common_DAG3, ]


#merging the phenos into one df: 
phenos_new <- merge(participants, phenos, by="DAG3_sampleID")
phenos_new <- phenos_new[,-1]
phenos_new <- phenos_new[phenos_new$PSEUDOIDEXT %in% dag3_metaphlan_exposure$PSEUDOIDEXT,]

dag3_metaphlan_exposure <- dag3_metaphlan_exposure[dag3_metaphlan_exposure$PSEUDOIDEXT %in% phenos_new$PSEUDOIDEXT, ]

#age/sex plot 

age_sex_dag3 <- phenos_new[,c("PSEUDOIDEXT", "ANTHRO.AGE", "ANTHRO.Sex")]
age_sex_dag3$ANTHRO.Sex <- as.factor(age_sex_dag3$ANTHRO.Sex)


png("age_sex_distribution_DAG3.png")

ggplot(age_sex_dag3, aes(x= ANTHRO.AGE, fill = ANTHRO.Sex)) + 
        geom_bar(data = subset(age_sex_dag3, ANTHRO.Sex == "F")) + 
        geom_bar(data = subset(age_sex_dag3, ANTHRO.Sex == "M"), aes(y=..count..*(-1))) +
        scale_y_continuous(breaks = seq(-4300,4300,100), labels = abs(seq(-4300,4300,100))) +
        xlab("Population") + ylab("Age")+
        coord_flip()

dev.off()

#create age groups to plot better: 
age_sex_dag3 <- age_sex_dag3 %>%
        mutate(age_groups = case_when(
                ANTHRO.AGE %in% c(20:30) ~ "20-30",
                ANTHRO.AGE %in% c(31:40) ~"31-40",
                ANTHRO.AGE %in% c(41:50) ~ "41-50",
                ANTHRO.AGE %in% c(51:60) ~ "51-60",
                ANTHRO.AGE %in% c(61:70) ~ "61-70", 
                ANTHRO.AGE %in% c(71:80) ~ "71-80",
                ANTHRO.AGE %in% c(81:100) ~ ">80"
        ))

age_sex_dag3$age_groups <- as.factor(age_sex_dag3$age_groups)

age_M_dag3 <- age_sex_dag3[age_sex_dag3$ANTHRO.Sex == "M",]
age_F_dag3 <- age_sex_dag3[age_sex_dag3$ANTHRO.Sex == "F",]


#plot age-sex distribution for age group: 
png("age_sex_distribution2_dag3.png")

gg <- ggplot(age_sex_dag3, aes(x= age_groups))

gg.male <- gg +
        geom_bar(data = subset(age_sex_dag3, ANTHRO.Sex == "M"), 
                 aes(y=..count../sum(..count..), fill = age_groups)) +
        scale_y_continuous('', labels = scales::percent) +
        theme(legend.position = 'none',
              axis.title.y = element_blank(),
              plot.title = element_text(size = 11.5),
              plot.margin = unit(c(0.1,0.2,0.1,-.1), "cm"),
              axis.ticks.y = element_blank(),
              axis.text.y = theme_bw()$axis.text.y, 
              panel.background = element_blank()) + 
        ggtitle("Male") + 
        coord_flip() 


gg.female <- gg +
        geom_bar(data = subset(age_sex_dag3, ANTHRO.Sex == "F"), 
                 aes(y=..count../sum(..count..), fill = age_groups)) +
        scale_y_continuous('', labels = scales::percent,
                           trans = 'reverse') +
        theme(legend.position = 'none',
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.background = element_blank(),
              plot.title = element_text(size = 11.5),
              plot.margin = unit(c(0.1,0,0.1,0.05), "cm"))+
        ggtitle("Female") + 
        coord_flip() +
        ylab ("Age")

grid.arrange (gg.female, gg.male, 
              widths =c(0.5,0.5),
              ncol=2)
dev.off()


#Checking work-related questionnaires, to get time of employment per occupational code: 
work_dag3 <- work2[work2$PSEUDOIDEXT %in% dag3_metaphlan_exposure$PSEUDOIDEXT,]
work_dag3_backup <- work_dag3
work_dag3$WORK8 <- as.factor(work_dag3$WORK8)

questionnaires <- c("Baseline assessment (1A)", "Second assessment (2A)")
work_1a_2a<- work_dag3[work_dag3$ENCOUNTERCODE %in% questionnaires,]

work_2a<- work_1a_2a[work_1a_2a$ENCOUNTERCODE == "Second assessment (2A)",]
work_1a<- work_1a_2a[work_1a_2a$ENCOUNTERCODE == "Baseline assessment (1A)",]

work_1a[sapply(work_1a, is.character)] <- lapply(work_1a[sapply(work_1a, is.character)],
                                                 as.factor)

work_2a[sapply(work_2a, is.character)] <- lapply(work_2a[sapply(work_2a, is.character)],
                                                 as.factor)

work_1a_clean <- work_1a[, !sapply(work_1a, function(col) nlevels(col) ==1)]
work_2a_clean <- work_2a[, !sapply(work_2a, function(col) nlevels(col) ==1)]

#check what is up with the unemployed/not working people: 
##Here, I managed to trace them back to studying or on leave by their work. I did not add this to this code,
## because it's rather just exploratory code and relatively long.
work_1a_unemp <- work_1a_clean[work_1a_clean$WORK8 == " ", ]

#continue the analysis by crossing work-related questions:
questions <- c("PSEUDOIDEXT", "WORK8", "WORK10", "WORK11", "WORK2CA")

work_1a_clean <- work_1a_clean[, which(names(work_1a_clean) %in% questions)]
work_2a_clean <- work_2a_clean[, which(names(work_2a_clean) %in% questions)]

dag3_exposure_work <- merge(dag3_metaphlan_exposure, work_1a_clean, by= "PSEUDOIDEXT")
levels(dag3_exposure_work$WORK8) [levels(dag3_exposure_work$WORK8) == " "] <- "Unemployed"
dag3_exposure_work$WORK8 <- factor(dag3_exposure_work$WORK8, levels = c("<1 year", "1 to 5 years", "6 to 10 years", ">10 years", "Unemployed"))

png("work8_freq.png")

ggplot(dag3_exposure_work, aes(x= WORK8)) +
        geom_histogram( stat = "count", binwidth = 100, boundary =0.5, fill = "#69b3a2", color = "#e9ecef", alpha=0.9) +
        labs(y = "Frequency", x = "\nTime spent on job") +
        theme(panel.background = element_blank(),
              panel.grid.major.y = element_line(size = 0.5, linetype = 'solid', colour = "grey"),
              axis.ticks = element_line(size = 0.5, colour = "grey"), 
              axis.ticks.length = unit(.3, "cm"),
              axis.title.x = element_text(size = (15)),
              axis.text.x = element_text( size = (10)),
              axis.text.y = element_text(size = (10)),
              axis.title.y = element_text( size =(15))
        )

dev.off()

#Cross time of employment with pesticide exposure levels to get cumulative exposure: 
dag3_exp_work <- dag3_exposure_work

dag3_exp_work <- dag3_exp_work %>% 
        mutate(pest_all.2 = case_when(
                pest_all %in% c("0") ~ 0, 
                pest_all %in% c("1") ~ 1, 
                pest_all %in% c("2") ~ 4,),
               pest_herbi.2 = case_when(
                       pest_herbi %in% c("0") ~ 0, 
                       pest_herbi %in% c("1") ~ 1, 
                       pest_herbi %in% c("2") ~ 4,), 
               pest_insec.2 = case_when(
                       pest_insec %in% c("0") ~ 0, 
                       pest_insec %in% c("1") ~ 1, 
                       pest_insec %in% c("2") ~ 4,),
               pest_fungi.2 = case_when(
                       pest_fungi %in% c("0") ~ 0, 
                       pest_fungi %in% c("1") ~ 1, 
                       pest_fungi %in% c("2") ~ 4,), 
               work_to_mult = case_when(
                       WORK8 %in% c(">10 years") ~ 10, 
                       WORK8 %in% c("6 to 10 years") ~ 8,
                       WORK8 %in% c("1 to 5 years") ~ 2.5,
                       WORK8 %in% c("<1 year") ~ 1,
                       WORK8 %in% c("Unemployed") ~ 1,
               )) 


dag3_exp_work$pest_all_cumulative <- dag3_exp_work$pest_all.2 * dag3_exp_work$work_to_mult
dag3_exp_work$pest_herbi_cumulative <- dag3_exp_work$pest_herbi.2 * dag3_exp_work$work_to_mult
dag3_exp_work$pest_insec_cumulative <- dag3_exp_work$pest_insec.2 * dag3_exp_work$work_to_mult
dag3_exp_work$pest_fungi_cumulative <- dag3_exp_work$pest_fungi.2 * dag3_exp_work$work_to_mult

###Do inverse rank transformation of the cumulative values: 
invrank= function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}

dag3_exp_work[, c("pest_all_cumulative", "pest_herbi_cumulative", "pest_insec_cumulative", "pest_fungi_cumulative")] <- apply(dag3_exp_work[,c("pest_all_cumulative", "pest_herbi_cumulative", "pest_insec_cumulative", "pest_fungi_cumulative")],2, invrank )

dag3_exp_work_final<-dag3_exp_work

#prepare dataframe for alpha diversity calculations
to_keep <- c("PSEUDOIDEXT", "pest_all", "pest_herbi", "pest_insec", "pest_fungi", 
             "pest_all_cumulative", "pest_herbi_cumulative", "pest_insec_cumulative", "pest_fungi_cumulative")
dag3_exp_work_final <- dag3_exp_work_final[to_keep]

to_factor <- c("pest_all", "pest_herbi", "pest_insec", "pest_fungi", 
             "pest_all_cumulative", "pest_herbi_cumulative", "pest_insec_cumulative", "pest_fungi_cumulative")

dag3_exp_work_final[to_factor] <- lapply(dag3_exp_work_final[to_factor], factor)


at_least_dag32 <-dag3_metaphlan_exposure[apply(dag3_metaphlan_exposure[7:10], 1, function(x) any(x>0)),]

##### for double check: age & sex distribution among exposed people
age_sex_dag3_exposed <- age_sex_dag3[age_sex_dag3$PSEUDOIDEXT %in% at_least_dag32$PSEUDOIDEXT, ]


###alpha diverstity calculation: 

library(vegan)

DAG3_metaphlan_exp <- DAG3_metaphlan_final[DAG3_metaphlan_final$PSEUDOIDEXT %in% dag3_metaphlan_exposure$PSEUDOIDEXT,]
rownames(DAG3_metaphlan_exp) <- DAG3_metaphlan_exp$PSEUDOIDEXT
DAG3_metaphlan_exp$PSEUDOIDEXT <- NULL

#use only spp level: 
DAG3_metaphlan_exp_spp <- DAG3_metaphlan_exp[,grep("t__", colnames(DAG3_metaphlan_exp), invert = T)]
DAG3_metaphlan_exp_spp <-DAG3_metaphlan_exp_spp[, grep("s__", colnames(DAG3_metaphlan_exp_spp))]
colnames(DAG3_metaphlan_exp_spp) = sub(".*[.]s__", "", colnames(DAG3_metaphlan_exp_spp))
chr_conversion <- colnames(DAG3_metaphlan_exp_spp)
DAG3_metaphlan_exp_spp[, chr_conversion] <- lapply(chr_conversion, function(x) as.numeric(DAG3_metaphlan_exp_spp[[x]]))


DAG3_aD_spp_all <- diversity(DAG3_metaphlan_exp_spp, index = "shannon")

#add alpha diversity to main df 
dag3_ad_spp_012 <- stack(DAG3_aD_spp_all)
dag3_ad_spp_012 <- dag3_ad_spp_012 %>% 
        rename(
                alpha.div.shannon = values, 
                PSEUDOIDEXT = ind
        )

#prep for running linear regressions between alpha diversity and pesticide exposure
dag3_exp_work_final <- merge(dag3_exp_work_final, dag3_ad_spp_012, by = "PSEUDOIDEXT")
dag3_exp_work_final <- merge(dag3_exp_work_final, age_sex_dag3, by = "PSEUDOIDEXT")
phenos_other <- phenos_new[,c("META.BATCH", "META.POOP.COLLECTION_SEASON", "PSEUDOIDEXT", "META.DNA.postclean.reads")]
dag3_exp_work_final <- merge(dag3_exp_work_final, phenos_other, by= "PSEUDOIDEXT")

dag3_exp_work_final2 <- dag3_exp_work_final
dag3_exp_work_final2[,6:9] <- apply(dag3_exp_work_final2[, 6:9], 2, function(x) as.numeric(as.character(x)))

##Here, I ran the linear regressions one by one and separately with/without covariates to compare the results and see whether
## the model is good (checking the effect estimates of each covariate). 
#alpha-pesticides lm
lm_alpha_all_dag3 <- lm(alpha.div.shannon ~  pest_all, data = dag3_exp_work_final)
lm_alpha_all_cum_dag3 <- lm(alpha.div.shannon ~  pest_all_cumulative, data = dag3_exp_work_final2)

lm_alpha_all_dag3_2 <- lm(alpha.div.shannon ~ ANTHRO.Sex + ANTHRO.AGE + META.POOP.COLLECTION_SEASON +pest_all, data = dag3_exp_work_final)
lm_alpha_all_cum_dag3_2 <- lm(alpha.div.shannon ~ ANTHRO.Sex + ANTHRO.AGE + META.POOP.COLLECTION_SEASON +pest_all_cumulative, data = dag3_exp_work_final2)

#alpha-herbicides lm
lm_alpha_herbi_dag3 <- lm(alpha.div.shannon ~  pest_herbi, data = dag3_exp_work_final)
lm_alpha_herbi_cum_dag3 <- lm(alpha.div.shannon ~  pest_herbi_cumulative, data = dag3_exp_work_final2)

lm_alpha_herbi_dag3_2 <- lm(alpha.div.shannon ~ ANTHRO.Sex + ANTHRO.AGE + META.POOP.COLLECTION_SEASON +pest_herbi, data = dag3_exp_work_final)
lm_alpha_herbi_cum_dag3_2 <- lm(alpha.div.shannon ~ ANTHRO.Sex + ANTHRO.AGE + META.POOP.COLLECTION_SEASON +pest_herbi_cumulative, data = dag3_exp_work_final2)

#alpha-insecticides lm
lm_alpha_insec_dag3 <- lm(alpha.div.shannon ~  pest_insec, data = dag3_exp_work_final)
lm_alpha_insec_cum_dag3 <- lm(alpha.div.shannon ~  pest_insec_cumulative, data = dag3_exp_work_final2)

lm_alpha_insec_dag3_2 <- lm(alpha.div.shannon ~ ANTHRO.Sex + ANTHRO.AGE + META.POOP.COLLECTION_SEASON +pest_insec, data = dag3_exp_work_final)
lm_alpha_insec_cum_dag3_2 <- lm(alpha.div.shannon ~ ANTHRO.Sex + ANTHRO.AGE + META.POOP.COLLECTION_SEASON +pest_insec_cumulative, data = dag3_exp_work_final2)

#alpha-fungicides lm
lm_alpha_fungi_dag3 <- lm(alpha.div.shannon ~  pest_fungi, data = dag3_exp_work_final)
lm_alpha_fungi_cum_dag3 <- lm(alpha.div.shannon ~  pest_fungi_cumulative, data = dag3_exp_work_final2)

lm_alpha_fungi_dag3_2 <- lm(alpha.div.shannon ~ ANTHRO.Sex + ANTHRO.AGE + META.POOP.COLLECTION_SEASON +pest_fungi, data = dag3_exp_work_final)
lm_alpha_fungi_cum_dag3_2 <- lm(alpha.div.shannon ~ ANTHRO.Sex + ANTHRO.AGE + META.POOP.COLLECTION_SEASON +pest_fungi_cumulative, data = dag3_exp_work_final2)

#run kruskal-wallis and pairwise test to check the results of the lm (here I only did it with one example)

kruskal.test(alpha.div.shannon ~ pest_all, data= dag3_exp_work_final)
pairwise.wilcox.test(dag3_exp_work_final$alpha.div.shannon, dag3_exp_work_final$pest_all, p.adjust.method = "BH")

#plot alpha-diversity: 
library(viridis)
library(ggpubr)
library(gridExtra)


pdf("alpha_div_shannon_redo.pdf", useDingbats = F)

comp_1 <- list(c("0", "1"), c("0","2"), c( "1", "2"))
labs <- c("No", "Low", "High")

p1 <- ggplot(dag3_exp_work_final , aes(x= pest_all, y= alpha.div.shannon, fill = pest_all)) + 
        geom_jitter(color = "grey", size= 0.4, alpha= 0.3)+
        geom_violin() + 
        scale_color_gradientn(colors= rainbow (3))+ 
        labs(y = "Shannon diversity", x = "\nPesticide exposure") +
        theme(  panel.background = element_blank(),
                axis.ticks = element_line(size = 0.5, colour = "grey"), 
                axis.ticks.length = unit(.3, "cm"),
                axis.title.x = element_text(size = (15)),
                axis.text.x = element_text( size = (10)),
                axis.text.y = element_text(size = (10)),
                axis.title.y = element_text( size =(15)),  
                legend.position = "none"
        )+
        stat_compare_means(comparisons = comp_1) + 
        scale_x_discrete(labels= labs) +
        stat_compare_means(label.y = 5.0)


p2 <- ggplot(dag3_exp_work_final, aes(x= pest_herbi, y= alpha.div.shannon, fill = pest_herbi)) + 
        geom_jitter(color = "grey", size= 0.4, alpha= 0.3)+
        geom_violin() + 
        scale_color_gradientn(colors= rainbow (3))+ 
        labs(y = "Shannon diversity", x = "\nHerbicide exposure") +
        theme(  panel.background = element_blank(),
                axis.ticks = element_line(size = 0.5, colour = "grey"), 
                axis.ticks.length = unit(.3, "cm"),
                axis.title.x = element_text(size = (15)),
                axis.text.x = element_text( size = (10)),
                axis.text.y = element_text(size = (10)),
                axis.title.y = element_text( size =(15)),  
                legend.position = "none"
        )+
        stat_compare_means(comparisons = comp_1) + 
        scale_x_discrete(labels= labs) +
        stat_compare_means(label.y = 5.0)


p3<-  ggplot(dag3_exp_work_final, aes(x= pest_insec, y= alpha.div.shannon, fill = pest_insec)) + 
        geom_jitter(color = "grey", size= 0.4, alpha= 0.3)+
        geom_violin() + 
        scale_color_gradientn(colors= rainbow (3))+ 
        labs(y = "Shannon diversity", x = "\nInsecticide exposure") +
        theme(  panel.background = element_blank(),
                axis.ticks = element_line(size = 0.5, colour = "grey"), 
                axis.ticks.length = unit(.3, "cm"),
                axis.title.x = element_text(size = (15)),
                axis.text.x = element_text( size = (10)),
                axis.text.y = element_text(size = (10)),
                axis.title.y = element_text( size =(15)),  
                legend.position = "none"
        )+
        stat_compare_means(comparisons = comp_1) + 
        scale_x_discrete(labels= labs) +
        stat_compare_means(label.y = 5.0)


p4<-  ggplot(dag3_exp_work_final, aes(x= pest_fungi, y= alpha.div.shannon, fill = pest_fungi)) + 
        geom_jitter(color = "grey", size= 0.4, alpha= 0.3)+
        geom_violin() + 
        scale_color_gradientn(colors= rainbow (3))+ 
        labs(y = "Shannon diversity", x = "\nFungicide exposure") +
        theme(  panel.background = element_blank(),
                axis.ticks = element_line(size = 0.5, colour = "grey"), 
                axis.ticks.length = unit(.3, "cm"),
                axis.title.x = element_text(size = (15)),
                axis.text.x = element_text( size = (10)),
                axis.text.y = element_text(size = (10)),
                axis.title.y = element_text( size =(15)),  
                legend.position = "none"
        )+
        stat_compare_means(comparisons = comp_1) + 
        scale_x_discrete(labels= labs) +
        stat_compare_means(label.y = 5.0)

ggarrange(p1, p2, p3, p4, 
          labels = c("a", "b", "c", "d"), 
          ncol=2, nrow = 2)

dev.off()

####beta diversity: 
library(phyloseq)

dist_BC <- vegdist(DAG3_metaphlan_exp_spp, method = "bray")

bray_test = adonis(dist_BC ~ dag3_exp_work_final$ANTHRO.Sex + dag3_exp_work_final$ANTHRO.AGE + dag3_exp_work_final$pest_all)
bray_test = adonis(dist_BC ~ dag3_exp_work_final$ANTHRO.Sex + dag3_exp_work_final$ANTHRO.AGE + dag3_exp_work_final$pest_fungi, permutations = 1000)

#statistical tests on particular bacterial features
library(foreach)

dag3_taxa_transf <- list(dag3_taxa_transf)
dag3_taxa_transf <- do.call(rbind.data.frame, dag3_taxa_transf)

DAG3_metaphlan_exp_spp2 <- dag3_taxa_transf[,grep("t__", colnames(dag3_taxa_transf), invert = T)] #366 taxa 
DAG3_metaphlan_exp_spp2 <-DAG3_metaphlan_exp_spp2[, grep("s__", colnames(DAG3_metaphlan_exp_spp2))]
colnames(DAG3_metaphlan_exp_spp2) = sub(".*[.]s__", "", colnames(DAG3_metaphlan_exp_spp2))
chr_conversion2 <- colnames(DAG3_metaphlan_exp_spp2)
DAG3_metaphlan_exp_spp2 <- as.data.frame(DAG3_metaphlan_exp_spp2)
DAG3_metaphlan_exp_spp2[, chr_conversion2] <- lapply(chr_conversion2, function(x) as.numeric(DAG3_metaphlan_exp_spp2[[x]]))


DAG3_metaphlan_spp_5percent <- DAG3_metaphlan_exp_spp2[, colSums(DAG3_metaphlan_exp_spp2>0) > 360 ] #139 spp are present in more than 5% of the samples 
DAG3_metaphlan_spp_filtered <- DAG3_metaphlan_exp_spp2[, colSums(DAG3_metaphlan_exp_spp2>0) >1 ] #530 spp are present in more than 1 sample

dag3_exp_work_final2<- dag3_exp_work_final
dag3_exp_work_final2 [, c(6:9)] <- sapply(dag3_exp_work_final2 [,c(6:9)], as.character)
dag3_exp_work_final2 [, c(6:9)] <- sapply(dag3_exp_work_final2 [,c(6:9)], as.numeric)

DAG3_metaphlan_spp_5percent2 <- subset(DAG3_metaphlan_spp_5percent, rownames(DAG3_metaphlan_spp_5percent) %in% dag3_exp_work_final2$PSEUDOIDEXT)

##Linear regressions between spp present in more than 5 percent of the participants and the exposure levels to each pesticide class
## after that, I ran an anova test between the linear regression with and without the pesticide exposure variable. 
## FDR corrected the p-val for those: 
results = data.frame()
results_anova = data.frame()

for(i in 1:ncol(DAG3_metaphlan_spp_5percent2)){
        for(j in 2:9) {
                lm0 = lm(DAG3_metaphlan_spp_5percent2[,i] ~ dag3_exp_work_final2$ANTHRO.AGE+ dag3_exp_work_final2$ANTHRO.Sex + dag3_exp_work_final2$META.POOP.COLLECTION_SEASON + dag3_exp_work_final2$META.DNA.postclean.reads)
                lm1 = lm(DAG3_metaphlan_spp_5percent2[,i] ~ dag3_exp_work_final2$ANTHRO.AGE+ dag3_exp_work_final2$ANTHRO.Sex + dag3_exp_work_final2$META.POOP.COLLECTION_SEASON + dag3_exp_work_final2$META.DNA.postclean.reads + dag3_exp_work_final2[,j])
                anova1 = anova(lm1, lm0)
                summary1 = summary(lm1)
                summary1$coefficients[8:nrow(summary1$coefficients),1:4]
                oneReport = data.frame(bac = colnames(DAG3_metaphlan_spp_5percent2)[i],
                                       pheno = colnames(dag3_exp_work_final2)[j],
                                       coef1 = summary1$coefficients[8,1],
                                       SE1 = summary1$coefficients[8,2],
                                       T1 = summary1$coefficients[8,3],
                                       P1 = summary1$coefficients[8,4], 
                                       coef2 = summary1$coefficients[nrow(summary1$coefficients),1],
                                       SE2 = summary1$coefficients[nrow(summary1$coefficients),2],
                                       T2 = summary1$coefficients[nrow(summary1$coefficients),3],
                                       P2 = summary1$coefficients[nrow(summary1$coefficients),4]
                )
                
                results_an = data.frame(bac = colnames(DAG3_metaphlan_spp_5percent2)[i],
                                        pheno = colnames(dag3_exp_work_final2)[j],
                                        P_anova = anova1$Pr[[2]])
                
                results = rbind(results,oneReport)
                results_anova = rbind(results_anova, results_an)
                
        }
}

results$FDR1 = p.adjust(results$P1 , method = "BH" )
results$FDR2 = p.adjust(results$P2 , method = "BH" )


results_psignf <- results[(results$P1 < 0.05 | results$P2 <0.05),]

write.table(results_psignf, "taxa_exp_lm_signf.txt", sep = "\t", row.names = T, col.names = T )

signif_bac <- unique(results_psignf$bac)

#double check: 
test_bacteroides0 <- lm(DAG3_metaphlan_spp_5percent2$Bacteroides_caccae  ~ dag3_exp_work_final2$ANTHRO.AGE+ dag3_exp_work_final2$ANTHRO.Sex + dag3_exp_work_final2$META.POOP.COLLECTION_SEASON + dag3_exp_work_final2$META.DNA.postclean.reads)
test_bacteroides1 <- lm(DAG3_metaphlan_spp_5percent2$Bacteroides_caccae  ~ dag3_exp_work_final2$ANTHRO.AGE+ dag3_exp_work_final2$ANTHRO.Sex + dag3_exp_work_final2$META.POOP.COLLECTION_SEASON + dag3_exp_work_final2$META.DNA.postclean.reads + dag3_exp_work_final2$pest_insec)
test_bac_anova <- anova(test_bacteroides1, test_bacteroides0)

test_bifido0 <- lm(DAG3_metaphlan_spp_5percent2$Bifidobacterium_bifidum  ~ dag3_exp_work_final2$ANTHRO.AGE+ dag3_exp_work_final2$ANTHRO.Sex + dag3_exp_work_final2$META.POOP.COLLECTION_SEASON + dag3_exp_work_final2$META.DNA.postclean.reads)
test_bifido1 <- lm(DAG3_metaphlan_spp_5percent2$Bifidobacterium_bifidum  ~ dag3_exp_work_final2$ANTHRO.AGE+ dag3_exp_work_final2$ANTHRO.Sex + dag3_exp_work_final2$META.POOP.COLLECTION_SEASON + dag3_exp_work_final2$META.DNA.postclean.reads + dag3_exp_work_final2$pest_herbi)
test_bif_anova <- anova(test_bifido1, test_bifido0)
summary_test <- summary(test_bif_anova)

pairwise.wilcox.test(DAG3_metaphlan_spp_5percent2$Alistipes_sp_An31A  , dag3_exp_work_final2$pest_herbi, p.adjust.method = "BH")
kruskal.test(DAG3_metaphlan_spp_5percent2$Alistipes_sp_An31A   ~ dag3_exp_work_final2$pest_herbi)



####PLOT
bac_to_plot <- DAG3_metaphlan_spp_5percent2[colnames(DAG3_metaphlan_spp_5percent2) %in% signif_bac]
bac_to_plot$PSEUDOIDEXT <- rownames(bac_to_plot)
bac_to_plot <- merge(dag3_exp_work_final2,bac_to_plot, by = "PSEUDOIDEXT")

#phenotype - exp correlations: 
phenos2 <- phenos_new
phenos2 <- phenos2 %>% mutate_if(is.character, as.factor)

phenos3 <- phenos2[grep("(PSEUDO|ANTHRO|EXP.DIET|MED|META).*", names(phenos2))]

#get binary phenotypes: 
is_binary <- function(x) (length(unique(x)) == 2)
binary_pheno <- sapply(phenos3, is_binary)
binary_pheno <- as.data.frame(binary_pheno)
binary_pheno$pheno = rownames(binary_pheno)
binary_pheno <- binary_pheno[binary_pheno$binary_pheno == "TRUE",]

binary_pheno <- phenos3[names(phenos3) %in% binary_pheno$pheno]

#run glm on binary traits and exposure

binary_pheno_results = data.frame()

for(i in 2:ncol(binary_pheno)){
        for(j in 2:9) {
                lm1 = glm(binary_pheno[,i] ~ dag3_exp_work_final2$ANTHRO.AGE+ dag3_exp_work_final2$ANTHRO.Sex + dag3_exp_work_final2$META.POOP.COLLECTION_SEASON + dag3_exp_work_final2$META.DNA.postclean.reads + dag3_exp_work_final2[,j], family = binomial())
                summary1 = summary(lm1)
                summary1$coefficients[nrow(summary1$coefficients),1:4]
                oneReport = data.frame(pheno = colnames(binary_pheno)[i],
                                       exp = colnames(dag3_exp_work_final2)[j],
                                      coef = summary1$coefficients[nrow(summary1$coefficients),1],
                                       SE = summary1$coefficients[nrow(summary1$coefficients),2],
                                       Z = summary1$coefficients[nrow(summary1$coefficients),3],
                                       P = summary1$coefficients[nrow(summary1$coefficients),4]
                )
                binary_pheno_results = rbind(binary_pheno_results,oneReport)
        }
}
binary_pheno_results$FDR = p.adjust(binary_pheno_results$P,method = "BH" )
binary_pheno_results_signif <- binary_pheno_results[binary_pheno_results$P < 0.05,]

signif_pheno_exp <- unique(binary_pheno_results_signif$pheno)


lm1 = glm(binary_pheno$MED.SURGERY.Appendectomy ~ dag3_exp_work_final2$ANTHRO.AGE+ dag3_exp_work_final2$ANTHRO.Sex + dag3_exp_work_final2$META.POOP.COLLECTION_SEASON + dag3_exp_work_final2$META.DNA.postclean.reads + dag3_exp_work_final2$pest_all, family = binomial())


###PATHWAYS: 

pathways <- read.csv(file.path("/groups/umcg-dag3/tmp01/DAG3_biobakery_v3_results_sorted/results_merged/DAG3_humann3_headersfixed_transposed_cleaned_normalized.csv"), header = T, stringsAsFactors = F)
pathways_sum <- read.table("pathway_cleaned_MG.txt", sep="\t", header = T)

pathways_sum2 <- pathways_sum[(pathways_sum$more_than_5 == "yes"),]
pathways_interest <- pathways_sum2$PWY

pathways2 <- pathways[,c(pathways_interest, "ID" )]

names(pathways2)[names(pathways2)== "ID"] <- "DAG3_sampleID"
pathways2 <- merge(participants,pathways2, by = "DAG3_sampleID")
pathways2$DAG3_sampleID <- NULL 

pathways2 <- subset(pathways2, pathways2$PSEUDOIDEXT %in% dag3_exp_work_final2$PSEUDOIDEXT)
pwys <- pathways2
row.names(pwys) <- pwys$PSEUDOIDEXT
pwys$PSEUDOIDEXT <- NULL

#transform data: 
pathways3 <- do_clr_externalWeighting(pwys, pwys)

#run linear reg:
dag3_exp_work_final2 <- subset(dag3_exp_work_final2, dag3_exp_work_final2$PSEUDOIDEXT %in% pathways2$PSEUDOIDEXT)

pwy_results = data.frame()

for(i in 1:ncol(pathways3)){
        for(j in 2:9) {
                lm1 = lm(pathways3[,i] ~ dag3_exp_work_final2$ANTHRO.AGE+ dag3_exp_work_final2$ANTHRO.Sex + dag3_exp_work_final2$META.POOP.COLLECTION_SEASON + dag3_exp_work_final2$META.DNA.postclean.reads + dag3_exp_work_final2[,j])
                summary1 = summary(lm1)
                summary1$coefficients[8:nrow(summary1$coefficients),1:4]
                oneReport = data.frame(PWY = colnames(pathways3)[i],
                                       exp = colnames(dag3_exp_work_final2)[j],
                                       coef1 = summary1$coefficients[8,1],
                                       SE1 = summary1$coefficients[8,2],
                                       T1 = summary1$coefficients[8,3],
                                       P1 = summary1$coefficients[8,4], 
                                       coef2 = summary1$coefficients[nrow(summary1$coefficients),1],
                                       SE2 = summary1$coefficients[nrow(summary1$coefficients),2],
                                       T2 = summary1$coefficients[nrow(summary1$coefficients),3],
                                       P2 = summary1$coefficients[nrow(summary1$coefficients),4]
                )
                pwy_results = rbind(pwy_results,oneReport)
        }
}
pwy_results$FDR1 = p.adjust(pwy_results$P1 , method = "BH" )
pwy_results$FDR2 = p.adjust(pwy_results$P2 , method = "BH" )

pwy_results_psignf <- pwy_results[(pwy_results$P1 < 0.05 | pwy_results$P2 <0.05),]

signif_pwy <- unique(pwy_results_psignf$PWY)

signif_pwy_sum <- subset(pathways_sum2, pathways_sum2$PWY %in% signif_pwy)

pwy_results_p2 <- merge(pwy_results_psignf, signif_pwy_sum, by = "PWY")
pwy_results_p2 <- pwy_results_p2[,-c(14:18)]
pwy_results_psignf[1:3,c(1:3,5,10)]
pwy_results_p2[1:3, c(1:3,5,10,16)]

"Taxa.with.the.pathway"

write.table(pwy_results_p2, "pwy_exp_lm_signf2.txt", sep = "\t", row.names = T, col.names = T )

signif_pwy <- unique(pwy_results_psignf$PWY)
