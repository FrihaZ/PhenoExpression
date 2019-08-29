
###############################################################################################################################
###############################################################################################################################

### Project: PhenoExpression ##################################################################################################
### Script: Expression_Atlas_Mice_Developmental_Stages.R ###################################################################################################
### Purpose: To compare gene expression of mice in different developmental stages using the Expression Atlas dataset ########################################################
### Author: Friha Zafar ####################################################################################################
### Date: 11/07/2019 ##########################################################################################################

###############################################################################################################################
###############################################################################################################################

## load packages ##############################################################################################################

library(dplyr);library(tidyr);library(stringr);library(readr);library(tidyverse);library(ggpubr);library(ggplot2)

################################################################################################################################
################################################################################################################################

# Load Files
EA.Mouse.Development<-read_tsv("D:/MSC RESEARCH PROJECT/GITHUB/data/E-ERAD-401-query-results.tpms.tsv",skip = 4)
View(EA.Mouse.Development)

mgi.ensembl.genemapping<-read.delim("D:/MSC RESEARCH PROJECT/GITHUB/data/mgi.ensembl.genemapping.txt")


mgi.phenotypes<-read.delim("D:/MSC RESEARCH PROJECT/GITHUB/data/mgi.phenotypes.txt")

mp.ancestor.nodes<- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/mp.ancestor.nodes.txt",delim = "\t")

mp.toplevels<- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/mp.toplevels.txt",delim = "\t") %>%
  rename(mp.ancestors=mp.term )  #renamed mp.term to mp.ancestors 

################################################################################################################################
################################################################################################################################

###############~~~~~~~~~~~~~~~~~~~~~~~~CLEAN DATA AND MAP GENE IDS~~~~~~~~~~~~#########

################################################################################################################################
################################################################################################################################

## Map the ENSEMBL gene ID to their MGI.ID


mgi.ensembl.genemapping$Gene.ID<-as.character(mgi.ensembl.genemapping$Gene.ID)

MGI.EA.Mouse.Development<-right_join(EA.Mouse.Development,mgi.ensembl.genemapping, by=c("Gene ID"="Gene.ID") )

View(MGI.EA.Mouse.Development)




#################################################################################################
########################################################################################################

### To link the phenotypes back to top level ###################################################################################

################################################################################################################################
################################################################################################################################

## Join files to one 

mp.toplevels<-left_join(mp.toplevels, mp.ancestor.nodes)   # the left_join keeps the data that is already in mp.toplevels, and add the mp.ancestors column

View(mp.toplevels)

mp.toplevels$mp.term<-as.character(mp.toplevels$mp.term )

################################################################################################################################
############################################################################################################################

## MAP THE MGI.ID to their MP terms and look at expression across the TS for each phenotype

mgi.phenotypes$MGI.ID<-as.character(mgi.phenotypes$MGI.ID)
MGI.EA.Mouse.Development$MGI.Accession.ID<-as.character(MGI.EA.Mouse.Development$MGI.Accession.ID)


MGI.MP.EA.Mouse.Development<-right_join(MGI.EA.Mouse.Development,mgi.phenotypes, by=c("MGI.Accession.ID" ="MGI.ID") )


mgi.phenotypes$MP.ID<-as.character(mgi.phenotypes$MP.ID)
MGI.MP.EA.Mouse.Development$MP.ID<-as.character(MGI.MP.EA.Mouse.Development$MP.ID)


MGI.TOP.MP.EA.Mouse.Development<-left_join(MGI.MP.EA.Mouse.Development,mp.toplevels, by=c( "MP.ID"= "mp.term" ) )


View(MGI.TOP.MP.EA.Mouse.Development)

##################################################################################

###EMBRYO PHENOTYPE


MGI.TOP.MP.EA.Mouse.Development_EMBRYO<-MGI.TOP.MP.EA.Mouse.Development%>%
  select_all()%>%
  filter(mp.description=="embryo phenotype")%>%
  distinct()

View(MGI.TOP.MP.EA.Mouse.Development_EMBRYO)




MGI.EMBRYO<-MGI.TOP.MP.EA.Mouse.Development_EMBRYO%>%
  gather("Development Stage","Expression", `4-somite stage` :`36-somite stage`)%>%
  select_all()%>%
  distinct()

View(MGI.EMBRYO)
# 

sort_MEAN.MGI.EMBRYO<-MGI.EMBRYO


sort_MEAN.MGI.EMBRYO$`Development Stage` <- factor(sort_MEAN.MGI.EMBRYO$`Development Stage`,
                                                                levels=rev(unique(sort_MEAN.MGI.EMBRYO$`Development Stage`)))

View(sort_MEAN.MGI.EMBRYO)


#write.csv(sort_MEAN.Pheno.Mouse.Development,'./Output_Files/Mice/EA_DMDDs DATASET/sort_MEAN.Pheno.Mouse.Development.csv')





sort_MEAN.MGI.EMBRYO_<-sort_MEAN.MGI.EMBRYO%>%
  select(`Development Stage` , Expression, MP.DESCRIPTION)%>%
  group_by(`Development Stage`, MP.DESCRIPTION)%>%
  mutate(Mean= mean(Expression, na.rm = TRUE), SD=sd(Expression, na.rm = TRUE),
         SE= sd(Expression, na.rm = TRUE)/sqrt(length(Expression)))
  distinct(MP.DESCRIPTION, Mean, SE, SD, `Development Stage`)

View(sort_MEAN.MGI.EMBRYO_)

sort_MEAN.MGI.EMBRYO_<-sort_MEAN.MGI.EMBRYO_%>%
  select_all()%>%
  drop_na() %>%
  distinct()



library(randomcoloR)
n <-402
palette_ <- distinctColorPalette(n)

MEAN.Pheno.Mouse.Development_P1<-ggplot(sort_MEAN.MGI.EMBRYO_, aes(x =sort_MEAN.MGI.EMBRYO_$`Development Stage`, 
                                                                               y =sort_MEAN.MGI.EMBRYO_$Mean,
                                                                               group=sort_MEAN.MGI.EMBRYO_$MP.DESCRIPTION,
                                                                               colour=sort_MEAN.MGI.EMBRYO_$MP.DESCRIPTION))%>%
  +labs(x= "Developmental Stage", y="Mean Gene Expression (TPM)", 
        subtitle= "Mean Gene Expression in Each Developmental Stage for 'Embryo Phenotype' Mammalian Phenotype Ontologies")%>%
  + scale_y_continuous(breaks = seq(0, 1000, by = 100))%>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10, angle = 45, hjust = 1, vjust = 1)) %>%
  +geom_line()%>%
  +geom_point(size = 4)%>%
  #+geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=1.3)%>%
  + scale_color_manual(values =palette_, aesthetics = "colour")%>%
  +guides(color=guide_legend(title="Top-Level MP-Term"))%>%
  +  xlim(rev(levels(sort_MEAN.MGI.EMBRYO_$`Development Stage`)))
#+coord_flip()

# +facet_wrap( . ~ Pheno.Mouse.Development$mp.description)


#MEAN.Pheno.Mouse.Development_P2



##!!!UNHASH TO SAVE PLOT!!! 

library(cowplot)
# 
save_plot(paste("./Plots/Mouse/EA_Developmental_Stages/Mean gene expression/MEAN.Pheno.Mouse.Development_P1.png"),
          MEAN.Pheno.Mouse.Development_P1 ,  base_width=900, base_height=900, units="mm") 
# #  



######
#THE PLOT



MEAN.Pheno.Mouse.Development_P2<-ggplot(sort_MEAN.MGI.EMBRYO_, aes(x =sort_MEAN.MGI.EMBRYO_$`Development Stage`, 
                                                                   y =sort_MEAN.MGI.EMBRYO_$Mean,
                                                                   group=sort_MEAN.MGI.EMBRYO_$MP.DESCRIPTION,
                                                                   colour=sort_MEAN.MGI.EMBRYO_$MP.DESCRIPTION))%>%
  +labs(x= "Developmental Stage", y="Mean Gene Expression (TPM)", 
        subtitle= "Mean Gene Expression in Each Developmental Stage for 'Embryo Phenotype' Mammalian Phenotype Ontologies")%>%
  + scale_y_continuous(breaks = seq(0, 1000, by = 100))%>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = 'none',
         #legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10, angle = 45, hjust = 1, vjust = 1)) %>%
  +geom_line()%>%
  +geom_point(size = 4)%>%
  #+geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=1.3)%>%
  + scale_color_manual(values =palette_, aesthetics = "colour")%>%
  +guides(color=guide_legend(title="Top-Level MP-Term"))%>%
  +  xlim(rev(levels(sort_MEAN.MGI.EMBRYO_$`Development Stage`)))
#+coord_flip()

# +facet_wrap( . ~ Pheno.Mouse.Development$mp.description)


#MEAN.Pheno.Mouse.Development_P2



##!!!UNHASH TO SAVE PLOT!!! 

library(cowplot)
# 
save_plot(paste("./Plots/Mouse/EA_Developmental_Stages/Mean gene expression/MEAN.Pheno.Mouse.Development_P2.png"),
          MEAN.Pheno.Mouse.Development_P2 ,  base_width=297, base_height=210, units="mm") 
# #  







###################################################################################################
##REMOVE MP.DESCRIPTION

MGI.TOP.MP.EA.Mouse.Development<-MGI.TOP.MP.EA.Mouse.Development%>%
  select(-MP.DESCRIPTION)%>%
  drop_na()%>%
  distinct()


#####################################################################################################
Pheno.Mouse.Development_<-MGI.TOP.MP.EA.Mouse.Development%>%
  select(-MP.ID)%>%
  filter(mp.description=="adipose tissue phenotype")%>%
  distinct()

str(Pheno.Mouse.Development_)


Pheno.Mouse.Development_<-Pheno.Mouse.Development_%>%
  gather("Development Stage","Expression", `4-somite stage` :`36-somite stage`)%>%
  select_all()%>%
  distinct()
# 

Pheno.Mouse.Development_$`Development Stage`<-as.vector(Pheno.Mouse.Development_$`Development Stage`)

Pheno.Mouse.Development__ordered<-Pheno.Mouse.Development_[order(Pheno.Mouse.Development_$`Development Stage`, decreasing = TRUE ),]
View(Pheno.Mouse.Development_)

sort_Pheno.Mouse.Development_<-Pheno.Mouse.Development_

library(gtools)
sort_Pheno.Mouse.Development_$`Development Stage`<-mixedsort(Pheno.Mouse.Development_$`Development Stage`)


View(sort_Pheno.Mouse.Development_)
sort_Pheno.Mouse.Development_$`Development Stage` <- factor(sort_Pheno.Mouse.Development_$`Development Stage`,
                                                           levels=rev(unique(sort_Pheno.Mouse.Development_$`Development Stage`)))


Pheno.Mouse.Development__try<-ggplot(sort_Pheno.Mouse.Development_, 
                                   aes(x =sort_Pheno.Mouse.Development_$`Development Stage`, 
                                       y =sort_Pheno.Mouse.Development_$Expression, 
                                       fill=sort_Pheno.Mouse.Development_$`Development Stage`))%>%
  +labs(x= "Development Stage", y="Gene Expression (TPM)", 
        subtitle= "Gene Expression in Each Developmental Stage for ")%>%
  +theme(axis.title=element_text(size=10,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),# linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10),
         legend.position = "none") %>%
  + scale_y_continuous(breaks = seq(0, 7000, by = 500))%>%
  +geom_boxplot(outlier.shape = 1)%>%
  +coord_flip()
#+facet_wrap( . ~ Pheno.Mouse.Development_$mp.description)

Pheno.Mouse.Development__try
















#########################################################################################

## Filter data for each phenotype 

list_phenotypes<-MGI.TOP.MP.EA.Mouse.Development%>%
  select(mp.description)%>%
  distinct()
View(list_phenotypes)

Pheno_DS<-function(phenotype_name){  
  Pheno.Mouse.Developments<-MGI.TOP.MP.EA.Mouse.Development%>%
    select(-MP.ID)%>%
    filter(mp.description==phenotype_name)%>%
    distinct()
  
  str(Pheno.Mouse.Developments)
  
  
  Pheno.Mouse.Developments<-Pheno.Mouse.Developments%>%
    gather("Development Stage","Expression", `4-somite stage` :`36-somite stage`)%>%
    select_all()%>%
    distinct()
  # 
  
  Pheno.Mouse.Developments$`Development Stage`<-as.vector(Pheno.Mouse.Developments$`Development Stage`)
  
  Pheno.Mouse.Developments_ordered<-Pheno.Mouse.Developments[order(Pheno.Mouse.Developments$`Development Stage`, decreasing = TRUE ),]
  View(Pheno.Mouse.Developments)
  
  sort_Pheno.Mouse.Developments<-Pheno.Mouse.Developments
  library(gtools)
  sort_Pheno.Mouse.Developments$`Development Stage`<-mixedsort(Pheno.Mouse.Developments$`Development Stage`)
  
  
  View(sort_Pheno.Mouse.Developments)
  sort_Pheno.Mouse.Developments$`Development Stage` <- factor(sort_Pheno.Mouse.Developments$`Development Stage`,
                                                             levels=rev(unique(sort_Pheno.Mouse.Developments$`Development Stage`)))
  
  
  Pheno.Mouse.Developments_P1<-ggplot(sort_Pheno.Mouse.Developments, 
                                     aes(x =sort_Pheno.Mouse.Developments$`Development Stage`, 
                                                                       y =sort_Pheno.Mouse.Developments$Expression, 
                                         fill=sort_Pheno.Mouse.Developments$`Development Stage`))%>%
    +labs(x= "Development Stage", y="Gene Expression (TPM)", 
          subtitle= paste0("Gene Expression in Each Developmental Stage for ", phenotype_name))%>%
    +theme(axis.title=element_text(size=10,face="bold"),
           plot.caption=element_text(face = "italic", size=12, hjust = 0),
           text = element_text(size=10),# linetype = "solid"),
           axis.text.y= element_text(size=10),
           axis.text.x= element_text(size=10),
           legend.position = "none") %>%
    + scale_y_continuous(breaks = seq(0, 50000, by = 1000))%>%
    +geom_boxplot(outlier.shape = 1)%>%
    +coord_flip()
      #+facet_wrap( . ~ Pheno.Mouse.Developments$mp.description)
  
  
  c<-gsub(" ", "_", phenotype_name, fixed=TRUE)
  c<-gsub("/", "_", c, fixed=TRUE)
  c<-gsub("-", "_", c, fixed=TRUE)
  #    
  print(c)
  # 
  library(cowplot)
  #    
  save_plot(paste("./Plots/Mouse/EA_Developmental_Stages/Numb_GENES_EACH_PHENOTYPE/", c,".png"),
            Pheno.Mouse.Developments_P1 ,base_height= 5.5 ,base_aspect_ratio = 2) 
  #    
  
  
}

# 
# Pheno_DS("liver/biliary system phenotype")
# 
# Pheno_DS("embryo phenotype")
# 
# Pheno_DS("growth/size/body region phenotype")
# 
# Pheno_DS("nervous system phenotype")
# 
# Pheno_DS("hematopoietic system phenotype")
# 
# Pheno_DS("homeostasis/metabolism phenotype")
# 
# Pheno_DS("vision/eye phenotype")
# 
# Pheno_DS("integument phenotype")
# 
# Pheno_DS("immune system phenotype")
# 
# Pheno_DS("mortality/aging")
# 
# Pheno_DS("cardiovascular system phenotype")
# 
# Pheno_DS("limbs/digits/tail phenotype")
# 
# Pheno_DS("craniofacial phenotype")
# 
# Pheno_DS("muscle phenotype")
# 
# Pheno_DS("endocrine/exocrine gland phenotype")
# 
# Pheno_DS("digestive/alimentary phenotype")
# 
# Pheno_DS("hearing/vestibular/ear phenotype")
# 
# Pheno_DS("skeleton phenotype")
# 
# Pheno_DS("respiratory system phenotype")
# 
# Pheno_DS("behavior/neurological phenotype")
# 
# Pheno_DS("taste/olfaction phenotype")
# 
# Pheno_DS("normal phenotype")
# 
# Pheno_DS("pigmentation phenotype")
# 
# Pheno_DS("reproductive system phenotype")
# 
# Pheno_DS("neoplasm")
# 
# Pheno_DS("renal/urinary system phenotype")
# 
# Pheno_DS("adipose tissue phenotype")
# 
# 
# # # 
# # # # !!!!!!! UNHASH TO USE THE FOR-LOOP !!!!!!!!!!
# # # 
# # for(p in 1:nrow(list_phenotypes)) {
# #     pp <- list_phenotypes[p,1]     #FIRST COLUMN
# #     print(pp) # TO CHECK IF THE CORRECT ROW IS BEING TAKEN
# #   
# #     c= pp 
# # #    
# # #    
# # #    
# #     Pheno.Mouse_P1<-Pheno_DS(phenotype_name=pp)
# #      
# # #    
# # #     # REMOVE SPACES, /, -, SO THAT THE FILE CAN BE SAVED 
# # #    
# #     c<-gsub(" ", "_", pp, fixed=TRUE)
# #     c<-gsub("/", "_", c, fixed=TRUE)
# #     c<-gsub("-", "_", c, fixed=TRUE)
# # #    
# #     print(c)
# # # 
# # # # SAVE PLOT----CHANGE LOCATION!
# #     library(cowplot)
# # #    
# #     save_plot(paste("./Plots/Mouse/EA_Developmental_Stages/Numb_GENES_EACH_PHENOTYPE/", c,".pdf"),
# #               Pheno.Mouse_P1 ,base_height= 5.5 ,base_aspect_ratio = 2) 
# # #    
# # #    
# #  }

##########################################
# mean for each phenotype at each stage


Pheno.Mouse.Development<-MGI.TOP.MP.EA.Mouse.Development%>%
  gather("Development Stage","Expression", `4-somite stage` :`36-somite stage`)%>%
  select_all()%>%
  distinct()
# 

sort_MEAN.Pheno.Mouse.Development<-Pheno.Mouse.Development

#library(gtools)
#sort_MEAN.Pheno.Mouse.Development$`Development Stage`<-mixedsort(sort_MEAN.Pheno.Mouse.Development$`Development Stage`)


sort_MEAN.Pheno.Mouse.Development$`Development Stage` <- factor(sort_MEAN.Pheno.Mouse.Development$`Development Stage`,
                                                                levels=rev(unique(sort_MEAN.Pheno.Mouse.Development$`Development Stage`)))

View(sort_MEAN.Pheno.Mouse.Development)


#write.csv(sort_MEAN.Pheno.Mouse.Development,'./Output_Files/Mice/EA_DMDDs DATASET/sort_MEAN.Pheno.Mouse.Development.csv')





sort_MEAN.Pheno.Mouse.Development<-sort_MEAN.Pheno.Mouse.Development%>%
  select(`Development Stage`, Expression, mp.description)%>%
  group_by(`Development Stage`, mp.description)%>%
  mutate(Mean= mean(Expression))%>%
  distinct()




View(sort_MEAN.Pheno.Mouse.Development)



# CREATE A PLOT FOR EACH PHENOTYPE WITH THE GENE EXPRESSIONS FOR EACH TS LEVEL


library(randomcoloR)
n <-70
palette <- distinctColorPalette(n)
# 

MEAN.Pheno.Mouse.Development_P1<-ggplot(sort_MEAN.Pheno.Mouse.Development, aes(x =sort_MEAN.Pheno.Mouse.Development$mp.description , 
                                                      y =sort_MEAN.Pheno.Mouse.Development$Mean,
                                                      group=sort_MEAN.Pheno.Mouse.Development$`Development Stage`,
                                                      color=sort_MEAN.Pheno.Mouse.Development$`Development Stage` ))%>%
  +labs(x= "Top-Level MP-Term", y="Mean Gene Expression (TPM)", 
        subtitle= "Mean Gene Expression in Each Developmental Stage for Each Mammalian Phenotype (MP) Top-Level")%>%
  + scale_y_continuous(breaks = seq(0, 200, by = 10))%>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10,angle = 45, hjust = 1, vjust = 1)) %>%
  +geom_point(size=5)%>%
  + scale_color_manual(values =palette, aesthetics = "colour")%>%
  +guides(color=guide_legend(title="Developmental Stage: "))%>%
  +geom_line()
  #+coord_flip()
  
 # +facet_wrap( . ~ Pheno.Mouse.Development$mp.description)
  #+geom_line(color="black")%>%


MEAN.Pheno.Mouse.Development_P1

# # 
# # 
 library(cowplot)
# #     
save_plot(paste("./Plots/Mouse/EA_Developmental_Stages/Mean gene expression/MEAN_PHENO_EXPRESSION.png"),
            MEAN.Pheno.Mouse.Development_P1 ,base_height= 10 ,base_aspect_ratio = 2) 
# #     
#     








##############################################################################################################


# For each developmental stage
library(randomcoloR)
number <-60
number <- distinctColorPalette(number)
# 
# 
# library(RColorBrewer)
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# pie(rep(1,n), col=sample(col_vector, n))



MEAN.Pheno.Mouse.Development_P2<-ggplot(sort_MEAN.Pheno.Mouse.Development, aes(x =sort_MEAN.Pheno.Mouse.Development$`Development Stage`, 
                                                                          y =sort_MEAN.Pheno.Mouse.Development$Mean,
                                                                          group=sort_MEAN.Pheno.Mouse.Development$mp.description,
                                                                          colour=sort_MEAN.Pheno.Mouse.Development$mp.description))%>%
  +labs(x= "Developmental Stage", y="Mean Gene Expression (TPM)", 
        subtitle= "Mean Gene Expression in Each Developmental Stage for Each Mammalian Phenotype (MP) Top-Level")%>%
  + scale_y_continuous(breaks = seq(0, 200, by = 10))%>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10, angle = 45, hjust = 1, vjust = 1)) %>%
  +geom_line()%>%
  +geom_point(size = 4)%>%
  + scale_color_manual(values =number, aesthetics = "colour")%>%
  +guides(color=guide_legend(title="Top-Level MP-Term"))%>%
  +  xlim(rev(levels(sort_MEAN.Pheno.Mouse.Development$`Development Stage`)))
  #+coord_flip()

# +facet_wrap( . ~ Pheno.Mouse.Development$mp.description)


MEAN.Pheno.Mouse.Development_P2

# 
# # 
# library(cowplot)
# # 
# save_plot(paste("./Plots/Mouse/EA_Developmental_Stages/Mean gene expression/MEAN_PHENO_EXPRESSION2.png"),
#           MEAN.Pheno.Mouse.Development_P2 ,base_height= 10 ,base_aspect_ratio = 2) 
# # 
# 














############################################################################################################################

