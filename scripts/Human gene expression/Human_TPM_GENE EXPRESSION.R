
###############################################################################################################################
###############################################################################################################################

### Project: PhenoExpression ##################################################################################################
### Script: Human_TPM_GENE EXPRESSION.R ###################################################################################################
### Purpose: To visualise the TPM data to get an overview ###############################################################################
### Author: Friha Zafar ####################################################################################################
### Date: 01/07/2019 ##########################################################################################################

###############################################################################################################################
###############################################################################################################################

## load packages ##############################################################################################################

library(dplyr);library(tidyr);library(stringr);library(readr);library(tidyverse);library(ggpubr)

################################################################################################################################
################################################################################################################################

## import files ################################################################################################################


Human_Genes<-read.csv('D:/MSC RESEARCH PROJECT//Output_Files/Human/human_genes.csv')

View(Human_Genes )
################################################################################################################

Human_Genes_cleaned<-Human_Genes%>%
  select_all()%>%
  gather("GTEX.tissues","Expression", Adipose...Subcutaneous:Whole.Blood)%>%
  select(HGNC.ID,GTEX.tissues,Expression)%>%
  arrange(desc(Expression))%>%
  distinct()
View(Human_Genes_cleaned)


Human_Genes_cleaned$GTEX.tissues<-gsub("...", "-", Human_Genes_cleaned$GTEX.tissues, fixed=TRUE)



library(gtools)
# Human_Genes_cleaned$GTEX.tissues<-mixedsort(Human_Genes_cleaned$GTEX.tissues)


Human_Genes_cleaned$GTEX.tissues<- factor(Human_Genes_cleaned$GTEX.tissues,
                                                                levels=rev(unique(Human_Genes_cleaned$GTEX.tissues)))


View(Human_Genes_cleaned)


#################################################################################################
######################################################################################################

# Gene expression for each Tissue
#########################################################################################
## PLOT ALL


library(ggplot2)
Human_Genes_cleaned_Plot_ALL<-ggplot(Human_Genes_cleaned, 
                                            aes(x =Human_Genes_cleaned$GTEX.tissues,
                                                y =Human_Genes_cleaned$Expression, 
                                                fill=Human_Genes_cleaned$GTEX.tissues))%>%
  +labs(x= "GTEX Tissue", y="Gene Expression (TPM)", 
        subtitle= "Human Gene Expression in Each Tissue")%>%
  +theme_classic()%>%
  +scale_y_continuous(breaks = seq(0, 246700, by = 10000))%>%
  +theme(axis.title=element_text(size=10,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),# linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10, angle = 45, hjust = 1, vjust = 1),
         legend.position = "none") %>%
  # +geom_point(shape = 1, color="red") %>%
  +stat_boxplot(geom = "errorbar")%>%
  +geom_boxplot(position=position_dodge(1),horizontal=TRUE,
                axes=FALSE,
                outlier.color = "red",
                outlier.fill = NA,
                outlier.shape = 2)%>%
  +coord_cartesian(ylim = c(0, 246700),expand = TRUE)
  # +coord_flip()

#Human_Genes_cleaned_Plot_ALL



##!!!UNHASH TO SAVE PLOT!!! 
# 
# library(cowplot)
# # 
# save_plot("./Plots/Human/TPM/Human_Genes_cleaned_Plot_ALL.jpeg",
#           Human_Genes_cleaned_Plot_ALL ,base_height= 8 ,base_aspect_ratio = 3) 


#############################################################################################################################
###########################################################################################################################

## Enlarged plot

library(ggplot2)
Human_Genes_cleaned_Plot_ALL_zoom<-ggplot(Human_Genes_cleaned, 
                                          aes(x =Human_Genes_cleaned$GTEX.tissues,
                                              y =Human_Genes_cleaned$Expression, 
                                         fill=Human_Genes_cleaned$GTEX.tissues))%>%
  +labs(x= "GTEX Tissue", y="Gene Expression (TPM)", 
        subtitle= "Human Gene Expression in Each Tissue")%>%
  +theme_classic()%>%
  +scale_y_continuous(breaks = seq(0, 246700 , by = 100))%>%
  +theme(axis.title=element_text(size=10,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),# linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10, angle = 45, hjust = 1, vjust = 1),
         legend.position = "none") %>%
  # +geom_point(shape = 1, color="red") %>%
  +stat_boxplot(geom = "errorbar")%>%
  +geom_boxplot(position=position_dodge(1),horizontal=TRUE,
                axes=FALSE,
                outlier.color = "red",
                outlier.fill = NA,
                outlier.shape = NA)%>%
  +coord_cartesian(ylim = c(0, 100),expand = TRUE)
# +coord_flip()

#Human_Genes_cleaned_Plot_ALL_zoom



##!!!UNHASH TO SAVE PLOT!!! 

# library(cowplot) 
# save_plot("./Plots/Human/TPM/Human_Genes_cleaned_Plot_ALL_zoom.jpeg",
#           Human_Genes_cleaned_Plot_ALL_zoom ,base_height= 8 ,base_aspect_ratio = 3) 

####################################################################################################
###########################################################################################################


##REMOVE OUTLIERS 


Human_Genes_cleaned_no_outliers<-Human_Genes_cleaned
outlier<-boxplot.stats(Human_Genes_cleaned_no_outliers$Expression)$out


Human_Genes_cleaned_no_outliers <- Human_Genes_cleaned_no_outliers[!Human_Genes_cleaned_no_outliers$Expression %in% outlier, ]

View(Human_Genes_cleaned_no_outliers)





library(ggplot2)
Human_Genes_cleaned_Plot_no_ouliers<-ggplot(Human_Genes_cleaned_no_outliers, 
                                 aes(x =Human_Genes_cleaned_no_outliers$GTEX.tissues, 
                                     y =Human_Genes_cleaned_no_outliers$Expression, 
                                     fill=Human_Genes_cleaned_no_outliers$GTEX.tissues))%>%
  +labs(x= "GTEX Tissue", y="Gene Expression (TPM)", 
        subtitle= "Human Gene Expression in Each Tissue")%>%
  +theme_classic()%>%
  +scale_y_continuous(breaks = seq(0, 60, by = 2))%>%
  +theme(axis.title=element_text(size=10,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),# linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10),
         legend.position = "none") %>%
  +geom_point(shape = 1, color="black") %>%
  +stat_boxplot(geom = "errorbar")%>%
  +geom_boxplot(outlier.color = "black",
                outlier.fill = "red",
                outlier.shape = 5) 
#+facet_wrap( . ~ Pheno.Mouse.Developments$mp.description)

#Human_Genes_cleaned_Plot_no_ouliers
 


##!!!UNHASH TO SAVE PLOT!!! 

# library(cowplot)
# #     
# save_plot("./Plots/Human/TPM/Human_Genes_cleaned_Plot_no_ouliers.jpeg",
#      Human_Genes_cleaned_Plot_no_ouliers ,base_height= 8 ,base_aspect_ratio = 3) 

#########################################################################################################
###########################################################################################################



