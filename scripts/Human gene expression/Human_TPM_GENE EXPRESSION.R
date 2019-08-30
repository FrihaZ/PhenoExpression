
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

View(Human_Genes)
################################################################################################################

Human_Genes_cleaned<-Human_Genes%>%
  select_all()%>%
  gather("GTEX.tissues","Expression", Adipose...Subcutaneous:Whole.Blood)%>%
  select(HGNC.ID,GTEX.tissues,Expression)%>%
  arrange(desc(Expression))%>%
  mutate(Log.10= log10(Expression))%>%
  distinct()



Human_Genes_cleaned$GTEX.tissues<-gsub("...", "-", Human_Genes_cleaned$GTEX.tissues, fixed=TRUE)



library(gtools)
Human_Genes_cleaned$GTEX.tissues<-mixedsort(Human_Genes_cleaned$GTEX.tissues)


Human_Genes_cleaned$GTEX.tissues<- factor(Human_Genes_cleaned$GTEX.tissues,
                                                                levels=rev(unique(Human_Genes_cleaned$GTEX.tissues)))


View(Human_Genes_cleaned)

filter1_Human_Genes_cleaned<-Human_Genes_cleaned[ 1:18895,]   

filter2_Human_Genes_cleaned<-Human_Genes_cleaned[ 18896:661299,]


filter3_Human_Genes_cleaned<-Human_Genes_cleaned[ 661300:906932,]

filter4_Human_Genes_cleaned<-Human_Genes_cleaned[ 906933:1001405,]


View(filter2_Human_Genes_cleaned)

##########################################################################################

Human_Genes_cleaned_p<-ggplot(Human_Genes_cleaned, aes(x=Log.10))%>%
  +geom_histogram(aes(y=..density..), colour="black", fill="white")%>%
  +geom_density(alpha=.2, fill="#FF6666") %>%
  + labs(x= "Gene Expression (Log10(TPM))", y="Density", 
         subtitle="The Density of Genes with Similar Gene Expression (Log10(TPM))")%>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=8),
         axis.text.x= element_text(size=10),
         panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  + scale_x_continuous(breaks = seq(-10, 10, by = 0.25))%>%
  + scale_y_continuous(breaks = seq(0, 1, by = 0.1))
  
Human_Genes_cleaned_p
# 

hist(Human_Genes_cleaned$Expression, breaks = 5)

# library(cowplot)
# 
# 
# save_plot("./Plots/Human/TPM/Human_Genes_cleaned_p.PDF",
#           Human_Genes_cleaned_p ,base_height= 7 ,base_aspect_ratio = 3) 
# 

# # 
# save_plot("./Plots/Human/TPM/Human_Genes_cleaned_p.jpeg",
#            Human_Genes_cleaned_p ,base_height= 7 ,base_aspect_ratio = 3) 
#  
# 


#################################################################################################
######################################################################################################

# Gene expression for each Tissue

#filter1_Human_Genes_cleaned

library(ggplot2)
Human_Genes_cleaned_Plot<-ggplot(filter1_Human_Genes_cleaned, 
                                    aes(x =filter1_Human_Genes_cleaned$GTEX.tissues, 
                                        y =filter1_Human_Genes_cleaned$Expression, 
                                        fill=filter1_Human_Genes_cleaned$GTEX.tissues))%>%
  +labs(x= "GTEX Tissue", y="Gene Expression (TPM)", 
        subtitle= "Human Gene Expression in Each Tissue")%>%
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
                outlier.shape = 5) %>%
  +coord_flip()
#+facet_wrap( . ~ Pheno.Mouse.Developments$mp.description)

Human_Genes_cleaned_Plot


library(cowplot)
#    
save_plot("./Plots/Human/TPM/Human_Genes_cleaned_Plot1.jpeg",
          Human_Genes_cleaned_Plot ,base_height= 8 ,base_aspect_ratio = 3) 
 




#######################################################################################
##############################################################################


library(ggplot2)
Human_Genes_cleaned_Plot2<-ggplot(filter2_Human_Genes_cleaned, 
                                 aes(x =filter2_Human_Genes_cleaned$GTEX.tissues, 
                                     y =filter2_Human_Genes_cleaned$Expression, 
                                     fill=filter2_Human_Genes_cleaned$GTEX.tissues))%>%
  +labs(x= "GTEX Tissue", y="Gene Expression (TPM)", 
        subtitle= "Human Gene Expression in Each Tissue")%>%
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
                outlier.shape = 5) %>%
  +coord_flip()
#+facet_wrap( . ~ Pheno.Mouse.Developments$mp.description)

Human_Genes_cleaned_Plot2


library(cowplot)
#    
save_plot("./Plots/Human/TPM/Human_Genes_cleaned_Plot2.jpeg",
          Human_Genes_cleaned_Plot2 ,base_height= 8 ,base_aspect_ratio = 3) 

#####################################################################################################


#filter3_Human_Genes_cleaned

library(ggplot2)
Human_Genes_cleaned_Plot3<-ggplot(filter3_Human_Genes_cleaned, 
                                 aes(x =filter3_Human_Genes_cleaned$GTEX.tissues, 
                                     y =filter3_Human_Genes_cleaned$Expression, 
                                     fill=filter3_Human_Genes_cleaned$GTEX.tissues))%>%
  +labs(x= "GTEX Tissue", y="Gene Expression (TPM)", 
        subtitle= "Human Gene Expression in Each Tissue")%>%
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
                outlier.shape = 5) %>%
  +coord_flip()
#+facet_wrap( . ~ Pheno.Mouse.Developments$mp.description)

Human_Genes_cleaned_Plot3


library(cowplot)
#    
save_plot("./Plots/Human/TPM/Human_Genes_cleaned_Plot3.jpeg",
          Human_Genes_cleaned_Plot3 ,base_height= 8 ,base_aspect_ratio = 3) 


########################################################################################


#4



#filter4_Human_Genes_cleaned

library(ggplot2)
Human_Genes_cleaned_Plot4<-ggplot(filter4_Human_Genes_cleaned, 
                                 aes(x =filter4_Human_Genes_cleaned$GTEX.tissues, 
                                     y =filter4_Human_Genes_cleaned$Expression, 
                                     fill=filter4_Human_Genes_cleaned$GTEX.tissues))%>%
  +labs(x= "GTEX Tissue", y="Gene Expression (TPM)", 
        subtitle= "Human Gene Expression in Each Tissue")%>%
  +theme(axis.title=element_text(size=10,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),# linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=5),
         legend.position = "none") %>%
  +geom_point(shape = 1, color="black") %>%
  +stat_boxplot(geom = "errorbar")%>%
  +geom_boxplot(outlier.color = "black",
                outlier.fill = "red",
                outlier.shape = 5) %>%
  + scale_y_continuous(breaks = seq(0, 1, by = 0.001))%>%

  +coord_flip()
#+facet_wrap( . ~ Pheno.Mouse.Developments$mp.description)

Human_Genes_cleaned_Plot4


library(cowplot)
#    
save_plot("./Plots/Human/TPM/Human_Genes_cleaned_Plot4.jpeg",
          Human_Genes_cleaned_Plot4 ,base_height= 3 ,base_aspect_ratio = 3) 





