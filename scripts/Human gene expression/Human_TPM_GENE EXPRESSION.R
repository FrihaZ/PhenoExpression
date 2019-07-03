
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


################################################################################################################

Human_Genes_cleaned<-Human_Genes%>%
  select_all()%>%
  gather("GTEX.tissues","Expression", Adipose...Subcutaneous:Whole.Blood)%>%
  select(HGNC.ID,GTEX.tissues,Expression)%>%
  arrange(desc(Expression))%>%
  mutate(Log.10= log10(Expression))%>%
  distinct()

View(Human_Genes_cleaned)

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
# library(cowplot)
# 
# 
# save_plot("./Plots/Human/TPM/Human_Genes_cleaned_p.PDF",
#           Human_Genes_cleaned_p ,base_height= 7 ,base_aspect_ratio = 3) 
# 

