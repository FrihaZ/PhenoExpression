
###############################################################################################################################
###############################################################################################################################

### Project: PhenoExpression ##################################################################################################
### Script: Human_Gene_Venn_Diagrams.R ###################################################################################################
### Purpose: To create venn diagrams to show statisticallt significant genes taht may overlap betwene thresholds ###############################################################################
### Author: Friha Zafar ####################################################################################################
### Date: 18/06/2019 ##########################################################################################################

###############################################################################################################################
###############################################################################################################################

## load packages ##############################################################################################################

library(dplyr);library(tidyr);library(stringr);library(readr);library(tidyverse);library(ggpubr)

################################################################################################################################
################################################################################################################################

## import files ################################################################################################################
PVALUE_GENE_Greater0<-read_csv("D:/MSC RESEARCH PROJECT/Output_Files/Human/Gene Pvalues/PVALUE_GENE_Greater0.csv")

PVALUE_GENE_0.1<-read_csv("D:/MSC RESEARCH PROJECT/Output_Files/Human/Gene Pvalues/PVALUE_GENE_0.1.csv")

PVALUE_GENE_1<-read_csv("D:/MSC RESEARCH PROJECT/Output_Files/Human/Gene Pvalues/PVALUE_GENE_1.csv")

############################################################################################################################
##################################################################################################################################
## CLEAN DATA

## Gene Expression >0 

PVALUE_GENE_Greater0<-PVALUE_GENE_Greater0%>%
  select(HGNC.ID, P.VALUE)%>%
  distinct()
View(PVALUE_GENE_Greater0)



## Gene Expression >0.1

PVALUE_GENE_0.1<-PVALUE_GENE_0.1%>%
  select(HGNC.ID, P.VALUE)%>%
  distinct()
#View(PVALUE_GENE_0.1)



## Gene Expression >1

PVALUE_GENE_1<-PVALUE_GENE_1%>%
  select(HGNC.ID, P.VALUE)%>%
  distinct()
View(PVALUE_GENE_1)

####################################################################################################################
#########################################################################################################################

## >0 TPM 

PVALUE_GENE_Greater0$Adjusted.Pvalue<-p.adjust(PVALUE_GENE_Greater0$P.VALUE, method="BH")

View(PVALUE_GENE_Greater0)


## >0.1 TPM 

PVALUE_GENE_0.1$Adjusted.Pvalue<-p.adjust(PVALUE_GENE_0.1$P.VALUE, method="BH")

View(PVALUE_GENE_0.1)


## >1 TPM
PVALUE_GENE_1$Adjusted.Pvalue<-p.adjust(PVALUE_GENE_1$P.VALUE, method="BH")

View(PVALUE_GENE_1)


##############################################################################################
##############################################################################################

#Filter for significant genes (<0.05)

## >0

ADJ.SIGN.PVALUE.Greater0<- PVALUE_GENE_Greater0 %>%
  select(HGNC.ID, Adjusted.Pvalue)%>%
  filter(Adjusted.Pvalue <= 0.05)%>%
  distinct(HGNC.ID, Adjusted.Pvalue)

View(ADJ.SIGN.PVALUE.Greater0)

#write.csv(ADJ.SIGN.PVALUE.Greater0,'./Output_Files/Human/Adjusted_Significant_PValues/ADJ.SIGN.PVALUE.Greater0.csv')



## >0.1 TPM

ADJ.SIGN.PVALUE_0.1<- PVALUE_GENE_0.1 %>%
  select(HGNC.ID, Adjusted.Pvalue)%>%
  filter(Adjusted.Pvalue <= 0.05)%>%
  distinct(HGNC.ID, Adjusted.Pvalue)

View(ADJ.SIGN.PVALUE_0.1)

#write.csv(ADJ.SIGN.PVALUE_0.1,'./Output_Files/Human/Adjusted_Significant_PValues/ADJ.SIGN.PVALUE_0.1.csv')



## >0.1 TPM

ADJ.SIGN.PVALUE_1<- PVALUE_GENE_1 %>%
  select(HGNC.ID, Adjusted.Pvalue)%>%
  filter(Adjusted.Pvalue <= 0.05)%>%
  distinct(HGNC.ID, Adjusted.Pvalue)

View(ADJ.SIGN.PVALUE_1)

#write.csv(ADJ.SIGN.PVALUE_1,'./Output_Files/Human/Adjusted_Significant_PValues/ADJ.SIGN.PVALUE_1.csv')






######################################################################################################################
###################################################################################################################################

## CREATE DATAFRAME THAT HAS NAMES OF THE GENES FOR EACH THRESHOLD


## make variables for each of the thresholds

SIGN.PVALUE.Greater0<- ADJ.SIGN.PVALUE.Greater0$HGNC.ID
SIGN.PVALUE_0.1<- ADJ.SIGN.PVALUE_0.1$HGNC.ID
SIGN.PVALUE_1<- ADJ.SIGN.PVALUE_1$HGNC.ID

library(qpcR)

All.Sign.Genes<- qpcR:::cbind.na(SIGN.PVALUE.Greater0,
                          SIGN.PVALUE_0.1,
                          SIGN.PVALUE_1)

View(All.Sign.Genes)


#write.csv(All.Sign.Genes,'./Output_Files/Human/Adjusted_Significant_PValues/ALL.ADJ.SIGN.Genes.csv')




library(VennDiagram)
library(grid)


vp <- venn.diagram(list("TPM >/= 0.1 TPM"= SIGN.PVALUE_0.1, 
                       "TPM >/= 1 TPM"=SIGN.PVALUE_1),
                   fill = c("seagreen3", "orchid3"),
                   alpha = 0.5, filename = NULL, cex = 1, 
                   counts.col = "red", euler= TRUE,cat.cex = 1,
                   scaled= FALSE, margin = 0.05, cat.pos = c(-0.2, 0.2));
grid.draw(vp)

library(grDevices)

#save plot

pdf(file="./Plots/Human/Statistical_Sign_Genes/Venn_Diagram.pdf")
grid.draw(vp)
dev.off()



