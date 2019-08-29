
###############################################################################################################################
###############################################################################################################################

### Project: PhenoExpression ##################################################################################################
### Script: clean_df_FUNCTION.R ###################################################################################################
### Purpose: To clean the human gene expression (GTEX) dataset  ########################################################
### Author: Friha Zafar ####################################################################################################
### Date: 18/06/2019 ##########################################################################################################

###############################################################################################################################
###############################################################################################################################

## load packages ##############################################################################################################

library(dplyr);library(tidyr);library(stringr);library(readr);library(tidyverse);library(ggpubr)

################################################################################################################################
################################################################################################################################

# Load Files

GTEX.Tissues.to.HPO<-read_csv("D:/MSC RESEARCH PROJECT/Output_Files/GTEX.Tissues.to.HPO.csv")

hp.mp.toplevelmapping<-read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/hp-mp-toplevelmapping.txt",delim = "\t")

mgi.phenotypes<-read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/mgi.phenotypes.txt",delim = "\t")
############################################################################################################################################################
################################################################################################################################################

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~Clean the Human data to have only the genes that are orthologs of mice ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###########################

## Function to change the df into long format

clean_df<-function(df){
  
  pheno_cleaned <- df %>% 
    gather("GTEX.tissues","Expression", `Adipose - Subcutaneous`:`Whole Blood`)%>%
    select(HGNC.ID,hpo.ancestors,hpo.description,GTEX.tissues,Expression)%>%
    distinct()
  
  # merge the tissues with their hpo.descp using the gtex file 
  
  
  # USE THE GTEX TISSUES TO GET THE HPO.DESCP-TISSUE RELATIONSHIP CORRECTED
  
  GTEX.Tissues.to.tissues<-GTEX.Tissues.to.HPO%>%
    select_all()%>%
    drop_na()%>%
    select(GTEX.tissues)
  
  ## There are some tissues that are found in different HPO.superclass.description tehrefore have been duplicated in the GTEX.Tissues.to.tissues
  
  # Merge the GTEX.Tissues.to.tissues and human_genes_TPM_greater0_pheno_count to get the duplicated tissues
  
  pheno_cleaned_Tissues<-merge(GTEX.Tissues.to.tissues, pheno_cleaned, 
                               by = c("GTEX.tissues", "GTEX.tissues"), all.x = TRUE)
  
  #View(pheno_cleaned_Tissues)
  
  
  # The human_genes_TPM_greater0_pheno_count and GTEX.Tissues.to.HPO joind to get the HPO.superclass.description for each tissue
  
  
  ## HERE I KEPT THE TWO COLUMNS WITH THE HPO.SUPERCLASS THAT IS LINKED TO THE GTEX TISSUE, AND THE HPO.DESC WHICH IS THE PHENOTYPE THE GENE IS ASSOCIATED WITH
  
  
  pheno_cleaned_Tissues<-inner_join(GTEX.Tissues.to.HPO,pheno_cleaned_Tissues)
  pheno_cleaned_Tissues<-pheno_cleaned_Tissues%>%
    select_all()%>%
    drop_na()%>%    # the NAs are dropped
    distinct(GTEX.tissues,hpo.description,HGNC.ID,HPO.superclass.description,Expression ) # duplicates are removed that are not needed
  
  
  
  # Genes that are Linked to HPO.Description
  
  linked.hpo.desc<-pheno_cleaned_Tissues%>%
    select(hpo.description,HGNC.ID)%>%
    distinct()
  
  
  
  #clean the pheno_cleaned_Tissues by removing the hpo.description that are linked to the gene
  
  pheno_cleaned_Tissues_<-pheno_cleaned_Tissues%>%
    select(GTEX.tissues,HGNC.ID,HPO.superclass.description,Expression)%>%
    distinct()
  return(pheno_cleaned_Tissues_)
  
}




##############################################################################################################################
############################################################################################################################

## Use function to clean the dfs
