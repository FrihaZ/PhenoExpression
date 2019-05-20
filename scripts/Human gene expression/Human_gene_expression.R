
################################################################################################################################
###############################################################################################################################

### Project: PhenoExpression ##################################################################################################
### Script: Human_gene_expression.R ###################################################################################################
### Purpose: To clean the Human gene expression data and compare to mice data ###############################################################################
### Author: Friha Zafar ####################################################################################################
### Date: 01/04/2019 ##########################################################################################################

###############################################################################################################################
###############################################################################################################################

## load packages ##############################################################################################################

library(dplyr);library(tidyr);library(stringr);library(readr);library(tidyverse);library(ggpubr)

################################################################################################################################
################################################################################################################################

## import files ################################################################################################################

human_genes <- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/GTEX.TPM.gz",delim = "\t", skip = 2)

gene_protein <- read_delim("D:/MSC RESEARCH PROJECT/Gene-symbol-checker/gene_with_protein_product.txt",delim = "\t")
source(file="D:/MSC RESEARCH PROJECT/Gene-symbol-checker/hgnc_symbol_checker.R")


hpo.ancestor.nodes <- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/hpo.ancestor.nodes.txt",delim = "\t")

hpo.toplevels.phenotypic.abnormalities.only <- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/hpo.toplevels.phenotypic.abnormalities.only.txt",delim = "\t")

#!!!!!!!!!MANUALLY REMOVE THE <TAB> AND INSERT REAL TABS IN THE HPO_phenotypes.txt FILE'S COLUMN NAMES!!!

HPO_phenotypes <- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/HPO_phenotypes.txt",delim = "\t")


GTEX.Tissues.to.HPO <- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/GTEX.Tissues.to.HPO.txt",delim = "\t")



###############################################################################################################################################

human_genes_filtered <-human_genes %>%
  select(Description)%>% # only the gene names are selected as these will be taken to find the hgnc.id using the hgnc.checker. 
  distinct(Description) %>%
  pull(Description)           # to make it into a vector for the hgnc.checker to work


##View(human_genes_filtered)   




################################################################################################################################
################################################################################################################################

#gene_protein <- read_delim("D:/MSC RESEARCH PROJECT/Gene-symbol-checker/gene_with_protein_product.txt",delim = "\t")


##View(gene_protein)


# Retrieve hgnc.checker function


#source(file="D:/MSC RESEARCH PROJECT/Gene-symbol-checker/hgnc_symbol_checker.R")



approved_hgnc.id <- hgnc.checker(human_genes_filtered, gene_protein)

##View(approved_hgnc.id)



################################################################################################################################
################################################################################################################################

# Remove rows that are empty: filter by  '-' and 'Ambiguous.Symbol'

approved_hgnc.id <- approved_hgnc.id %>% 
  filter(HGNC.ID  != "-") %>%  #Remove ambiguous expression data
  filter(Type != "Ambiguous.Symbol") 

##View(approved_hgnc.id)


################################################################################################################################
################################################################################################################################

# Convert gene expression values (TPM) to binary variable. Three different thresholds:
# TPM > 0 (yes,no)
# TPM > 0.1 (yes, no)
# TPM> 1 (yes, no)



## Join  the approved_hgnc.id with the human_genes table to get the gene expressions 

#str(approved_hgnc.id)

human_genes <-left_join(human_genes, approved_hgnc.id, by=c("Description" = "Gene.Symbol"))  

human_genes <-human_genes%>%
  select(HGNC.ID,gene_id, Description, Type,`Adipose - Subcutaneous`:`Whole Blood`)%>% 
  drop_na() %>%  #removed NAs that for the genes that didn't have a hgncd.id 
  select(-gene_id)%>%   #removed the gene_id 
  distinct()

options("scipen"=100)

#View(human_genes)



# Convert gene expression values (TPM) to binary variable. Three different thresholds:



#human_genes<-format(human_genes, scientific=F)


human_genes_TPM_None<- human_genes

human_genes_TPM_None[4:56 ][ human_genes_TPM_None[4:56 ] > 0 ] <- "No"

human_genes_TPM_None[4:56 ][ human_genes_TPM_None[4:56 ] <= 0 ] <- "Yes"

#View(human_genes_TPM_None)

#ACTRT2


# TPM > 0 (yes,no)
##View(human_genes)

human_genes_TPM_greater0<- human_genes


human_genes_TPM_greater0[4:56 ][ human_genes_TPM_greater0[4:56 ] > 0 ] <- "Yes"

human_genes_TPM_greater0[4:56 ][ human_genes_TPM_greater0[4:56 ] <= 0 ] <- "No"

##View(human_genes_TPM_greater0)


# TPM > 0.1 (yes, no)

human_genes_TPM_0.1<- human_genes

human_genes_TPM_0.1[4:56 ][ human_genes_TPM_0.1[4:56 ] >= 0.1 ] <- "Yes"
human_genes_TPM_0.1[4:56 ][ human_genes_TPM_0.1[4:56 ] < 0.1 ] <- "No"




#View(human_genes_TPM_0.1)

# TPM> 1 (yes, no)

human_genes_TPM_1<- human_genes

human_genes_TPM_1[4:56 ][ human_genes_TPM_1[4:56 ] >= 1 ] <- "Yes"
human_genes_TPM_1[4:56 ][ human_genes_TPM_1[4:56 ] < 1 ] <- "No"



#View(human_genes_TPM_1)





################################################################################################################################
################################################################################################################################


# Convert gene symbols in HPO_phenotypes to hgnc.id


###MANUALLY REMOVE THE <TAB> AND INSERT REAL TABS IN THE HPO_phenotypes.txt FILE'S COLUMN NAMES!!!

#HPO_phenotypes <- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/HPO_phenotypes.txt",delim = "\t")

# The dataframe is turned into a vector and only the `entrez-gene-symbol`is used 
HPO_phenotypes_filtered <-HPO_phenotypes %>%
  select(`entrez-gene-symbol`)%>%
  distinct(`entrez-gene-symbol`) %>%
  pull(`entrez-gene-symbol`)           # to make it into a vector


#View(HPO_phenotypes_filtered)  



# Use the hgnc.checker function

HPO_phenotypes_hgnc.id <- hgnc.checker(HPO_phenotypes_filtered, gene_protein)



#View(HPO_phenotypes_hgnc.id)



################################################################################################################################
################################################################################################################################

# Remove rows that are empty: filter by  '-' and 'Notfound.ProteinCoding.Symbol'

HPO_phenotypes_hgnc.id <- HPO_phenotypes_hgnc.id %>% 
  filter(HGNC.ID  != "-") %>%  #Remove ambiguous expression data
  filter(Type != "Notfound.ProteinCoding.Symbol" | Type != "Ambiguous.Symbol") 

#View(HPO_phenotypes_hgnc.id)


################################################################################################################################
################################################################################################################################


## Join  the HPO_phenotypes_hgnc.id with the HPO_phenotypes table 


HPO_phenotypes <-left_join(HPO_phenotypes, HPO_phenotypes_hgnc.id, by=c(`entrez-gene-symbol` = "Gene.Symbol"))  

#View(HPO_phenotypes)


################################################################################################################################
################################################################################################################################


# Map the HPO-Term-ID to the top level nodes using hpo.ancestor.nodes.txt file

#hpo.ancestor.nodes <- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/hpo.ancestor.nodes.txt",delim = "\t")

##View(hpo.ancestor.nodes)
#hpo.toplevels.phenotypic.abnormalities.only <- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/hpo.toplevels.phenotypic.abnormalities.only.txt",delim = "\t")

##View(hpo.toplevels.phenotypic.abnormalities.only)


hpo.ancestor.nodes <- left_join(hpo.ancestor.nodes, hpo.toplevels.phenotypic.abnormalities.only, by= c("hpo.ancestors"="hpo.term"))


##View(hpo.ancestor.nodes)
hpo.ancestor.nodes<-hpo.ancestor.nodes%>%
  select(hpo.term, hpo.ancestors, hpo.description)%>%
  drop_na()%>%
  distinct()


#View(hpo.ancestor.nodes)

# Join the two tables 


HPO_phenotypes <-left_join(HPO_phenotypes, hpo.ancestor.nodes, by=c("HPO-Term-ID" = "hpo.term"))  



#REMOVED NAs
HPO_phenotypes<-HPO_phenotypes %>%
  select_all()%>%
  drop_na()%>%       
  distinct()

#View(HPO_phenotypes)

#####################################################################################################################

# Map genes with gene expression to the HPO annotations (top level)


#####################################################################################################################


## To all disease associated genes, irrespective of threshold

all_disease_genes_phenotype <-left_join(human_genes, HPO_phenotypes, by=c("HGNC.ID"= "HGNC.ID"))



all_disease_genes_pheno <-all_disease_genes_phenotype%>%
  select_all()%>%
  drop_na() %>%  #removed NAs that for the genes that didn't have a hgncd.id 
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`)%>%   #removed the entrez-gene-id etc
  distinct()


#View(all_disease_genes_pheno)




#################################################################


# Remove all genes that are disease related, and keep all those that are not disease related. 

all_genes_pheno_No_disease<-all_disease_genes_phenotype[!(all_disease_genes_phenotype$HGNC.ID %in% all_disease_genes_pheno$HGNC.ID),]

all_genes_pheno_No_disease<-all_genes_pheno_No_disease%>%
  select_all()%>%
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`,
         -`HPO-Term-ID`,-`hpo.ancestors`, -`hpo.description` )%>%   #removed the entrez-gene-id etc
  distinct()

#View(all_genes_pheno_No_disease)




######################################################################################################################
##########################################################################################################



#######################################################################################
## NO EXPRESSION

# DISEASE RELATED GENES

## E-P+

human_genes_TPM_None_phenotype<- left_join(human_genes_TPM_None, HPO_phenotypes, by=c("HGNC.ID"= "HGNC.ID"))



human_genes_TPM_None_pheno <-human_genes_TPM_None_phenotype%>%
  select_all()%>%
  drop_na() %>%  #removed NAs that for the genes that didn't have a hgncd.id 
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`)%>%   #removed the entrez-gene-id etc
  distinct()


#View(human_genes_TPM_None_pheno)




## E- P-
# Remove all genes that are disease related, and keep all those that are not disease related. 

human_genes_TPM_None_pheno_No_disease<-human_genes_TPM_None_phenotype[!(human_genes_TPM_None_phenotype$HGNC.ID %in% human_genes_TPM_None_pheno$HGNC.ID),]

human_genes_TPM_None_pheno_No_disease<-human_genes_TPM_None_pheno_No_disease%>%
  select_all()%>%
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`,
         -`HPO-Term-ID`,-`hpo.ancestors`, -`hpo.description` )%>%   #removed the entrez-gene-id etc
  distinct()

#View(human_genes_TPM_None_pheno_No_disease)









#human_genes_TPM_None


#######################################################################################
# DISEASE RELATED GENES



human_genes_TPM_greater0_phenotype<- left_join(human_genes_TPM_greater0, HPO_phenotypes, by=c("HGNC.ID"= "HGNC.ID"))



human_genes_TPM_greater0_pheno <-human_genes_TPM_greater0_phenotype%>%
  select_all()%>%
  drop_na() %>%  #removed NAs that for the genes that didn't have a hgncd.id 
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`)%>%   #removed the entrez-gene-id etc
  distinct()


#View(human_genes_TPM_greater0_pheno)





# Remove all genes that are disease related, and keep all those that are not disease related. 

human_genes_greater0_pheno_No_disease<-human_genes_TPM_greater0_phenotype[!(human_genes_TPM_greater0_phenotype$HGNC.ID %in% human_genes_TPM_greater0_pheno$HGNC.ID),]

human_genes_greater0_pheno_No_disease<-human_genes_greater0_pheno_No_disease%>%
  select_all()%>%
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`,
         -`HPO-Term-ID`,-`hpo.ancestors`, -`hpo.description` )%>%   #removed the entrez-gene-id etc
  distinct()

#View(human_genes_greater0_pheno_No_disease)


########################################################################################################
# Map genes with gene expression to the HPO annotations (top level)

human_genes_TPM_0.1_phenotype<- left_join(human_genes_TPM_0.1, HPO_phenotypes, by=c("HGNC.ID"= "HGNC.ID"))


# DISEASE RELATED GENES
human_genes_TPM_0.1_pheno <-human_genes_TPM_0.1_phenotype%>%
  select_all()%>%
  drop_na() %>%  #removed NAs that for the genes that didn't have a hgncd.id 
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`)%>%   #removed the entrez-gene-id etc
  distinct()


#View(human_genes_TPM_0.1_pheno)



# Remove all genes that are disease related, and keep all those that are not disease related. 

human_genes_TPM_0.1_pheno_No_disease<-human_genes_TPM_0.1_phenotype[!(human_genes_TPM_0.1_phenotype$HGNC.ID %in% human_genes_TPM_0.1_pheno$HGNC.ID),]

human_genes_TPM_0.1_pheno_No_disease<-human_genes_TPM_0.1_pheno_No_disease%>%
  select_all()%>%
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`,
         -`HPO-Term-ID`,-`hpo.ancestors`, -`hpo.description` )%>%   #removed the entrez-gene-id etc
  distinct()

#View(human_genes_TPM_0.1_pheno_No_disease)



#########################################################################################################

# DISEASE RELATED GENES

human_genes_TPM_1_phenotype<- left_join(human_genes_TPM_1, HPO_phenotypes, by=c("HGNC.ID"= "HGNC.ID"))

#View(human_genes_TPM_1_phenotype)

human_genes_TPM_1_pheno <-human_genes_TPM_1_phenotype%>%
  select_all()%>%
  drop_na() %>%  #removed NAs that for the genes that didn't have a hgncd.id 
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`)%>%   #removed the entrez-gene-id etc
  distinct()


#View(human_genes_TPM_1_pheno)

str(human_genes_TPM_1_pheno)


# Remove all genes that are disease related, and keep all those that are not disease related. 

human_genes_TPM_1_pheno_No_disease<-human_genes_TPM_1_phenotype[!(human_genes_TPM_1_phenotype$HGNC.ID %in% human_genes_TPM_1_pheno$HGNC.ID),]

human_genes_TPM_1_pheno_No_disease<-human_genes_TPM_1_pheno_No_disease%>%
  select_all()%>%
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`,
         -`HPO-Term-ID`,-`hpo.ancestors`, -`hpo.description` )%>%   #removed the entrez-gene-id etc
  distinct()

#View(human_genes_TPM_1_pheno_No_disease)




#########################################################################################################
#########################################################################################################

## ALL NON-DISEASE GENES!

# Remove the phenotype description ,hpo.ancestors HPO.Term and count the number of Yes for each row (gene) for 
# each threshold


#########################################################################################################
#########################################################################################################





# For expressions =0, calculate the number of Yes and No

human_genes_TPM_None_pheno_No_disease$NUMBER_OF_YES <- rowSums(human_genes_TPM_None_pheno_No_disease == "Yes")
human_genes_TPM_None_pheno_No_disease$NUMBER_OF_NO <- rowSums(human_genes_TPM_None_pheno_No_disease == "No")

human_genes_TPM_None_pheno_No_disease<-human_genes_TPM_None_pheno_No_disease %>%
  select_all()%>%
  distinct()

#View(human_genes_TPM_None_pheno_No_disease)




#########################################################################################################
#########################################################################################################

## COUNT THE NUMBER OF NON-DISEASE GENES EXPRESSED FOR EACH TISSUE



human_genes_TPM_None_tissue_ND<-human_genes_TPM_None_pheno_No_disease%>%
  select_all()%>%
  select(-NUMBER_OF_NO, -NUMBER_OF_YES, -`HGNC.ID`, -`Description`, -`Type.x`)

human_genes_TPM_None_tissue_ND <- colSums(human_genes_TPM_None_tissue_ND == "Yes")

human_genes_TPM_None_tissue_ND<-data.frame(human_genes_TPM_None_tissue_ND)

human_genes_TPM_None_tissue_ND<-setNames(cbind(rownames(human_genes_TPM_None_tissue_ND), human_genes_TPM_None_tissue_ND,
                                               row.names = NULL), 
                                             c("Tissue", "Number_of_Yes"))
#View(human_genes_TPM_None_tissue_ND)




#####~~~~~~!!!!!!!!!!!!!!IMPORTANT KEEP FOR DISSERTATION!!!!!!!!!!!~~~~~~##############

# PLOT THE NUMBER OF GENES THAT ARE EXPRESSED IN EACH TISSUE, where gene expression is >0 and not disease specific

library(ggplot2)


#dev.new(width=60, height=30) # to open in a new window

p_None_tissue_ND<-ggplot(data=human_genes_TPM_None_tissue_ND,
                             aes(x=reorder(human_genes_TPM_None_tissue_ND$Tissue,
                                           -human_genes_TPM_None_tissue_ND$Number_of_Yes),
                                 y=human_genes_TPM_None_tissue_ND$Number_of_Yes,
                                 fill=human_genes_TPM_None_tissue_ND$Number_of_Yes))

p_None_tissue_ND <-p_None_tissue_ND %>%
  
  + labs(x= "Tissue", y="Number of Genes",         
         subtitle= "Genes Not Associated with Disease")%>%
  
  + scale_fill_gradient(low = "blue", 
                        high = "red")%>%
  + scale_y_continuous(breaks = seq(0, 14000, by = 1000))%>%
  
  + guides(fill=guide_legend(title="Number of Genes Expressed:  ")) %>%
  + theme(axis.title=element_text(size=14,face="bold"),
          plot.caption=element_text(face = "italic", size=12, hjust = 0),
          text = element_text(size=10),
          legend.position = "bottom",legend.direction = "horizontal",
          legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
          axis.text.y= element_text(size=8),
          axis.text.x= element_text(size=10),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +geom_bar(stat = "identity")%>%    # to create a stacked barchart
  +geom_text(aes(label=paste0(round(human_genes_TPM_None_tissue_ND$Number_of_Yes,1))),size = 3, position = position_stack(vjust = 0.5), colour=c("white"), fontface='bold') %>%
  
  + coord_flip()%>%
  + ggsave(filename=paste("./Plots/Human/Tissue/","p_None_tissue_ND.png",sep=" "), limitsize = TRUE)


p_None_tissue_ND









#########################################################################################################


# For expressions >0, calculate the number of Yes and No

human_genes_greater0_pheno_No_disease$NUMBER_OF_YES <- rowSums(human_genes_greater0_pheno_No_disease == "Yes")
human_genes_greater0_pheno_No_disease$NUMBER_OF_NO <- rowSums(human_genes_greater0_pheno_No_disease == "No")

human_genes_greater0_pheno_No_disease<-human_genes_greater0_pheno_No_disease %>%
  select_all()%>%
  distinct()

#View(human_genes_greater0_pheno_No_disease)




#########################################################################################################
#########################################################################################################

## COUNT THE NUMBER OF NON-DISEASE GENES EXPRESSED FOR EACH TISSUE

human_genes_TPM_greater0_tissue_ND<-human_genes_greater0_pheno_No_disease%>%
  select_all()%>%
  select(-NUMBER_OF_NO, -NUMBER_OF_YES, -`HGNC.ID`, -`Description`, -`Type.x`)

human_genes_TPM_greater0_tissue_ND <- colSums(human_genes_TPM_greater0_tissue_ND == "Yes")

human_genes_TPM_greater0_tissue_ND<-data.frame(human_genes_TPM_greater0_tissue_ND)

human_genes_TPM_greater0_tissue_ND<-setNames(cbind(rownames(human_genes_TPM_greater0_tissue_ND), human_genes_TPM_greater0_tissue_ND, row.names = NULL), 
         c("Tissue", "Number_of_Yes"))
#View(human_genes_TPM_greater0_tissue_ND)



#####~~~~~~!!!!!!!!!!!!!!IMPORTANT KEEP FOR DISSERTATION!!!!!!!!!!!~~~~~~##############

# PLOT THE NUMBER OF GENES THAT ARE EXPRESSED IN EACH TISSUE, where gene expression is >0 and not disease specific

library(ggplot2)


#dev.new(width=60, height=30) # to open in a new window

p_greater0_tissue_ND<-ggplot(data=human_genes_TPM_greater0_tissue_ND,
                aes(x=reorder(human_genes_TPM_greater0_tissue_ND$Tissue,
                              -human_genes_TPM_greater0_tissue_ND$Number_of_Yes),
                y=human_genes_TPM_greater0_tissue_ND$Number_of_Yes,
                fill=human_genes_TPM_greater0_tissue_ND$Number_of_Yes))

p_greater0_tissue_ND <-p_greater0_tissue_ND %>%
  
  + labs(x= "Tissue", y="Number of Genes",         
         subtitle= "Genes Not Associated with Disease")%>%

  + scale_fill_gradient(low = "blue", 
                        high = "red")%>%
  + scale_y_continuous(breaks = seq(0, 14000, by = 1000))%>%
  
  + guides(fill=guide_legend(title="Number of Genes Expressed:  ")) %>%
  + theme(axis.title=element_text(size=14,face="bold"),
          plot.caption=element_text(face = "italic", size=12, hjust = 0),
          text = element_text(size=10),
          legend.position = "bottom",legend.direction = "horizontal",
          legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
          axis.text.y= element_text(size=8),
          axis.text.x= element_text(size=10),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +geom_bar(stat = "identity")%>%    # to create a stacked barchart
  +geom_text(aes(label=paste0(round(human_genes_TPM_greater0_tissue_ND$Number_of_Yes,1))),size = 3, position = position_stack(vjust = 0.5), colour=c("white"), fontface='bold') %>%
  
  + coord_flip()%>%
  + ggsave(filename=paste("./Plots/Human/Tissue/","p_greater0_tissue_ND.png",sep=" "), limitsize = TRUE)












#########################################################################################################

# For expression greater than 0.1



human_genes_TPM_0.1_pheno_No_disease$NUMBER_OF_YES <- rowSums(human_genes_TPM_0.1_pheno_No_disease == "Yes")
human_genes_TPM_0.1_pheno_No_disease$NUMBER_OF_NO <- rowSums(human_genes_TPM_0.1_pheno_No_disease == "No")

human_genes_TPM_0.1_pheno_No_disease<-human_genes_TPM_0.1_pheno_No_disease %>%
  select_all()%>%
  distinct()

#View(human_genes_TPM_0.1_pheno_No_disease)



## COUNT THE NUMBER OF GENES EXPRESSED FOR EACH TISSUE

human_genes_TPM_0.1_tissue_ND<-human_genes_TPM_0.1_pheno_No_disease%>%
  select_all()%>%
  select(-NUMBER_OF_NO, -NUMBER_OF_YES, -`HGNC.ID`, -`Description`, -`Type.x`)

human_genes_TPM_0.1_tissue_ND <- colSums(human_genes_TPM_0.1_tissue_ND == "Yes")


#View(human_genes_TPM_0.1_tissue_ND)


human_genes_TPM_0.1_tissue_ND<-data.frame(human_genes_TPM_0.1_tissue_ND)

human_genes_TPM_0.1_tissue_ND<-setNames(cbind(rownames(human_genes_TPM_0.1_tissue_ND), 
                                              human_genes_TPM_0.1_tissue_ND, row.names = NULL), 
                                          c("Tissue", "Number_of_Yes"))
#View(human_genes_TPM_0.1_tissue_ND)



#####~~~~~~!!!!!!!!!!!!!!IMPORTANT KEEP FOR DISSERTATION!!!!!!!!!!!~~~~~~##############

# PLOT THE NUMBER OF GENES THAT ARE EXPRESSED IN EACH TISSUE, where gene expression is >/=0.1 and not disease specific

library(ggplot2)


#dev.new(width=60, height=30) # to open in a new window

p_TPM_0.1_tissue_ND<-ggplot(data=human_genes_TPM_0.1_tissue_ND,
                           aes(x=reorder(human_genes_TPM_0.1_tissue_ND$Tissue, 
                                          -human_genes_TPM_0.1_tissue_ND$Number_of_Yes),
                                y=human_genes_TPM_0.1_tissue_ND$Number_of_Yes, 
                                fill=human_genes_TPM_0.1_tissue_ND$Number_of_Yes))
p_TPM_0.1_tissue_ND<-p_TPM_0.1_tissue_ND %>%
  
  
  + labs(x= "Tissue", y="Number of Genes", 
         subtitle= "Genes Not Associated with Disease")%>%
  + scale_fill_gradient(low = "blue", 
                        high = "red")%>% 
  + scale_y_continuous(breaks = seq(0, 14000, by = 1000))%>%
  
  + guides(fill=guide_legend(title="Number of Genes Expressed:  ")) %>%
  + theme(axis.title=element_text(size=14,face="bold"),
          plot.caption=element_text(face = "italic", size=12, hjust = 0),
          text = element_text(size=10),
          legend.position = "bottom",legend.direction = "horizontal",
          legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
          axis.text.y= element_text(size=8),
          axis.text.x= element_text(size=10),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +geom_bar(stat = "identity")%>%    # to create a stacked barchart
  +geom_text(aes(label=paste0(round(human_genes_TPM_0.1_tissue_ND$Number_of_Yes,1))),size = 3, position = position_stack(vjust = 0.5), colour=c("white"), fontface='bold') %>%
  
  + coord_flip()%>%
  + ggsave(filename=paste("./Plots/Human/Tissue/","p_TPM_0.1_tissue_ND.png",sep=" "), limitsize = TRUE)



#########################################################################################################

# For expression greater than 1



human_genes_TPM_1_pheno_No_disease$NUMBER_OF_YES <- rowSums(human_genes_TPM_1_pheno_No_disease == "Yes")
human_genes_TPM_1_pheno_No_disease$NUMBER_OF_NO <- rowSums(human_genes_TPM_1_pheno_No_disease == "No")

human_genes_TPM_1_pheno_No_disease<-human_genes_TPM_1_pheno_No_disease %>%
  select_all()%>%
  distinct()

#View(human_genes_TPM_1_pheno_No_disease)
## COUNT THE NUMBER OF GENES EXPRESSED FOR EACH TISSUE

human_genes_TPM_1_tissue_ND<-human_genes_TPM_1_pheno_No_disease%>%
  select_all()%>%
  select(-NUMBER_OF_NO, -NUMBER_OF_YES, -`HGNC.ID`, -`Description`, -`Type.x`)

human_genes_TPM_1_tissue_ND <- colSums(human_genes_TPM_1_tissue_ND == "Yes")


#View(human_genes_TPM_1_tissue_ND)



human_genes_TPM_1_tissue_ND<-data.frame(human_genes_TPM_1_tissue_ND)

human_genes_TPM_1_tissue_ND<-setNames(cbind(rownames(human_genes_TPM_1_tissue_ND), 
                                         human_genes_TPM_1_tissue_ND, row.names = NULL), 
                                     c("Tissue", "Number_of_Yes"))
#View(human_genes_TPM_1_tissue_ND)



#####~~~~~~!!!!!!!!!!!!!!IMPORTANT KEEP FOR DISSERTATION!!!!!!!!!!!~~~~~~##############

# PLOT THE NUMBER OF GENES THAT ARE EXPRESSED IN EACH TISSUE, where gene expression is >/=0.1 and not disease specific

library(ggplot2)


#dev.new(width=60, height=30) # to open in a new window

p_TPM_1_tissue_ND<-ggplot(data=human_genes_TPM_1_tissue_ND,
                          aes(x=reorder(human_genes_TPM_1_tissue_ND$Tissue, -human_genes_TPM_1_tissue_ND$Number_of_Yes),
                               y=human_genes_TPM_1_tissue_ND$Number_of_Yes,
                               fill=human_genes_TPM_1_tissue_ND$Number_of_Yes))
p_TPM_1_tissue_ND<-p_TPM_1_tissue_ND %>%
  
  + labs(x= "Tissue", y="Number ofGenes",
         subtitle= "Genes Not Associated with Disease")%>%
  
  + scale_fill_gradient(low = "blue", 
                        high = "red")%>% 
  + scale_y_continuous(breaks = seq(0, 14000, by = 1000))%>%
  
  + guides(fill=guide_legend(title="Number of Genes Expressed:  ")) %>%
  + theme(axis.title=element_text(size=14,face="bold"),
          plot.caption=element_text(face = "italic", size=12, hjust = 0),
          text = element_text(size=10),
          legend.position = "bottom",legend.direction = "horizontal",
          legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
          axis.text.y= element_text(size = 8),
          axis.text.x= element_text(size=10),
          
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +geom_bar(stat = "identity")%>%    # to create a stacked barchart
  +geom_text(aes(label=paste0(round(human_genes_TPM_1_tissue_ND$Number_of_Yes,1))),size = 3, position = position_stack(vjust = 0.5), colour=c("white"), fontface='bold') %>%
  + coord_flip()%>%
  + ggsave(filename=paste("./Plots/Human/Tissue/","p_TPM_1_tissue_ND.png",sep=" "), limitsize = TRUE)




#########################################################################################################
#########################################################################################################

## ALL DISEASE GENES!

# Remove the phenotype description ,hpo.ancestors HPO.Term and count the number of Yes for each row (gene) for 
# each threshold


#########################################################################################################



#Threshold =0 for all disease associated genes

human_genes_TPM_None_pheno_disease<-human_genes_TPM_None_pheno%>%
  select_all()%>%
  select(-`HPO-Term-ID` , -`hpo.ancestors`,-`hpo.description`)%>%
  distinct()



#View(human_genes_TPM_None_pheno_disease)

# The number of Tissues that each gene is expressed in

human_genes_TPM_None_pheno_disease$NUMBER_OF_YES <- rowSums(human_genes_TPM_None_pheno_disease == "Yes")
human_genes_TPM_None_pheno_disease$NUMBER_OF_NO <- rowSums(human_genes_TPM_None_pheno_disease == "No")

human_genes_TPM_None_pheno_disease<-human_genes_TPM_None_pheno_disease %>%
  select_all()%>%
  distinct()

#View(human_genes_TPM_None_pheno_disease)



## COUNT THE NUMBER OF GENES EXPRESSED FOR EACH TISSUE

human_genes_TPM_None_pheno_disease_tissue_D<-human_genes_TPM_None_pheno_disease%>%
  select_all()%>%
  select(-NUMBER_OF_NO, -NUMBER_OF_YES, -`HGNC.ID`, -`Description`, -`Type.x`)

human_genes_TPM_None_pheno_disease_tissue_D <- colSums(human_genes_TPM_None_pheno_disease_tissue_D == "Yes")


#View(human_genes_TPM_None_pheno_disease_tissue_D)


human_genes_TPM_None_pheno_disease_tissue_D<-data.frame(human_genes_TPM_None_pheno_disease_tissue_D)

human_genes_TPM_None_pheno_disease_tissue_D<-setNames(cbind(rownames(human_genes_TPM_None_pheno_disease_tissue_D), 
                                                            human_genes_TPM_None_pheno_disease_tissue_D, row.names = NULL), 
                                                          c("Tissue", "Number_of_Yes"))
#View(human_genes_TPM_None_pheno_disease_tissue_D)





#####~~~~~~!!!!!!!!!!!!!!IMPORTANT KEEP FOR DISSERTATION!!!!!!!!!!!~~~~~~##############

# PLOT THE NUMBER OF GENES THAT ARE EXPRESSED IN EACH TISSUE, where gene expression is >0 and not disease specific

library(ggplot2)


#dev.new(width=60, height=30) # to open in a new window

p_TPM_None_tissue_D<-ggplot(data=human_genes_TPM_None_pheno_disease_tissue_D,
                                aes(x=reorder(human_genes_TPM_None_pheno_disease_tissue_D$Tissue, 
                                              -human_genes_TPM_None_pheno_disease_tissue_D$Number_of_Yes),
                                    y=human_genes_TPM_None_pheno_disease_tissue_D$Number_of_Yes, 
                                    fill=human_genes_TPM_None_pheno_disease_tissue_D$Number_of_Yes))
p_TPM_None_tissue_D<-p_TPM_None_tissue_D %>%
  + labs(x= "Tissue", y="Number of Genes", 
         subtitle= "Genes Associated with Disease")%>%
  + scale_fill_gradient(low = "blue", 
                        high = "red")%>% 
  + scale_y_continuous(breaks = seq(0, 14000, by = 1000))%>%
  
  + guides(fill=guide_legend(title="Number of Genes Expressed:  ")) %>%
  + theme(axis.title=element_text(size=14,face="bold"),
          plot.caption=element_text(face = "italic", size=12, hjust = 0),
          text = element_text(size=10),
          legend.position = "bottom",legend.direction = "horizontal",
          legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
          axis.text.y= element_text(size=8),
          axis.text.x= element_text(size=10),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +geom_bar(stat = "identity")%>%    # to create a stacked barchart
  +geom_text(aes(label=paste0(round(human_genes_TPM_None_pheno_disease_tissue_D$Number_of_Yes,1))),size = 3, position = position_stack(vjust = 0.5), colour=c("white"), fontface='bold') %>%
  
  + coord_flip()%>%
  + ggsave(filename=paste("./Plots/Human/Tissue/","p_TPM_None_tissue_D.png",sep=" "), limitsize = TRUE)



###################################################################################################################


#Threshold >0 for all disease associated genes

human_genes_TPM_greater0_pheno_disease<-human_genes_TPM_greater0_pheno%>%
  select_all()%>%
  select(-`HPO-Term-ID` , -`hpo.ancestors`,-`hpo.description`)%>%
  distinct()



##View(human_genes_TPM_greater0_pheno_disease)

# The number of Tissues that each gene is expressed in

human_genes_TPM_greater0_pheno_disease$NUMBER_OF_YES <- rowSums(human_genes_TPM_greater0_pheno_disease == "Yes")
human_genes_TPM_greater0_pheno_disease$NUMBER_OF_NO <- rowSums(human_genes_TPM_greater0_pheno_disease == "No")

human_genes_TPM_greater0_pheno_disease<-human_genes_TPM_greater0_pheno_disease %>%
  select_all()%>%
  distinct()

##View(human_genes_TPM_greater0_pheno_disease)



## COUNT THE NUMBER OF GENES EXPRESSED FOR EACH TISSUE

human_genes_TPM_greater0_pheno_disease_tissue_D<-human_genes_TPM_greater0_pheno_disease%>%
  select_all()%>%
  select(-NUMBER_OF_NO, -NUMBER_OF_YES, -`HGNC.ID`, -`Description`, -`Type.x`)

human_genes_TPM_greater0_pheno_disease_tissue_D <- colSums(human_genes_TPM_greater0_pheno_disease_tissue_D == "Yes")


##View(human_genes_TPM_greater0_pheno_disease_tissue_D)


human_genes_TPM_greater0_pheno_disease_tissue_D<-data.frame(human_genes_TPM_greater0_pheno_disease_tissue_D)

human_genes_TPM_greater0_pheno_disease_tissue_D<-setNames(cbind(rownames(human_genes_TPM_greater0_pheno_disease_tissue_D), 
                                                                human_genes_TPM_greater0_pheno_disease_tissue_D, row.names = NULL), 
                                                     c("Tissue", "Number_of_Yes"))
##View(human_genes_TPM_greater0_pheno_disease_tissue_D)





#####~~~~~~!!!!!!!!!!!!!!IMPORTANT KEEP FOR DISSERTATION!!!!!!!!!!!~~~~~~##############

# PLOT THE NUMBER OF GENES THAT ARE EXPRESSED IN EACH TISSUE, where gene expression is >0 and not disease specific

library(ggplot2)


#dev.new(width=60, height=30) # to open in a new window

p_TPM_greater0_tissue_D<-ggplot(data=human_genes_TPM_greater0_pheno_disease_tissue_D,
                            aes(x=reorder(human_genes_TPM_greater0_pheno_disease_tissue_D$Tissue, 
                                                             -human_genes_TPM_greater0_pheno_disease_tissue_D$Number_of_Yes),
                                 y=human_genes_TPM_greater0_pheno_disease_tissue_D$Number_of_Yes, 
                                 fill=human_genes_TPM_greater0_pheno_disease_tissue_D$Number_of_Yes))
p_TPM_greater0_tissue_D<-p_TPM_greater0_tissue_D %>%
  + labs(x= "Tissue", y="Number of Genes", 
         subtitle= "Genes Associated with Disease")%>%
  + scale_fill_gradient(low = "blue", 
                        high = "red")%>% 
  + scale_y_continuous(breaks = seq(0, 14000, by = 1000))%>%
  
  + guides(fill=guide_legend(title="Number of Genes Expressed:  ")) %>%
  + theme(axis.title=element_text(size=14,face="bold"),
          plot.caption=element_text(face = "italic", size=12, hjust = 0),
          text = element_text(size=10),
          legend.position = "bottom",legend.direction = "horizontal",
          legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
          axis.text.y= element_text(size=8),
          axis.text.x= element_text(size=10),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +geom_bar(stat = "identity")%>%    # to create a stacked barchart
  +geom_text(aes(label=paste0(round(human_genes_TPM_greater0_pheno_disease_tissue_D$Number_of_Yes,1))),size = 3, position = position_stack(vjust = 0.5), colour=c("white"), fontface='bold') %>%
  
  + coord_flip()%>%
  + ggsave(filename=paste("./Plots/Human/Tissue/","p_TPM_greater0_tissue_D.png",sep=" "), limitsize = TRUE)



###################################################################################################################


#Threshold >0.1 for all disease associated genes

human_genes_TPM_0.1_pheno_disease<-human_genes_TPM_0.1_pheno%>%
  select_all()%>%
  select(-`HPO-Term-ID` , -`hpo.ancestors`,-`hpo.description`)%>%
  distinct()



##View(human_genes_TPM_0.1_pheno_disease)

#str(human_genes_TPM_0.1_pheno_disease)


# The number of Tissues that each gene is expressed in

human_genes_TPM_0.1_pheno_disease$NUMBER_OF_YES <- rowSums(human_genes_TPM_0.1_pheno_disease == "Yes")
human_genes_TPM_0.1_pheno_disease$NUMBER_OF_NO <- rowSums(human_genes_TPM_0.1_pheno_disease == "No")

human_genes_TPM_0.1_pheno_disease<-human_genes_TPM_0.1_pheno_disease %>%
  select_all()%>%
  distinct()

##View(human_genes_TPM_0.1_pheno_disease)



## COUNT THE NUMBER OF GENES EXPRESSED FOR EACH TISSUE

human_genes_TPM_0.1_pheno_disease_tissue_D<-human_genes_TPM_0.1_pheno_disease%>%
  select_all()%>%
  select(-NUMBER_OF_NO, -NUMBER_OF_YES, -`HGNC.ID`, -`Description`, -`Type.x`)

human_genes_TPM_0.1_pheno_disease_tissue_D <- colSums(human_genes_TPM_0.1_pheno_disease_tissue_D == "Yes")


##View(human_genes_TPM_0.1_pheno_disease_tissue_D)


human_genes_TPM_0.1_pheno_disease_tissue_D<-data.frame(human_genes_TPM_0.1_pheno_disease_tissue_D)

human_genes_TPM_0.1_pheno_disease_tissue_D<-setNames(cbind(rownames(human_genes_TPM_0.1_pheno_disease_tissue_D), 
                                                           human_genes_TPM_0.1_pheno_disease_tissue_D, row.names = NULL), 
                                        c("Tissue", "Number_of_Yes"))
##View(human_genes_TPM_0.1_pheno_disease_tissue_D)





#####~~~~~~!!!!!!!!!!!!!!IMPORTANT KEEP FOR DISSERTATION!!!!!!!!!!!~~~~~~##############

# PLOT THE NUMBER OF GENES THAT ARE EXPRESSED IN EACH TISSUE, where gene expression is >/=0.1 and not disease specific

# NEW GRAPH
library(ggplot2)


#dev.new(width=60, height=30) # to open in a new window

p_TPM_0.1_tissue_D<-ggplot(data=human_genes_TPM_0.1_pheno_disease_tissue_D
                            , aes(x=reorder(human_genes_TPM_0.1_pheno_disease_tissue_D$Tissue, 
                                            -human_genes_TPM_0.1_pheno_disease_tissue_D$Number_of_Yes),
                                  y=human_genes_TPM_0.1_pheno_disease_tissue_D$Number_of_Yes, 
                                            
                                  fill=human_genes_TPM_0.1_pheno_disease_tissue_D$Number_of_Yes))
p_TPM_0.1_tissue_D<-p_TPM_0.1_tissue_D %>%
  
  
  + labs(x= "Tissue", y="Number of Genes", 
         subtitle= "Genes Associated with Disease")%>%
  + scale_fill_gradient(low = "blue", 
                        high = "red")%>% 
  + scale_y_continuous(breaks = seq(0, 14000, by = 1000))%>%
  
  + guides(fill=guide_legend(title="Number of Genes Expressed:  ")) %>%

    + theme(axis.title=element_text(size=14,face="bold"),
            plot.caption=element_text(face = "italic", size=12, hjust = 0),
            text = element_text(size=10),
            legend.position = "bottom",legend.direction = "horizontal",
          legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
          axis.text.y= element_text(size=8),
          axis.text.x= element_text(size=10),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +geom_bar(stat = "identity")%>%    # to create a stacked barchart
  +geom_text(aes(label=paste0(round(human_genes_TPM_0.1_pheno_disease_tissue_D$Number_of_Yes,1))),
             size = 3, position = position_stack(vjust = 0.5), colour=c("white"), fontface='bold') %>%
  
  + coord_flip()%>%
  + ggsave(filename=paste("./Plots/Human/Tissue/","p_TPM_0.1_tissue_D.png",sep=" "), limitsize = TRUE)

#######################################################################################################################################################


#Threshold >/=1 for all disease associated genes

human_genes_TPM_1_pheno_disease<-human_genes_TPM_1_pheno%>%
  select_all()%>%
  select(-`HPO-Term-ID` , -`hpo.ancestors`,-`hpo.description`)%>%
  distinct()



##View(human_genes_TPM_1_pheno_disease)


# The number of Tissues that each gene is expressed in

human_genes_TPM_1_pheno_disease$NUMBER_OF_YES <- rowSums(human_genes_TPM_1_pheno_disease == "Yes")
human_genes_TPM_1_pheno_disease$NUMBER_OF_NO <- rowSums(human_genes_TPM_1_pheno_disease == "No")

human_genes_TPM_1_pheno_disease<-human_genes_TPM_1_pheno_disease %>%
  select_all()%>%
  distinct()

##View(human_genes_TPM_1_pheno_disease)



## COUNT THE NUMBER OF GENES EXPRESSED FOR EACH TISSUE

human_genes_TPM_1_pheno_disease_tissue_D<-human_genes_TPM_1_pheno_disease%>%
  select_all()%>%
  select(-NUMBER_OF_NO, -NUMBER_OF_YES, -`HGNC.ID`, -`Description`, -`Type.x`)

human_genes_TPM_1_pheno_disease_tissue_D <- colSums(human_genes_TPM_1_pheno_disease_tissue_D == "Yes")


##View(human_genes_TPM_1_pheno_disease_tissue_D)


human_genes_TPM_1_pheno_disease_tissue_D<-data.frame(human_genes_TPM_1_pheno_disease_tissue_D)

human_genes_TPM_1_pheno_disease_tissue_D<-setNames(cbind(rownames(human_genes_TPM_1_pheno_disease_tissue_D), 
                                                                human_genes_TPM_1_pheno_disease_tissue_D, row.names = NULL), 
                                                          c("Tissue", "Number_of_Yes"))
##View(human_genes_TPM_1_pheno_disease_tissue_D)



human_genes_TPM_1_pheno_disease_tissue_D$Tissue <- reorder(human_genes_TPM_1_pheno_disease_tissue_D$Tissue
                                                    , human_genes_TPM_1_pheno_disease_tissue_D$Number_of_Yes)

##View(human_genes_TPM_1_pheno_disease_tissue_D)


#####~~~~~~!!!!!!!!!!!!!!IMPORTANT KEEP FOR DISSERTATION!!!!!!!!!!!~~~~~~##############

# PLOT THE NUMBER OF GENES THAT ARE EXPRESSED IN EACH TISSUE, where gene expression is >/=0.1 and not disease specific

library(ggplot2)

#dev.new(width=60, height=30) # to open in a new window





p_TPM_1_tissue_D<-ggplot(data=human_genes_TPM_1_pheno_disease_tissue_D
                           , aes(x=reorder(human_genes_TPM_1_pheno_disease_tissue_D$Tissue, 
                                           -human_genes_TPM_1_pheno_disease_tissue_D$Number_of_Yes),
                                 y=human_genes_TPM_1_pheno_disease_tissue_D$Number_of_Yes, 
                                 
                                 fill=human_genes_TPM_1_pheno_disease_tissue_D$Number_of_Yes))
p_TPM_1_tissue_D<-p_TPM_1_tissue_D %>%
  
  
  + labs(x= "Tissue", y="Number of Genes", 
         subtitle= "Genes Associated with Disease")%>%
  + scale_fill_gradient(low = "blue", 
                        high = "red")%>% 
  + scale_y_continuous(breaks = seq(0, 14000, by = 1000))%>%
  
  + guides(fill=guide_legend(title="Number of Genes Expressed:  ")) %>%
  + theme(axis.title=element_text(size=14,face="bold"),
          plot.caption=element_text(face = "italic", size=12, hjust = 0),
          text = element_text(size=10),
          legend.position = "bottom",legend.direction = "horizontal",
          legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
          axis.text.y= element_text(size=8),
          axis.text.x= element_text(size=10),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +geom_bar(stat = "identity")%>%    # to create a stacked barchart
  +geom_text(aes(label=paste0(round(human_genes_TPM_1_pheno_disease_tissue_D$Number_of_Yes,1))),size = 3, position = position_stack(vjust = 0.5), colour=c("white"), fontface='bold') %>%
  
  + coord_flip()%>%
  + ggsave(filename=paste("./Plots/Human/Tissue/","p_TPM_1_tissue_D.png",sep=" "), limitsize = TRUE)



#####################################################################################################################################################
########################################################################################################################################################

# SHOW GRAPHS IN A REPRESENTABLE WAY

library(gridExtra)
library(grid)

# NO EXPRESSION :=0

dev.new(width=60, height=30) # to open in a new window


grid.arrange(
  p_TPM_None_tissue_D,
  p_None_tissue_ND,
  nrow = 1,
  top = textGrob(
    "Gene Expression =0 (TPM)",
    gp = gpar(fontface = 2, fontsize = 18),
    hjust = 0.5
  ),
  bottom = textGrob(
    "",
    gp = gpar(fontface = 3, fontsize = 9),
    hjust = 0.5  )
)

########################################################################################################################################################
# >0


dev.new(width=60, height=30) # to open in a new window


grid.arrange(
  p_TPM_greater0_tissue_D,
  p_greater0_tissue_ND,
  nrow = 1,
  top = textGrob(
    "Gene Expression >0 (TPM)",
    gp = gpar(fontface = 2, fontsize = 18),
    hjust = 0.5
  ),
  bottom = textGrob(
    "",
    gp = gpar(fontface = 3, fontsize = 9),
    hjust = 0.5  )
)

########################################################################################################################################################
# 0.1



dev.new(width=60, height=30) # to open in a new window


grid.arrange(
  p_TPM_0.1_tissue_D ,
  p_TPM_0.1_tissue_ND,
  nrow = 1,
  top = textGrob(
    "Gene Expression >/=0.1 (TPM)",
    gp = gpar(fontface = 2, fontsize = 18),
    hjust = 0.5
  ),
  bottom = textGrob(
    "",
    gp = gpar(fontface = 3, fontsize = 9),
    hjust = 0.5  )
)


#####################################################################################################################################################

# >1


dev.new(width=60, height=30) # to open in a new window


grid.arrange(
  p_TPM_0.1_tissue_D ,
  p_TPM_0.1_tissue_ND,
  nrow = 1,
  top = textGrob(
    "Gene Expression >/=1 (TPM)",
    gp = gpar(fontface = 2, fontsize = 18),
    hjust = 0.5
  ),
  bottom = textGrob(
    "",
    gp = gpar(fontface = 3, fontsize = 9),
    hjust = 0.5  )
)



#####################################################################################################################################################
########################################################################################################################################################


# Merge the HGNC.ID dataframes for each threshold to get a disease and non-disease box-plot

##View(human_genes_TPM_None_pheno_disease)

HGNC.ID_DISEASE_None<- human_genes_TPM_None_pheno_disease%>%
  select(NUMBER_OF_YES)

HGNC.ID_NO_DISEASE_None<- human_genes_TPM_None_pheno_No_disease%>%
  select(NUMBER_OF_YES)

names(HGNC.ID_DISEASE_None)[names(HGNC.ID_DISEASE_None) == 'NUMBER_OF_YES'] <- 'VALUE'
names(HGNC.ID_NO_DISEASE_None)[names(HGNC.ID_NO_DISEASE_None) == 'NUMBER_OF_YES'] <- 'VALUE'





##View(HGNC.ID_DISEASE)
##View(HGNC.ID_NO_DISEASE)

HGNC.ID_DISEASE_None["GROUP"] <- "Disease Associated" # That creates the new column named "MY_NEW_COLUMN" 
HGNC.ID_NO_DISEASE_None["GROUP"] <- "Not Disease Associated" # That creates the new column named "MY_NEW_COLUMN" filled with "NA"


HGNC.ID_DISEASE_NON_DI_None<- rbind(HGNC.ID_DISEASE_None, 
                                    HGNC.ID_NO_DISEASE_None)
# # Make Violin plot

#View(HGNC.ID_DISEASE_NON_DI_None)
#table(HGNC.ID_DISEASE_NON_DI_greater0$GROUP)

dev.new(width=60, height=30) # to open in a new window

p_None_violin<-ggplot(HGNC.ID_DISEASE_NON_DI_None, 
                          aes(x=HGNC.ID_DISEASE_NON_DI_None$GROUP, 
                              y=HGNC.ID_DISEASE_NON_DI_None$VALUE, fill=GROUP))+ geom_violin(trim=FALSE)
p_None_violin<-p_None_violin%>%
  + scale_y_continuous(breaks = seq(-25, 60, by = 20))%>%
  
  + labs(x= "Disease/Non-Disease Associated Genes", y="Number of Tissues", 
         title= "The Number of Tissues in which Disease and Non-Disease Associated Genes are Expressed at a Threshold of = 0 TPM ") %>%
  +stat_summary(fun.data="mean_sdl", 
                geom="pointrange", width=0.2, colour="black" )%>%
  +stat_compare_means(method ="wilcox.test",label.y = 80,label.x = 1.5,paired = FALSE, 
                      aes(label = paste0(..method.., "\n", "p =", ..p.format..)))    # Add global p-valu   # Add global p-valu


p_None_violin

#res <- wilcox.test(VALUE ~ GROUP, data = HGNC.ID_DISEASE_NON_DI_greater0,
#                   exact = FALSE)



HGNC.ID_None_MEAN_SD <- aggregate(VALUE~ GROUP, HGNC.ID_DISEASE_NON_DI_None, function(x) c(mean = mean(x), sd = sd(x)))

HGNC.ID_None_MEAN_SD




################################################################################################################################
## >0

# Merge the HGNC.ID dataframes for each threshold to get a disease and non-disease box-plot

##View(human_genes_TPM_greater0_pheno_disease)

HGNC.ID_DISEASE<- human_genes_TPM_greater0_pheno_disease%>%
  select(NUMBER_OF_YES)

HGNC.ID_NO_DISEASE<- human_genes_greater0_pheno_No_disease%>%
  select(NUMBER_OF_YES)

names(HGNC.ID_DISEASE)[names(HGNC.ID_DISEASE) == 'NUMBER_OF_YES'] <- 'VALUE'
names(HGNC.ID_NO_DISEASE)[names(HGNC.ID_NO_DISEASE) == 'NUMBER_OF_YES'] <- 'VALUE'


                                  


##View(HGNC.ID_DISEASE)
##View(HGNC.ID_NO_DISEASE)

HGNC.ID_DISEASE["GROUP"] <- "Disease Associated" # That creates the new column named "MY_NEW_COLUMN" 
HGNC.ID_NO_DISEASE["GROUP"] <- "Not Disease Associated" # That creates the new column named "MY_NEW_COLUMN" filled with "NA"


HGNC.ID_DISEASE_NON_DI_greater0<- rbind(HGNC.ID_DISEASE, 
                                    HGNC.ID_NO_DISEASE)
# # Make Violin plot

##View(HGNC.ID_DISEASE_NON_DI_greater0)
#table(HGNC.ID_DISEASE_NON_DI_greater0$GROUP)

dev.new(width=60, height=30) # to open in a new window

p_greater0_violin<-ggplot(HGNC.ID_DISEASE_NON_DI_greater0, 
                          aes(x=HGNC.ID_DISEASE_NON_DI_greater0$GROUP, 
                              y=HGNC.ID_DISEASE_NON_DI_greater0$VALUE, fill=GROUP))+ geom_violin(trim=FALSE)
p_greater0_violin<-p_greater0_violin%>%
  + scale_y_continuous(breaks = seq(-25, 80, by = 20))%>%
  
  + labs(x= "Disease/Non-Disease Associated Genes", y="Number of Tissues", 
         title= "The Number of Tissues in which Disease and Non-Disease Associated Genes are Expressed at a Threshold of >0 ") %>%
  +stat_summary(fun.data="mean_sdl", 
                geom="pointrange", width=0.2, colour="black" )%>%
  +stat_compare_means(method ="wilcox.test",label.y = 80,label.x = 1.5,paired = FALSE, 
                      aes(label = paste0(..method.., "\n", "p =", ..p.format..)))    # Add global p-valu   # Add global p-valu


p_greater0_violin

#res <- wilcox.test(VALUE ~ GROUP, data = HGNC.ID_DISEASE_NON_DI_greater0,
#                   exact = FALSE)



HGNC.ID_0_MEAN_SD <- aggregate(VALUE~ GROUP, HGNC.ID_DISEASE_NON_DI_greater0, function(x) c(mean = mean(x), sd = sd(x)))

HGNC.ID_0_MEAN_SD

 
 
 
################################################################################################################################
## 0.1


# Merge the HGNC.ID dataframes fro each threshold to get a disease and non-disease box-plot

HGNC.ID_DISEASE_0.1_d<- human_genes_TPM_0.1_pheno_disease%>%
  select(NUMBER_OF_YES)

HGNC.ID_NO_DISEASE_0.1_nd<- human_genes_TPM_0.1_pheno_No_disease%>%
  select(NUMBER_OF_YES)

names(HGNC.ID_DISEASE_0.1_d)[names(HGNC.ID_DISEASE_0.1_d) == 'NUMBER_OF_YES'] <- 'VALUE'
names(HGNC.ID_NO_DISEASE_0.1_nd)[names(HGNC.ID_NO_DISEASE_0.1_nd) == 'NUMBER_OF_YES'] <- 'VALUE'



# #View(HGNC.ID_DISEASE_0.1_d)
# #View(HGNC.ID_NO_DISEASE_0.1_nd)

HGNC.ID_DISEASE_0.1_d["GROUP"] <- "Disease Associated" # That creates the new column named "MY_NEW_COLUMN" filled with "NA"
HGNC.ID_NO_DISEASE_0.1_nd["GROUP"] <- "Not Disease Associated" # That creates the new column named "MY_NEW_COLUMN" filled with "NA"


HGNC.ID_DISEASE_NON_DI_0.1<- rbind(HGNC.ID_DISEASE_0.1_d, 
                                   HGNC.ID_NO_DISEASE_0.1_nd)
# # Make Violin plot

# #View(HGNC.ID_DISEASE_NON_DI_0.1)

dev.new(width=60, height=30) # to open in a new window
 
p_0.1_violin<- ggplot(HGNC.ID_DISEASE_NON_DI_0.1, 
                      aes(x=HGNC.ID_DISEASE_NON_DI_0.1$GROUP, 
                          y=HGNC.ID_DISEASE_NON_DI_0.1$VALUE, fill= GROUP))+ geom_violin(trim=FALSE)
p_0.1_violin<-p_0.1_violin%>%
  + scale_y_continuous(breaks = seq(-25, 80, by = 20))%>%
  + labs(x= "Disease/Non-Disease Associated Genes", y="Number of Tissues", 
         title= "The Number of Tissues in which Disease and Non-Disease Associated Genes are Expressed at a Threshold of >/=0.1 ") %>%
  +stat_summary(fun.data="mean_sdl", 
                geom="pointrange", width=0.2, colour="black" )%>%
  +stat_compare_means(method ="wilcox.test",label.y = 80,label.x = 1.5,paired = FALSE, 
                      aes(label = paste0(..method.., "\n", "p =", ..p.format..)))    # Add global p-valu   # Add global p-valu


p_0.1_violin


### CALCULATE THE MEAN AND THE SD FOR EACH GROUP

HGNC.ID_0.1_MEAN_SD <- aggregate(VALUE~ GROUP, HGNC.ID_DISEASE_NON_DI_0.1, function(x) c(mean = mean(x), sd = sd(x)))

HGNC.ID_0.1_MEAN_SD





#Add the mean/ sd into a table and onto the graph

################################################################################################################################
## >1


# Merge the HGNC.ID dataframes fro each threshold to get a disease and non-disease box-plot

HGNC.ID_DISEASE_1_d<- human_genes_TPM_1_pheno_disease%>%
  select(NUMBER_OF_YES, HGNC.ID)%>%
  distinct()

HGNC.ID_NO_DISEASE_1_nd<- human_genes_TPM_1_pheno_No_disease%>%
  select(NUMBER_OF_YES, HGNC.ID)%>%
  distinct()

names(HGNC.ID_DISEASE_1_d)[names(HGNC.ID_DISEASE_1_d) == 'NUMBER_OF_YES'] <- 'VALUE'
names(HGNC.ID_NO_DISEASE_1_nd)[names(HGNC.ID_NO_DISEASE_1_nd) == 'NUMBER_OF_YES'] <- 'VALUE'


# 
# #View(HGNC.ID_DISEASE_1_d)
# #View(HGNC.ID_NO_DISEASE_1_nd)

HGNC.ID_DISEASE_1_d["GROUP"] <- "Disease Associated" # That creates the new column named "MY_NEW_COLUMN" 
HGNC.ID_NO_DISEASE_1_nd["GROUP"] <- "Not Disease Associated" # That creates the new column named "MY_NEW_COLUMN"


HGNC.ID_DISEASE_NON_DI_1<- rbind(HGNC.ID_DISEASE_1_d, 
                                 HGNC.ID_NO_DISEASE_1_nd)
# # Make Violin plot

# #View(HGNC.ID_DISEASE_NON_DI_1)




dev.new(width=60, height=30) # to open in a new window

p_1_violin__<-ggplot(HGNC.ID_DISEASE_NON_DI_1, 
                     aes(x=HGNC.ID_DISEASE_NON_DI_1$GROUP, 
                         y=HGNC.ID_DISEASE_NON_DI_1$VALUE, fill=GROUP))+ geom_violin(trim=FALSE)
p_1_violin__<-p_1_violin__%>%
  + scale_y_continuous(breaks = seq(-25, 80, by = 20))%>%
  + labs(x= "Disease/Non-Disease Associated Genes", y="Number of Tissues", 
         title= "The Number of Tissues in which Disease and Non-Disease Associated Genes are Expressed at a Threshold of >/=1 ") %>%
  +stat_summary(fun.data="mean_sdl", 
                geom="pointrange", width=0.2, colour="black" )%>%
  +stat_compare_means(method ="wilcox.test",label.y = 80,label.x = 1.5,paired = FALSE, 
                     aes(label = paste0(..method.., "\n", "p =", ..p.format..)))    # Add global p-valu   # Add global p-valu



p_1_violin__


### CALCULATE THE MEAN AND THE SD FOR EACH GROUP

HGNC.ID_1_MEAN_SD <- aggregate(VALUE~ GROUP, HGNC.ID_DISEASE_NON_DI_1, function(x) c(mean = mean(x), sd = sd(x)))

HGNC.ID_1_MEAN_SD


#####################################################################################################################################################
########################################################################################################################################################

#FOR EACH HPO.DESCRIPTION COUNT THE NUMBER OF GENES EXPRESSED FOR ALL GENES WITH HPO.DESC WHICH ARE ALL DISEASE ASSOCIATED

#all_disease_genes_pheno 

# #View(all_disease_genes_hpo_desc)

# TO REMOVED EVERYTHING ELSE ASIDE FROM THE GENE'S HGNC.ID AND THE HPO.DESCRIPTION

all_disease_genes_hpo_desc<- all_disease_genes_pheno%>%
  select(hpo.description,HGNC.ID)%>%
  distinct()

all_disease_genes_hpo_desc<-all_disease_genes_hpo_desc%>%
  select_all()%>%
  dplyr::group_by(hpo.description,HGNC.ID) %>%
  dplyr::summarise(n=n()) %>%
  distinct()


names(all_disease_genes_hpo_desc)[names(all_disease_genes_hpo_desc)=="n"] <- "Frequency"

all_disease_genes_hpo_desc_aggr<-aggregate(Frequency~ hpo.description, data = all_disease_genes_hpo_desc, sum)

# #View(all_disease_genes_hpo_desc_aggr)

dev.new(width=60, height=30) # to open in a new window

p_all_disease_genes_hpo_desc<-ggplot(data=all_disease_genes_hpo_desc_aggr
                         , aes(x=reorder(all_disease_genes_hpo_desc_aggr$hpo.description,- all_disease_genes_hpo_desc_aggr$Frequency),
                               y=all_disease_genes_hpo_desc_aggr$Frequency, 
                               fill=all_disease_genes_hpo_desc_aggr$Frequency))

p_all_disease_genes_hpo_desc<-p_all_disease_genes_hpo_desc %>%
  
  + labs(x= "HPO.Description", y="Number of Genes", 
         title= "The Number of Genes Expressed in Different HPO.Descriptions ")%>%
  + scale_fill_gradient(low = "blue", 
                        high = "red")%>% 
  + scale_y_continuous(breaks = seq(0, 14000, by = 1000))%>%
  + theme(legend.title = element_blank(),
          legend.position = "none",
          axis.title=element_text(size=14,face="bold"),
          plot.caption=element_text(face = "italic", size=12, hjust = 0),
          text = element_text(size=10),
          axis.text.x= element_text(size=10),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +geom_bar(stat = "identity")%>%    # to create a stacked barchart
  +geom_text(show.legend = FALSE,aes(label=paste0(round(all_disease_genes_hpo_desc_aggr$Frequency,1))),size = 5, 
             position = position_stack(vjust = 0.5), colour=c("black"), fontface='bold') %>%
  
  + coord_flip()
  + ggsave(filename=paste("./Plots/Human/HPO.DESCRIPTION/p_all_disease_genes_hpo_desc.png",sep=" "), limitsize = TRUE)


p_all_disease_genes_hpo_desc







#######################################################################################################################################################
#########################################################################################################################################################

## Load the 'GTEX.Tissues.to.HPO.txt' file so that the HPO.Descriptions for each threshold can be plotted against the tissues



#GTEX.Tissues.to.HPO <- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/GTEX.Tissues.to.HPO.txt",delim = "\t")

#REMOVED EMPTY COLUMNS BY SELECTING THOSE THAT WERE NEEDED
GTEX.Tissues.to.HPO <-GTEX.Tissues.to.HPO %>%
  select(GTEX.tissues, HPO.id, HPO.description,
        HPO.superclass.id, HPO.superclass.description, HPO.superclass.physiological.system)%>%
  distinct()


######################################################################################################################################


###########~~~~~~~~~~~~~~~~~~~~~~~~~~~~CREATE A FUNCTION~~~~~~~~~~~~~~~~~~~~~~~~~~~################
# GTEX.Tissues= GTEX.Tissues.to.HPO 
# human_genes= human_genes_TPM_greater0_pheno
# Abnormality = "Abnormality of the nervous system"
# HPO.ID= "HP:0000707" 
# System
#Threshold
GTEX.to.human_genes <- function(GTEX.Tissues ,human_genes, Abnormality,
                       HPO.ID, System, Threshold,System_Threshold){
  library(dplyr)
  
  GTEX.ABN.<-GTEX.Tissues %>%
    select_all()%>%
    filter(HPO.superclass.description== Abnormality | HPO.superclass.id == HPO.ID)%>%
    distinct()
  
  ##View(GTEX.ABN.)
    
  ## Change the dataframe so that the tissues are all in rows with their expressions
  
  human_genes_cleaned <- human_genes %>% 
    gather(Tissues,Expression, `Adipose - Subcutaneous`:`Whole Blood`)
  
  
  human_genes_cleaned_ABN<-human_genes_cleaned %>%
    select(HGNC.ID,hpo.description,Tissues, Expression)%>%
    filter( hpo.description==Abnormality)%>%
    filter(Tissues %in% GTEX.ABN.$GTEX.tissues)%>%
    
    distinct()
  
  
  ##View(human_genes_cleaned_ABN)
  
  
  ############ COUNT THE NUMBER OF YES PER TISSUE ########
  
  # 
  human_genes_cleaned_ABN_COUNT<-human_genes_cleaned_ABN%>%
    select_all()%>%
    dplyr::group_by(Tissues,Expression) %>%
    dplyr::summarise(n=n()) %>%
    distinct()
  
  names(human_genes_cleaned_ABN_COUNT)[names(human_genes_cleaned_ABN_COUNT)=="n"] <- "Frequency"
#  #View(human_genes_cleaned_ABN_COUNT)
  
  
  
  
  library(ggplot2)
  GTEX.human_genes_cleaned_ABN_COUNT_p<-ggplot(data=human_genes_cleaned_ABN_COUNT
                                                 , aes(x= human_genes_cleaned_ABN_COUNT$Tissues,
                                                       y=human_genes_cleaned_ABN_COUNT$Frequency,
                                                       fill = Expression))
  
  
  GTEX.human_genes_cleaned_ABN_COUNT_p<-GTEX.human_genes_cleaned_ABN_COUNT_p %>%
    
    + labs(x= "Tissue", y="Number of Genes", 
           title= paste('The Number of Human Genes Expressed \nin',System, 'at a Threshold of: \n', Threshold, ' TPM' ),
           tag = paste())%>%

    + guides(fill=guide_legend(title="Gene Expression:  ")) %>%
    + scale_y_continuous(breaks = seq(0, 14000, by = 250))%>%
    
    + theme(axis.title=element_text(size=14,face="bold"),
            plot.caption=element_text(face = "italic", size=12, hjust = 0),
            text = element_text(size=10),
            legend.position = "bottom",legend.direction = "horizontal",
            legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
            axis.text.x = element_text(size = 7),
            axis.text.y = element_text(size = 12),
            panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
    +scale_fill_manual("Gene Expression:  ", values = c("No" = "cyan1", "Yes"="hotpink1"))%>%
    
    
    + geom_bar(stat = "identity")%>%    # to create a stacked barchart
    +geom_text(aes(label=paste0(round(human_genes_cleaned_ABN_COUNT$Frequency,1))),size = 4, 
               position = position_stack(vjust = 0.5), colour=c("black"), fontface='bold') %>%
    
    + coord_flip()%>%
    + ggsave(filename=paste("./Plots/Human/System/", System_Threshold,".png",sep=" "),
             width = 10, height = 10, dpi = 150, units = "in", device='png')
  
  }



#######################################################################################################################################


All_systems<-GTEX.Tissues.to.HPO%>%
  select(HPO.superclass.description,HPO.superclass.id,HPO.superclass.physiological.system)%>%
  filter(HPO.superclass.description %in% all_disease_genes_pheno$hpo.description |
           HPO.superclass.id %in% all_disease_genes_pheno$hpo.ancestors )%>%
  distinct()

#View(All_systems)

df_threshold<-function(df, threshold, threshold_2){
  GTEX.to.human_genes(GTEX.Tissues=GTEX.Tissues.to.HPO,
                      human_genes=df  ,
                      Abnormality= "Abnormality of connective tissue", HPO.ID= "HP:0003549",
                      System="connective tissue", Threshold=threshold, System_Threshold= paste("connective_tissue","_",threshold_2)) 
  
  GTEX.to.human_genes(GTEX.Tissues=GTEX.Tissues.to.HPO,
                      human_genes=df  ,
                      Abnormality= "Abnormality of the endocrine system", HPO.ID= "HP:0000818",
                      System="endocrine system", Threshold=threshold, System_Threshold= paste("endocrine_system","_",threshold_2)) 
  
  GTEX.to.human_genes(GTEX.Tissues=GTEX.Tissues.to.HPO,
                      human_genes=df  ,
                      Abnormality= "Abnormality of the cardiovascular system", HPO.ID= "HP:0001626",
                      System="cardiovascular system", Threshold=threshold, System_Threshold= paste("cardiovascular_system","_",threshold_2)) 
  
  GTEX.to.human_genes(GTEX.Tissues=GTEX.Tissues.to.HPO,
                      human_genes=df  ,
                      Abnormality= "Abnormality of the genitourinary system", HPO.ID= "HP:0000119",
                      System="genitourinary system", Threshold=threshold, System_Threshold= paste("genitourinary_system","_",threshold_2)) 
  
  
  GTEX.to.human_genes(GTEX.Tissues=GTEX.Tissues.to.HPO,
                      human_genes=df  ,
                      Abnormality= "Abnormality of the nervous system", HPO.ID= "HP:0000707",
                      System="nervous system", Threshold=threshold, System_Threshold= paste("nervous_system","_",threshold_2)) 
  
  
  GTEX.to.human_genes(GTEX.Tissues=GTEX.Tissues.to.HPO,
                      human_genes=df  ,
                      Abnormality= "Abnormality of the breast", HPO.ID= "HP:0000769",
                      System="breast", Threshold=threshold, System_Threshold= paste("breast","_",threshold_2)) 
  
  
  GTEX.to.human_genes(GTEX.Tissues=GTEX.Tissues.to.HPO,
                      human_genes=df  ,
                      Abnormality= "Abnormality of the digestive system", HPO.ID= "HP:0025031",
                      System="digestive system", Threshold=threshold, System_Threshold= paste("digestive_system","_",threshold_2)) 
  

  GTEX.to.human_genes(GTEX.Tissues=GTEX.Tissues.to.HPO,
                      human_genes=df  ,
                      Abnormality= "Abnormality of the respiratory system", HPO.ID= "HP:0002086",
                      System="respiratory system", Threshold=threshold, System_Threshold= paste("respiratory_system","_",threshold_2)) 
  
  GTEX.to.human_genes(GTEX.Tissues=GTEX.Tissues.to.HPO,
                      human_genes=df  ,
                      Abnormality= "Abnormality of head or neck", HPO.ID= "HP:0000152",
                      System="head or neck", Threshold=threshold, System_Threshold= paste("head_or_neck","_",threshold_2)) 
  
  GTEX.to.human_genes(GTEX.Tissues=GTEX.Tissues.to.HPO,
                      human_genes=df  ,
                      Abnormality= "Abnormality of the skeletal system", HPO.ID= "HP:0000924",
                      System="skeletal system", Threshold=threshold, System_Threshold= paste("skeletal_system","_",threshold_2)) 
  
  
  GTEX.to.human_genes(GTEX.Tissues=GTEX.Tissues.to.HPO,
                      human_genes=df  ,
                      Abnormality= "Abnormality of the musculature", HPO.ID= "HP:0003011",
                      System="musculature", Threshold=threshold, System_Threshold= paste("musculature","_",threshold_2))
  
  
  GTEX.to.human_genes(GTEX.Tissues=GTEX.Tissues.to.HPO,
                      human_genes=df  ,
                      Abnormality= "Abnormality of the integument", HPO.ID= "HP:0001574",
                      System="integument", Threshold=threshold, System_Threshold= paste("integument","_",threshold_2)) 
  
  
  GTEX.to.human_genes(GTEX.Tissues=GTEX.Tissues.to.HPO,
                      human_genes=df  ,
                      Abnormality= "Abnormality of the immune system", HPO.ID= "HP:0002715",
                      System="immune system", Threshold=threshold, System_Threshold= paste("immune_system","_",threshold_2)) 
  
  

  GTEX.to.human_genes(GTEX.Tissues=GTEX.Tissues.to.HPO,
                      human_genes=df  ,
                      Abnormality= "Abnormality of blood and blood-forming tissues", HPO.ID= "HP:0001871",
                      System="blood and blood-forming tissues", Threshold=threshold, System_Threshold= paste("blood_and_blood-forming_tissues","_",threshold_2)) 
  
  
  GTEX.to.human_genes(GTEX.Tissues=GTEX.Tissues.to.HPO,
                      human_genes=df  ,
                      Abnormality= "Abnormality of metabolism/homeostasis", HPO.ID= "HP:0001939",
                      System="metabolism/homeostasis", Threshold=threshold, System_Threshold= paste("metabolism_homeostasis","_",threshold_2)) 
  
}

# 
# df_threshold(df=human_genes_TPM_None_pheno, threshold="= 0", threshold_2="NONE")
# 
# 
# df_threshold(df=human_genes_TPM_greater0_pheno, threshold="> 0", threshold_2="0")
# 
# df_threshold(df=human_genes_TPM_0.1_pheno, threshold=">/= 0.1", threshold_2= "0.1")
# 
# df_threshold(df=human_genes_TPM_1_pheno , threshold=">/= 1", threshold_2="1")
# 
# 
# 
# 
# 




###################################################################################################################################################################

###################################################################################################################################################################


####~~~~~~~~ TO COMPARE EACH TISSUE'S NUMBER OF GENES EXPRESSED AGAINST ALL OTHER TISSUES THAT ARE NOT IN THE SAME SYSTEM ~~~~~~~~~~######
#

# ###~~~~~ THRESHOLD >0 
# 
# ##~~ Change the dataframe so that the tissues are all in rows with their expressions
# 
# 
# human_genes_TPM_greater0_pheno_cleaned <- human_genes_TPM_greater0_pheno %>% 
#   gather("GTEX.tissues","Expression", `Adipose - Subcutaneous`:`Whole Blood`)%>%
#   select(`HPO-Term-ID`,hpo.ancestors,hpo.description,GTEX.tissues,Expression)%>%
#   distinct()
# 
# 
# human_genes_TPM_greater0_pheno_count<-human_genes_TPM_greater0_pheno_cleaned%>%
#   dplyr::group_by(GTEX.tissues,Expression) %>%
#   dplyr::summarise(n=n()) 
# 
# GTEX.Tissues.to.tissues<-GTEX.Tissues.to.HPO%>%
#   select_all()%>%
#   drop_na()%>%
#   select(GTEX.tissues)
# 
# 
# #View(GTEX.Tissues.to.tissues)
# ## There are some tissues that are found in different HPO.superclass.description tehrefore have been duplicated in the GTEX.Tissues.to.tissues
# 
# # Merge the GTEX.Tissues.to.tissues and human_genes_TPM_greater0_pheno_count to get the duplicated tissues
# 
# human_genes_TPM_greater0_pheno_Tissues<-merge(GTEX.Tissues.to.tissues, human_genes_TPM_greater0_pheno_count, by = c("GTEX.tissues", "GTEX.tissues"), all.x = TRUE)
# 
# 
# # The human_genes_TPM_greater0_pheno_count and GTEX.Tissues.to.HPO joind to get the HPO.superclass.description for each tissue
# 
# human_genes_TPM_greater0_pheno_Tissues<-inner_join(GTEX.Tissues.to.HPO,human_genes_TPM_greater0_pheno_Tissues)
# human_genes_TPM_greater0_pheno_Tissues<-human_genes_TPM_greater0_pheno_Tissues%>%
#   select_all()%>%
#   drop_na()%>%    # the NAs are dropped
#   distinct(GTEX.tissues,HPO.id, HPO.superclass.id,HPO.superclass.description,Expression,n ) # duplicates are removed that are not needed
# 
# ##View(human_genes_TPM_greater0_pheno_Tissues)
# 
# 
# #View(human_genes_TPM_greater0_pheno_Tissues)
# 
# human_genes_TPM_greater0_pheno_Tissues$Expression <- as.factor(human_genes_TPM_greater0_pheno_Tissues$Expression)
# human_genes_TPM_greater0_pheno_Tissues$GTEX.tissues <- as.factor(human_genes_TPM_greater0_pheno_Tissues$GTEX.tissues)
# human_genes_TPM_greater0_pheno_Tissues$HPO.superclass.description <- as.factor(human_genes_TPM_greater0_pheno_Tissues$HPO.superclass.description)
# 
# human_genes_TPM_greater0_pheno_Tissues$HPO.id <- as.factor(human_genes_TPM_greater0_pheno_Tissues$HPO.id)
# human_genes_TPM_greater0_pheno_Tissues$HPO.superclass.id <- as.factor(human_genes_TPM_greater0_pheno_Tissues$HPO.superclass.id)
# 
# #View(human_genes_TPM_greater0_pheno_Tissues)
# # 
# # human_genes_TPM_greater0_pheno_cleaned_wide<-human_genes_TPM_greater0_pheno_Tissues %>%
# #   spread(Expression,n)%>%
# #   distinct()
# 
# 
# # #View(human_genes_TPM_greater0_pheno_cleaned_wide)
# 
# ####~~~~~~~~~~~~~~~~~ To do a Fischer test for every Tissue vs Tissues not in the same System~~~~~~~
# 
# 
# 
# 
# greater0_pheno_fischer<-human_genes_TPM_greater0_pheno_Tissues%>%
#   select(GTEX.tissues,HPO.superclass.description,Expression,n)%>%
#   distinct()
# 
# 
# # 
# # #View(greater0_pheno_fischer)
# # #View(human_genes_TPM_greater0_tissue_ND)
# # 
# 
# 
# greater0_pheno_fischer_Adipose_Subcutaneous<-greater0_pheno_fischer%>%
#   select(GTEX.tissues,HPO.superclass.description,Expression,n)%>%
#   filter(HPO.superclass.description !="Abnormality of connective tissue")%>%
#   filter(GTEX.tissues != "Adipose - Subcutaneous")%>%)%>%
#   distinct(GTEX.tissues,Expression,n)
# 
# # row.names(greater0_pheno_fischer_Adipose_Subcutaneous) <- 1:nrow(greater0_pheno_fischer_Adipose_Subcutaneous)
# # #greater0_pheno_fischer_Adipose_Subcutaneous<-t(greater0_pheno_fischer_Adipose_Subcutaneous)
# 
# 
# 
# #changge them into factors
# greater0_pheno_fischer_Adipose_Subcutaneous$GTEX.tissues <- as.factor(greater0_pheno_fischer_Adipose_Subcutaneous$GTEX.tissues)
# greater0_pheno_fischer_Adipose_Subcutaneous$Expression <- as.factor(greater0_pheno_fischer_Adipose_Subcutaneous$Expression)
# str(greater0_pheno_fischer_Adipose_Subcutaneous)
# 
# 
# 
# 
# greater0_pheno_fischer_Adipose_Subc<-greater0_pheno_fischer%>%
#   select(GTEX.tissues,HPO.superclass.description,Expression,n)%>%
#   filter(HPO.superclass.description =="Abnormality of connective tissue")%>%
#   filter(GTEX.tissues == "Adipose - Subcutaneous")%>%
#   distinct(GTEX.tissues,Expression,n)%>%
#   spread(Expression,n)%>%
#   distinct(GTEX.tissues,No, Yes)
# #View(greater0_pheno_fischer_Adipose_Subc)
# 
# 
# # superclass="Abnormality of connective tissue", 
# # Tissue="Adipose - Subcutaneous"
# 
# ## superclass="Abnormality of the genitourinary system", 
# # Tissue="Kidney - Cortex"
# 
# #All tissues without the Adipose_Subcutaneous
# 
# all_tissues<-greater0_pheno_fischer_Adipose_Subcutaneous %>%
#   spread(Expression,n)%>%
#   distinct()
# 
# #View(all_tissues)
# 
# # # trun the factors into numerics
# # all_tissues$No <- as.numeric(as.character(all_tissues$No))
# # 
# # all_tissues$Yes <- as.numeric(as.character(all_tissues$Yes))
# 
# 
# # Sum all of the No and Yes rows to get total 
# 
# all_tissues_<-all_tissues%>%
#   select_all()%>%
#   summarise_at(c("No","Yes"),sum) %>% #sums the columns
#   distinct(No, Yes)%>%
#   mutate(GTEX.tissues="Total")%>%
#   distinct(GTEX.tissues,No,Yes)
#   
#   
# #View(all_tissues_)
# 
# # add name of the row
# #all_tissues_$GTEX.tissues<- "Total"
# 
# #Add column name
# #colnames(all_tissues_) <- c("No","Yes", "GTEX.tissues")
# ##View(all_tissues_)
# 
# # Bind the two datfarames to create a matrix
# total_tissue<-rbind(greater0_pheno_fischer_Adipose_Subc,all_tissues_)
# 
# 
# #define the rownames 
# library(magrittr)
# total_tissue<-total_tissue %>%
#   set_rownames(.$GTEX.tissues) %>% 
#   select(-GTEX.tissues)
# #View(total_tissue)
# 
# #rownames(total_tissue) <- total_tissue$GTEX.tissues
# # 
# # # remove GTEX.tissues as it was an extra column
# # total_tissue<-total_tissue%>%
# #   select(-GTEX.tissues)
# 
# #View(total_tissue)
# 
# 
# # perform Fisher test
# 
# pval<-fisher.test(total_tissue,alternative='two.sided')  # THIS WORKED
# pval
# 
# pvalue<-pval$p.value
# 
# # create a matrix with all tissue and pvalue
# p.value.df<-data.frame("Tissue"= c("Adipose - Subcutaneous"),
#               "HPO.DESCRIPTION"= c("Abnormality of connective tissue"),
#               p.value=pvalue)
# 
# #View(p.value.df)
# 
#  
# # 
# # superclass="Abnormality of connective tissue", 
# # Tissue="Adipose - Subcutaneous"
# 
# 

####~~~~~~~~~~~~~CREATE A DATAFRAME WITH ALL TISSUES AND SUPERCLASSES (HPO.DESC/TOPLEVEL)~~~~~~~~~~~~~~~~~~~~########
all_disease_genes_pheno_cleaned <- all_disease_genes_pheno %>% 
  gather("GTEX.tissues","Expression", `Adipose - Subcutaneous`:`Whole Blood`)%>%
  select(HGNC.ID,`HPO-Term-ID`,hpo.ancestors,hpo.description,GTEX.tissues,Expression)%>%
  distinct()
#View(all_disease_genes_pheno_cleaned)

GTEX.Tissues.to.tissues<-GTEX.Tissues.to.HPO%>%
  select_all()%>%
  drop_na()%>%
  select(GTEX.tissues)

## There are some tissues that are found in different HPO.superclass.description tehrefore have been duplicated in the GTEX.Tissues.to.tissues

# Merge the GTEX.Tissues.to.tissues and human_genes_TPM_greater0_pheno_count to get the duplicated tissues

all_disease_genes_pheno_Tissues<-merge(GTEX.Tissues.to.tissues, all_disease_genes_pheno_cleaned, by = c("GTEX.tissues", "GTEX.tissues"), all.x = TRUE)


# The human_genes_TPM_greater0_pheno_count and GTEX.Tissues.to.HPO joind to get the HPO.superclass.description for each tissue

all_disease_genes_pheno_Tissues<-inner_join(GTEX.Tissues.to.HPO,all_disease_genes_pheno_Tissues)
all_disease_genes_pheno_Tissues<-all_disease_genes_pheno_Tissues%>%
  select_all()%>%
  drop_na()%>%    # the NAs are dropped
  distinct(GTEX.tissues,HGNC.ID,HPO.superclass.description) # duplicates are removed that are not needed

#View(all_disease_genes_pheno_Tissues)






#######~~~~~~~~~~~~~~~TRY 2~~~~~~~~~~~~~~~~~~~#


#THRESHOLD DATAFRAME WITH THE GENE EXPRESSIONS

#View(human_genes_TPM_greater0_pheno)


pheno_cleaned <- human_genes_TPM_greater0_pheno %>% 
  gather("GTEX.tissues","Expression", `Adipose - Subcutaneous`:`Whole Blood`)%>%
  select(HGNC.ID,hpo.ancestors,hpo.description,GTEX.tissues,Expression)%>%
  distinct()
View(pheno_cleaned)

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


######!!!!!!!!!!!! HERE I KEPT THE TWO COLUMNS WITH THE HPO.SUPERCLASS THAT IS LINKED TO THE GTEX TISSUE, AND THE HPO.DESC WHICH IS THE PHENOTYPE THE GENE IS ASSOCIATED WITH


pheno_cleaned_Tissues<-inner_join(GTEX.Tissues.to.HPO,pheno_cleaned_Tissues)
pheno_cleaned_Tissues<-pheno_cleaned_Tissues%>%
  select_all()%>%
  drop_na()%>%    # the NAs are dropped
  distinct(GTEX.tissues,hpo.description,HGNC.ID,HPO.superclass.description,Expression ) # duplicates are removed that are not needed

View(pheno_cleaned_Tissues)

###

# COUNT THE NUMBER OF GENES EXPRESSED

pheno_cleaned_Tissues_COUNT<-pheno_cleaned_Tissues%>%
  dplyr::group_by(HPO.superclass.description,Expression) %>%  # group by gene here 
  dplyr::mutate(n=n()) %>%     #mutate keeps the other columns
  mutate(freq= n /sum(as.numeric(n)))%>%
  distinct(HGNC.ID,HPO.superclass.description,Expression, n, freq ) # duplicates are removed that are not needed

  
#View(pheno_cleaned_Tissues_COUNT)

# 
# pheno_cleaned_Tissues_COUNT_hpo.super<-merge(pheno_cleaned_Tissues, pheno_cleaned_Tissues_COUNT, 
#                              by = c("HGNC.ID", "HGNC.ID"), all.x = TRUE)


pheno_cleaned_Tissues_COUNT_hpo.super<-pheno_cleaned_Tissues_COUNT_hpo.super%>%
  select(HGNC.ID, HPO.superclass.description,Expression, freq,n )%>% # duplicates are removed that are not needed
  drop_na()%>%
  distinct()
  

#change the name of EXPRESSION.Y to Expression
names(pheno_cleaned_Tissues_COUNT_hpo.super)[names(pheno_cleaned_Tissues_COUNT_hpo.super) == 'Expression.y'] <- 'Expression'


# CHANGE THEM INTO FACTORS
pheno_cleaned_Tissues_COUNT_hpo.super$Expression <- as.factor(pheno_cleaned_Tissues_COUNT_hpo.super$Expression)
pheno_cleaned_Tissues_COUNT_hpo.super$HPO.superclass.description <- as.factor(pheno_cleaned_Tissues_COUNT_hpo.super$HPO.superclass.description)
pheno_cleaned_Tissues_COUNT_hpo.super$HGNC.ID <- as.factor(pheno_cleaned_Tissues_COUNT_hpo.super$HGNC.ID)

#View(pheno_cleaned_Tissues_COUNT_hpo.super)

####~~~~~~~~~~~~~~~~~ To do a Fischer test for every Tissue vs Tissues not in the same System~~~~~~~




greater0_pheno_fischer<-pheno_cleaned_Tissues_COUNT_hpo.super%>%
  select(HGNC.ID,GTEX.tissues, HPO.superclass.description,Expression,n, freq)%>%
  distinct()

greater0_pheno_fischer_NO_TISSUE<-greater0_pheno_fischer%>%
  select(HGNC.ID,GTEX.tissues,HPO.superclass.description,Expression,n,freq)%>%
  filter(HPO.superclass.description != "Abnormality of the genitourinary system" )%>%
  filter(GTEX.tissues != "Testis")%>%
  distinct(HGNC.ID,GTEX.tissues,Expression,n,freq)

#View(greater0_pheno_fischer_NO_TISSUE)
# superclass="Abnormality of the genitourinary system", 
#Tissue="Testis")
#
#change them into factors
greater0_pheno_fischer_NO_TISSUE$HGNC.ID <- as.factor(greater0_pheno_fischer_NO_TISSUE$HGNC.ID)
greater0_pheno_fischer_NO_TISSUE$Expression <- as.factor(greater0_pheno_fischer_NO_TISSUE$Expression)




greater0_pheno_fischer_tissue<-greater0_pheno_fischer%>%
  select(GTEX.tissues,HPO.superclass.description,Expression,n)%>%
  filter(HPO.superclass.description =="Abnormality of the genitourinary system" )%>%
  filter(GTEX.tissues == "Testis" )%>%
  distinct(GTEX.tissues,Expression,n)%>%
  spread(Expression,n)%>%
  distinct(GTEX.tissues,No, Yes)


#All tissues without the Adipose_Subcutaneous

All_tissues<-greater0_pheno_fischer_NO_TISSUE %>%
  spread(Expression,n)%>%
  distinct()

##View(all_tissues)


# Sum all of the No and Yes rows to get total 

All_tissues_<-All_tissues%>%
  select_all()%>%
  summarise_at(c("No","Yes"),sum) %>% #sums the columns
  distinct(No, Yes)%>%
  mutate(GTEX.tissues="Total")%>%
  distinct(GTEX.tissues,No,Yes)


##View(all_tissues_)


# Bind the two datfarames to create a matrix
Total_Tissues<-rbind(greater0_pheno_fischer_tissue,All_tissues_)


#define the rownames 
library(magrittr)
Total_Tissues<-Total_Tissues %>%
  set_rownames(.$GTEX.tissues) %>% 
  select(-GTEX.tissues)


#View(Total_Tissues)


# perform Fisher test

PVals<-fisher.test(Total_Tissues,alternative='two.sided')  # THIS WORKED

PVALUES<-PVals$p.value

PVALUES









#################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~THE FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#################################



TISSUE_FISHER_TEST<-function(df,superclass, Tissue){  
  
# For each threhsold the number of genes being expressed is counted 
  
  #change df to  long format using gather
  # keep gene id here
  pheno_cleaned <- df %>% 
    gather("GTEX.tissues","Expression", `Adipose - Subcutaneous`:`Whole Blood`)%>%
    select(GTEX.tissues,`HPO-Term-ID`,hpo.ancestors,hpo.description,GTEX.tissues,Expression)%>%
    distinct()
  #View(pheno_cleaned)
  
  
  
  pheno_cleaned_count<-pheno_cleaned%>%
    dplyr::group_by(GTEX.tissues,Expression) %>%  # group by gene here 
    dplyr::summarise(n=n()) 
  
  #View(pheno_cleaned_count)
  
  
  
  GTEX.Tissues.to.tissues<-GTEX.Tissues.to.HPO%>%
    select_all()%>%
    drop_na()%>%
    select(GTEX.tissues)
  
  ## There are some tissues that are found in different HPO.superclass.description tehrefore have been duplicated in the GTEX.Tissues.to.tissues
  
  # Merge the GTEX.Tissues.to.tissues and human_genes_TPM_greater0_pheno_count to get the duplicated tissues
  
  pheno_cleaned_Tissues<-merge(GTEX.Tissues.to.tissues, pheno_cleaned_count, 
                               by = c("GTEX.tissues", "GTEX.tissues"), all.x = TRUE)
  
  
  # The human_genes_TPM_greater0_pheno_count and GTEX.Tissues.to.HPO joind to get the HPO.superclass.description for each tissue
  
  pheno_cleaned_Tissues<-inner_join(GTEX.Tissues.to.HPO,pheno_cleaned_Tissues)
  pheno_cleaned_Tissues<-pheno_cleaned_Tissues%>%
    select_all()%>%
    drop_na()%>%    # the NAs are dropped
    distinct(GTEX.tissues,HPO.id, HPO.superclass.id,HPO.superclass.description,Expression,n ) # duplicates are removed that are not needed
  
  #View(pheno_cleaned_Tissues)
  
  
 # #View(human_genes_TPM_greater0_pheno_Tissues)
  
  pheno_cleaned_Tissues$Expression <- as.factor(pheno_cleaned_Tissues$Expression)
  pheno_cleaned_Tissues$GTEX.tissues <- as.factor(pheno_cleaned_Tissues$GTEX.tissues)
  pheno_cleaned_Tissues$HPO.superclass.description <- as.factor(pheno_cleaned_Tissues$HPO.superclass.description)
  
  pheno_cleaned_Tissues$HPO.id <- as.factor(pheno_cleaned_Tissues$HPO.id)
  pheno_cleaned_Tissues$HPO.superclass.id <- as.factor(pheno_cleaned_Tissues$HPO.superclass.id)
  

  ####~~~~~~~~~~~~~~~~~ To do a Fischer test for every Tissue vs Tissues not in the same System~~~~~~~
  
  
  
  
  greater0_pheno_fischer<-pheno_cleaned_Tissues%>%
    select(GTEX.tissues,HPO.superclass.description,Expression,n)%>%
    distinct()

  greater0_pheno_fischer_NO_TISSUE<-greater0_pheno_fischer%>%
    select(GTEX.tissues,HPO.superclass.description,Expression,n)%>%
    filter(HPO.superclass.description != superclass )%>%
    filter(GTEX.tissues != Tissue)%>%
    distinct(GTEX.tissues,Expression,n)

  
  
  #change them into factors
  greater0_pheno_fischer_NO_TISSUE$GTEX.tissues <- as.factor(greater0_pheno_fischer_NO_TISSUE$GTEX.tissues)
  greater0_pheno_fischer_NO_TISSUE$Expression <- as.factor(greater0_pheno_fischer_NO_TISSUE$Expression)

  
  
  
  greater0_pheno_fischer_tissue<-greater0_pheno_fischer%>%
    select(GTEX.tissues,HPO.superclass.description,Expression,n)%>%
    filter(HPO.superclass.description ==superclass)%>%
    filter(GTEX.tissues == Tissue )%>%
    distinct(GTEX.tissues,Expression,n)%>%
    spread(Expression,n)%>%
    distinct(GTEX.tissues,No, Yes)

  
  #All tissues without the Adipose_Subcutaneous
  
  All_tissues<-greater0_pheno_fischer_NO_TISSUE %>%
    spread(Expression,n)%>%
    distinct()
  
  ##View(all_tissues)
  
  
  # Sum all of the No and Yes rows to get total 
  
  All_tissues_<-All_tissues%>%
    select_all()%>%
    summarise_at(c("No","Yes"),sum) %>% #sums the columns
    distinct(No, Yes)%>%
    mutate(GTEX.tissues="Total")%>%
    distinct(GTEX.tissues,No,Yes)
  
  
  ##View(all_tissues_)
  

  # Bind the two datfarames to create a matrix
  Total_Tissues<-rbind(greater0_pheno_fischer_tissue,All_tissues_)
  
  
  #define the rownames 
  library(magrittr)
  Total_Tissues<-Total_Tissues %>%
    set_rownames(.$GTEX.tissues) %>% 
    select(-GTEX.tissues)
 
  
  #View(Total_Tissues)
  
  
  # perform Fisher test
  
  PVals<-fisher.test(Total_Tissues,alternative='two.sided')  # THIS WORKED

  PVALUES<-PVals$p.value

  return(PVALUES)
   
  }


# #TEST FUNCTION
TISSUE_FISHER_TEST(df=human_genes_TPM_greater0_pheno,
                    superclass="Abnormality of the genitourinary system", 
                    Tissue="Testis")
# superclass= "Abnormality of connective tissue"
# Tissue="Adipose - Subcutaneous"




#Group by tissue and gene 


###########~~~~~~~~~~~~~~~~USE OF FOR-LOOP ~~~~~#

#THE FOR LOOP ITERATES THROUGH THE all_disease_genes_pheno_Tissues DATAFRAME WHICH HAS ALL TISSUES AND HPO.DESC(TOPLEVEL)

# CREATE A data frame with two columns and 55 rows, WHICH HAS THE TISSUE NAMES, HPO.SUPERCLASS AND AN EMPTY COLUMN


# NO GENE EXPRESSION 
PVALUE_Tissues_None <- data.frame(all_disease_genes_pheno_Tissues$GTEX.tissues,
                                      all_disease_genes_pheno_Tissues$HPO.superclass.description,
                                      P.VALUE=vector(length=55)) 


#View(human_genes_TPM_greater0_pheno)
for(i in 1:nrow(all_disease_genes_pheno_Tissues)) {
  ii <- all_disease_genes_pheno_Tissues[i,1]     #FIRST COLUMN
  jj <-all_disease_genes_pheno_Tissues[i,2]      #SECOND COLUMN
  
  print(ii)                                     # TO CHECK IF THE CORRECT ROW IS BEING TAKEN
  print(jj)
  
        
  ## THE OUTCOME OF THE FUNCTION IS SAVED AS PVAL WHICH IS THE FISHER TEST PVALUE 
  PVAL<-TISSUE_FISHER_TEST(df=human_genes_TPM_None_pheno, # THE FUNCTION NEEDS THE DF THAT HAS THE GENE GENE EXPRESSIONS AT SPECIFIC THRESHOLDS
  
                           superclass= paste(jj), 
                           Tissue=paste(ii))
  
  print(PVAL)                  # THE PVALUE IS PRINTED
  
  PVALUE_Tissues_None$P.VALUE[i]<-paste(PVAL)              # THE PVALUE IS ADDED INTO THE DF CREATED ABOVE
  
}

# THE DATAFRAME IS #ViewED
#View(PVALUE_Tissues_None)



# > 0
PVALUE_Tissues_GREATER0 <- data.frame(all_disease_genes_pheno_Tissues$GTEX.tissues,
                                      all_disease_genes_pheno_Tissues$HPO.superclass.description,
                                      P.VALUE=vector(length=55)) 

for(i in 1:nrow(all_disease_genes_pheno_Tissues)) {
  ii <- all_disease_genes_pheno_Tissues[i,1]     #FIRST COLUMN
  jj <-all_disease_genes_pheno_Tissues[i,2]      #SECOND COLUMN
 
  print(ii)                                     # TO CHECK IF THE CORRECT ROW IS BEING TAKEN
  print(jj)
  
           
  
  ## THE OUTCOME OF THE FUNCTION IS SAVED AS PVAL WHICH IS THE FISHER TEST PVALUE 
  PVAL<-TISSUE_FISHER_TEST(df=human_genes_TPM_greater0_pheno, # THE FUNCTION NEEDS THE DF THAT HAS THE GENE GENE EXPRESSIONS AT SPECIFIC THRESHOLDS
                           superclass= paste(jj), 
                           Tissue=paste(ii))
  
  print(PVAL)                  # THE PVALUE IS PRINTED

  PVALUE_Tissues_GREATER0$P.VALUE[i]<-paste(PVAL)              # THE PVALUE IS ADDED INTO THE DF CREATED ABOVE

}

# THE DATAFRAME IS #ViewED
#View(PVALUE_Tissues_GREATER0)



# > 0.1
PVALUE_Tissues_0.1 <- data.frame(all_disease_genes_pheno_Tissues$GTEX.tissues,
                                      all_disease_genes_pheno_Tissues$HPO.superclass.description,
                                      P.VALUE=vector(length=55)) 

for(i in 1:nrow(all_disease_genes_pheno_Tissues)) {
  ii <- all_disease_genes_pheno_Tissues[i,1]     #FIRST COLUMN
  jj <-all_disease_genes_pheno_Tissues[i,2]      #SECOND COLUMN
  
  print(ii)                                     # TO CHECK IF THE CORRECT ROW IS BEING TAKEN
  print(jj)
  
            
  
  ## THE OUTCOME OF THE FUNCTION IS SAVED AS PVAL WHICH IS THE FISHER TEST PVALUE 
  PVAL<-TISSUE_FISHER_TEST(df=human_genes_TPM_0.1_pheno, # THE FUNCTION NEED THE DF THAT HAS TEH GENE GENE EXPRESSIONS AT SPECIFIC THRESHOLDS
                           superclass= paste(jj), 
                           Tissue=paste(ii))
  
  print(PVAL)                  # THE PVALUE IS PRINTED
  
  PVALUE_Tissues_0.1$P.VALUE[i]<-paste(PVAL)              # THE PVALUE IS ADDED INTO THE DF CREATED ABOVE
  
}

# THE DATAFRAME IS #ViewED
#View(PVALUE_Tissues_0.1)

# > 1
PVALUE_Tissues_1 <- data.frame(all_disease_genes_pheno_Tissues$GTEX.tissues,
                                 all_disease_genes_pheno_Tissues$HPO.superclass.description,
                                 P.VALUE=vector(length=55)) 

for(i in 1:nrow(all_disease_genes_pheno_Tissues)) {
  ii <- all_disease_genes_pheno_Tissues[i,1]     #FIRST COLUMN
  jj <-all_disease_genes_pheno_Tissues[i,2]      #SECOND COLUMN
  
  print(ii)                                     # TO CHECK IF THE CORRECT ROW IS BEING TAKEN
  print(jj)
        
  
  ## THE OUTCOME OF THE FUNCTION IS SAVED AS PVAL WHICH IS THE FISHER TEST PVALUE 
  PVAL<-TISSUE_FISHER_TEST(df=human_genes_TPM_1_pheno,# THE FUNCTION NEED THE DF THAT HAS TEH GENE GENE EXPRESSIONS AT SPECIFIC THRESHOLDS
                           superclass= paste(jj), 
                           Tissue=paste(ii))
  
  print(PVAL)                  # THE PVALUE IS PRINTED
  
  PVALUE_Tissues_1$P.VALUE[i]<-paste(PVAL)              # THE PVALUE IS ADDED INTO THE DF CREATED ABOVE
  
}

# THE DATAFRAME IS #ViewED
#View(PVALUE_Tissues_1)









#######################~~~~~~~~~~~~~~~~~~~~~~~~~Gene Level Analysis ~~~~~~~~~~~~~~~~~~~~~~~~###################


