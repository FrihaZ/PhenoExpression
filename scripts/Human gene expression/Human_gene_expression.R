
###############################################################################################################################
###############################################################################################################################

### Project: PhenoExpression ##################################################################################################
### Script: Human_gene_expression.R ###################################################################################################
### Purpose: To clean the Human gene expression data and compare to mice data ###############################################################################
### Author: Friha Zafar ####################################################################################################
### Date: 01/04/2019 ##########################################################################################################

###############################################################################################################################
###############################################################################################################################

## load packages ##############################################################################################################

library(dplyr);library(tidyr);library(stringr);library(readr);library(tidyverse);

################################################################################################################################
################################################################################################################################

## import files ################################################################################################################

human_genes <- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/GTEX.TPM.gz",delim = "\t", skip = 2)


human_genes_filtered <-human_genes %>%
  select(Description)%>%
  distinct(Description) %>%
  pull(Description)           # to make it into a vector


View(human_genes)  




################################################################################################################################
################################################################################################################################

gene_protein <- read_delim("D:/MSC RESEARCH PROJECT/Gene-symbol-checker/gene_with_protein_product.txt",delim = "\t")


View(gene_protein)


# Retrieve hgnc.checker function


source(file="D:/MSC RESEARCH PROJECT/Gene-symbol-checker/hgnc_symbol_checker.R")



approved_hgnc.id <- hgnc.checker(human_genes_filtered, gene_protein)

View(approved_hgnc.id)



################################################################################################################################
################################################################################################################################

# Remove rows that are empty: filter by  '-' and 'Ambiguous.Symbol'

approved_hgnc.id <- approved_hgnc.id %>% 
  filter(HGNC.ID  != "-") %>%  #Remove ambiguous expression data
  filter(Type != "Ambiguous.Symbol") 

View(approved_hgnc.id)


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

View(human_genes)



# Convert gene expression values (TPM) to binary variable. Three different thresholds:
# TPM > 0 (yes,no)

#str(human_genes)

human_genes_TPM_greater0<- human_genes


human_genes_TPM_greater0[4:56 ][ human_genes_TPM_greater0[4:56 ] > 0 ] <- "Yes"

human_genes_TPM_greater0[4:56 ][ human_genes_TPM_greater0[4:56 ] <= 0 ] <- "No"

View(human_genes_TPM_greater0)


# TPM > 0.1 (yes, no)

human_genes_TPM_0.1<- human_genes

human_genes_TPM_0.1[4:56 ][ human_genes_TPM_0.1[4:56 ] >= 0.1 ] <- "Yes"
human_genes_TPM_0.1[4:56 ][ human_genes_TPM_0.1[4:56 ] < 0.1 ] <- "No"




View(human_genes_TPM_0.1)

# TPM> 1 (yes, no)

human_genes_TPM_1<- human_genes

human_genes_TPM_1[4:56 ][ human_genes_TPM_1[4:56 ] >= 1 ] <- "Yes"
human_genes_TPM_1[4:56 ][ human_genes_TPM_1[4:56 ] < 1 ] <- "No"



View(human_genes_TPM_1)





################################################################################################################################
################################################################################################################################


# Convert gene symbols in HPO_phenotypes to hgnc.id


###MANUALLY REMOVE THE <TAB> AND INSERT REAL TABS IN THE HPO_phenotypes.txt FILE'S COLUMN NAMES!!!

HPO_phenotypes <- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/HPO_phenotypes.txt",delim = "\t")

# The dataframe is turned into a vector and only the `entrez-gene-symbol`is used 
HPO_phenotypes_filtered <-HPO_phenotypes %>%
  select(`entrez-gene-symbol`)%>%
  distinct(`entrez-gene-symbol`) %>%
  pull(`entrez-gene-symbol`)           # to make it into a vector


View(HPO_phenotypes_filtered)  



# Use the hgnc.checker function

HPO_phenotypes_hgnc.id <- hgnc.checker(HPO_phenotypes_filtered, gene_protein)



View(HPO_phenotypes_hgnc.id)



################################################################################################################################
################################################################################################################################

# Remove rows that are empty: filter by  '-' and 'Notfound.ProteinCoding.Symbol'

HPO_phenotypes_hgnc.id <- HPO_phenotypes_hgnc.id %>% 
  filter(HGNC.ID  != "-") %>%  #Remove ambiguous expression data
  filter(Type != "Notfound.ProteinCoding.Symbol" | Type != "Ambiguous.Symbol") 

View(HPO_phenotypes_hgnc.id)


################################################################################################################################
################################################################################################################################


## Join  the HPO_phenotypes_hgnc.id with the HPO_phenotypes table 


HPO_phenotypes <-left_join(HPO_phenotypes, HPO_phenotypes_hgnc.id, by=c(`entrez-gene-symbol` = "Gene.Symbol"))  

View(HPO_phenotypes)


################################################################################################################################
################################################################################################################################


# Map the HPO-Term-ID to the top level nodes using hpo.ancestor.nodes.txt file

hpo.ancestor.nodes <- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/hpo.ancestor.nodes.txt",delim = "\t")

#View(hpo.ancestor.nodes)
hpo.toplevels.phenotypic.abnormalities.only <- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/hpo.toplevels.phenotypic.abnormalities.only.txt",delim = "\t")

#View(hpo.toplevels.phenotypic.abnormalities.only)


hpo.ancestor.nodes <- left_join(hpo.ancestor.nodes, hpo.toplevels.phenotypic.abnormalities.only, by= c("hpo.ancestors"="hpo.term"))


#View(hpo.ancestor.nodes)
hpo.ancestor.nodes<-hpo.ancestor.nodes%>%
  select(hpo.term, hpo.ancestors, hpo.description)%>%
  drop_na()%>%
  distinct()


View(hpo.ancestor.nodes)

# Join the two tables 


HPO_phenotypes <-left_join(HPO_phenotypes, hpo.ancestor.nodes, by=c("HPO-Term-ID" = "hpo.term"))  



#REMOVED NAs
HPO_phenotypes<-HPO_phenotypes %>%
  select_all()%>%
  drop_na()%>%       
  distinct()

View(HPO_phenotypes)

#####################################################################################################################

# Map genes with gene expression to the HPO annotations (top level)


human_genes_TPM_greater0_phenotype<- left_join(human_genes_TPM_greater0, HPO_phenotypes, by=c("HGNC.ID"= "HGNC.ID"))



human_genes_TPM_greater0_pheno <-human_genes_TPM_greater0_phenotype%>%
  select_all()%>%
  drop_na() %>%  #removed NAs that for the genes that didn't have a hgncd.id 
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`)%>%   #removed the entrez-gene-id etc
  distinct()


View(human_genes_TPM_greater0_pheno)





# Remove all genes that are disease related, and keep all those that are not disease related. 

human_genes_greater0_pheno_No_disease<-human_genes_TPM_greater0_phenotype[!(human_genes_TPM_greater0_phenotype$HGNC.ID %in% human_genes_TPM_greater0_pheno$HGNC.ID),]

human_genes_greater0_pheno_No_disease<-human_genes_greater0_pheno_No_disease%>%
  select_all()%>%
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`,
         -`HPO-Term-ID`,-`hpo.ancestors`, -`hpo.description` )%>%   #removed the entrez-gene-id etc
  distinct()

View(human_genes_greater0_pheno_No_disease)


########################################################################################################

human_genes_TPM_0.1_phenotype<- left_join(human_genes_TPM_0.1, HPO_phenotypes, by=c("HGNC.ID"= "HGNC.ID"))


# DISEASE RELATED GENES
human_genes_TPM_0.1_pheno <-human_genes_TPM_0.1_phenotype%>%
  select_all()%>%
  drop_na() %>%  #removed NAs that for the genes that didn't have a hgncd.id 
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`)%>%   #removed the entrez-gene-id etc
  distinct()


View(human_genes_TPM_0.1_pheno)



# Remove all genes that are disease related, and keep all those that are not disease related. 

human_genes_TPM_0.1_pheno_No_disease<-human_genes_TPM_0.1_phenotype[!(human_genes_TPM_0.1_phenotype$HGNC.ID %in% human_genes_TPM_0.1_pheno$HGNC.ID),]

human_genes_TPM_0.1_pheno_No_disease<-human_genes_TPM_0.1_pheno_No_disease%>%
  select_all()%>%
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`,
         -`HPO-Term-ID`,-`hpo.ancestors`, -`hpo.description` )%>%   #removed the entrez-gene-id etc
  distinct()

View(human_genes_TPM_0.1_pheno_No_disease)



#########################################################################################################

# DISEASE RELATED GENES

human_genes_TPM_1_phenotype<- left_join(human_genes_TPM_1, HPO_phenotypes, by=c("HGNC.ID"= "HGNC.ID"))

View(human_genes_TPM_1_phenotype)

human_genes_TPM_1_pheno <-human_genes_TPM_1_phenotype%>%
  select_all()%>%
  drop_na() %>%  #removed NAs that for the genes that didn't have a hgncd.id 
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`)%>%   #removed the entrez-gene-id etc
  distinct()


View(human_genes_TPM_1_pheno)

str(human_genes_TPM_1_pheno)


# Remove all genes that are disease related, and keep all those that are not disease related. 

human_genes_TPM_1_pheno_No_disease<-human_genes_TPM_1_phenotype[!(human_genes_TPM_1_phenotype$HGNC.ID %in% human_genes_TPM_1_pheno$HGNC.ID),]

human_genes_TPM_1_pheno_No_disease<-human_genes_TPM_1_pheno_No_disease%>%
  select_all()%>%
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`,
         -`HPO-Term-ID`,-`hpo.ancestors`, -`hpo.description` )%>%   #removed the entrez-gene-id etc
  distinct()

View(human_genes_TPM_1_pheno_No_disease)




#########################################################################################################
#########################################################################################################

## ALL NON-DISEASE GENES!

# Remove the phenotype description ,hpo.ancestors HPO.Term and count the number of Yes for each row (gene) for 
# each threshold


#########################################################################################################

# For expressions >0, calculate the number of Yes and No

human_genes_greater0_pheno_No_disease$NUMBER_OF_YES <- rowSums(human_genes_greater0_pheno_No_disease == "Yes")
human_genes_greater0_pheno_No_disease$NUMBER_OF_NO <- rowSums(human_genes_greater0_pheno_No_disease == "No")

human_genes_greater0_pheno_No_disease<-human_genes_greater0_pheno_No_disease %>%
  select_all()%>%
  distinct()

View(human_genes_greater0_pheno_No_disease)




#####~~~~~~!!!!!!!!!!!!!!IMPORTANT KEEP FOR DISSERTATION!!!!!!!!!!!~~~~~~##############

# PLOT THE NUMBER OF GENES THAT ARE EXPRESSED IN EACH TISSUE, where gene expression is >0 and not disease specific

library(Hmisc)
library(ggplot2)


human_genes_greater0_pheno_No_disease_top30 <- filter(human_genes_greater0_pheno_No_disease, human_genes_greater0_pheno_No_disease$NUMBER_OF_YES >=30)

dev.new(width=60, height=30) # to open in a new window

p_greater0_genes_ND<-ggplot(data=human_genes_greater0_pheno_No_disease_top30
                             , aes(x=human_genes_greater0_pheno_No_disease_top30$HGNC.ID,
                                   y=human_genes_greater0_pheno_No_disease_top30$NUMBER_OF_YES,
                                   fill=human_genes_greater0_pheno_No_disease_top30$NUMBER_OF_YES))
p_greater0_genes_ND %>%
  
  + labs(x= "HGNC.ID", y="Number of Tissues the Gene is Expressed in", title= "The Number of Genes Expressed in Different Tissues.",
         tag="D",
         caption = "The data shows the number of tissues in which the genes (HGNC.ID) are expressed in at threshold of >0 in different human tissues, where the genes are not known to be associated with disease.") %>%
  + scale_fill_gradient(low = "blue", 
                        high = "red")%>%
  + guides(fill=guide_legend(title="Number of Genes Expressed:  ")) %>%
  + theme(plot.caption=element_text(face = "italic", size=12, hjust = 0),
          text = element_text(size=12),
          legend.position = "bottom",legend.direction = "horizontal",
          legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +geom_bar(stat = "identity")    # to create a stacked barchart












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
View(human_genes_TPM_greater0_tissue_ND)



#####~~~~~~!!!!!!!!!!!!!!IMPORTANT KEEP FOR DISSERTATION!!!!!!!!!!!~~~~~~##############

# PLOT THE NUMBER OF GENES THAT ARE EXPRESSED IN EACH TISSUE, where gene expression is >0 and not disease specific

library(Hmisc)
library(ggplot2)


dev.new(width=60, height=30) # to open in a new window

p_greater0_tissue_ND<-ggplot(data=human_genes_TPM_greater0_tissue_ND
          , aes(x=reorder(human_genes_TPM_greater0_tissue_ND$Tissue, -human_genes_TPM_greater0_tissue_ND$Number_of_Yes),
                y=reorder(human_genes_TPM_greater0_tissue_ND$Number_of_Yes, -human_genes_TPM_greater0_tissue_ND$Tissue),
                fill=human_genes_TPM_greater0_tissue_ND$Number_of_Yes))
p_greater0_tissue_ND %>%
  
  + labs(x= "Tissue", y="Number of Non-Disease Genes Expressed at a Threshold >0", title= "The Number of Genes Expressed in Different Tissues.",
         tag="A",
         caption = "The data shows the number of genes expressed at a threshold of >0 in different human tissues, where the genes are not known to be associated with disease.") %>%
  + scale_fill_gradient(low = "blue", 
                        high = "red")%>%
  + guides(fill=guide_legend(title="Number of Genes Expressed:  ")) %>%
  + theme(plot.caption=element_text(face = "italic", size=12, hjust = 0),
          text = element_text(size=12),
          legend.position = "bottom",legend.direction = "horizontal",
          legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +geom_bar(stat = "identity")%>%    # to create a stacked barchart
  + coord_flip()






# # Make Violin plot
# 
# dev.new(width=60, height=30) # to open in a new window
# 
# p_greater0_tissue_violin<-ggplot(data=human_genes_TPM_greater0_tissue_30
#                                  , aes(x=human_genes_TPM_greater0_tissue_30$Tissue,
#                                        y=human_genes_TPM_greater0_tissue_30$Number_of_Yes, fill=human_genes_TPM_greater0_tissue_30$Number_of_Yes))
# 
# p_greater0_tissue_violin%>%
#   + geom_violin(trim=TRUE) %>%
#   + geom_boxplot(width=20)%>%
#   
#   + labs(x= "Tissue", y="Number of Genes Expressed at a Threshold >0 (30)", title= "The Number of Genes Expressed in Different Tissues.",
#          tag="A",
#          caption = "The data shows the number of genes expressed at a threshold of >0 in different human tissues.(30)") %>%
#   + scale_fill_gradient(low = "white", 
#                         high = "blue")%>%
#   + guides(fill=guide_legend(title="Number of Genes Expressed:  ")) %>%
#   + theme(plot.caption=element_text(face = "italic", size=12, hjust = 0),
#           text = element_text(size=12),
#           legend.position = "bottom",legend.direction = "horizontal",
#           legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
#           axis.text.x = element_text(angle = 90, hjust = 1),
#           panel.border = element_rect(colour = "black", fill=NA, size=1))










#########################################################################################################

# For expression greater than 0.1



human_genes_TPM_0.1_pheno_No_disease$NUMBER_OF_YES <- rowSums(human_genes_TPM_0.1_pheno_No_disease == "Yes")
human_genes_TPM_0.1_pheno_No_disease$NUMBER_OF_NO <- rowSums(human_genes_TPM_0.1_pheno_No_disease == "No")

human_genes_TPM_0.1_pheno_No_disease<-human_genes_TPM_0.1_pheno_No_disease %>%
  select_all()%>%
  distinct()

View(human_genes_TPM_0.1_pheno_No_disease)



## COUNT THE NUMBER OF GENES EXPRESSED FOR EACH TISSUE

human_genes_TPM_0.1_tissue_ND<-human_genes_TPM_0.1_pheno_No_disease%>%
  select_all()%>%
  select(-NUMBER_OF_NO, -NUMBER_OF_YES, -`HGNC.ID`, -`Description`, -`Type.x`)

human_genes_TPM_0.1_tissue_ND <- colSums(human_genes_TPM_0.1_tissue_ND == "Yes")


View(human_genes_TPM_0.1_tissue_ND)


human_genes_TPM_0.1_tissue_ND<-data.frame(human_genes_TPM_0.1_tissue_ND)

human_genes_TPM_0.1_tissue_ND<-setNames(cbind(rownames(human_genes_TPM_0.1_tissue_ND), 
                                              human_genes_TPM_0.1_tissue_ND, row.names = NULL), 
                                          c("Tissue", "Number_of_Yes"))
View(human_genes_TPM_0.1_tissue_ND)



#####~~~~~~!!!!!!!!!!!!!!IMPORTANT KEEP FOR DISSERTATION!!!!!!!!!!!~~~~~~##############

# PLOT THE NUMBER OF GENES THAT ARE EXPRESSED IN EACH TISSUE, where gene expression is >/=0.1 and not disease specific

library(ggplot2)


dev.new(width=60, height=30) # to open in a new window

p_TPM_0.1_tissue_ND<-ggplot(data=human_genes_TPM_0.1_tissue_ND
                          , aes(x=reorder(human_genes_TPM_0.1_tissue_ND$Tissue, 
                                          -human_genes_TPM_0.1_tissue_ND$Number_of_Yes),
                                y=reorder(human_genes_TPM_0.1_tissue_ND$Number_of_Yes, 
                                          -human_genes_TPM_0.1_tissue_ND$Tissue),
                                fill=human_genes_TPM_0.1_tissue_ND$Number_of_Yes))
p_TPM_0.1_tissue_ND %>%
  
  
  + labs(x= "Tissue", y="Number of Non-Disease Genes Expressed at a Threshold >/=0.1", title= "The Number of Genes Expressed in Different Tissues.",
         tag="B",
         caption = "The data shows the number of genes expressed at a threshold of >/=0.1 in different human tissues, where the genes are not known to be associated with disease.") %>%
  + scale_fill_gradient(low = "blue", 
                        high = "red")%>% 
  + guides(fill=guide_legend(title="Number of Genes Expressed:  ")) %>%
  + theme(plot.caption=element_text(face = "italic", size=12, hjust = 0),
          text = element_text(size=12),
          legend.position = "bottom",legend.direction = "horizontal",
          legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +geom_bar(stat = "identity")%>%    # to create a stacked barchart
  + coord_flip()











#########################################################################################################

# For expression greater than 1



human_genes_TPM_1_pheno_No_disease$NUMBER_OF_YES <- rowSums(human_genes_TPM_1_pheno_No_disease == "Yes")
human_genes_TPM_1_pheno_No_disease$NUMBER_OF_NO <- rowSums(human_genes_TPM_1_pheno_No_disease == "No")

human_genes_TPM_1_pheno_No_disease<-human_genes_TPM_1_pheno_No_disease %>%
  select_all()%>%
  distinct()

View(human_genes_TPM_1_pheno_No_disease)
## COUNT THE NUMBER OF GENES EXPRESSED FOR EACH TISSUE

human_genes_TPM_1_tissue_ND<-human_genes_TPM_1_pheno_No_disease%>%
  select_all()%>%
  select(-NUMBER_OF_NO, -NUMBER_OF_YES, -`HGNC.ID`, -`Description`, -`Type.x`)

human_genes_TPM_1_tissue_ND <- colSums(human_genes_TPM_1_tissue_ND == "Yes")


View(human_genes_TPM_1_tissue_ND)



human_genes_TPM_1_tissue_ND<-data.frame(human_genes_TPM_1_tissue_ND)

human_genes_TPM_1_tissue_ND<-setNames(cbind(rownames(human_genes_TPM_1_tissue_ND), 
                                         human_genes_TPM_1_tissue_ND, row.names = NULL), 
                                     c("Tissue", "Number_of_Yes"))
View(human_genes_TPM_1_tissue_ND)



#####~~~~~~!!!!!!!!!!!!!!IMPORTANT KEEP FOR DISSERTATION!!!!!!!!!!!~~~~~~##############

# PLOT THE NUMBER OF GENES THAT ARE EXPRESSED IN EACH TISSUE, where gene expression is >/=0.1 and not disease specific

library(ggplot2)


dev.new(width=60, height=30) # to open in a new window

p_TPM_1_tissue_ND<-ggplot(data=human_genes_TPM_1_tissue_ND
                         , aes(x=reorder(human_genes_TPM_1_tissue_ND$Tissue, -human_genes_TPM_1_tissue_ND$Number_of_Yes),
                               y=reorder(human_genes_TPM_1_tissue_ND$Number_of_Yes, -human_genes_TPM_1_tissue_ND$Tissue),
                               fill=human_genes_TPM_1_tissue_ND$Number_of_Yes))
p_TPM_1_tissue_ND %>%
  
  + labs(x= "Tissue", y="Number of Non-Disease Genes Expressed at a Threshold >/= 1", title= "The Number of Genes Expressed in Different Tissues.",
         tag="C",
         caption = "The data shows the number of genes expressed at a threshold of >/= 1 in different human tissues, where the genes are not known to be associated with disease.") %>%
  + scale_fill_gradient(low = "blue", 
                        high = "red")%>% 
  + guides(fill=guide_legend(title="Number of Genes Expressed:  ")) %>%
  + theme(plot.caption=element_text(face = "italic", size=12, hjust = 0),
          text = element_text(size=12),
          legend.position = "bottom",legend.direction = "horizontal",
          legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +geom_bar(stat = "identity")%>%    # to create a stacked barchart
  + coord_flip()





















#   distinct(`Adipose - Subcutaneous`, `Adipose - Visceral (Omentum)`,
#            `Adrenal Gland`,`Artery - Aorta`,`Artery - Coronary`,`Artery - Tibial`,
#            `Bladder`,`Brain - Amygdala`,`Brain - Anterior cingulate cortex (BA24)`,
#            `Brain - Caudate (basal ganglia)`,`Brain - Cerebellar Hemisphere`,
#            `Brain - Cerebellum`,`Brain - Cortex`,`Brain - Frontal Cortex (BA9)`,
#            `Brain - Hippocampus`,`Brain - Hypothalamus`,`Brain - Nucleus accumbens (basal ganglia)`,
#            `Brain - Putamen (basal ganglia)`,`Brain - Spinal cord (cervical c-1)`,
#            `Brain - Substantia nigra`,`Breast - Mammary Tissue`,`Cells - EBV-transformed lymphocytes`,
#            `Cells - Transformed fibroblasts`,`Cervix - Ectocervix`,`Cervix - Endocervix`, `Colon - Sigmoid`,
#            `Colon - Transverse`,`Esophagus - Gastroesophageal Junction`,`Esophagus - Mucosa`,
#            `Esophagus - Muscularis`,`Fallopian Tube`,`Heart - Atrial Appendage`,`Heart - Left Ventricle`,
#            `Kidney - Cortex`,`Liver`,`Lung`,`Minor Salivary Gland`,`Muscle - Skeletal`,`Nerve - Tibial`,
#            `Ovary`,`Pancreas`,`Pituitary`,`Prostate`,`Skin - Not Sun Exposed (Suprapubic)`,
#            `Skin - Sun Exposed (Lower leg)`, `Small Intestine - Terminal Ileum`,`Spleen`,
#            `Stomach`,`Testis`,`Thyroid`,`Uterus`,`Vagina`,`Whole Blood`)


