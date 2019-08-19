
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

library(dplyr);library(tidyr);library(stringr);library(readr);library(tidyverse);library(ggpubr)

################################################################################################################################
################################################################################################################################

## import files ################################################################################################################

human_genes <- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/GTEX.TPM.gz",delim = "\t", skip = 2)

gene_protein <- read_delim("D:/MSC RESEARCH PROJECT/Gene-symbol-checker/gene_with_protein_product.txt",delim = "\t")

#Retrieve hgnc.checker function 

source(file="D:/MSC RESEARCH PROJECT/Gene-symbol-checker/hgnc_symbol_checker.R")


hpo.ancestor.nodes <- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/hpo.ancestor.nodes.txt",delim = "\t")

hpo.toplevels.phenotypic.abnormalities.only <- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/hpo.toplevels.phenotypic.abnormalities.only.txt",delim = "\t")

#!!!!!!!!!MANUALLY REMOVE THE <TAB> AND INSERT REAL TABS IN THE HPO_phenotypes.txt FILE'S COLUMN NAMES!!!

HPO_phenotypes <- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/HPO_phenotypes.txt",delim = "\t")


GTEX.Tissues.to.HPO <- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/GTEX.Tissues.to.HPO.txt",delim = "\t")

#######################################################################################################

###############################################################################################################################################

#########~~~~~~~~~~~~~~~~~~~CLEAN GTEX HUMAN GENE EXPRESSION DF ~~~~~~~~~~~~######################################

human_genes_filtered <-human_genes %>%
  select(Description)%>% # only the gene names are selected as these will be taken to find the hgnc.id using the hgnc.checker. 
  distinct(Description) %>%
  pull(Description)           # to make it into a vector for the hgnc.checker to work


View(human_genes_filtered)   


## Use the hgnc.checker function to vonvert gene names TO hgnc.ids

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


## Join  the approved_hgnc.id with the human_genes table to get the gene expressions 


human_genes <-left_join(human_genes, approved_hgnc.id, by=c("Description" = "Gene.Symbol"))  

##########################################################################################
## KEEP THE ENTREZ IDS  FOR GENE ENRICHMENT ANALYSIS

human_genes_with_GENE_ID <-human_genes%>%
  select(HGNC.ID)%>% 
  drop_na() %>%  #removed NAs that for the genes that didn't have a hgncd.id 
  distinct()

# open gene_protein

gene_protein_ENTREZ<-gene_protein%>%
  select(entrez_id, hgnc_id)%>%
  distinct()


#View(gene_protein_ENTREZ)

names(gene_protein_ENTREZ)[names(gene_protein_ENTREZ) == 'hgnc_id'] <- 'HGNC.ID'

human_genes_with_GENE_ID_entrez<-right_join(gene_protein_ENTREZ, human_genes_with_GENE_ID)


#View(human_genes_with_GENE_ID_entrez)


#!!!!! UNHASH TO SAVE FILE!!!!!!!

# write.csv(human_genes_with_GENE_ID_entrez,'./Output_Files/Human/Expression_1/human_genes_with_GENE_ID_entrez.csv')


##################################################################################################################################
##########################################################################################################


######~~~~~~~~ Clean human_genes by removing the gene_id ~~~~~~~~~~~~~~#########################

human_genes <-human_genes%>%
  select(HGNC.ID,gene_id, Description, Type,`Adipose - Subcutaneous`:`Whole Blood`)%>% 
  drop_na() %>%  #removed NAs that for the genes that didn't have a hgncd.id 
  select(-gene_id)%>%   #removed the gene_id 
  distinct()
options("scipen"=100)

#!!!!! UNHASH TO SAVE FILE!!!!!!!

# write.csv(human_genes,'./Output_Files/Human/human_genes.csv')


##################################################################################################################################
###############################################################################################################

####~~~~~~~~~~~~~ Convert gene expression values (TPM) to binary variable. Three different thresholds: ~~~~~~~~~~##################

# TPM > 0 (yes,no)
# TPM > 0.1 (yes, no)
# TPM> 1 (yes, no)

## GENES WITH NO GENE EXPRESSION (TPM =0) ###########################################################
human_genes_TPM_None<- human_genes

human_genes_TPM_None[4:56 ][ human_genes_TPM_None[4:56 ] > 0 ] <- "No"

human_genes_TPM_None[4:56 ][ human_genes_TPM_None[4:56 ] <= 0 ] <- "Yes"

View(human_genes_TPM_None)

#!!!!! UNHASH TO SAVE FILE!!!!!!!


# write.csv(human_genes_TPM_None,'./Output_Files/Human/No_Expression/human_genes_TPM_None.csv')

#####################################################################################################

# GENES WITHTPM > 0 (yes,no)##################################################################

human_genes_TPM_greater0<- human_genes


human_genes_TPM_greater0[4:56 ][ human_genes_TPM_greater0[4:56 ] > 0 ] <- "Yes"

human_genes_TPM_greater0[4:56 ][ human_genes_TPM_greater0[4:56 ] <= 0 ] <- "No"

#!!!!! UNHASH TO SAVE FILE!!!!!!!

# write.csv(human_genes_TPM_greater0,'./Output_Files/Human/Expression_greater_0/human_genes_TPM_greater0.csv')


##View(human_genes_TPM_greater0)

###############################################################################################################


# GENES WITH TPM > 0.1 (yes, no) ########################################################################

human_genes_TPM_0.1<- human_genes

human_genes_TPM_0.1[4:56 ][ human_genes_TPM_0.1[4:56 ] >= 0.1 ] <- "Yes"
human_genes_TPM_0.1[4:56 ][ human_genes_TPM_0.1[4:56 ] < 0.1 ] <- "No"


#!!!!! UNHASH TO SAVE FILE!!!!!!!

# write.csv(human_genes_TPM_0.1,'./Output_Files/Human/Expression_0.1/human_genes_TPM_0.1.csv')
#View(human_genes_TPM_0.1)

#########################################################################################


# GENES WITH TPM > 1 (yes, no) ###############################################################

human_genes_TPM_1<- human_genes

human_genes_TPM_1[4:56 ][ human_genes_TPM_1[4:56 ] >= 1 ] <- "Yes"
human_genes_TPM_1[4:56 ][ human_genes_TPM_1[4:56 ] < 1 ] <- "No"


#!!!!! UNHASH TO SAVE FILE!!!!!!!

# write.csv(human_genes_TPM_1,'./Output_Files/Human/Expression_1/human_genes_TPM_1.csv')

#View(human_genes_TPM_1)





################################################################################################################################
################################################################################################################################


# Convert gene symbols in HPO_phenotypes to hgnc.id


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


hpo.ancestor.nodes <- left_join(hpo.ancestor.nodes, hpo.toplevels.phenotypic.abnormalities.only, by= c("hpo.ancestors"="hpo.term"))


##View(hpo.ancestor.nodes)
hpo.ancestor.nodes<-hpo.ancestor.nodes%>%
  select(hpo.term, hpo.ancestors, hpo.description)%>%
  drop_na()%>%
  distinct()


#View(hpo.ancestor.nodes)

#################################################################################################
####################################################################################################

# Join the two DFs 


HPO_phenotypes <-left_join(HPO_phenotypes, hpo.ancestor.nodes, by=c("HPO-Term-ID" = "hpo.term"))  



#REMOVED NAs
HPO_phenotypes<-HPO_phenotypes %>%
  select_all()%>%
  drop_na()%>%       
  distinct()

#View(HPO_phenotypes)

#####################################################################################################
#####################################################################################################################

####~~~~~~~~~~~~~~~~~ Map genes with gene expression to the HPO annotations (top level)~~~~~~~~~~~~~~~~~~##############




## To all disease associated genes, irrespective of threshold

all_disease_genes_phenotype <-left_join(human_genes, HPO_phenotypes, by=c("HGNC.ID"= "HGNC.ID"))



all_disease_genes_pheno <-all_disease_genes_phenotype%>%
  select_all()%>%
  drop_na() %>%  #removed NAs that for the genes that didn't have a hgncd.id 
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`)%>%   #removed the entrez-gene-id etc
  distinct()



#!!!!! UNHASH TO SAVE FILE!!!!!!!

# write.csv(all_disease_genes_pheno,'./Output_Files/Human/all_disease_genes_pheno.csv')


# count the number of genes that are disease phenotype linked
all_disease_genes_pheno_count<-all_disease_genes_pheno%>%
  select(HGNC.ID)%>%
  distinct()


#!!!!! UNHASH TO SAVE FILE!!!!!!!

# write.csv(all_disease_genes_pheno_count,'./Output_Files/Human/all_disease_genes_pheno_count.csv')


#################################################################


# Remove all genes that are disease related, and keep all those that are not disease related. 

all_genes_pheno_No_disease<-all_disease_genes_phenotype[!(all_disease_genes_phenotype$HGNC.ID %in% all_disease_genes_pheno$HGNC.ID),]

all_genes_pheno_No_disease<-all_genes_pheno_No_disease%>%
  select_all()%>%
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`,
         -`HPO-Term-ID`,-`hpo.ancestors`, -`hpo.description` )%>%   #removed the entrez-gene-id etc
  distinct()

#View(all_genes_pheno_No_disease)

#!!!!! UNHASH TO SAVE FILE!!!!!!!

# write.csv(all_genes_pheno_No_disease,'./Output_Files/Human/all_genes_pheno_No_disease.csv')


# count the number of genes that are NOT disease phenotype linked
all_genes_pheno_No_disease_count<-all_genes_pheno_No_disease%>%
  select(HGNC.ID)%>%
  distinct()
#View(all_genes_pheno_No_disease_count)


# write.csv(all_genes_pheno_No_disease_count,'./Output_Files/Human/all_genes_pheno_No_disease_count.csv')


#nrow(all_genes_pheno_No_disease_count)

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

#!!!!! UNHASH TO SAVE FILE!!!!!!!


# write.csv(human_genes_TPM_None_pheno,
#           './Output_Files/Human/No_Expression/human_genes_TPM_None_pheno.csv')


######################################################################################################

## EXPRESSION:- PHENOTYPE:-
# Remove all genes that are disease related, and keep all those that are not disease related. 

human_genes_TPM_None_pheno_No_disease<-human_genes_TPM_None_phenotype[!(human_genes_TPM_None_phenotype$HGNC.ID %in% human_genes_TPM_None_pheno$HGNC.ID),]

human_genes_TPM_None_pheno_No_disease<-human_genes_TPM_None_pheno_No_disease%>%
  select_all()%>%
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`,
         -`HPO-Term-ID`,-`hpo.ancestors`, -`hpo.description` )%>%   #removed the entrez-gene-id etc
  distinct()

#View(human_genes_TPM_None_pheno_No_disease)

#!!!!! UNHASH TO SAVE FILE!!!!!!!

 
# write.csv(human_genes_TPM_None_pheno_No_disease,
#           './Output_Files/Human/No_Expression/human_genes_TPM_None_pheno_No_disease.csv')



#human_genes_TPM_None


#######################################################################################

# DISEASE RELATED GENES
#EXPRESSION:+ PHENOTYPE:+


human_genes_TPM_greater0_phenotype<- left_join(human_genes_TPM_greater0, HPO_phenotypes, by=c("HGNC.ID"= "HGNC.ID"))



human_genes_TPM_greater0_pheno <-human_genes_TPM_greater0_phenotype%>%
  select_all()%>%
  drop_na() %>%  #removed NAs that for the genes that didn't have a hgncd.id 
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`)%>%   #removed the entrez-gene-id etc
  distinct()


#View(human_genes_TPM_greater0_pheno)

#!!!!! UNHASH TO SAVE FILE!!!!!!!

 
# write.csv(human_genes_TPM_greater0_pheno,
#           './Output_Files/Human/Expression_greater_0/human_genes_TPM_greater0_pheno.csv')



# Remove all genes that are disease related, and keep all those that are not disease related. 

human_genes_greater0_pheno_No_disease<-human_genes_TPM_greater0_phenotype[!(human_genes_TPM_greater0_phenotype$HGNC.ID %in% human_genes_TPM_greater0_pheno$HGNC.ID),]

human_genes_greater0_pheno_No_disease<-human_genes_greater0_pheno_No_disease%>%
  select_all()%>%
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`,
         -`HPO-Term-ID`,-`hpo.ancestors`, -`hpo.description` )%>%   #removed the entrez-gene-id etc
  distinct()

#View(human_genes_greater0_pheno_No_disease)

#!!!!! UNHASH TO SAVE FILE!!!!!!!


# write.csv(human_genes_greater0_pheno_No_disease,
#           './Output_Files/Human/Expression_greater_0/human_genes_greater0_pheno_No_disease.csv')



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

#!!!!! UNHASH TO SAVE FILE!!!!!!!

# write.csv(human_genes_TPM_0.1_pheno,
#           './Output_Files/Human/Expression_0.1/human_genes_TPM_0.1_pheno.csv')
# 

# Remove all genes that are disease related, and keep all those that are not disease related. 

human_genes_TPM_0.1_pheno_No_disease<-human_genes_TPM_0.1_phenotype[!(human_genes_TPM_0.1_phenotype$HGNC.ID %in% human_genes_TPM_0.1_pheno$HGNC.ID),]

human_genes_TPM_0.1_pheno_No_disease<-human_genes_TPM_0.1_pheno_No_disease%>%
  select_all()%>%
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`,
         -`HPO-Term-ID`,-`hpo.ancestors`, -`hpo.description` )%>%   #removed the entrez-gene-id etc
  distinct()

#View(human_genes_TPM_0.1_pheno_No_disease)

#!!!!! UNHASH TO SAVE FILE!!!!!!!

# write.csv(human_genes_TPM_0.1_pheno_No_disease,
#           './Output_Files/Human/Expression_0.1/human_genes_TPM_0.1_pheno_No_disease.csv')



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

#!!!! UNHASH TO SAVE FILE !!!!!!!

# write.csv(human_genes_TPM_1_pheno,
#           './Output_Files/Human/Expression_1/human_genes_TPM_1_pheno.csv')


# Remove all genes that are disease related, and keep all those that are not disease related. 

human_genes_TPM_1_pheno_No_disease<-human_genes_TPM_1_phenotype[!(human_genes_TPM_1_phenotype$HGNC.ID %in% human_genes_TPM_1_pheno$HGNC.ID),]

human_genes_TPM_1_pheno_No_disease<-human_genes_TPM_1_pheno_No_disease%>%
  select_all()%>%
  select(-`#Format: entrez-gene-id`,-`entrez-gene-symbol` ,-`Type.y`, -`HPO-Term-Name`,
         -`HPO-Term-ID`,-`hpo.ancestors`, -`hpo.description` )%>%   #removed the entrez-gene-id etc
  distinct()

#View(human_genes_TPM_1_pheno_No_disease)

#!!!!! UNHASH TO SAVE FILE!!!!!!!
 
# write.csv(human_genes_TPM_1_pheno_No_disease,
#           './Output_Files/Human/Expression_1/human_genes_TPM_1_pheno_No_disease.csv')



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

# PLOT THE NUMBER OF GENES THAT ARE EXPRESSED IN EACH TISSUE, where gene expression is =0
#and not disease associated

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
  
  + coord_flip()

p_None_tissue_ND


###!!! UNHASH TO SAVE PLOT!!!!!
# 
# library(cowplot)
# save_plot("./Plots/Human/Tissue/p_None_tissue_ND.jpeg",
#           p_None_tissue_ND ,base_height= 10 ,base_aspect_ratio = 2) 








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
  
  + coord_flip()



###!!! UNHASH TO SAVE PLOT!!!!!
# 
# library(cowplot)
# save_plot("./Plots/Human/Tissue/p_greater0_tissue_ND.png",
#           p_greater0_tissue_ND ,base_height= 10 ,base_aspect_ratio = 2) 












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
  
  + coord_flip()


###!!! UNHASH TO SAVE PLOT!!!!!
# 
# library(cowplot)
# save_plot("./Plots/Human/Tissue/p_TPM_0.1_tissue_ND.png",
#           p_TPM_0.1_tissue_ND ,base_height= 10 ,base_aspect_ratio = 2) 




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
  + coord_flip()


###!!! UNHASH TO SAVE PLOT!!!!!
# 
# library(cowplot)
# save_plot("./Plots/Human/Tissue/p_TPM_1_tissue_ND.png",
#           p_TPM_1_tissue_ND ,base_height= 10 ,base_aspect_ratio = 2) 




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
  
  + coord_flip()


###!!! UNHASH TO SAVE PLOT!!!!!
# 
# library(cowplot)
# save_plot("./Plots/Human/Tissue/p_TPM_None_tissue_D.png",
#           p_TPM_None_tissue_D ,base_height= 10 ,base_aspect_ratio = 2) 



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
  
  + coord_flip()


###!!! UNHASH TO SAVE PLOT!!!!!
# 
# library(cowplot)
# save_plot("./Plots/Human/Tissue/p_TPM_greater0_tissue_D.png",
#           p_TPM_greater0_tissue_D ,base_height= 10 ,base_aspect_ratio = 2) 

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
  
  + coord_flip()


###!!! UNHASH TO SAVE PLOT!!!!!
# 
# library(cowplot)
# save_plot("./Plots/Human/Tissue/p_TPM_0.1_tissue_D.png",
#           p_TPM_0.1_tissue_D ,base_height= 10 ,base_aspect_ratio = 2) 


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
  
  + coord_flip()


###!!! UNHASH TO SAVE PLOT!!!!!
# 
# library(cowplot)
# save_plot("./Plots/Human/Tissue/p_TPM_1_tissue_D.png",
#           p_TPM_1_tissue_D ,base_height= 10 ,base_aspect_ratio = 2) 

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
         title= "Gene Expression = 0 TPM ") %>%
  + theme(axis.title=element_text(size=14,face="bold"),
          plot.caption=element_text(face = "italic", size=12, hjust = 0),
          text = element_text(size=12),
          axis.text.y= element_text(size=12),
          axis.text.x= element_text(size=12))%>%
  +stat_summary(fun.data="mean_sdl", 
                geom="pointrange", width=0.2, colour="black" )%>%
  +stat_compare_means(method ="wilcox.test",label.y = 80,label.x = 1.5,paired = FALSE)    # Add global p-valu   # Add global p-valu


p_None_violin

# 
# library(cowplot)
# save_plot("./Plots/Human/Disease_Non-Disease/p_None_violin.jpeg",
#           p_None_violin ,base_height= 10 ,base_aspect_ratio = 2) 
#  


##CALCULATE THE MEAN AND SD

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
         title= "Gene Expression >0 TPM ") %>%
  
  + theme(axis.title=element_text(size=14,face="bold"),
          plot.caption=element_text(face = "italic", size=12, hjust = 0),
          text = element_text(size=12),
          axis.text.y= element_text(size=12),
          axis.text.x= element_text(size=12))%>%
  
  +stat_summary(fun.data="mean_sdl", 
                geom="pointrange", width=0.2, colour="black" )%>%
  +stat_compare_means(method ="wilcox.test",label.y = 80,label.x = 1.5,paired = FALSE)    # Add global p-valu   # Add global p-valu


p_greater0_violin

### !!!!!!!!!!!!UNHASH TO SAVE PLOT !!!!!!!!
# 
# library(cowplot)
# save_plot("./Plots/Human/Disease_Non-Disease/p_greater0_violin.jpeg",
#           p_greater0_violin ,base_height= 10 ,base_aspect_ratio = 2) 



## CALCULATE THE MEAN AND SD

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


dev.new(width=60, height=30) # to open in a new window
 
p_0.1_violin<- ggplot(HGNC.ID_DISEASE_NON_DI_0.1, 
                      aes(x=HGNC.ID_DISEASE_NON_DI_0.1$GROUP, 
                          y=HGNC.ID_DISEASE_NON_DI_0.1$VALUE, fill= GROUP))+ geom_violin(trim=FALSE)
p_0.1_violin<-p_0.1_violin%>%
  + scale_y_continuous(breaks = seq(-25, 80, by = 20))%>%
  + labs(x= "Disease/Non-Disease Associated Genes", y="Number of Tissues", 
         title= "Gene Expression >/= 0.1 TPM ") %>%
  + theme(axis.title=element_text(size=14,face="bold"),
          plot.caption=element_text(face = "italic", size=12, hjust = 0),
          text = element_text(size=12),
          axis.text.y= element_text(size=12),
          axis.text.x= element_text(size=12))%>%
  +stat_summary(fun.data="mean_sdl", 
                geom="pointrange", width=0.2, colour="black")%>%
  +stat_compare_means(method ="wilcox.test",label.y = 80,label.x = 1.5,paired = FALSE)    # Add global p-valu   # Add global p-valu


p_0.1_violin

# 
# library(cowplot)
# save_plot("./Plots/Human/Disease_Non-Disease/p_0.1_violin.jpeg",
#           p_0.1_violin ,base_height= 10 ,base_aspect_ratio = 2) 
#  

### CALCULATE THE MEAN AND THE SD FOR EACH GROUP

HGNC.ID_0.1_MEAN_SD <- aggregate(VALUE~ GROUP, HGNC.ID_DISEASE_NON_DI_0.1, function(x) c(mean = mean(x), sd = sd(x)))

HGNC.ID_0.1_MEAN_SD



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


HGNC.ID_DISEASE_1_d["GROUP"] <- "Disease Associated" # That creates the new column named "MY_NEW_COLUMN" 
HGNC.ID_NO_DISEASE_1_nd["GROUP"] <- "Not Disease Associated" # That creates the new column named "MY_NEW_COLUMN"


HGNC.ID_DISEASE_NON_DI_1<- rbind(HGNC.ID_DISEASE_1_d, 
                                 HGNC.ID_NO_DISEASE_1_nd)




# # Make Violin plot





dev.new(width=60, height=30) # to open in a new window

p_1_violin__<-ggplot(HGNC.ID_DISEASE_NON_DI_1, 
                     aes(x=HGNC.ID_DISEASE_NON_DI_1$GROUP, 
                         y=HGNC.ID_DISEASE_NON_DI_1$VALUE, fill=GROUP))+ geom_violin(trim=FALSE)
p_1_violin__<-p_1_violin__%>%
  + scale_y_continuous(breaks = seq(-25, 80, by = 20))%>%
  + labs(x= "Disease/Non-Disease Associated Genes", y="Number of Tissues", 
         title= "Gene Expression >/= 1 TPM ") %>%
  + theme(axis.title=element_text(size=14,face="bold"),
          plot.caption=element_text(face = "italic", size=12, hjust = 0),
          text = element_text(size=12),
          axis.text.y= element_text(size=12),
          axis.text.x= element_text(size=12))%>%
  +stat_summary(fun.data="mean_sdl", 
                geom="pointrange", width=0.2, colour="black" )%>%
  +stat_compare_means(method ="wilcox.test",label.y = 80,label.x = 1.5,paired = FALSE)    # Add global p-valu   # Add global p-valu



p_1_violin__


##!!!!!UNHASH TO SAVE PLOT !!!!!!!
# 
# library(cowplot)
# save_plot("./Plots/Human/Disease_Non-Disease/p_1_violin__.jpeg",
#           p_1_violin__ ,base_height= 10 ,base_aspect_ratio = 2) 
#  


### CALCULATE THE MEAN AND THE SD FOR EACH GROUP

HGNC.ID_1_MEAN_SD <- aggregate(VALUE~ GROUP, HGNC.ID_DISEASE_NON_DI_1, function(x) c(mean = mean(x), sd = sd(x)))

HGNC.ID_1_MEAN_SD


#####################################################################################################################################################
########################################################################################################################################################

#FOR EACH HPO.DESCRIPTION COUNT THE NUMBER OF GENES EXPRESSED FOR ALL GENES WITH HPO.DESC WHICH ARE ALL DISEASE ASSOCIATED



# TO REMOVE EVERYTHING ELSE ASIDE FROM THE GENE'S HGNC.ID AND THE HPO.DESCRIPTION

all_disease_genes_pheno_expression<-all_disease_genes_pheno


## ALL DISEASE GENES THAT HAVE AN EXPRESSION >/= 0.1

all_disease_genes_pheno_expression[4:56 ][ all_disease_genes_pheno_expression[4:56 ] >= 0.1 ] <- "Yes"

all_disease_genes_pheno_expression[4:56 ][ all_disease_genes_pheno_expression[4:56 ] < 0.1] <- "No"


                            
all_disease_genes_pheno_expression_expr <- all_disease_genes_pheno_expression %>% 
  gather("GTEX.tissues","Expression", `Adipose - Subcutaneous`:`Whole Blood`)%>%
  select(HGNC.ID,`HPO-Term-ID`,hpo.ancestors,hpo.description,GTEX.tissues,Expression)%>%
  distinct()

View(all_disease_genes_pheno_expression_expr)

all_disease_genes_hpo_desc<- all_disease_genes_pheno_expression_expr%>%
  select(hpo.description,HGNC.ID, Expression)%>%
  filter(Expression == "Yes")%>%
  distinct()

View(all_disease_genes_hpo_desc)
all_disease_genes_hpo_desc<-all_disease_genes_hpo_desc%>%
  select_all()%>%
  dplyr::group_by(hpo.description,HGNC.ID) %>%
  dplyr::summarise(n=n()) %>%
  distinct()

#count the number of disease genes with hpo.des
all_disease_genes_hpo_desc_COUNT<-ungroup(all_disease_genes_hpo_desc)%>%
  select(HGNC.ID )%>%
  distinct()

names(all_disease_genes_hpo_desc)[names(all_disease_genes_hpo_desc)=="n"] <- "Frequency"

all_disease_genes_hpo_desc_aggr<-aggregate(Frequency~ hpo.description, data = all_disease_genes_hpo_desc, sum)

View(all_disease_genes_hpo_desc_aggr)

## PLOT THE DISEASE ASSOCIATED GENES FOR ALL GENES WITH TPM >0.1


p_all_disease_genes_hpo_desc<-ggplot(data=all_disease_genes_hpo_desc_aggr
                         , aes(x=reorder(all_disease_genes_hpo_desc_aggr$hpo.description,- all_disease_genes_hpo_desc_aggr$Frequency),
                               y=all_disease_genes_hpo_desc_aggr$Frequency, 
                               fill=all_disease_genes_hpo_desc_aggr$Frequency))

p_all_disease_genes_hpo_desc<-p_all_disease_genes_hpo_desc %>%
  
  + labs(x= "Top-Level HPO Term", y="Number of Genes", 
         title= "Number of Genes with a Phenotype to the Top-Level HPO Term")%>%
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


#p_all_disease_genes_hpo_desc



##!!!!!UNHASH TO SAVE PLOT!!!!!!!!!!!
# 
# library(cowplot)
# save_plot("./Plots/Human/HPO.DESCRIPTION/p_all_disease_genes_hpo_desc.jpeg",
#           p_all_disease_genes_hpo_desc ,base_height= 10 ,base_aspect_ratio = 2) 
#  






#######################################################################################################################################################
#########################################################################################################################################################

## Load the 'GTEX.Tissues.to.HPO.txt' file so that the HPO.Descriptions 
#for each threshold can be plotted against the tissues



#REMOVED EMPTY COLUMNS BY SELECTING THOSE THAT WERE NEEDED
GTEX.Tissues.to.HPO <-GTEX.Tissues.to.HPO %>%
  select(GTEX.tissues, HPO.id, HPO.description,
        HPO.superclass.id, HPO.superclass.description, HPO.superclass.physiological.system)%>%
  distinct()


#!!!!! UNHASH TO SAVE FILE!!!!!!!

# write.csv(GTEX.Tissues.to.HPO,"./Output_Files/GTEX.Tissues.to.HPO.csv")

######################################################################################################################################


###########~~~~~~~~~~~~~~A FUNCTION FOR LINKING THE GTEX TISSUES TO THEIR HPO DESCRIPTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~################




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
    
    + coord_flip()
  
  # 
  # library(cowplot)
  # save_plot(paste("./Plots/Human/System/", System_Threshold,".jpeg"),
  #           GTEX.human_genes_cleaned_ABN_COUNT_p ,base_height= 10 ,base_aspect_ratio = 2) 
  # 
  
  }

# TEST


# GTEX.to.human_genes(GTEX.Tissues=GTEX.Tissues.to.HPO,
#                     human_genes=human_genes_TPM_1_pheno  ,
#                     Abnormality= "Abnormality of the nervous system", HPO.ID= "HP:0000707",
#                     System="nervous system", Threshold=">/= 1", System_Threshold= paste("nervous_system","_","1")) 
# # df_threshold(df=human_genes_TPM_1_pheno , threshold=">/= 1", threshold_2="1")


#######################################################################################################################################


All_systems<-GTEX.Tissues.to.HPO%>%
  select(HPO.superclass.description,HPO.superclass.id,HPO.superclass.physiological.system)%>%
  filter(HPO.superclass.description %in% all_disease_genes_pheno$hpo.description |
           HPO.superclass.id %in% all_disease_genes_pheno$hpo.ancestors )%>%
  distinct()

#View(All_systems)

##### A FUNCTION TO RUN AFTER THE ABOVE FUNCTION ################

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
# # 
# df_threshold(df=human_genes_TPM_None_pheno, threshold="= 0", threshold_2="NONE")
# # 
# # 
# df_threshold(df=human_genes_TPM_greater0_pheno, threshold="> 0", threshold_2="0")
# # 
# df_threshold(df=human_genes_TPM_0.1_pheno, threshold=">/= 0.1", threshold_2= "0.1")
# # 
# df_threshold(df=human_genes_TPM_1_pheno , threshold=">/= 1", threshold_2="1")
# 
# 
# 
# 
# 




###################################################################################################################################################################

###################~~~~~~~~CREATE A FUNCTION THAT COUNTS THE NUMBER OF TISSUES THAT A GENE IS BEING EXPRESSED IN ~~~~~~~~~~~~~~~~~~~~~#############
##################~~~~~~~~~DEPENDING ON WHETHER THE TISSUE IS ASSOCIATED TO THE GENE'S HPO PHENOTYPE OR NOT ~~~~~~~~~~~~~~~~~~~~~###############

#######################################################################################################################################

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
################################################################################################

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


#RETRIEVE THE TOTAL NUMBER OF GENES THAT HAVE A LINKED AND NOT LINKED PHENOTYPE TO GET A LIST OF GENES FOR THE FOR LOOP

#### This is to create a df of genes that can be used to iterate through the 'For loop'
Pheno_cleaned <- all_disease_genes_pheno %>% 
  gather("GTEX.tissues","Expression", `Adipose - Subcutaneous`:`Whole Blood`)%>%
  select(HGNC.ID,hpo.ancestors,hpo.description,GTEX.tissues,Expression)%>%
  distinct()
View(Pheno_cleaned)

# merge the tissues with their hpo.descp using the gtex file 


# USE THE GTEX TISSUES TO GET THE HPO.DESCP-TISSUE RELATIONSHIP CORRECTED

GTEX.Tissues.to.tissues<-GTEX.Tissues.to.HPO%>%
  select_all()%>%
  drop_na()%>%
  select(GTEX.tissues)

## There are some tissues that are found in different HPO.superclass.description tehrefore have been duplicated in the GTEX.Tissues.to.tissues

# Merge the GTEX.Tissues.to.tissues and human_genes_TPM_greater0_pheno_count to get the duplicated tissues

Pheno_Cleaned_Tissues<-merge(GTEX.Tissues.to.tissues, Pheno_cleaned, 
                             by = c("GTEX.tissues", "GTEX.tissues"), all.x = TRUE)

#View(pheno_cleaned_Tissues)


# The human_genes_TPM_greater0_pheno_count and GTEX.Tissues.to.HPO joind to get the HPO.superclass.description for each tissue



Pheno_Cleaned_Tissues_GTEX<-inner_join(GTEX.Tissues.to.HPO,Pheno_Cleaned_Tissues)
Pheno_Cleaned_Tissues_GTEX<-Pheno_Cleaned_Tissues_GTEX%>%
  select_all()%>%
  drop_na()%>%    # the NAs are dropped
  distinct(GTEX.tissues,hpo.description,HGNC.ID,HPO.superclass.description,Expression ) # duplicates are removed that are not needed

View(Pheno_Cleaned_Tissues_GTEX)

#########################################################################


### Genes that are Linked to HPO.Description

Linked.HPO.desc<-Pheno_Cleaned_Tissues_GTEX%>%
  select(hpo.description,HGNC.ID)%>%
  distinct()


View(Linked.HPO.desc)

# Remove HPO.Description to get the genes that ARE NOT HPO.DESC LINKED AND HAVE GENE EXPRESSION IN OTHER TISSUES

Pheno_Cleaned_Tissues_GTEX<-Pheno_Cleaned_Tissues_GTEX%>%
  select(GTEX.tissues,HGNC.ID,HPO.superclass.description,Expression)%>%
  distinct()

View(Pheno_Cleaned_Tissues_GTEX)


# Keep rows that do NOT have the same "HPO.superclass.description"="hpo.description", to find the phenos that are not linked
NonLinked_PHENO<-anti_join(Pheno_Cleaned_Tissues_GTEX,Linked.HPO.desc, by=c("HGNC.ID"="HGNC.ID","HPO.superclass.description"="hpo.description")) # keep rows with matching ID



# Keep rows that ARE the same (NOTE THERE MAY BE SOME GENES THAT DID NOT HAVE AN OBSERVED PHENOTYPE BUT HAD A LINKED HPO.DESC THEREFORE THEY MAY HAVE NAS)
Linked_PHENO<-right_join(Pheno_Cleaned_Tissues_GTEX,Linked.HPO.desc, by=c("HGNC.ID"="HGNC.ID","HPO.superclass.description"="hpo.description")) # keep rows with matching ID

# Remove NAs, genes that have linked  expression but no tissues or expression linked 

Linked_PHENO<-Linked_PHENO%>%
  select_all()%>%
  drop_na()%>%
  distinct()



View(Linked_PHENO)


# KEEP GENES THAT HAVE BOTH LINKED AND NON-LINKED 

NonLinked_PHENO<-semi_join(NonLinked_PHENO,Linked_PHENO, by=c("HGNC.ID"="HGNC.ID")) # keep rows with matching ID

View(NonLinked_PHENO)
# COUNT the number of tissues 

#change expression to factor

NonLinked_PHENO$Expression <- as.factor(NonLinked_PHENO$Expression)

Linked_PHENO$Expression <- as.factor(Linked_PHENO$Expression)


#Remove GTEX.TISSUEs To prevent duplication of  the gene

Linked_PHENO_HGNC<-ungroup(Linked_PHENO)%>%
  select(HGNC.ID, -GTEX.tissues)%>%
  dplyr::distinct()

all_genes2<-right_join(Linked_PHENO_HGNC, NonLinked_PHENO)
all_genes2<-all_genes2%>%
  select(HGNC.ID)%>%
  distinct()

View(all_genes2)    # GENES THAT HAVE HPO PHENOTYPE AND EXPRESSION IN HPO PHENOTYPE(LINKED AND UNLINKED)





####~~~~~~~~~~~~~~~~~ To do a Fischer test for every Tissue vs Tissues not in the same System~~~~~~~

#################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~THE GENE FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#################################



GENE_FISHER_TEST<-function(df, Gene){  


  
  #THRESHOLD DATAFRAME WITH THE GENE EXPRESSIONS
  
  #View(human_genes_TPM_greater0_pheno)
  
  
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
  

  # Keep rows that do NOT have the same "HPO.superclass.description"="hpo.description", to find the phenos
  #that are not linked
  
  nonlinked_PHENO<-anti_join(pheno_cleaned_Tissues_,linked.hpo.desc, by=c("HGNC.ID"="HGNC.ID","HPO.superclass.description"="hpo.description")) # keep rows with matching ID
  
  
  # Keep rows that ARE the same (NOTE THERE MAY BE SOME GENES THAT DID NOT HAVE AN OBSERVED PHENOTYPE 
  #BUT HAD A LINKED HPO.DESC THEREFORE THEY MAY HAVE NAS)
  linked_PHENO<-right_join(pheno_cleaned_Tissues_,linked.hpo.desc, by=c("HGNC.ID"="HGNC.ID","HPO.superclass.description"="hpo.description")) # keep rows with matching ID
  
  # Remove NAs, genes that have linked  expression but no tissues or expression linked 
  
  linked_PHENO<-linked_PHENO%>%
    select_all()%>%
    drop_na()%>%
    distinct()
  
  #View(linked_PHENO)
  
  
  # KEEP GENES THAT HAVE BOTH LINKED AND NON-LINKED 
  
  nonlinked_PHENO<-semi_join(nonlinked_PHENO,linked_PHENO, by=c("HGNC.ID"="HGNC.ID")) # keep rows with matching ID
  
  # CHANGE EXPRESSION TO FACTOR SO THAT THE .drop=FALSE can be used for genes that only have Yes or No
  nonlinked_PHENO$Expression <- as.factor(nonlinked_PHENO$Expression)
  
  linked_PHENO$Expression <- as.factor(linked_PHENO$Expression)

  
  linked_PHENO_COUNT2<-linked_PHENO%>%
    dplyr::select(GTEX.tissues,HGNC.ID, Expression)%>%
    dplyr::distinct()
  
  linked_PHENO_COUNT3<-linked_PHENO_COUNT2%>%
    dplyr::select_all()%>%
    dplyr::group_by(HGNC.ID, GTEX.tissues,Expression,.drop = FALSE) %>%  # group by gene here 
    dplyr::summarise(Count=n()) %>%     #mutate keeps the other columns
    drop_na()%>%
    dplyr::distinct(HGNC.ID,GTEX.tissues,Expression,Count ) 
  
  #View(linked_PHENO_COUNT3)
  
  nonlinked_PHENO_COUNT2<-nonlinked_PHENO%>%
    dplyr::select(GTEX.tissues,HGNC.ID, Expression)%>%
    dplyr::distinct()
  
#  View(nonlinked_PHENO_COUNT2)
  
  nonlinked_PHENO_COUNT3<-nonlinked_PHENO_COUNT2%>%
    dplyr::select_all()%>%
    dplyr::group_by(HGNC.ID, GTEX.tissues,Expression,.drop = FALSE) %>%  # group by gene here 
    dplyr::summarise(Count=n()) %>%     #mutate keeps the other columns
    drop_na()%>%
    dplyr::distinct(HGNC.ID,GTEX.tissues,Expression,Count ) 
  

  
  All__genes<-ungroup(nonlinked_PHENO_COUNT3)%>%   #gtex.tissue was being grouped with the hgnc.id so i removed the grouping
    select(HGNC.ID, -GTEX.tissues)%>%
    distinct()

  
  All__genes2<-ungroup(linked_PHENO_COUNT3)%>%
    select(HGNC.ID, -GTEX.tissues)%>%
    distinct()
  
  all__genes<-right_join(All__genes, All__genes2)
  
  #CHANGE THE DF SO THAT WE HAVE YES AND NO AS COLUMNS
  
  nonlinked_PHENO_COUNT_spread<-nonlinked_PHENO_COUNT3 %>%
    spread(Expression,Count)%>%
    dplyr::distinct()

  linked_PHENO_Count_spread<-linked_PHENO_COUNT3 %>%
    spread(Expression,Count)%>%
    dplyr::distinct()
  

  #Filter gene
  
  linked_PHENO_Count_GENESPECIFIC<-linked_PHENO_Count_spread%>%
    dplyr::select_all()%>%
    dplyr::filter(HGNC.ID == Gene )%>%
    dplyr::distinct() # duplicates are removed that are not needed
  
    

  nonlinked_PHENO_COUNT_GENESPECIFIC<-nonlinked_PHENO_COUNT_spread%>%
    dplyr::select_all()%>%
    dplyr::filter(HGNC.ID == Gene )%>%
    dplyr::distinct() # duplicates are removed that are not needed
  
  
  # 
  nonlinked_PHENO_COUNT_GENESPECIFIC_SUM<-nonlinked_PHENO_COUNT_GENESPECIFIC%>%
    group_by(HGNC.ID)%>%
    summarise_at(c("No","Yes"),sum) %>% #sums the columns
    mutate(Phenotype="Not Linked Phenotype Observed")%>%
    distinct()
  

  
  
  linked_PHENO_Count_GENESPECIFIC_SUM<-linked_PHENO_Count_GENESPECIFIC%>%
    group_by(HGNC.ID)%>%
    summarise_at(c("No","Yes"),sum) %>% #sums the columns
    mutate(Phenotype="Linked Phenotype Observed")%>%
    distinct()
  
  

  
  # Bind the two datfarames to create a matrix
  Linked_NON_LINKED_PHENO<-rbind(linked_PHENO_Count_GENESPECIFIC_SUM,nonlinked_PHENO_COUNT_GENESPECIFIC_SUM)
  
  #View(Linked_NON_LINKED_PHENO)
  
  #define the rownames 
  library(magrittr)
  Linked_NON_LINKED_PHENO<-Linked_NON_LINKED_PHENO %>%
    set_rownames(.$Phenotype) %>% 
    select(-HGNC.ID, -Phenotype)
  
  
  View(Linked_NON_LINKED_PHENO)
  
  
  # perform Fisher test
  
  GENEPVals<-fisher.test(Linked_NON_LINKED_PHENO,alternative = "two.sided" )  # THIS WORKED
  
  GENEPVALUES<-GENEPVals$p.value
  
  return(GENEPVALUES)
  
  
  }


Gene_FISHER_TEST(df=human_genes_TPM_1_pheno, Gene="HGNC:1382")  
  

###########~~~~~~~~~~~~~~~~USE OF FOR-LOOP ~~~~~#

#THE FOR LOOP ITERATES THROUGH THE all_disease_genes_pheno_Tissues DATAFRAME WHICH HAS ALL TISSUES AND HPO.DESC(TOPLEVEL)

# CREATE A data frame with two columns and 55 rows, WHICH HAS THE TISSUE NAMES, HPO.SUPERCLASS AND AN EMPTY COLUMN
View(all_genes2)



#View(human_genes_TPM_greater0_pheno)

## !!!!!!!!!!!!!!!!!!!!!!!!! LOOPS HAVE BEEN '#' HASHED OUT PLS REMOVE '#' !!!!!


## TPM >0 TPM (Yes is for genes that have a gene expression of >0):
 

# NO GENE EXPRESSION 
# PVALUE_GENE_Greater0 <- data.frame(all_genes2,
#                                   P.VALUE=vector(length=3678)) 

# for(i in 1:nrow(all_genes2)) {
#   ii <- all_genes2[i,1]     #FIRST COLUMN
#   print(ii)                                     # TO CHECK IF THE CORRECT ROW IS BEING TAKEN
# 
#   
#   ## THE OUTCOME OF THE FUNCTION IS SAVED AS PVAL WHICH IS THE FISHER TEST PVALUE 
#   
#   PVAL<-Gene_FISHER_TEST(df=human_genes_TPM_greater0_pheno, Gene= paste(ii)) 
# 
#   print(PVAL)                  # THE PVALUE IS PRINTED
#   
#   PVALUE_GENE_None$P.VALUE[i]<-paste(PVAL)              # THE PVALUE IS ADDED INTO THE DF CREATED ABOVE
#   
#   }
# 
# write.csv(PVALUE_GENE_Greater0,'./Output_Files/Human/PVALUE_GENE_Greater0.csv')
#   




## No gene expression (Yes is for genes that have a gene expression of 0)

# 
# PVALUE_GENE_None <- data.frame(all_genes2,
#                                P.VALUE=vector(length=3678)) 


#View(human_genes_TPM_greater0_pheno)
# 
# for(i in 1:nrow(all_genes2)) {
#   ii <- all_genes2[i,1]     #FIRST COLUMN
#   print(ii)                                     # TO CHECK IF THE CORRECT ROW IS BEING TAKEN
#   
#   
#   ## THE OUTCOME OF THE FUNCTION IS SAVED AS PVAL WHICH IS THE FISHER TEST PVALUE 
#   
#   PVALN<-Gene_FISHER_TEST(df=human_genes_TPM_None_pheno, Gene= paste(ii)) 
#   
#   print(PVALN)                  # THE PVALUE IS PRINTED
#   
#   PVALUE_GENE_None$P.VALUE[i]<-paste(PVALN)              # THE PVALUE IS ADDED INTO THE DF CREATED ABOVE
#   
# }
# 
# write.csv(PVALUE_GENE_None,'./Output_Files/Human/PVALUE_GENE_None.csv')



## For gene expression >0.1

# 
# PVALUE_GENE_0.1 <- data.frame(all_genes2,
#                               P.VALUE=vector(length=3678)) 
# 
# for(i in 1:nrow(all_genes2)) {
#   ii <- all_genes2[i,1]     #FIRST COLUMN
#   print(ii)                                     # TO CHECK IF THE CORRECT ROW IS BEING TAKEN
#   
#   
#   ## THE OUTCOME OF THE FUNCTION IS SAVED AS PVAL WHICH IS THE FISHER TEST PVALUE 
#   
#   PVAL<-Gene_FISHER_TEST(df=human_genes_TPM_0.1_pheno, Gene= paste(ii)) 
#   
#   print(PVAL)                  # THE PVALUE IS PRINTED
#   
#   PVALUE_GENE_0.1$P.VALUE[i]<-paste(PVAL)              # THE PVALUE IS ADDED INTO THE DF CREATED ABOVE
#   
# }
# 
# write.csv(PVALUE_GENE_0.1,'./Output_Files/Human/PVALUE_GENE_0.1.csv')


## For gene expression >1

# 
# PVALUE_GENE_1 <- data.frame(all_genes2,
#                             P.VALUE=vector(length=3678)) 
# 
# 
# #View(human_genes_TPM_greater0_pheno)
# for(i in 1:nrow(all_genes2)) {
#   ii <- all_genes2[i,1]     #FIRST COLUMN
#   print(ii)                                     # TO CHECK IF THE CORRECT ROW IS BEING TAKEN
#   
#   
#   ## THE OUTCOME OF THE FUNCTION IS SAVED AS PVAL WHICH IS THE FISHER TEST PVALUE 
#   
#   PVAL<-Gene_FISHER_TEST(df=human_genes_TPM_1_pheno, Gene= paste(ii)) 
#   
#   print(PVAL)                  # THE PVALUE IS PRINTED
#   
#   PVALUE_GENE_1$P.VALUE[i]<-paste(PVAL)              # THE PVALUE IS ADDED INTO THE DF CREATED ABOVE
#   
# }
# 
# write.csv(PVALUE_GENE_1,'./Output_Files/Human/PVALUE_GENE_1.csv')











#######~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ADJUST PVALUE FOR ALL THRESHOLDS~~~~~~~~~~~~~~~~~~~~~~~~~~~#################


## Gene Expression >0
PVALUE_GENE_Greater0<-read_csv("D:/MSC RESEARCH PROJECT/Output_Files/PVALUE_GENE_Greater0.csv")

PVALUE_GENE_Greater0<-PVALUE_GENE_Greater0%>%
  select(HGNC.ID, P.VALUE)%>%
  distinct()
#View(PVALUE_GENE_Greater0)

## The Benjamini & Hochberg (1995) "BH" was used 

PVALUE_GENE_Greater0$Adjusted.Pvalue<-p.adjust(PVALUE_GENE_Greater0$P.VALUE, method="BH")

View(PVALUE_GENE_Greater0)


## Gene Expression =0


PVALUE_GENE_None<-read_csv("D:/MSC RESEARCH PROJECT/Output_Files/PVALUE_GENE_None.csv")

PVALUE_GENE_None<-PVALUE_GENE_None%>%
  select(HGNC.ID, P.VALUE)%>%
  distinct()
#View(PVALUE_GENE_None)

## The Benjamini & Hochberg (1995) "BH" was used 

PVALUE_GENE_None$Adjusted.Pvalue<-p.adjust(PVALUE_GENE_None$P.VALUE, method="BH")

View(PVALUE_GENE_None)



## Gene Expression >0.1


PVALUE_GENE_0.1<-read_csv("D:/MSC RESEARCH PROJECT/Output_Files/PVALUE_GENE_0.1.csv")

PVALUE_GENE_0.1<-PVALUE_GENE_0.1%>%
  select(HGNC.ID, P.VALUE)%>%
  distinct()
#View(PVALUE_GENE_0.1)

## The Benjamini & Hochberg (1995) "BH" was used 

PVALUE_GENE_0.1$Adjusted.Pvalue<-p.adjust(PVALUE_GENE_0.1$P.VALUE, method="BH")

View(PVALUE_GENE_0.1)




## Gene Expression >1


PVALUE_GENE_1<-read_csv("D:/MSC RESEARCH PROJECT/Output_Files/PVALUE_GENE_1.csv")

PVALUE_GENE_1<-PVALUE_GENE_1%>%
  select(HGNC.ID, P.VALUE)%>%
  distinct()
View(PVALUE_GENE_1)

## The Benjamini & Hochberg (1995) "BH" was used 

PVALUE_GENE_1$Adjusted.Pvalue<-p.adjust(PVALUE_GENE_1$P.VALUE, method="BH")

View(PVALUE_GENE_1)




##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PREPARE DATAFRAMES FOR GEA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~############


## Get teh Entrez Gene IDs for each gene (From the first dataframe human_gene)

human_genes_with_GENE_ID_<-human_genes_with_GENE_ID%>%
  select(HGNC.ID, gene_id)%>%
  distinct()

View(human_genes_with_GENE_ID_)


## >0 

PVALUE_GENE_Greater0_entrez<-right_join(human_genes_with_GENE_ID_, PVALUE_GENE_Greater0, by=c("HGNC.ID" ="HGNC.ID" ) )


PVALUE_GENE_Greater0_entrez<-PVALUE_GENE_Greater0_entrez%>%
  select(gene_id, Adjusted.Pvalue)%>%
  distinct()

View(PVALUE_GENE_Greater0_entrez)


# Filter <0.05

Greater0.test.set<- PVALUE_GENE_Greater0_entrez %>%
  select(gene_id, Adjusted.Pvalue)%>%
  filter(Adjusted.Pvalue <= 0.05)%>%
  distinct(gene_id, Adjusted.Pvalue)

View(Greater0.test.set)

## No gene expression

PVALUE_GENE_None_entrez<-right_join(human_genes_with_GENE_ID_, PVALUE_GENE_None, by=c("HGNC.ID" ="HGNC.ID" ) )


PVALUE_GENE_None_entrez<-PVALUE_GENE_None_entrez%>%
  select(gene_id, Adjusted.Pvalue)%>%
  distinct()
View(PVALUE_GENE_None_entrez)


# Filter <0.05

None.test.set<- PVALUE_GENE_None_entrez %>%
  select(gene_id, Adjusted.Pvalue)%>%
  filter(Adjusted.Pvalue <= 0.05)%>%
  distinct(gene_id, Adjusted.Pvalue)

View(None.test.set)

## >=/= 0.1 expression

PVALUE_GENE_0.1_entrez<-right_join(human_genes_with_GENE_ID_, PVALUE_GENE_0.1, by=c("HGNC.ID" ="HGNC.ID" ) )


PVALUE_GENE_0.1_entrez<-PVALUE_GENE_0.1_entrez%>%
  select(gene_id, Adjusted.Pvalue)%>%
  distinct()
View(PVALUE_GENE_None_entrez)


# Filter <0.05

test.set_0.1<- PVALUE_GENE_0.1_entrez %>%
  select(gene_id, Adjusted.Pvalue)%>%
  filter(Adjusted.Pvalue <= 0.05)%>%
  distinct(gene_id, Adjusted.Pvalue)

View(test.set_0.1)



## >=/= 1 expression

PVALUE_GENE_1_entrez<-right_join(human_genes_with_GENE_ID_, PVALUE_GENE_1, by=c("HGNC.ID" ="HGNC.ID" ) )


PVALUE_GENE_1_entrez<-PVALUE_GENE_1_entrez%>%
  select(gene_id, Adjusted.Pvalue)%>%
  distinct()
# 
# # remove "." from the entrez gene ids
# 
# t<-t(data.frame(strsplit(PVALUE_GENE_1_entrez$gene_id, ".", fixed = TRUE)))
# View(t)
# t<-data.frame(t)
# PVALUE_GENE_1_entrez$GENE_SPLIT<-t$X1
# view(PVALUE_GENE_1_entrez)
# 

# Keep only the gene ids that don't have a "." with the pvalue


PVALUE_GENE_1_entrez<-PVALUE_GENE_1_entrez%>%
  select(GENE_SPLIT, Adjusted.Pvalue)%>%
  distinct()


# Make the dataframe into a vector

PVALUE_GENE_1_entrez<-as.vector(PVALUE_GENE_1_entrez)

View(PVALUE_GENE_1_entrez)

# Filter <0.05


test.set_1<- PVALUE_GENE_1_entrez %>%
  select(GENE_SPLIT, Adjusted.Pvalue)%>%
  filter(Adjusted.Pvalue <= 0.05)%>%
  distinct(GENE_SPLIT, Adjusted.Pvalue)



# Make the dataframe into a vector


test.set_1<-as.vector(test.set_1)

View(test.set_1)



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GENE ENRICHMENT ANALYSIS ~~~~~~~~~~~~~~~~~~~~~~~~~#####

## GO annotations and enrichment ################################################################################################

############# LOAD PACKAGES ###########################################################
# !!!!!!!!!!!!!!!!!!!!! REMOVE THE "#" TO RUN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

# source("https://bioconductor.org/biocLite.R")
# biocLite("org.Hs.eg.db")
# # 
# source("https://bioconductor.org/biocLite.R") 
# biocLite("GOstats") 

# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# # The following initializes usage of Bioc devel
# BiocManager::install(version='devel')
# 
# BiocManager::install("BiocGenerics")
# 
# 
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("org.Hs.eg.db")
# 
# 
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GOstats")
# 
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("multtest")
# 
# 
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GO.db")
# 
# library(org.Hs.eg.db);library(GOstats);library(multtest); library(GO.db)
# # 
# 

##############################################################################################


###### FUNCTION FOR THE GEA ###################



#!!!!!!!!!! UNHASH TO RUN THE FUNCTION !!!!!!!!!!!!!!!!

#go.enrichment <- function(test.set,reference.set,ontology.type ="BP",label=NULL) {
   
   
 
#   GO.param <- new("GOHyperGParams",
#                   
#                   geneIds = test.set,
#                   
#                   universeGeneIds = reference.set,
#                   
#                   ontology = ontology.type,
#                   
#                   annotation = "org.Hs.eg.db",
#                   
#                   testDirection = "over",
#                   
#                   pvalueCutoff = 1,
#                   
#                   conditional = T)
#   
#   
#   
#   GO.sign.results <- as.data.frame(summary(hyperGTest(GO.param))) %>%
#     
#     mutate (Pvalue.BH = p.adjust(Pvalue,method="BH")) %>%
#     
#     filter(Pvalue.BH < 0.05) %>%
#     
#     dplyr::select(1,3:7,2,8)
#   
#   
#   
#   
#   
#   if(missing(label)) {
#     
#     
#     
#     return(GO.sign.results)
#     
#     
#     
#   } else {
#     
#     
#     
#     return(GO.sign.results %>%
#              
#              mutate(Category = label))
#     
#     
#     
#   }
#   
#   
#   
# }



#######~~~~~~~~~~~~~~~~~RUN GEA FUNCTION FOR >0.1 & >1 TPM FOR BP, CC, MF ~~~~~~~~~~~~~~~~##############

##########!!!!!!!!!!!!!!!!!!!!!!UNHASH TO RUN CODE !!!!!!!!!!!!!!!!!!!!!############################



#################################### TPM >1 #########################################################################################



################################ Biological Pathway ################################################################################

#go.enrichment_BP_1<-go.enrichment(test.set =test.set_1 ,reference.set=PVALUE_GENE_1_entrez,ontology.type ="BP",label=NULL)

#write.csv(go.enrichment_BP_1,'./Output_Files/Human/go.enrichment_BP_1.csv')

##########################################################################################



############################### Cellular Component #################################################

#go.enrichment_CC_1<-go.enrichment(test.set =test.set_1 ,reference.set=PVALUE_GENE_1_entrez,ontology.type ="CC",label=NULL)

#write.csv(go.enrichment_CC_1,'./Output_Files/Human/go.enrichment_CC_1.csv')

#############################################################################################



############################## Molecular Function ##################################################

#go.enrichment_MF_1<-go.enrichment(test.set =test.set_1 ,reference.set=PVALUE_GENE_1_entrez,ontology.type ="MF",label=NULL)

#write.csv(go.enrichment_MF_1,'./Output_Files/Human/go.enrichment_MF_1.csv')

################################################################################################################




# VIEW DFS FOR >1TPM RUNS ##########################################################

# View(go.enrichment_BP_1)
# View(go.enrichment_CC_1)
# View(go.enrichment_MF_1)

############################################################################################## 


###################################### TPM >0.1 #################################################


##################################### Biological Pathway ############################################

#go.enrichment_BP_0.1<-go.enrichment(test.set =test.set_0.1 ,reference.set=PVALUE_GENE_0.1_entrez,ontology.type ="BP",label=NULL)

#write.csv(go.enrichment_BP_0.1,'./Output_Files/Human/go.enrichment_BP_0.1.csv')

#########################################################################################################



#################################### Cellular Component #########################################################

#go.enrichment_CC_0.1<-go.enrichment(test.set =test.set_0.1 ,reference.set=PVALUE_GENE_0.1_entrez,ontology.type ="CC",label=NULL)

#write.csv(go.enrichment_CC_0.1,'./Output_Files/Human/go.enrichment_CC_0.1.csv')

#########################################################################################################################



###################################### Molecular Function #######################################################################

#go.enrichment_MF_0.1<-go.enrichment(test.set =test.set_0.1 ,reference.set=PVALUE_GENE_0.1_entrez,ontology.type ="MF",label=NULL)

#write.csv(go.enrichment_MF_0.1,'./Output_Files/Human/go.enrichment_MF_0.1.csv')

######################################################################################################################################

## VIEW DFS FOD >0.1 RUNS 


# View(go.enrichment_BP_0.1)
# View(go.enrichment_CC_0.1)
# View(go.enrichment_MF_0.1)
