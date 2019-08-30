###############################################################################################################################
###############################################################################################################################

### Project: PhenoExpression ##################################################################################################
### Script: GXD_Tidy_Data.R ###################################################################################################
### Purpose: Tidy gene expression data from GXD ###############################################################################
### Author: Pilar Cacheiro, Friha Zafar ####################################################################################################
### Date: 04/03/2019, 08/03/2019, 22/03/2019 ##########################################################################################################

###############################################################################################################################
###############################################################################################################################

## load packages ##############################################################################################################

library(dplyr);library(tidyr);library(readr);

################################################################################################################################
################################################################################################################################

## import files ################################################################################################################

## import and tidy gene expression for TS28 mutant lines


gxd.ts28 <- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/MGIgeneExpressionQuery_20190304_mutant.txt.gz",delim = "\t") %>%
  select(MGI.Gene.ID ,Structure,Detected) %>% # select columns of interest 
  # we will try to always use MGI IDs for mouse genes / HGNC ids for human genes (these are stable identifiers)
  distinct(MGI.Gene.ID, Structure,Detected) %>% # retain unique rows
  arrange(Detected) %>% ## sort by this variable
  group_by(MGI.Gene.ID,Structure) %>% # group data by gene and tissue %>%
  summarise(Detected_all = paste0(Detected,collapse="|")) %>% ## collpase all the expression data availabe for a given
  # gene and tissue  in one cell
  mutate(Detected_all = ifelse(Detected_all =="No","No",
                               ifelse(Detected_all == "Yes","Yes","Ambiguous_NotSpecified_Contradictory"))) %>%
  # recode variable based on values 
  rename(MGI.ID = MGI.Gene.ID, TS28.mutant.expression.detected = Detected_all) %>% ## rename variables
  select(MGI.ID,Structure,TS28.mutant.expression.detected) %>% # select columns of interest
  distinct(MGI.ID,Structure,TS28.mutant.expression.detected) %>% # check again for duplicated rows
  mutate(Gene.Anatomy = paste0(MGI.ID,"-",Structure)) # create a new variable (Gene - Anatomy (Structure) term)


## count the number of KO genes 
gxd.ts28_count<-gxd.ts28%>%
  select(MGI.ID)%>%
  drop_na()%>%
  distinct()
nrow(gxd.ts28_count)

#################################################################################################################
#############################################################################################################################

##import MGI gene - phenotype data

mgi.genepheno <-  read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/MGI_GenePheno.rpt",delim="\t",
                             col_names = c("Allelic.Composition","Allele.Symbols","Allele.ID",
                                           "Genetic.Background","MP.ID","PubMed.ID","MGI.ID","MGI.Genotype")) %>%
  # import file, add column names
  select(MGI.ID,MP.ID) %>% # select columns of interest
  distinct(MGI.ID, MP.ID) # unique gene-phenotype associations

## import phenotype ontology terms to anatomy terms mapping file

mapping <- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/MP_EMAPA.rpt",delim="\t",col_names = c("MP.ID","MP.Term","EMAPA.ID","EMAPA.Term"))


################################################################################################################################
################################################################################################################################

## map mgi mouse phentoypes to tissues

mgi.genepheno.tissue <- mgi.genepheno %>%
  inner_join(mapping,by = c("MP.ID" = "MP.ID")) %>% # join two dataframe by MP.ID (inner_join: only MP IDS with tissue mapping)
  ## a given mp term can map to different tissues
  mutate(Gene.Anatomy = paste0(MGI.ID,"-",EMAPA.Term)) # create a new variable (Gene - Anatomy term)


## map mgi mouse phentoypes and the corresponding tissues to gene expression results


mgi.genepheno.tissue.ts28 <- mgi.genepheno.tissue %>%
  left_join(gxd.ts28, by=c("Gene.Anatomy"  ="Gene.Anatomy")) %>% # join two dataframe by Gene.Anatomy
  # left_join: keep all the data in the first dataframe
  rename(MGI.ID = MGI.ID.x) %>% # rename variable
  select(MGI.ID,MP.ID,MP.Term,EMAPA.ID,EMAPA.Term,Gene.Anatomy,TS28.mutant.expression.detected) %>%
  # select columns of interest
  distinct(MGI.ID,MP.ID,MP.Term,EMAPA.ID,EMAPA.Term,Gene.Anatomy,TS28.mutant.expression.detected) # unique rows



################################################################################################################################
################################################################################################################################








#####~~~~~~~~~~~~~~~~Wildtype Gene Expression Data Analysis~~~~~~~~~~~~ ###################################################################################


################################################################################################################################
################################################################################################################################

## Importing file ################################################################################################################################
## Cleaning the data ###########################################################################################################

gxd_wt.ts28<- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/MGIgeneExpressionQuery_20190304_wildtype.txt",delim = "\t")  %>%
  select(MGI.Gene.ID,Structure,Detected) %>% # select columns of interest
  # we will try to always use MGI IDs for mouse genes / HGNC ids for human genes (these are stable identifiers)
  distinct(MGI.Gene.ID, Structure,Detected) %>% # retain unique rows
  arrange(Detected) %>% ## sort by this variable
  group_by(MGI.Gene.ID,Structure) %>% # group data by gene and tissue %>%
  summarise(Detected_all = paste0(Detected,collapse="|")) %>% ## collpase all the expression data availabe for a given
  # gene and tissue  in one cell
  mutate(Detected_all = ifelse(Detected_all =="No","No",
                               ifelse(Detected_all == "Yes","Yes","Ambiguous_NotSpecified_Contradictory"))) %>%
  # recode variable based on values
  rename(MGI.ID = MGI.Gene.ID, TS28.wt.expression.detected = Detected_all) %>% ## rename variables
  select(MGI.ID,Structure,TS28.wt.expression.detected) %>% # select columns of interest
  distinct(MGI.ID,Structure,TS28.wt.expression.detected) %>% # check again for duplicated rows
  mutate(Gene.Anatomy = paste0(MGI.ID,"-",Structure)) # create a new variable (Gene - Anatomy (Structure) term)



gxd_wt.ts28_count<-gxd_wt.ts28%>%
  select(MGI.ID)%>%
  drop_na()%>%
  distinct()
nrow(gxd_wt.ts28_count)
################################################################################################################################
################################################################################################################################

## map mgi mouse phentoypes and the corresponding tissues to gene expression results of WT line


mgi.genepheno.tissue.ts28_wt <- mgi.genepheno.tissue %>%
  left_join(gxd_wt.ts28, by=c("Gene.Anatomy"  ="Gene.Anatomy")) %>% # join two dataframe by Gene.Anatomy
  # left_join: keep all the data in the first dataframe
  rename(MGI.ID = MGI.ID.x) %>% # rename variable
  select(MGI.ID,MP.ID,MP.Term,EMAPA.ID,EMAPA.Term,Gene.Anatomy,TS28.wt.expression.detected) %>%
  # select columns of interest
  distinct(MGI.ID,MP.ID,MP.Term,EMAPA.ID,EMAPA.Term,Gene.Anatomy,TS28.wt.expression.detected) # unique rows

################################################################################################################################
################################################################################################################################

### To remove 'NA', and 'Ambiguous_NotSpecified_Contradictory' data
mgi.genepheno.tissue.ts28_wt_mutant<- inner_join(mgi.genepheno.tissue.ts28, mgi.genepheno.tissue.ts28_wt) %>%
  select(MGI.ID, MP.ID, MP.Term, EMAPA.ID, EMAPA.Term, Gene.Anatomy, TS28.mutant.expression.detected, TS28.wt.expression.detected) %>%
  drop_na() %>%  #dropped all rows with 'NA'
  filter(TS28.mutant.expression.detected  != "Ambiguous_NotSpecified_Contradictory") %>%  #Remove ambiguous expression data
  filter(TS28.wt.expression.detected != "Ambiguous_NotSpecified_Contradictory") 
  
  
View(mgi.genepheno.tissue.ts28_wt_mutant)

################################################################################################################################
################################################################################################################################

### To link the phenotypes back to top level ###################################################################################

## Open files 

mp.ancestor.nodes<- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/mp.ancestor.nodes.txt",delim = "\t")
#View(mp.ancestor.nodes)

mp.toplevels<- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/mp.toplevels.txt",delim = "\t") %>%
  rename(mp.ancestors=mp.term )  #renamed mp.term to mp.ancestors 
################################################################################################################################
################################################################################################################################

## Join dfs to one 

mp.toplevels<-left_join(mp.toplevels, mp.ancestor.nodes)   # the left_join keeps the data that is already in mp.toplevels, and add the mp.ancestors column

#View(mp.toplevels)

################################################################################################################################
################################################################################################################################

## Link this back to the mgi.genepheno.tissue.ts28_wt_mutant

mgi.genepheno.tissue.ts28_wt_mutant <- left_join(mgi.genepheno.tissue.ts28_wt_mutant, mp.toplevels, by=c("MP.ID" = "mp.term") )
#View(mgi.genepheno.tissue.ts28_wt_mutant)


################################################################################################################################
################################################################################################################################

## To create a Constigency table for the top ancestor node and the expression detection

library(plyr)
library(dplyr)

## Join the two expression columns together 

mgi.genepheno.tissue.ts28_wt_mutant<-mgi.genepheno.tissue.ts28_wt_mutant%>%   # joined mutant and WT expressions
  mutate(Mutant.WT.Expression = paste0(TS28.mutant.expression.detected,"-",TS28.wt.expression.detected)) %>%
  dplyr::distinct(MGI.ID, MP.ID, MP.Term, mp.description,EMAPA.ID, EMAPA.Term, Gene.Anatomy,mp.ancestors, TS28.mutant.expression.detected, TS28.wt.expression.detected,Mutant.WT.Expression)  # To remove duplicates
  


View(mgi.genepheno.tissue.ts28_wt_mutant)


#!!!!!!!!UNHASH TO SAVE FILE !!!!!!!!!!!!!!!!


# write.csv(mgi.genepheno.tissue.ts28_wt_mutant,
#           './Output_Files/Mice/mgi.genepheno.tissue.ts28_wt_mutant.csv')

#####################################################################################################################################
#####################################################################################################################################

### COUNT WT-KO GENES
Count_WT_KO<-mgi.genepheno.tissue.ts28_wt_mutant%>%
  select(MGI.ID)%>%
  drop_na()%>%
  distinct()

nrow(Count_WT_KO)

#######################################################################################################################
###############################################################################################################################


# Constingency table with gene names and count the number of No-No, No-Yes, Yes-Yes, Yes-No

Express_freq2<-ftable(mgi.genepheno.tissue.ts28_wt_mutant$MP.ID,
                      mgi.genepheno.tissue.ts28_wt_mutant$mp.description, 
                     mgi.genepheno.tissue.ts28_wt_mutant$Mutant.WT.Expression)

#View(Express_freq2)

######################################################################################################################
#####################################################################################################################################

#  Count the number of No-No, No-Yes, Yes-Yes, Yes-No


Express_freq_Perc<-mgi.genepheno.tissue.ts28_wt_mutant %>%
  dplyr::group_by(mp.description,Mutant.WT.Expression) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(Percentage = n/sum(n)*100)


names(Express_freq_Perc)[names(Express_freq_Perc)=="n"] <- "Frequency"

View(Express_freq_Perc)

#####################################################################################################################################
#####################################################################################################################################

# A graphical representation of the data

library(ggplot2)
library(ggthemes)


p_Express_freq_Perc<-ggplot(data=Express_freq_Perc
          , aes(x= Express_freq_Perc$mp.description
                ,y=Express_freq_Perc$Percentage, fill=Express_freq_Perc$Mutant.WT.Expression))

p_Express_freq_Perc <-p_Express_freq_Perc%>%
  + labs(x= "Top-Level Phenotype", y="Percentage (%) of the Expression Profile for Each Phenotype", 
         title= "The Gene Expression Profiles for Mutant and Wild-Type Mice.")%>%

  + guides(fill=guide_legend(title="Keys (Mutant-WT):  ")) %>%
  + theme(axis.title=element_text(size=14,face="bold"),
          plot.caption=element_text(face = "italic", size=12, hjust = 0),
          text = element_text(size=12),
          legend.position = "bottom",legend.direction = "horizontal",
          legend.background = element_rect(fill="lightyellow", size=2, linetype = "solid"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +scale_fill_manual("Keys (Mutant-WT):  ", values = c("No-No" = "red1", "No-Yes" = "darkorange1", "Yes-No" = "darkturquoise", "Yes-Yes"="darkorchid1"))%>%
  +geom_bar(stat = "identity")%>%    # to create a stacked barchart
  +geom_text(aes(label=paste0(round(Express_freq_Perc$Percentage,0),"%")),size = 3, position = position_stack(vjust = 0.5), colour=c("white"), fontface='bold') %>%
  + coord_flip()

#p_Express_freq_Perc

# 
# ###!!!!!!!!!!!!UNHASH TO SAVE PLOT !!!!!!!!!!!!!!!!!!!!
#  
# # library(cowplot)
# save_plot("./Plots/Mouse/MGI GXD/All_phenotypes/Each_Phenotype.jpeg",
#             p_Express_freq_Perc ,base_height= 5.5 ,base_aspect_ratio = 2) 



################################################################################################################################
################################################################################################################################


# Discordant graph

Express_freq_Discordant <-filter(mgi.genepheno.tissue.ts28_wt_mutant,  Mutant.WT.Expression == "Yes-No" | Mutant.WT.Expression == "No-Yes") 

Express_freq_Discordant_PERC<-Express_freq_Discordant %>%
  dplyr::group_by(mp.description,Mutant.WT.Expression) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(Percentage = n/sum(n)*100)

names(Express_freq_Discordant_PERC)[names(Express_freq_Discordant_PERC)=="n"] <- "Frequency"


View(Express_freq_Discordant_PERC)

################################################################################################################################
################################################################################################################################

## Plot the Discordant Graph

p_Discordant<-ggplot(data=Express_freq_Discordant_PERC
          , aes(x= Express_freq_Discordant_PERC$mp.description
                ,y=Express_freq_Discordant_PERC$Percentage, fill=Express_freq_Discordant_PERC$Mutant.WT.Expression))
p_Discordant<-p_Discordant %>%
  
  + labs(x= "Top-Level Phenotype", y="Percentage(%) of the Expression Profile", 
         title= "The Discordant Gene Expression Profiles for Mutant and Wild-Type Mice.") %>%
  
  + guides(fill=guide_legend(title="Keys (Mutant-WT):  ")) %>%
  + theme(axis.title=element_text(size=14,face="bold"),
          plot.caption=element_text(face = "italic", size=12, hjust = 0),
          text = element_text(size=10),
          legend.position = "bottom",legend.direction = "horizontal",
          legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +scale_fill_manual("Key (Mutant-WT):  ", values = c("No-No" = "red1", "No-Yes" = "darkorange1", "Yes-No" = "darkturquoise", "Yes-Yes"="darkorchid1"))%>%
  + geom_bar(stat = "identity")%>%    # to create a stacked barchart
  
  +geom_text(aes(label=paste0(round(Express_freq_Discordant_PERC$Percentage,0),"%")),size = 3, position = position_stack(vjust = 0.5), colour=c("white"), fontface='bold') %>%
  
  + coord_flip()

#p_Discordant

###!!!!!!!!!!!!UNHASH TO SAVE PLOT !!!!!!!!!!!!!!!!!!!!
# 
# library(cowplot)
# save_plot("./Plots/Mouse/MGI GXD/All_phenotypes/p_Discordant.jpeg",
#          p_Discordant ,base_height= 5.5 ,base_aspect_ratio = 2) 
# 

#################################################################################################################################
################################################################################################################################

# Concordant graph


Express_freq_Concordant <-filter(mgi.genepheno.tissue.ts28_wt_mutant, Mutant.WT.Expression == "Yes-Yes" | Mutant.WT.Expression == "No-No") 

Express_freq_Concordant_PERC<-Express_freq_Concordant %>%
  dplyr::group_by(mp.description,Mutant.WT.Expression) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(Percentage = n/sum(n)*100)

names(Express_freq_Concordant_PERC)[names(Express_freq_Concordant_PERC)=="n"] <- "Frequency"


View(Express_freq_Concordant_PERC)

#!!!!!!!!!!!UNHASH TO SAVE FILE !!!!!!!!!!!!

#write.csv(Express_freq_Perc,'./Output_Files/go.Express_freq_Perc_Mice.csv')


p_Concordant<-ggplot(data=Express_freq_Concordant_PERC
                     , aes(x= Express_freq_Concordant_PERC$mp.description
                           ,y=Express_freq_Concordant_PERC$Percentage , fill=Express_freq_Concordant_PERC$Mutant.WT.Expression))
p_Concordant<-p_Concordant %>%
  
  + labs(x= "Top-Level Phenotype", y="Frequency of the Expression Profile", 
         title= "The Concordant Gene Expression Profiles for Mutant and Wild-Type Mice.")%>%
  
  + guides(fill=guide_legend(title="Keys (Mutant-WT):  ")) %>%
  + theme(axis.title=element_text(size=14,face="bold"),
          plot.caption=element_text(face = "italic", size=12, hjust = 0),
          text = element_text(size=10),
          legend.position = "bottom",legend.direction = "horizontal",
          legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +scale_fill_manual("Key (Mutant-WT):  ", values = c("No-No" = "red1", "No-Yes" = "darkorange1", "Yes-No" = "darkturquoise", "Yes-Yes"="darkorchid1"))%>%
  + geom_bar(stat = "identity")%>%    # to create a stacked barchart
  
  +geom_text(aes(label=paste0(round(Express_freq_Concordant_PERC$Percentage,0),"%")),size = 3, position = position_stack(vjust = 0.5), colour=c("white"), fontface='bold') %>%
  
  
  + coord_flip()



###!!!!!!!!!!!!UNHASH TO SAVE PLOT !!!!!!!!!!!!!!!!!!!!
#  
# library(cowplot)
# save_plot("./Plots/Mouse/MGI GXD/All_phenotypes/p_Concordant.jpeg",
#            p_Concordant ,base_height= 5.5 ,base_aspect_ratio = 2) 

################################################################################################################################
################################################################################################################################
 
##~~~~~~~~~~~~~~~~~~~~ DIFFERENTIAL GENE EXPRESSION BETWEEN KO AND WT FOR EACH PHENOTYPE ~~~~~~~~~~############################################


#####~~~~~~~~~~~~~~~~~~~~~CREATE FUNCTION TO CALCULATE THE FISHER TEST FOR ALL TISSUES~~~~~~~~~~~~~~~~~~#######

## A list of all mp.Descriptions

list_mouse_tissue<-mgi.genepheno.tissue.ts28_wt_mutant%>%
  select(mp.description)%>%
  distinct()

View(list_mouse_tissue)

##################################################################################################################################
##################################################################################################################################


########~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tissue_Gene_Mouse FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~################

Tissue_Gene_Mouse<-function(pheno.tissue, mgi.mouse.df){
  
  mgi.mouse.df$TS28.mutant.expression.detected <- as.factor(mgi.mouse.df$TS28.mutant.expression.detected)
  
  mgi.mouse.df$TS28.wt.expression.detected <- as.factor(mgi.mouse.df$TS28.wt.expression.detected)
  
  Tissue_pheno<-mgi.mouse.df%>%
    select(mp.description, MGI.ID, TS28.mutant.expression.detected, TS28.wt.expression.detected)%>%
    filter(mp.description==pheno.tissue)%>%
    distinct()
  
  Tissue_pheno_mutant<-Tissue_pheno%>%
    select_all()%>%
    dplyr::group_by(mp.description,TS28.mutant.expression.detected, .drop = FALSE ) %>%
    dplyr::summarise(n=n()) %>%
    distinct()
  
  Tissue_pheno_WT<-Tissue_pheno%>%
    select_all()%>%
    dplyr::group_by(mp.description,TS28.wt.expression.detected,.drop = FALSE  ) %>%
    dplyr::summarise(n=n()) %>%
    distinct()
  
  
  # CHANGE COLUMN NAME SO THAT THEY CAN BE MERGED
  
  names(Tissue_pheno_mutant)[names(Tissue_pheno_mutant)=="n"] <- "Mutant"
  
  names(Tissue_pheno_mutant)[names(Tissue_pheno_mutant)=="TS28.mutant.expression.detected"] <- "Gene Expression"
  
  
  
  names(Tissue_pheno_WT)[names(Tissue_pheno_WT)=="n"] <- "WT"
  
  names(Tissue_pheno_WT)[names(Tissue_pheno_WT)=="TS28.wt.expression.detected"] <- "Gene Expression"
  
  
  ## MERGE BOTH
  Tissue_pheno_WT_muta<-merge.data.frame(Tissue_pheno_mutant,Tissue_pheno_WT)
  
  #CLEAN the dataframe
  
  Tissue_pheno_WT_muta<-Tissue_pheno_WT_muta%>%
    select_all()%>%
    distinct(`Gene Expression`, Mutant, WT)
  # 
   #library(magrittr)
  Tissue_pheno_WT_muta_matrix<-Tissue_pheno_WT_muta %>%
     select(-`Gene Expression`)
  # 
  Tissue_pheno_WT_muta_matrix<-t(as.matrix(Tissue_pheno_WT_muta_matrix))
   
  View(Tissue_pheno_WT_muta_matrix)
  
   Tissue_pheno_WT_muta_fisher<-fisher.test(as.matrix(Tissue_pheno_WT_muta_matrix),alternative = "two.sided"  )
   p_values<-Tissue_pheno_WT_muta_fisher$p.value
   
  return(p_values)
}

#####################################################################################################

# TESTING FUNCTION# 
# Tissue_Gene_Mouse(pheno.tissue="cellular phenotype", 
#                           mgi.mouse.df=mgi.genepheno.tissue.ts28_wt_mutant)
# #  

##############~~~~~~~~~RUN FUNCTION ~~~~~~~~~~~~~~~~~~##################################

#### FOR LOOP TO ITERATE THROUGH A COLUMN OF ALL MP.DESCRIPTIONS
PVALUE_Tissues <- data.frame(list_mouse_tissue ,
                             P.VALUE=vector(length=25))


for(i in 1:nrow(list_mouse_tissue)) {
  ii <- list_mouse_tissue[i,1]     #FIRST COLUMN
  
    
  print(ii)
 
  ### THE OUTCOME OF THE FUNCTION IS SAVED AS PVAL WHICH IS THE FISHER TEST PVALUE 
  PVAL<-Tissue_Gene_Mouse(pheno.tissue=paste(ii), 
                           mgi.mouse.df=mgi.genepheno.tissue.ts28_wt_mutant)
  
  print(PVAL)                  # THE PVALUE IS PRINTED
 
  
  PVALUE_Tissues$P.VALUE[i]<-paste(PVAL)              # THE PVALUE IS ADDED INTO THE DF CREATED ABOVE
  
  
}


PVALUE_Tissues$Adjusted.P.value<-p.adjust(PVALUE_Tissues$P.VALUE, method="BH")


View(PVALUE_Tissues)


#write.csv(PVALUE_Tissues,'./Output_Files/Mice/PVALUE_Tissues_Physiologicalsystems.csv')


###############################################################################################################
###############################################################################################################
## FILTER TO ONLY HAVE SIGNIFICANT PHYSIOLOGICAL SYSTEMS

PVALUE_Tissues_sign<-PVALUE_Tissues%>%
  select_all()%>%
  filter(Adjusted.P.value<0.05)%>%
  distinct()

View(PVALUE_Tissues_sign)


## !!!!!!UNHASH TO SAVE FILE !!!!!!!!!!!!!!!!!!!!!!!!!!

# write.csv(PVALUE_Tissues_sign,'./Output_Files/Mice/PVALUE_Tissues_SIGNIFICANT.csv')


###############################################################################################################
##############################################################################################################




#########~~~~~~~~~~~~~~~DIFFERENTIAL GENE EXPRESSION BETWEEN SPECIFIC PHENOTYPES ~~~~~~~~~######################


################################################################################################################################
################################################################################################################################

#All phenotypes are stored into a list 

all_pheno<- mgi.genepheno.tissue.ts28_wt_mutant %>%
  select(mp.description)%>%
  distinct(mp.description)

all_pheno

################################################################################################################################
################################################################################################################################

# Phenotype specific
Phenotypes<-data.frame( mgi.genepheno.tissue.ts28_wt_mutant$MP.Term,
                        mgi.genepheno.tissue.ts28_wt_mutant$mp.description,
                        mgi.genepheno.tissue.ts28_wt_mutant$Mutant.WT.Expression )


names(Phenotypes)[names(Phenotypes)=="mgi.genepheno.tissue.ts28_wt_mutant.mp.description"] <- "mp.description"
names(Phenotypes)[names(Phenotypes)=="mgi.genepheno.tissue.ts28_wt_mutant.MP.Term"] <- "MP.Term"
names(Phenotypes)[names(Phenotypes)=="mgi.genepheno.tissue.ts28_wt_mutant.Mutant.WT.Expression"] <- "Mutant.WT.Expression"

View(Phenotypes)

##################################################################################################################################
#Function to create a plot for each phenotype.

pheno<- function(phenoname){
    

# Retrieve all phenotype names, match  'phenoname' to the mp.description  
  
  x_table<-Phenotypes%>%
    select(mp.description,MP.Term,Mutant.WT.Expression)%>%
    filter(mp.description== phenoname )%>%
    
    dplyr::group_by(MP.Term,Mutant.WT.Expression) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::mutate(Percentage = n/sum(n)*100) %>%
    
    distinct(MP.Term,Mutant.WT.Expression, Percentage, n)
  
  names(x_table)[names(x_table)=="n"] <- "Frequency"
  
  

  # A graphical representation of the data
  
  library(ggplot2)
  
  phenoname_p<-ggplot(data=x_table
                   , aes(x= x_table$MP.Term
                         ,y=x_table$Percentage, fill=x_table$Mutant.WT.Expression))

  phenoname_2<-  gsub("/", "-", phenoname , fixed = TRUE)
  
  phenoname_p<-phenoname_p %>%
    
    + labs(x= "Phenotype", y="Frequency of the Expression Profile", title= paste('The Gene Expression Profiles for Mutant and Wild-Type Mice for the Top-Level MP-Term: ', phenoname),
           tag= paste(phenoname),
           caption = paste("The data shows the gene expression profiles found in", phenoname ,"in either mutant or wild-type mice and the percentage (%) of frequency.")) %>%
    
    + guides(fill=guide_legend(title="Keys (Mutant-WT):  ")) %>%
    + theme(axis.title=element_text(size=14,face="bold"),
            plot.caption=element_text(face = "italic", size=12, hjust = 0),
            text = element_text(size=10),
            legend.position = "bottom",legend.direction = "horizontal",
            legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10),
            panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
    +scale_fill_manual("Keys (Mutant-WT):  ", values = c("No-No" = "red1", "No-Yes" = "darkorange1", "Yes-No" = "darkturquoise", "Yes-Yes"="darkorchid1"))%>%
    + geom_bar(stat = "identity")%>%    # to create a stacked barchart
    
    
    +geom_text(aes(label=paste0(round(x_table$Percentage,1),"%")),size = 3, position = position_stack(vjust = 0.5), colour=c("white"), fontface='bold') %>%
    
    + coord_flip()
  
  
  library(cowplot)
   save_plot(paste("./Plots/Mouse/MGI GXD/Phenotypes/",phenoname_2,".jpeg"),
                   phenoname_p ,base_height= 10 ,base_aspect_ratio = 2) 
  

}    

################################################################################################################################
################################################################################################################################

## TEST FUNCTION 

#pheno("vision/eye phenotype")

# RUN THE FUNCTION FOR ALL PHENOTYPES


apply(all_pheno,1,pheno)

###########################################################################################################
############################################################################################################################

## DIFFERENTIAL GENE EXPRESSION FOR EACH TISSUE VS MGI.ID

# create graphs for each tissue (EMAPA.TERM) with the gene id --MGI-ID


Tissue_freq_Perc<-mgi.genepheno.tissue.ts28_wt_mutant %>%
  dplyr::group_by(MGI.ID,EMAPA.Term,Mutant.WT.Expression) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(Percentage = n/sum(n)*100)


names(Tissue_freq_Perc)[names(Tissue_freq_Perc)=="n"] <- "Frequency"

View(Tissue_freq_Perc)



###############################################################################################################################
################################################################################################################################

# All tissues

all_tissues<- mgi.genepheno.tissue.ts28_wt_mutant %>%
  select(EMAPA.Term)%>%
  distinct(EMAPA.Term)



View(all_tissues)

#Each tissue and their expression percentage


tissue_table<-Tissue_freq_Perc%>%
  select(MGI.ID,EMAPA.Term,Mutant.WT.Expression)%>%
  filter(EMAPA.Term== "brain" )%>%
  dplyr::group_by(MGI.ID,Mutant.WT.Expression) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(Percentage = n/sum(n)*100) %>%
  
  distinct(MGI.ID,Mutant.WT.Expression, Percentage, n)

names(tissue_table)[names(tissue_table)=="n"] <- "Frequency"

View(tissue_table)





###############################################################################################################################
################################################################################################################################
# # 
#Function to create a plot for each phenotype.


Tissue<- function(tissuename){
  
  
  tissue_table<-Tissue_freq_Perc%>%
    select(MGI.ID,EMAPA.Term,Mutant.WT.Expression)%>%
    filter(EMAPA.Term== tissuename)%>%
    
    dplyr::group_by(MGI.ID,Mutant.WT.Expression) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::mutate(Percentage = n/sum(n)*100) %>%
    
    distinct(MGI.ID,Mutant.WT.Expression, Percentage, n)
  
  names(tissue_table)[names(tissue_table)=="n"] <- "Frequency"
  
  View(tissue_table)
  
  library(ggplot2)
  library(ggthemes)
  
  
  #dev.new(width=40, height=20) # to open in a new window
  
  tissue_p <-ggplot(data=tissue_table
                    , aes(x= tissue_table$MGI.ID
                          ,y=tissue_table$Percentage, fill=tissue_table$Mutant.WT.Expression))
  
  
  tissue_p <-tissue_p%>%
    
    + labs(x= "Gene ID", y="Frequency of the Expression Profile", title= paste('The Gene Expression Profiles for Mutant and Wild-Type Mice: ', tissuename),
           tag= paste(tissuename),
           caption = paste("The data shows the gene expression profiles found in", tissuename ,"in either mutant or wild-type mice and the percentage (%) of frequency.")) %>%
    
    + guides(fill=guide_legend(title="Keys (Mutant-WT):  ")) %>%
    + theme(axis.title=element_text(size=14,face="bold"),
            plot.caption=element_text(face = "italic", size=12, hjust = 0),
            text = element_text(size=10),
            legend.position = "bottom",legend.direction = "horizontal",
            legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10),
            panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
    +scale_fill_manual("Key (Mutant-WT):  ", values = c("No-No" = "coral2", "No-Yes" = "slateblue1", "Yes-No" = "lightgreen", "Yes-Yes"="plum1"))%>%
    + geom_bar(stat = "identity")%>%    # to create a stacked barchart
    
    
    +geom_text(aes(label=paste0(round(tissue_table$Percentage,1),"%")),size = 3, position = position_stack(vjust = 0.5), colour=c("white"), fontface='bold') %>%
    
    + coord_flip() 
  
  library(cowplot)
  save_plot(paste("./Plots/Mouse/MGI GXD/Tissue/",tissuename,".jpeg"),
            tissue_p ,base_height= 10 ,base_aspect_ratio = 2) 
  
  
}   




###############################################################################################################################
################################################################################################################################


## TEST FUNCTION 

#Tissue("brain")

### Loop through all tissues using the function


apply(all_tissues,1,Tissue)

###############################################################################################################################
################################################################################################################################
















################################################################################################################################
################################################################################################################################

### Add the frequency percentage of tissues (EMAPA.TERM) with the mp. descriptions and MGID


# Create table with the EMAPA.term (tissue), MGI.ID (gene), mp.description (top level) and gene expression

expression_gene<-mgi.genepheno.tissue.ts28_wt_mutant %>%
  select(EMAPA.Term,MGI.ID,mp.description, Mutant.WT.Expression)%>%
  distinct(EMAPA.Term,MGI.ID,mp.description, Mutant.WT.Expression)            #making sure there are no duplicates
  

View(expression_gene)


# Calculating the frequencies and percentage of gene expression for tissues

expression_genes2<- expression_gene %>%
  dplyr::group_by(EMAPA.Term,Mutant.WT.Expression) %>%           
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(Percentage = n/sum(n)*100)
  

# To join the two tables

expression_genes2 <- left_join(expression_genes2, expression_gene, by=c("EMAPA.Term" = "EMAPA.Term", "Mutant.WT.Expression" = "Mutant.WT.Expression") )

expression_genes2 %>%
  select(EMAPA.Term,MGI.ID,mp.description, Mutant.WT.Expression, n, Percentage)%>%
  distinct(EMAPA.Term,MGI.ID,mp.description, Mutant.WT.Expression, n, Percentage)     #making sure there are no duplicates
  
View(expression_genes2)



################################################################################################################################
################################################################################################################################


