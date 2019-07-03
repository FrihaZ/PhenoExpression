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

library(dplyr);library(tidyr);library(stringr);library(readr); library(ggplot2); library(plyr)

################################################################################################################################
################################################################################################################################

#!!!! PLEASE CREATE THESE FILES INSIDE THE DIRECTORY THAT YOU WILL BE SAVING THIS SCRIPT INTO: !!!!!!!!!!!!
#   ./Plots/Tissue/
#   ./Plots/Phenotypes/

## import files ################################################################################################################

## import and tidy gene expression for TS28 mutant lines

gxd.ts28 <- read_delim("./data/MGIgeneExpressionQuery_20190304_mutant.txt.gz",delim = "\t") %>%
  select(MGI.Gene.ID,Structure,Detected) %>% # select columns of interest
  # we will try to always use MGI IDs for mouse genes / HGNC ids for human genes (these are stable identifiers)
  distinct(MGI.Gene.ID, Structure,Detected) %>% # retain unique rows
  arrange(Detected) %>% ## sort by this variable
  group_by(MGI.Gene.ID,Structure) %>% # group data by gene and tissue %>%
  summarise(Detected_all = paste0(Detected,collapse="|")) %>% ## collpase all the expression data availabe for a given
  # gene and tissue  in one cell
  mutate(Detected_all = ifelse(Detected_all =="No","No",
                              ifelse(Detected_all == "Yes","Yes","Ambiguous_NotSpecified_Contradictory"))) %>%
  rename(MGI.ID = MGI.Gene.ID, TS28.mutant.expression.detected = Detected_all) %>% ## rename variables
  select(MGI.ID,Structure,TS28.mutant.expression.detected) %>% # select columns of interest
  distinct(MGI.ID,Structure,TS28.mutant.expression.detected) %>% # check again for duplicated rows
  mutate(Gene.Anatomy = paste0(MGI.ID,"-",Structure)) # create a new variable (Gene - Anatomy (Structure) term)


## import MGI gene - phenotype data

mgi.genepheno <-  read_delim("./data/MGI_GenePheno.rpt",delim="\t",
                             col_names = c("Allelic.Composition","Allele.Symbols","Allele.ID",
                                           "Genetic.Background","MP.ID","PubMed.ID","MGI.ID","MGI.Genotype")) %>%
  # import file, add column names
  select(MGI.ID,MP.ID) %>% # select columns of interest
  distinct(MGI.ID, MP.ID) # unique gene-phenotype associations

## import phenotype ontology terms to anatomy terms mapping file

mapping <- read_delim("./data/MP_EMAPA.rpt",delim="\t",col_names = c("MP.ID","MP.Term","EMAPA.ID","EMAPA.Term"))


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

##### Wildtype Gene Expression Data Analysis ###################################################################################


################################################################################################################################
################################################################################################################################

## Importing file
## Cleaning the data ###########################################################################################################

gxd_wt.ts28<- read_delim("./data/MGIgeneExpressionQuery_20190304_wildtype.txt",delim = "\t")  %>%
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
View(mp.ancestor.nodes)

mp.toplevels<- read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/mp.toplevels.txt",delim = "\t") %>%
  rename(mp.ancestors=mp.term )  #renamed mp.term to mp.ancestors
################################################################################################################################
################################################################################################################################

## Join files to one

mp.toplevels<-left_join(mp.toplevels, mp.ancestor.nodes)   # the left_join keeps the data that is already in mp.toplevels, and add the mp.ancestors column

View(mp.toplevels)

################################################################################################################################
################################################################################################################################

## Link this back to the mgi.genepheno.tissue.ts28_wt_mutant

mgi.genepheno.tissue.ts28_wt_mutant <- left_join(mgi.genepheno.tissue.ts28_wt_mutant, mp.toplevels, by=c("MP.ID" = "mp.term") )
View(mgi.genepheno.tissue.ts28_wt_mutant)

################################################################################################################################
################################################################################################################################
library(ggplot2)
library(plyr)
library(dplyr)

## To create a Constigency table for the top ancestor node and the expression detection

## Join the two expression columns together 

mgi.genepheno.tissue.ts28_wt_mutant<-mgi.genepheno.tissue.ts28_wt_mutant%>%   # joined mutant and WT expressions
  mutate(Mutant.WT.Expression = paste0(TS28.mutant.expression.detected,"-",TS28.wt.expression.detected)) %>%
  dplyr::distinct(MGI.ID, MP.ID, MP.Term, mp.description,EMAPA.ID, EMAPA.Term, Gene.Anatomy,mp.ancestors, TS28.mutant.expression.detected, TS28.wt.expression.detected,Mutant.WT.Expression)  # To remove duplicates



View(mgi.genepheno.tissue.ts28_wt_mutant)



# Constingency table with gene names and count the number of No-No, No-Yes, Yes-Yes, Yes-No

Express_freq2<-ftable(mgi.genepheno.tissue.ts28_wt_mutant$MP.ID,
                      mgi.genepheno.tissue.ts28_wt_mutant$mp.description, 
                      mgi.genepheno.tissue.ts28_wt_mutant$Mutant.WT.Expression)

#Express_freq2



#  Count the number of No-No, No-Yes, Yes-Yes, Yes-No


Express_freq_Perc<-mgi.genepheno.tissue.ts28_wt_mutant %>%
  dplyr::group_by(mp.description,Mutant.WT.Expression) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(Percentage = n/sum(n)*100)


names(Express_freq_Perc)[names(Express_freq_Perc)=="n"] <- "Frequency"

View(Express_freq_Perc)

################################################################################################################################
################################################################################################################################



# A graphical representation of the data

library(ggplot2)
library(ggthemes)
#View(Express_freq_Perc)

#dev.new(width=60, height=30) # to open in a new window

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
  +geom_text(aes(label=paste0(round(Express_freq_Perc$Percentage,1),"%")),size = 2, position = position_stack(vjust = 0.5), colour=c("black"), fontface='bold') %>%
  + coord_flip()

p_Express_freq_Perc

library(cowplot)
save_plot("./Plots/Mouse/All_phenotypes/Each_Phenotype.jpeg",
          p_Express_freq_Perc ,base_height= 5.5 ,base_aspect_ratio = 2) 


#ggsave(filename=paste("./Plots/Phenotypes/Each_Phenotype.png",sep=" "), limitsize = TRUE)


+ theme(axis.title=element_text(size=14,face="bold"),
        plot.caption=element_text(face = "italic", size=12, hjust = 0),
        text = element_text(size=12),
        axis.text.y= element_text(size=12),
        axis.text.x= element_text(size=12))%>%


################################################################################################################################
################################################################################################################################


# Discordant graph

Express_freq_Discordant <-filter(mgi.genepheno.tissue.ts28_wt_mutant,  Mutant.WT.Expression == "Yes-No" | Mutant.WT.Expression == "No-Yes") 

Express_freq_Discordant_PERC<-Express_freq_Discordant %>%
  dplyr::group_by(mp.description,Mutant.WT.Expression) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(Percentage = n/sum(n)*100)

names(Express_freq_Discordant_PERC)[names(Express_freq_Discordant_PERC)=="n"] <- "Frequency"


#  select(n, mp.description, Mutant.WT.Expression ) %>%
#   mutate(Perc= round((Exp,"%")%>%
#   distinct(Perc, Freq, mp.description, Mutant.WT.Expression)# unique rows
# 

View(Express_freq_Discordant_PERC)



dev.new(width=40, height=20) # to open in a new window

p_Discordant<-ggplot(data=Express_freq_Discordant_PERC
          , aes(x= Express_freq_Discordant_PERC$mp.description
                ,y=Express_freq_Discordant_PERC$Percentage, fill=Express_freq_Discordant_PERC$Mutant.WT.Expression))
p_Discordant %>%
  
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
  
  +geom_text(aes(label=paste0(round(Express_freq_Discordant_PERC$Percentage,1),"%")),size = 3, position = position_stack(vjust = 0.5), colour=c("white"), fontface='bold') %>%
  
  + coord_flip()


# Concordant graph



Express_freq_Concordant <-filter(mgi.genepheno.tissue.ts28_wt_mutant, Mutant.WT.Expression == "Yes-Yes" | Mutant.WT.Expression == "No-No") 

Express_freq_Concordant_PERC<-Express_freq_Concordant %>%
  dplyr::group_by(mp.description,Mutant.WT.Expression) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(Percentage = n/sum(n)*100)

names(Express_freq_Concordant_PERC)[names(Express_freq_Concordant_PERC)=="n"] <- "Frequency"


View(Express_freq_Concordant_PERC)
write.csv(Express_freq_Perc,'./Output_Files/go.Express_freq_Perc_Mice.csv')



dev.new(width=40, height=20) # to open in a new window

p_Concordant<-ggplot(data=Express_freq_Concordant_PERC
                     , aes(x= Express_freq_Concordant_PERC$mp.description
                           ,y=Express_freq_Concordant_PERC$Percentage , fill=Express_freq_Concordant_PERC$Mutant.WT.Expression))
p_Concordant %>%
  
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
  
  +geom_text(aes(label=paste0(round(Express_freq_Concordant_PERC$Percentage,1),"%")),size = 3, position = position_stack(vjust = 0.5), colour=c("white"), fontface='bold') %>%
  
  
  + coord_flip()








################################################################################################################################
################################################################################################################################

########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~EACH PHYSIOLOGICAL SYSTEM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##############

################################################################################################################################
################################################################################################################################


# ## Look deeper into the Phenotype testing the function outside
# 
# mgi.genepheno.tissue.ts28_wt_mutant$TS28.mutant.expression.detected <- as.factor(mgi.genepheno.tissue.ts28_wt_mutant$TS28.mutant.expression.detected)
# 
# mgi.genepheno.tissue.ts28_wt_mutant$TS28.wt.expression.detected <- as.factor(mgi.genepheno.tissue.ts28_wt_mutant$TS28.wt.expression.detected)
# 
# Cellular_pheno<-mgi.genepheno.tissue.ts28_wt_mutant%>%
#   select(mp.description, MGI.ID, TS28.mutant.expression.detected, TS28.wt.expression.detected)%>%
#   filter(mp.description=="skeleton phenotype")%>%
#   distinct()
# 
# 
# Cellular_pheno_mutant<-Cellular_pheno%>%
#   select_all()%>%
#   dplyr::group_by(mp.description,TS28.mutant.expression.detected,.drop = FALSE) %>%
#   dplyr::summarise(n=n()) %>%
#   distinct()
# 
# View(Cellular_pheno_mutant)
# 
# Cellular_pheno_WT<-Cellular_pheno%>%
#   select_all()%>%
#   dplyr::group_by(mp.description,TS28.wt.expression.detected, .drop = FALSE ) %>%
#   dplyr::summarise(n=n()) %>%
#   distinct()
# View(Cellular_pheno_WT)
# 
# 
# # CHANGE COLUMN NAME SO THAT THEY CAN BE MERGED
# 
# names(Cellular_pheno_mutant)[names(Cellular_pheno_mutant)=="n"] <- "Mutant"
# 
# names(Cellular_pheno_mutant)[names(Cellular_pheno_mutant)=="TS28.mutant.expression.detected"] <- "Gene Expression"
# 
# 
# 
# names(Cellular_pheno_WT)[names(Cellular_pheno_WT)=="n"] <- "WT"
# 
# names(Cellular_pheno_WT)[names(Cellular_pheno_WT)=="TS28.wt.expression.detected"] <- "Gene Expression"
# 
# 
# ## MERGE BOTH
# Cellular_pheno_WT_muta<-merge.data.frame(Cellular_pheno_mutant,Cellular_pheno_WT)
# 
# #CLEAN the dataframe
# 
# Cellular_pheno_WT_mutant<-Cellular_pheno_WT_muta%>%
#   select_all()%>%
#   distinct(`Gene Expression`, Mutant, WT)
#   
# 
# View(Cellular_pheno_WT_muta)
# 
# 
# #library(magrittr)
# Cellular_pheno_WT_muta_matrix<-Cellular_pheno_WT_mutant %>%
# #  set_rownames(.$`Gene Expression`)%>%
#   select(-`Gene Expression`)
# Cellular_pheno_WT_muta_matrix<-as.matrix(Cellular_pheno_WT_muta_matrix)
# 
# View(Cellular_pheno_WT_muta_matrix)
# 
# #Fisher Test is conducted
# 
# if (nrow(Cellular_pheno_WT_muta_matrix)==2){
#   
#   print(fisher.test(as.matrix(Cellular_pheno_WT_muta_matrix)))
# } else {
#   return ("0")
# }



###############################################################################################################
###############################################################################################################

#####~~~~~~~~~~~~~~~~~~~~~CREATE FUNCTION TO CALCULATE THE FISHER TEST FOR ALL TISSUES~~~~~~~~~~~~~~~~~~#######

## A list of all mp.Descriptions

list_mouse_tissue<-mgi.genepheno.tissue.ts28_wt_mutant%>%
  select(mp.description)%>%
  distinct()

View(list_mouse_tissue)


#### ~~~~~~~FUNCTION~~~~~~~~~~~#

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
  Tissue_pheno_WT_muta_matrix<-as.matrix(Tissue_pheno_WT_muta_matrix)
   
  View(Tissue_pheno_WT_muta_matrix)
  
   Tissue_pheno_WT_muta_fisher<-fisher.test(as.matrix(Tissue_pheno_WT_muta_matrix))
   p_values<-Tissue_pheno_WT_muta_fisher$p.value
   
  return(p_values)
}
# TESTING FUNCTION
# 
# Tissue_Gene_Mouse(pheno.tissue="skeleton phenotype", 
#                          mgi.mouse.df=mgi.genepheno.tissue.ts28_wt_mutant)
#  


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

View(PVALUE_Tissues)

#write.csv(PVALUE_Tissues,'./Output_Files/Mice/PVALUE_Tissues_Physiologicalsystems.csv')


###############################################################################################################
###############################################################################################################
## FILTER TO ONLY HAVE SIGNIFICANT PHYSIOLOGICAL SYSTEMS

PVALUE_Tissues_sign<-PVALUE_Tissues%>%
  select_all()%>%
  filter(P.VALUE<0.05)%>%
  distinct()

View(PVALUE_Tissues_sign)
#write.csv(PVALUE_Tissues_sign,'./Output_Files/Mice/PVALUE_Tissues_SIGNIFICANT.csv')


###############################################################################################################
###############################################################################################################












###############################################################################################################
###############################################################################################################

## OLD CHI-SQUARED TEST

# Creating a variable for each Mutant and WT count

Yes_Yes <- sum(mgi.genepheno.tissue.ts28_wt_mutant$Mutant.WT.Expression == "Yes-Yes")
Yes_Yes

Yes_No <- sum(mgi.genepheno.tissue.ts28_wt_mutant$Mutant.WT.Expression == "Yes-No")
Yes_No

No_Yes <- sum(mgi.genepheno.tissue.ts28_wt_mutant$Mutant.WT.Expression == "No-Yes")
No_Yes

No_No <- sum(mgi.genepheno.tissue.ts28_wt_mutant$Mutant.WT.Expression == "No-No")
No_No

# Creating Matrix

WT_MUTANT = matrix( c(Yes_Yes,No_Yes,Yes_No, No_No ), # the data elements 
            nrow=2,              # number of rows 
            ncol=2,              # number of columns 
            byrow = TRUE)


row.names(WT_MUTANT) <- c("WT_Yes", "WT_No")   # row names

colnames(WT_MUTANT) <- c("Mutant_Yes", "Mutant_No")   #column names


View(WT_MUTANT)

################################################################################################################################
################################################################################################################################


# Carry out Chi-squared test

chisq<- chisq.test(WT_MUTANT)

chisq$expected
chisq$observed

chisq$p.value 

# = 0.002184877 if alpha value is 0.05 then the p-value is <alpha therefore the H0 can be rejected, and the H1 can be accepted. 
# The observed values did not occur by chance 



################################################################################################################################
################################################################################################################################

# Phenotype specific
Phenotypes<-data.frame( mgi.genepheno.tissue.ts28_wt_mutant$MP.Term,
                        mgi.genepheno.tissue.ts28_wt_mutant$mp.description,mgi.genepheno.tissue.ts28_wt_mutant$Mutant.WT.Expression )


names(Phenotypes)[names(Phenotypes)=="mgi.genepheno.tissue.ts28_wt_mutant.mp.description"] <- "mp.description"
names(Phenotypes)[names(Phenotypes)=="mgi.genepheno.tissue.ts28_wt_mutant.MP.Term"] <- "MP.Term"
names(Phenotypes)[names(Phenotypes)=="mgi.genepheno.tissue.ts28_wt_mutant.Mutant.WT.Expression"] <- "Mutant.WT.Expression"

View(Phenotypes)


################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################

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
  library(ggthemes)
  
  
  #dev.new(width=40, height=20) # to open in a new window
  
  phenoname_p<-ggplot(data=x_table
                   , aes(x= x_table$MP.Term
                         ,y=x_table$Percentage, fill=x_table$Mutant.WT.Expression))

  phenoname_2<-  gsub("/", "-", phenoname , fixed = TRUE)
  
  phenoname_p %>%
    
    + labs(x= "Top-Level Phenotype", y="Frequency of the Expression Profile", title= paste('The Gene Expression Profiles for Mutant and Wild-Type Mice: ', phenoname),
           tag= paste(phenoname),
           caption = paste("The data shows the gene expression profiles found in", phenoname ,"in either mutant or wild-type mice and the percentage (%) of frequency.")) %>%
    
    + guides(fill=guide_legend(title="Keys (Mutant-WT):  ")) %>%
    + theme(plot.caption=element_text(face = "italic", size=10, hjust = 0),
            text = element_text(size=10),
            legend.position = "bottom",legend.direction = "horizontal",
            legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
            axis.text.x = element_text(angle = 90, hjust = 1),
            panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
    +scale_fill_manual("Keys (Mutant-WT):  ", values = c("No-No" = "red1", "No-Yes" = "darkorange1", "Yes-No" = "darkturquoise", "Yes-Yes"="darkorchid1"))%>%
    + geom_bar(stat = "identity")%>%    # to create a stacked barchart
    
    
    +geom_text(aes(label=paste0(round(x_table$Percentage,1),"%")),size = 2, position = position_stack(vjust = 0.5), colour=c("black"), fontface='bold') %>%
    
    + coord_flip() %>%
    + ggsave(filename=paste("./Plots/Phenotypes/",phenoname_2,".png",sep=" "), limitsize = TRUE)
  
}    
    ################################################################################################################################
    ################################################################################################################################

# Use the function for each phenotype


#pheno("vision/eye phenotype")


#apply(all_pheno,1,pheno)






###############################################################################################################################
################################################################################################################################

# create graphs for each tissue (EMAPA.TERM) with the gene id --MGI-ID


Tissue_freq_Perc<-mgi.genepheno.tissue.ts28_wt_mutant %>%
  dplyr::group_by(EMAPA.ID,EMAPA.Term,Mutant.WT.Expression) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(Percentage = n/sum(n)*100)


names(Tissue_freq_Perc)[names(Tissue_freq_Perc)=="n"] <- "Frequency"

View(Tissue_freq_Perc)


#All tissues 

dev.new(width=40, height=20) # to open in a new window

Tissue_p<-ggplot(data=Tissue_freq_Perc
                     , aes(x= Tissue_freq_Perc$EMAPA.Term
                           ,y=Tissue_freq_Perc$Percentage , fill=Tissue_freq_Perc$Mutant.WT.Expression))
Tissue_p %>%
  
  + labs(x= "Tissue", y="Frequency of the Expression Profile", title= "The Gene Expression Profiles for Mutant and Wild-Type Mice Found in Different Tissues.",
         tag="E",
         caption = "The data shows the concordant gene expression profiles found in either mutant or wild-type mice and the percentage (%) of frequency.") %>%
  
  + guides(fill=guide_legend(title="Keys (Mutant-WT):  ")) %>%
  + theme(plot.caption=element_text(face = "italic", size=12, hjust = 0),
          text = element_text(size=12),
          legend.position = "bottom",legend.direction = "horizontal",
          legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +scale_fill_manual("Key (Mutant-WT):  ", values = c("No-No" = "coral2", "No-Yes" = "slateblue1", "Yes-No" = "lightgreen", "Yes-Yes"="plum1"))%>%
  + geom_bar(stat = "identity")%>%    # to create a stacked barchart
  
  +geom_text(aes(label=paste0(round(Tissue_freq_Perc$Percentage,1),"%")),size = 2, position = position_stack(vjust = 0.5), colour=c("black"), fontface='bold') %>%
  
  
  + coord_flip()




###############################################################################################################################
################################################################################################################################

# All tissues

all_tissues<- mgi.genepheno.tissue.ts28_wt_mutant %>%
  select(EMAPA.Term)%>%
  distinct(EMAPA.Term)

#all_tissues <- split(all_tissues, seq(nrow))


View(all_tissues)

#Each tissue and their expression percentage


tissue_table<-Tissue_freq_Perc%>%
  select(EMAPA.ID,EMAPA.Term,Mutant.WT.Expression)%>%
  filter(EMAPA.Term== "brain" )%>%
  
  dplyr::group_by(EMAPA.ID,Mutant.WT.Expression) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(Percentage = n/sum(n)*100) %>%
  
  distinct(EMAPA.ID,Mutant.WT.Expression, Percentage, n)

names(tissue_table)[names(tissue_table)=="n"] <- "Frequency"

View(tissue_table)



# Plot for each tissue:

dev.new(width=40, height=20) # to open in a new window


Tissue_plot<-ggplot(data=tissue_table
                    , aes(x= tissue_table$EMAPA.ID
                          ,y=tissue_table$Percentage, fill=tissue_table$Mutant.WT.Expression))


Tissue_plot %>%
  
  + labs(x= "Gene ID", y="Frequency of the Expression Profile", title= paste('The Gene Expression Profiles for Mutant and Wild-Type Mice: '),
         tag= "brain",  #paste(phenoname),
         caption = paste("The data shows the gene expression profiles found in brain in either mutant or wild-type mice and the percentage (%) of frequency.")) %>%
  
  + guides(fill=guide_legend(title="Keys (Mutant-WT):  ")) %>%
  + theme(plot.caption=element_text(face = "italic", size=10, hjust = 0),
          text = element_text(size=10),
          legend.position = "bottom",legend.direction = "horizontal",
          legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +scale_fill_manual("Key (Mutant-WT):  ", values = c("No-No" = "coral2", "No-Yes" = "slateblue1", "Yes-No" = "lightgreen", "Yes-Yes"="plum1"))%>%
  + geom_bar(stat = "identity")%>%    # to create a stacked barchart
  
  
  +geom_text(aes(label=paste0(round(tissue_table$Percentage,1),"%")),size = 2, position = position_stack(vjust = 0.5), colour=c("black"), fontface='bold') %>%
  
  + coord_flip() 
  
  #+ ggsave(filename=paste("./Plots/",phenoname_2,".png",sep=" "), limitsize = TRUE)

   


###############################################################################################################################
################################################################################################################################
# # 
#Function to create a plot for each phenotype.

Tissue<- function(tissuename){
  
  
  tissue_table<-Tissue_freq_Perc%>%
    select(EMAPA.ID,EMAPA.Term,Mutant.WT.Expression)%>%
    filter(EMAPA.Term== tissuename)%>%
    
    dplyr::group_by(EMAPA.ID,Mutant.WT.Expression) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::mutate(Percentage = n/sum(n)*100) %>%
    
    distinct(EMAPA.ID,Mutant.WT.Expression, Percentage, n)
  
  names(tissue_table)[names(tissue_table)=="n"] <- "Frequency"
  
  View(tissue_table)
  
  library(ggplot2)
  library(ggthemes)
  
  
  #dev.new(width=40, height=20) # to open in a new window
  
  tissue_p <-ggplot(data=tissue_table
                      , aes(x= tissue_table$EMAPA.ID
                            ,y=tissue_table$Percentage, fill=tissue_table$Mutant.WT.Expression))
  
  
  tissue_p %>%
    
    + labs(x= "Gene ID", y="Frequency of the Expression Profile", title= paste('The Gene Expression Profiles for Mutant and Wild-Type Mice: ', tissuename),
           tag= paste(tissuename),
           caption = paste("The data shows the gene expression profiles found in", tissuename ,"in either mutant or wild-type mice and the percentage (%) of frequency.")) %>%
    
    + guides(fill=guide_legend(title="Keys (Mutant-WT):  ")) %>%
    + theme(plot.caption=element_text(face = "italic", size=10, hjust = 0),
            text = element_text(size=10),
            legend.position = "bottom",legend.direction = "horizontal",
            legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
            axis.text.x = element_text(angle = 90, hjust = 1),
            panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
    +scale_fill_manual("Key (Mutant-WT):  ", values = c("No-No" = "coral2", "No-Yes" = "slateblue1", "Yes-No" = "lightgreen", "Yes-Yes"="plum1"))%>%
    + geom_bar(stat = "identity")%>%    # to create a stacked barchart
    
    
    +geom_text(aes(label=paste0(round(tissue_table$Percentage,1),"%")),size = 2, position = position_stack(vjust = 0.5), colour=c("black"), fontface='bold') %>%
    
    + coord_flip() %>%
    + ggsave(filename=paste("./Plots/Tissue/",tissuename,".png",sep=" "), limitsize = TRUE)
  
}   

Tissue("brain")




###############################################################################################################################
################################################################################################################################

#Loop through all tissues 


apply(all_tissues,1,Tissue)


###############################################################################################################################
################################################################################################################################



## Top level phenotype with Discordant expression at 100% and then <100% 


Discordant_expression <-filter(mgi.genepheno.tissue.ts28_wt_mutant,  Mutant.WT.Expression == "Yes-No" | Mutant.WT.Expression == "No-Yes") 


Discordant_expression_freq<-Discordant_expression %>%
  select(MP.Term,mp.description, Mutant.WT.Expression)%>%
  dplyr::group_by(mp.description,Mutant.WT.Expression) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(Percentage = n/sum(n)*100)

names(Discordant_expression_freq)[names(Discordant_expression_freq)=="n"] <- "Frequency"

View(Discordant_expression_freq)


# Filter for Full Discordance (100%)

Discordant_expression_hundred <-filter(Discordant_expression_freq,  Percentage == 100.00) 

View(Discordant_expression_hundred)



## Plot graph


library(ggplot2)
library(ggthemes)


dev.new(width=40, height=20) # to open in a new window


Discordant_expression_hundred_p <-ggplot(data=Discordant_expression_hundred
                  , aes(x= Discordant_expression_hundred$mp.description
                        ,y=Discordant_expression_hundred$Percentage, fill=Discordant_expression_hundred$Mutant.WT.Expression))


Discordant_expression_hundred_p %>%
  
  + labs(x= "Top Level Phenotype", y="Percentage (%) of the Expression Profile", title= paste('The Gene Expression Profiles for Mutant and Wild-Type Mice with 100% Discordance.'),
         tag= "F",
         caption = paste("The data shows the gene expression profiles found in either mutant or wild-type mice and the percentage (%) of frequency.")) %>%
  
  + guides(fill=guide_legend(title="Keys (Mutant-WT):  ")) %>%
  + theme(plot.caption=element_text(face = "italic", size=10, hjust = 0),
          text = element_text(size=10),
          legend.position = "bottom",legend.direction = "horizontal",
          legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +scale_fill_manual("Key (Mutant-WT):  ", values = c("No-No" = "coral2", "No-Yes" = "slateblue1", "Yes-No" = "lightgreen", "Yes-Yes"="plum1"))%>%
  + geom_bar(stat = "identity")%>%    # to create a stacked barchart
  
  
  +geom_text(aes(label=paste0(round(Discordant_expression_hundred$Percentage,1),"%")),size = 2, position = position_stack(vjust = 0.5), colour=c("black"), fontface='bold') %>%
  
  + coord_flip() 



################################################################################################################################
################################################################################################################################


## Mp.term instead of top level phenotype


Discordant_expression_mpterm <-filter(mgi.genepheno.tissue.ts28_wt_mutant,  Mutant.WT.Expression == "Yes-No" | Mutant.WT.Expression == "No-Yes") 


Discordant_expression_mpterm<-Discordant_expression_mpterm %>%
  select(MP.Term,mp.description, Mutant.WT.Expression)%>%
  dplyr::group_by(MP.Term,Mutant.WT.Expression) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(Percentage = n/sum(n)*100)

names(Discordant_expression_mpterm)[names(Discordant_expression_mpterm)=="n"] <- "Frequency"

View(Discordant_expression_mpterm)


# Full Discordance (100%)

Discordant_expression_hundred_mpterm <-filter(Discordant_expression_mpterm,  Percentage == 100.00) 

View(Discordant_expression_hundred_mpterm)



Discordant_expression_hundred_Mp.term_p <-ggplot(data=Discordant_expression_hundred_mpterm
                                         , aes(x= Discordant_expression_hundred_mpterm$MP.Term
                                               ,y=Discordant_expression_hundred_mpterm$Percentage, fill=Discordant_expression_hundred_mpterm$Mutant.WT.Expression))


Discordant_expression_hundred_Mp.term_p %>%
  
  + labs(x= "Top Level Phenotype", y="Percentage (%) of the Expression Profile", title= paste('The Gene Expression Profiles for Mutant and Wild-Type Mice.'),
         tag= "H",
         caption = paste("The data shows the gene expression profiles found in either mutant or wild-type mice and the percentage (%) of frequency.")) %>%
  
  + guides(fill=guide_legend(title="Keys (Mutant-WT):  ")) %>%
  + theme(plot.caption=element_text(face = "italic", size=10, hjust = 0),
          text = element_text(size=10),
          legend.position = "bottom",legend.direction = "horizontal",
          legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +scale_fill_manual("Key (Mutant-WT):  ", values = c("No-No" = "coral2", "No-Yes" = "slateblue1", "Yes-No" = "lightgreen", "Yes-Yes"="plum1"))%>%
  + geom_bar(stat = "identity")%>%    # to create a stacked barchart
  
  
  +geom_text(aes(label=paste0(round(Discordant_expression_hundred_mpterm$Percentage,1),"%")),size = 2, position = position_stack(vjust = 0.5), colour=c("black"), fontface='bold') %>%
  
  + coord_flip() 





################################################################################################################################
################################################################################################################################

# Less than 100% discordance using mp.description

Discordant_expression_less_hundred <-filter(Discordant_expression_freq,  Percentage < 100.00) 

View(Discordant_expression_less_hundred)


## Plot 

library(ggplot2)
library(ggthemes)


dev.new(width=60, height=30) # to open in a new window


Discordant_expression_less_hundred_p <-ggplot(data=Discordant_expression_less_hundred
                                         , aes(x= Discordant_expression_less_hundred$mp.description
                                               ,y=Discordant_expression_less_hundred$Percentage, fill=Discordant_expression_less_hundred$Mutant.WT.Expression))


Discordant_expression_less_hundred_p %>%
  
  + labs(x= "Top Level Phenotype", y="Percentage (%) of the Expression Profile", title= paste('The Gene Expression Profiles for Mutant and Wild-Type Mice with <100 Discordance.'),
         tag= "G",
         caption = paste("The data shows the gene expression profiles found in either mutant or wild-type mice and the percentage (%) of frequency.")) %>%
  
  + guides(fill=guide_legend(title="Keys (Mutant-WT):  ")) %>%
  + theme(plot.caption=element_text(face = "italic", size=10, hjust = 0),
          text = element_text(size=10),
          legend.position = "bottom",legend.direction = "horizontal",
          legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +scale_fill_manual("Key (Mutant-WT):  ", values = c("No-No" = "coral2", "No-Yes" = "slateblue1", "Yes-No" = "lightgreen", "Yes-Yes"="plum1"))%>%
  + geom_bar(stat = "identity")%>%    # to create a stacked barchart
  
  
  +geom_text(aes(label=paste0(round(Discordant_expression_less_hundred$Percentage,1),"%")),size = 2, position = position_stack(vjust = 0.5), colour=c("black"), fontface='bold') %>%
  
  + coord_flip() 



################################################################################################################################
################################################################################################################################


## Discordant Expression using the MP.TERM



Discordant_expression_mp.term <-filter(mgi.genepheno.tissue.ts28_wt_mutant,  Mutant.WT.Expression == "Yes-No" | Mutant.WT.Expression == "No-Yes") 


Discordant_expression_mp.term<-Discordant_expression_mp.term %>%
  select(MP.Term, Mutant.WT.Expression)%>%
  dplyr::group_by(MP.Term,Mutant.WT.Expression) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(Percentage = n/sum(n)*100)

names(Discordant_expression_mp.term)[names(Discordant_expression_mp.term)=="n"] <- "Frequency"

View(Discordant_expression_mp.term)
Discordant_expression_mp.term <-filter(Discordant_expression_mp.term,  Percentage < 100.00) 

View(Discordant_expression_mp.term)


####  NO MP.TERM WITH LESS THAN 100% DISCORDANCE FOUND



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



