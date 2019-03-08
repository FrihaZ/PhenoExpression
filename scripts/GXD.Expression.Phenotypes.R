###############################################################################################################################
###############################################################################################################################

### Project: PhenoExpression ##################################################################################################
### Script: GXD_Tidy_Data.R ###################################################################################################
### Purpose: Tidy gene expression data from GXD ###############################################################################
### Author: Pilar Cacheiro, Friha Zafar ####################################################################################################
### Date: 04/03/2019, 08/03/2019 ##########################################################################################################

###############################################################################################################################
###############################################################################################################################

## load packages ##############################################################################################################

library(dplyr);library(tidyr);library(stringr);library(readr); library(ggplot2); library(plyr)

################################################################################################################################
################################################################################################################################

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

## To create a Constigency table for the top ancestor node and the expression detection

library(ggplot2)
library(plyr)
library(dplyr)

## Join the two expression columns together

mgi.genepheno.tissue.ts28_wt_mutant<-mgi.genepheno.tissue.ts28_wt_mutant%>%
  mutate(Mutant.WT.Expression = paste0(TS28.mutant.expression.detected,"-",TS28.wt.expression.detected))

View(mgi.genepheno.tissue.ts28_wt_mutant)


#  Count the number of No-No, No-Yes, Yes-Yes, Yes-No


Express_freq<-ftable(mgi.genepheno.tissue.ts28_wt_mutant$MP.ID,
                    mgi.genepheno.tissue.ts28_wt_mutant$mp.description,
                    mgi.genepheno.tissue.ts28_wt_mutant$Mutant.WT.Expression)
Express_freq

################################################################################################################################
################################################################################################################################
