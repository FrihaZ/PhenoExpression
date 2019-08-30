
###############################################################################################################################
###############################################################################################################################

### Project: PhenoExpression ##################################################################################################
### Script: Human_vs_EA_DMDD_Mouse.R ###################################################################################################
### Purpose: To compare between the human and the mouse genes (Expression Atlas-DMDD DATASET) that may exhibit similar gene expressions ########################################################
### Author: Friha Zafar ####################################################################################################
### Date: 18/06/2019 ##########################################################################################################

###############################################################################################################################
###############################################################################################################################

## load packages ##############################################################################################################

library(dplyr);library(tidyr);library(readr);library(tidyverse);library(ggpubr)

################################################################################################################################
################################################################################################################################

# Load Files

human.mouse.orthologs<-read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/human.mouse.orthologs.txt",delim = "\t", skip = 1)

human_0.1_pheno<-read_csv("D:/MSC RESEARCH PROJECT/Output_Files/Human/Expression_0.1/human_genes_TPM_0.1_pheno.csv")

human_1_pheno<-read_csv("D:/MSC RESEARCH PROJECT/Output_Files/Human/Expression_1/human_genes_TPM_1_pheno.csv")

mgi.genepheno.tissue.ts28_wt_mutant<-read_csv("D:/MSC RESEARCH PROJECT/Output_Files/Mice/mgi.genepheno.tissue.ts28_wt_mutant.csv")

GTEX.Tissues.to.HPO<-read_csv("D:/MSC RESEARCH PROJECT/Output_Files/GTEX.Tissues.to.HPO.csv")

hp.mp.toplevelmapping<-read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/hp-mp-toplevelmapping.txt",delim = "\t")

mgi.phenotypes<-read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/mgi.phenotypes.txt",delim = "\t")


## load 'clean_df' function:

source(file="D:/MSC RESEARCH PROJECT/Scripts/Human+Mouse/clean_df_FUNCTION.R")




############################################################################################################################################################
################################################################################################################################################

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~Clean the Human data to have only the genes that are orthologs of mice ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###########################

## Use function to clean the dfs

human_0.1_pheno_cleaned <-clean_df(human_0.1_pheno)

names(human_0.1_pheno_cleaned)[names(human_0.1_pheno_cleaned)=="Expression"] <- "TPM >/= 0.1 TPM"

View(human_0.1_pheno_cleaned)

#############

human_1_pheno_cleaned <-clean_df(human_1_pheno)


names(human_1_pheno_cleaned)[names(human_1_pheno_cleaned)=="Expression"] <- "TPM >/= 1 TPM"


View(human_1_pheno_cleaned)

########################################################################################################
##################################################################################################################

## Join the two dfs

All_genes_expression<- right_join(human_0.1_pheno_cleaned,
                                  human_1_pheno_cleaned, by=c("HGNC.ID"="HGNC.ID", 
                                                                    "GTEX.tissues"="GTEX.tissues",
                                                                    "HPO.superclass.description"="HPO.superclass.description"))

View(All_genes_expression)


## !!!!!!!!!!!!! UNHASH TO SAVE FILE !!!!!!!!!!!!!!1


#write.csv(All_genes_expression,"D:/MSC RESEARCH PROJECT/Output_Files/Human/All_genes_expression.csv") 


##############################################################################################################################
############################################################################################################################

## Combine the Human expressions and the Mouse orthologs


Human_genes_expression_orthologs<-left_join(All_genes_expression,
                                             human.mouse.orthologs, by=c("HGNC.ID"="HGNC.ID"))

View(Human_genes_expression_orthologs)



##############################################################################################################################
############################################################################################################################

## Mouse Knockout expression cleaned

Mice.mutant.expression.clean<- mgi.genepheno.tissue.ts28_wt_mutant%>%
  select(MGI.ID, mp.description,TS28.mutant.expression.detected, MP.ID, mp.ancestors, MP.Term)%>%
  distinct()

#View(mgi.genepheno.tissue.ts28_wt_mutant)


##############################################################################################################################
############################################################################################################################

## Join mouse KO mouse expression with the phenotypes associated to the knockout gene 
## mgi.phenotypes: "MP.DESCRIPTION" ==  Mice.mutant.expression.clean: "MP.Term"

names(mgi.phenotypes)[names(mgi.phenotypes)=="MP.DESCRIPTION"] <- "MP.Term"


 
Mice.mutant.all.expression.phenotypes<-merge(Mice.mutant.expression.clean,
                                              mgi.phenotypes)

View(Mice.mutant.all.expression.phenotypes)










########################################################################################################
#######################################################################################################################

### Filter to only count genes with an expression

Mice.mutant.expression.phenotypes<-Mice.mutant.all.expression.phenotypes%>%
  select_all()%>%
  drop_na()%>%
  distinct()

## change column name
names(Mice.mutant.expression.phenotypes)[names(Mice.mutant.expression.phenotypes)=="mp.ancestors"] <- "mp.top.term"

#View(Mice.mutant.expression.phenotypes)

##############################################################################################################################
############################################################################################################################

## Join mouse KO mouse expression with the top level mapping of the hpo descriptions

hp.mp.toplevelmapping<-hp.mp.toplevelmapping%>%
  select(hp.top.term, hp.top.description, mp.top.term, mp.top.description )%>%
  distinct()

#View(hp.mp.toplevelmapping)


Mouse.Expression.pheno.hpo<-merge(Mice.mutant.expression.phenotypes, hp.mp.toplevelmapping)


names(Mouse.Expression.pheno.hpo)[names(Mouse.Expression.pheno.hpo)=="hp.top.description"] <- "HPO.superclass.description"


#### JOIN THE HGNC.ID/MGI.IF MAPPINGS TO MOUSE EXPRESSION

Mouse.Expression.pheno.hpo<-left_join(Mouse.Expression.pheno.hpo,
                                            human.mouse.orthologs, by=c("MGI.ID"="MGI.ID"))

View(Mouse.Expression.pheno.hpo)


##COUNT THE NUMBER OF MOUSE GENES 
countMOUSE<-Mouse.Expression.pheno.hpo%>%
  select(MGI.ID, TS28.mutant.expression.detected)%>%
  filter(TS28.mutant.expression.detected=="Yes")%>%
  distinct(MGI.ID)

nrow(countMOUSE)


##############################################################################################################################
############################################################################################################################

## Combine the Human expressions and the Mouse orthologs using the superclass and MG.ID

## clean
Human_genes_expression_orthologs<-Human_genes_expression_orthologs%>%
  select_all()%>%
  drop_na()%>%
  distinct()

#View(Human_genes_expression_orthologs)
## merge
Human.Mouse.Expression.pheno.hpo<-merge(Mouse.Expression.pheno.hpo,
                                        Human_genes_expression_orthologs)

## clean
Human.Mouse.Expression.pheno.hpo<-Human.Mouse.Expression.pheno.hpo%>%
  select_all()%>%
  drop_na()%>%
  distinct()
View(Human.Mouse.Expression.pheno.hpo)






##############################################################################################################################
                  
#~~~~~~~~~~~~~~~~~~~~ VISUALISATION OF THE DIFFERENCE BETWEEN HUMAN AND MICE GENE EXPRESSION~~~~~~~~~~~~~~~~#####

############################################################################################################################

## Count the number of Human Genes with Mice Orthologs



## 0.1
View(Human_genes_expression_orthologs)

Human_genes_expression_0.1_<- Human_genes_expression_orthologs%>%
  select(HPO.superclass.description,HGNC.ID,GTEX.tissues, `TPM >/= 0.1 TPM` )%>%
  filter(`TPM >/= 0.1 TPM`=="Yes")%>%
  distinct()

## COUNT THE NUMBER OF GENES EXPRESSED
counthuman_0.1<-Human_genes_expression_0.1_%>%
  select(HGNC.ID)%>%
  distinct()

nrow(counthuman_0.1)

###################
## Count the number of Human genes with Similar HPO/MP terms


Human_genes_expression_0.1.top<- Human.Mouse.Expression.pheno.hpo%>%
  select(HPO.superclass.description,HGNC.ID,GTEX.tissues, `TPM >/= 0.1 TPM` )%>%
  distinct()

## COUNT THE NUMBER OF GENES EXPRESSED
counthuman_0.1_top<-Human_genes_expression_0.1.top%>%
  select(HGNC.ID)%>%
  distinct()

nrow(counthuman_0.1_top)


## 1


Human_genes_expression_1_<- Human_genes_expression_orthologs%>%
  select(HPO.superclass.description,HGNC.ID, GTEX.tissues,`TPM >/= 1 TPM` )%>%
  filter(`TPM >/= 1 TPM`=="Yes")%>%
  distinct()



## COUNT THE NUMBER OF GENES EXPRESSED
counthuman_1<-Human_genes_expression_1_%>%
  select(HGNC.ID)%>%
  distinct()

nrow(counthuman_1)



########



################################################################################################################

## merge
Human_genes_expression_ALL<- merge(Human_genes_expression_0.1_, Human_genes_expression_1_)

#########################################################################################################

## clean
Human_genes_expression_ALL<-Human_genes_expression_ALL%>%
  select_all()%>%
  drop_na()%>%
  distinct()

#View(Human_genes_expression_ALL)



#####################################################################################################################
################################################################################################################################

# merge mouse and human  expression with HPO/MP Terms again to ensure there are not wrongly mapped


Human.Mouse.Expression.pheno.hpo_<- left_join(Human.Mouse.Expression.pheno.hpo,hp.mp.toplevelmapping,
                                    by=c("mp.description"="mp.top.description", "HPO.superclass.description"="hp.top.term"))
View(Human.Mouse.Expression.pheno.hpo_)

Human.Mouse.Expression.pheno.hpo_CLEANED<-Human.Mouse.Expression.pheno.hpo_%>%
  select(MGI.ID, HGNC.ID, mp.description,GTEX.tissues,HPO.superclass.description,
         MP.Term, TS28.mutant.expression.detected,`TPM >/= 0.1 TPM`, `TPM >/= 1 TPM` )%>%
  drop_na()%>%
  distinct()


View(Human.Mouse.Expression.pheno.hpo_CLEANED)
######################################################################################################################

# ##!!!!!!!!UNHASH TO SAVE FILE !!!!!!!
# write.csv(Human.Mouse.Expression.pheno.hpo_CLEANED,
#           './Output_Files/Human_vs_Mouse/Human.Mouse.Expression.pheno.hpo_CLEANED.csv')



#### READ CSV
# 
# Human.Mouse.Expression.pheno.hpo_CLEANED<-read.csv('./Output_Files/Human_vs_Mouse/Human.Mouse.Expression.pheno.hpo_CLEANED.csv')
#  
# View(Human.Mouse.Expression.pheno.hpo_CLEANED)
#####################################################################################################################
################################################################################################################################

# Visualise Mice.Human.expression

Mice.Human.expression.COUNT_<-Human.Mouse.Expression.pheno.hpo_CLEANED%>%
  select(HGNC.ID ,GTEX.tissues,MP.Term, HPO.superclass.description, mp.description,
         MGI.ID,`TPM >/= 1 TPM`, `TPM >/= 0.1 TPM`, TS28.mutant.expression.detected)%>%
  distinct() %>% # check again for duplicated rows
  dplyr::mutate(HPO.MP = paste0(HPO.superclass.description,"+",mp.description))%>%
  distinct(HPO.MP, HGNC.ID ,GTEX.tissues,MP.Term, MGI.ID,`TPM >/= 1 TPM`, `TPM >/= 0.1 TPM` , TS28.mutant.expression.detected)


names(Mice.Human.expression.COUNT_)[names(Mice.Human.expression.COUNT_)== "TPM >/= 0.1 TPM" ] <- "Human Genes with Expression >/= 0.1 TPM"

names(Mice.Human.expression.COUNT_)[names(Mice.Human.expression.COUNT_)=="TPM >/= 1 TPM"] <- "Human Genes with Expression >/= 1 TPM"

names(Mice.Human.expression.COUNT_)[names(Mice.Human.expression.COUNT_)=="TS28.mutant.expression.detected"] <- "KO Gene Expression"


View(Mice.Human.expression.COUNT_)



Mice.Human.expression.CLEANED <- Mice.Human.expression.COUNT_ %>% 
  gather("Expression Type","Gene Expression",`Human Genes with Expression >/= 0.1 TPM`,
         `Human Genes with Expression >/= 1 TPM`, `KO Gene Expression` )%>%
  select(MGI.ID, HGNC.ID,`Expression Type`,`Gene Expression`, HPO.MP, GTEX.tissues,MP.Term)%>%
  distinct(MGI.ID, HGNC.ID,GTEX.tissues,MP.Term,`Expression Type`,`Gene Expression`, HPO.MP)

View(Mice.Human.expression.CLEANED)



#################################################################################################



## Filter gene expression FOR NERVOUS SYSTEM 


Mice.Human.expression.NS<-Mice.Human.expression.CLEANED%>%
  select_all()%>%
  filter(HPO.MP=="Abnormality of the nervous system+nervous system phenotype")%>%
  distinct(HPO.MP, HGNC.ID,MGI.ID,GTEX.tissues, MP.Term,`Gene Expression`, `Expression Type` )
View(Mice.Human.expression.NS)



## Change yes to 1 and no to 2 (2 will be changed to 0 later, dcast does not recognise 0 )

Mice.Human.expression.NS[Mice.Human.expression.NS[, "Gene Expression"] == "Yes", "Gene Expression"] <- "1"

Mice.Human.expression.NS[Mice.Human.expression.NS[, "Gene Expression"] == "No", "Gene Expression"] <- "2"

Mice.Human.expression.NS[Mice.Human.expression.NS[, "Gene Expression"] == "No", "Gene Expression"] <- "2"


View(Mice.Human.expression.NS)


library(reshape2)

data_wide <- dcast(Mice.Human.expression.NS,  GTEX.tissues+MP.Term +HPO.MP+MGI.ID+ HGNC.ID ~ `Expression Type`, value.var="Gene Expression")


View(data_wide)


### CHANGE 2s BACK TO 0S 

data_wide[data_wide[, "Human Genes with Expression >/= 0.1 TPM"] == "2", "Human Genes with Expression >/= 0.1 TPM"] <- "0"

data_wide[data_wide[, "Human Genes with Expression >/= 1 TPM"] == "2", "Human Genes with Expression >/= 1 TPM"] <- "0"

data_wide[data_wide[, "KO Gene Expression"] == "2", "KO Gene Expression"] <- "0"

##CONVERT THE NUMBERS TO NEMRICS

data_wide$`Human Genes with Expression >/= 0.1 TPM`<-as.numeric(data_wide$`Human Genes with Expression >/= 0.1 TPM`)

data_wide$`Human Genes with Expression >/= 1 TPM`<-as.numeric(data_wide$`Human Genes with Expression >/= 1 TPM`)

data_wide$`KO Gene Expression`<-as.numeric(data_wide$`KO Gene Expression`)


View(data_wide)

### CREATE UPSETR PLOT

library('UpSetR')

upset(data_wide, 
      sets = c("KO Gene Expression",
               "Human Genes with Expression >/= 0.1 TPM", 
               "Human Genes with Expression >/= 1 TPM"), 
      order.by="degree", matrix.color="red", point.size=5,
      mainbar.y.label = "Number of Expressions Detected in Tissues 
      Associated with Nervous System Phenotype", sets.x.label = "Number of Expressions Detected for each Condition", 
      sets.bar.color=c("maroon","blue","orange"))



#####################################################################################################################
################################################################################################################################
### Cardiovascular



Mice.Human.expression.CD<-Mice.Human.expression.CLEANED%>%
  select_all()%>%
  filter(HPO.MP=="Abnormality of the cardiovascular system+cardiovascular system phenotype")%>%
  distinct(HPO.MP, HGNC.ID,MGI.ID,GTEX.tissues, MP.Term,`Gene Expression`, `Expression Type` )
View(Mice.Human.expression.CD)



## Change yes to 1 and no to 2 (2 will be changed to 0 later, dcast does not recognise 0 )

Mice.Human.expression.CD[Mice.Human.expression.CD[, "Gene Expression"] == "Yes", "Gene Expression"] <- "1"

Mice.Human.expression.CD[Mice.Human.expression.CD[, "Gene Expression"] == "No", "Gene Expression"] <- "2"

Mice.Human.expression.CD[Mice.Human.expression.CD[, "Gene Expression"] == "No", "Gene Expression"] <- "2"


View(Mice.Human.expression.CD)


library(reshape2)

data_wide_CD <- dcast(Mice.Human.expression.CD,  GTEX.tissues+MP.Term +HPO.MP+MGI.ID+ HGNC.ID ~ `Expression Type`, value.var="Gene Expression")


View(data_wide_CD)


### CHANGE 2s BACK TO 0S 

data_wide_CD[data_wide_CD[, "Human Genes with Expression >/= 0.1 TPM"] == "2", "Human Genes with Expression >/= 0.1 TPM"] <- "0"

data_wide_CD[data_wide_CD[, "Human Genes with Expression >/= 1 TPM"] == "2", "Human Genes with Expression >/= 1 TPM"] <- "0"

data_wide_CD[data_wide_CD[, "KO Gene Expression"] == "2", "KO Gene Expression"] <- "0"

##CONVERT THE NUMBERS TO NEMRICS

data_wide_CD$`Human Genes with Expression >/= 0.1 TPM`<-as.numeric(data_wide_CD$`Human Genes with Expression >/= 0.1 TPM`)

data_wide_CD$`Human Genes with Expression >/= 1 TPM`<-as.numeric(data_wide_CD$`Human Genes with Expression >/= 1 TPM`)

data_wide_CD$`KO Gene Expression`<-as.numeric(data_wide_CD$`KO Gene Expression`)


View(data_wide_CD)

### CREATE UPSETR PLOT

library('UpSetR')

upset(data_wide_CD, 
      sets = c("KO Gene Expression",
               "Human Genes with Expression >/= 0.1 TPM", 
               "Human Genes with Expression >/= 1 TPM"), 
      order.by="degree", matrix.color="red", point.size=5,
      mainbar.y.label = "Number of Expressions Detected in Tissues 
      Associated with Cardiovascular System Phenotype", sets.x.label = "Number of Expressions Detected for each Condition", 
      sets.bar.color=c("blue","orange", "maroon"))



