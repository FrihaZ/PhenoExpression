
###############################################################################################################################
###############################################################################################################################

### Project: PhenoExpression ##################################################################################################
### Script: Human_vs_Mouse.R ###################################################################################################
### Purpose: To compare between the human and the mouse genes that may exhibit similar gene expressions ########################################################
### Author: Friha Zafar ####################################################################################################
### Date: 18/06/2019 ##########################################################################################################

###############################################################################################################################
###############################################################################################################################

## load packages ##############################################################################################################

library(dplyr);library(tidyr);library(stringr);library(readr);library(tidyverse);library(ggpubr)

################################################################################################################################
################################################################################################################################

# Load Files

human.mouse.orthologs<-read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/human.mouse.orthologs.txt",delim = "\t", skip = 1)

human_greater0_pheno<-read_csv("D:/MSC RESEARCH PROJECT/Output_Files/Human/Expression_greater_0/human_genes_TPM_greater0_pheno.csv")

human_0.1_pheno<-read_csv("D:/MSC RESEARCH PROJECT/Output_Files/Human/Expression_0.1/human_genes_TPM_0.1_pheno.csv")

human_1_pheno<-read_csv("D:/MSC RESEARCH PROJECT/Output_Files/Human/Expression_1/human_genes_TPM_1_pheno.csv")

mgi.genepheno.tissue.ts28_wt_mutant<-read_csv("D:/MSC RESEARCH PROJECT/Output_Files/Mice/mgi.genepheno.tissue.ts28_wt_mutant.csv")

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

human_greater0_pheno_cleaned <-clean_df(human_greater0_pheno)

names(human_greater0_pheno_cleaned)[names(human_greater0_pheno_cleaned)=="Expression"] <- "TPM > 0 TPM"

View(human_greater0_pheno_cleaned)


human_0.1_pheno_cleaned <-clean_df(human_0.1_pheno)

names(human_0.1_pheno_cleaned)[names(human_0.1_pheno_cleaned)=="Expression"] <- "TPM >/= 0.1 TPM"

View(human_0.1_pheno_cleaned)


human_1_pheno_cleaned <-clean_df(human_1_pheno)

names(human_1_pheno_cleaned)[names(human_1_pheno_cleaned)=="Expression"] <- "TPM >/= 1 TPM"


str(human_1_pheno_cleaned)


All_genes_expression<- right_join(human_greater0_pheno_cleaned,
                                  human_0.1_pheno_cleaned, by=c("HGNC.ID"="HGNC.ID", 
                                                                    "GTEX.tissues"="GTEX.tissues",
                                                                    "HPO.superclass.description"="HPO.superclass.description"))


All_genes_expression<- right_join(All_genes_expression,
                                  human_1_pheno_cleaned, by=c("HGNC.ID"="HGNC.ID", 
                                                                     "GTEX.tissues"="GTEX.tissues",
                                                                     "HPO.superclass.description"="HPO.superclass.description"))

View(All_genes_expression)


write.csv(All_genes_expression,"D:/MSC RESEARCH PROJECT/Output_Files/Human/All_genes_expression.csv") 


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

View(mgi.genepheno.tissue.ts28_wt_mutant)


##############################################################################################################################
############################################################################################################################

## Join mouse KO mouse expression with the phenotypes associated to the knockout gene 
## mgi.phenotypes: "MP.DESCRIPTION" ==  Mice.mutant.expression.clean: "MP.Term"

names(mgi.phenotypes)[names(mgi.phenotypes)=="MP.DESCRIPTION"] <- "MP.Term"
# 
Mice.mutant.expression.phenotypes<-merge(Mice.mutant.expression.clean,
                                              mgi.phenotypes)



View(Mice.mutant.expression.phenotypes)

Mice.mutant.expression.phenotypes<-Mice.mutant.expression.phenotypes%>%
  select_all()%>%
  drop_na()%>%
  distinct()

## change column name
names(Mice.mutant.expression.phenotypes)[names(Mice.mutant.expression.phenotypes)=="mp.ancestors"] <- "mp.top.term"

View(Mice.mutant.expression.phenotypes)

##############################################################################################################################
############################################################################################################################

## Join mouse KO mouse expression with the top level mapping of the hpo descriptions

hp.mp.toplevelmapping<-hp.mp.toplevelmapping%>%
  select(hp.top.term, hp.top.description, mp.top.term, mp.top.description )%>%
  distinct()

View(hp.mp.toplevelmapping)


Mouse.Expression.pheno.hpo<-merge(Mice.mutant.expression.phenotypes, hp.mp.toplevelmapping)


names(Mouse.Expression.pheno.hpo)[names(Mouse.Expression.pheno.hpo)=="hp.top.description"] <- "HPO.superclass.description"

View(Mouse.Expression.pheno.hpo)

##############################################################################################################################
############################################################################################################################

## Combine the Human expressions and the Mouse orthologs using the superclass and MG.ID

## clean
Human_genes_expression_orthologs<-Human_genes_expression_orthologs%>%
  select_all()%>%
  drop_na()%>%
  distinct()

view(Human_genes_expression_orthologs)
## merge
Human.Mouse.Expression.pheno.hpo<-merge(Mouse.Expression.pheno.hpo, Human_genes_expression_orthologs)

## clean
Human.Mouse.Expression.pheno.hpo<-Human.Mouse.Expression.pheno.hpo%>%
  select_all()%>%
  drop_na()%>%
  distinct()

View(Human.Mouse.Expression.pheno.hpo)



##############################################################################################################################
                  
#~~~~~~~~~~~~~~~~~~~~ VISUALISATION OF THE DIFFERENCE BETWEEN HUMAN AND MICE GENE EXPRESSION~~~~~~~~~~~~~~~~#####

############################################################################################################################

## Make the human gene expression into counts


Human_genes_expression_COUNT_0<-Human_genes_expression_orthologs%>%
  select_all()%>%  
  dplyr::group_by(HPO.superclass.description,`TPM > 0 TPM` ) %>%
  dplyr::summarise(n=n()) %>%
  distinct()


Human_genes_expression_COUNT_0.1<-Human_genes_expression_orthologs%>%
  select_all()%>%  
  dplyr::group_by(HPO.superclass.description,`TPM >/= 0.1 TPM` ) %>%
  dplyr::summarise(n=n()) %>%
  distinct()


Human_genes_expression_COUNT_1<-Human_genes_expression_orthologs%>%
  select_all()%>%  
  dplyr::group_by(HPO.superclass.description,`TPM >/= 1 TPM` ) %>%
  dplyr::summarise(n=n()) %>%
  distinct()


# change column names
names(Human_genes_expression_COUNT_0)[names(Human_genes_expression_COUNT_0)=="n"] <- "Number of genes at > 0 TPM"
names(Human_genes_expression_COUNT_0)[names(Human_genes_expression_COUNT_0)=="TPM > 0 TPM"] <- "Expression"

names(Human_genes_expression_COUNT_0.1)[names(Human_genes_expression_COUNT_0.1)=="n"] <- "Number of genes at > 0.1 TPM"

names(Human_genes_expression_COUNT_0.1)[names(Human_genes_expression_COUNT_0.1)=="TPM >/= 0.1 TPM"] <- "Expression"

names(Human_genes_expression_COUNT_1)[names(Human_genes_expression_COUNT_1)=="n"] <- "Number of genes at > 1 TPM"

names(Human_genes_expression_COUNT_1)[names(Human_genes_expression_COUNT_1)=="TPM >/= 1 TPM"] <- "Expression"

## merge
Human_genes_expression_COUNT<- merge(Human_genes_expression_COUNT_0,Human_genes_expression_COUNT_0.1)

Human_genes_expression_COUNT<- merge(Human_genes_expression_COUNT,Human_genes_expression_COUNT_1)

## clean
Human_genes_expression_COUNT<-Human_genes_expression_COUNT%>%
  select_all()%>%
  drop_na()%>%
  distinct()

View(Human_genes_expression_COUNT)

#####################################################################################################################
################################################################################################################################

# Make the mouse gene expression into counts

View(Mice.mutant.expression.phenotypes)
Mouse.Expression.COUNT<- Mice.mutant.expression.phenotypes%>%
  select( TS28.mutant.expression.detected, mp.description, MGI.ID )%>%
  dplyr::group_by(mp.description,TS28.mutant.expression.detected) %>%
  dplyr::summarise(KO.MICE.EXPRESSION=n()) %>%
  distinct(TS28.mutant.expression.detected, mp.description, KO.MICE.EXPRESSION)

#change column name
names(Mouse.Expression.COUNT)[names(Mouse.Expression.COUNT)=="TS28.mutant.expression.detected"] <- "Expression"


View(Mouse.Expression.COUNT)


Mouse.Expression.COUNT.HPO<- merge(Mouse.Expression.COUNT,hp.mp.toplevelmapping,
                                   by.x = "mp.description", by.y="mp.top.description")
names(Mouse.Expression.COUNT.HPO)[names(Mouse.Expression.COUNT.HPO)=="hp.top.description"] <- "HPO.superclass.description"
View(Mouse.Expression.COUNT.HPO)


#####################################################################################################################
################################################################################################################################

# merge human and mouse
Mice.Human.expression.COUNT<- merge(Human_genes_expression_COUNT,Mouse.Expression.COUNT.HPO)

Mice.Human.expression.COUNT<-Mice.Human.expression.COUNT%>%
  select_all()%>%
  distinct()

View(Mice.Human.expression.COUNT)


#####################################################################################################################
################################################################################################################################

# Visualise Mice.Human.expression.COUNT
Mice.Human.expression.COUNT_<-Mice.Human.expression.COUNT%>%
  select_all()%>% # select columns of interest
  distinct() %>% # check again for duplicated rows
  mutate(HPO.MP = paste0(HPO.superclass.description,"+",mp.description))%>%
  distinct(HPO.MP, Expression, `Number of genes at > 0 TPM`, `Number of genes at > 0.1 TPM`,
           `Number of genes at > 1 TPM`, KO.MICE.EXPRESSION)


names(Mice.Human.expression.COUNT_)[names(Mice.Human.expression.COUNT_)=="Number of genes at > 0 TPM"] <- "Human Genes with Expression > 0 TPM"

names(Mice.Human.expression.COUNT_)[names(Mice.Human.expression.COUNT_)=="Number of genes at > 0.1 TPM"] <- "Human Genes with Expression >/= 0.1 TPM"

names(Mice.Human.expression.COUNT_)[names(Mice.Human.expression.COUNT_)=="Number of genes at > 1 TPM"] <- "Human Genes with Expression >/= 1 TPM"

names(Mice.Human.expression.COUNT_)[names(Mice.Human.expression.COUNT_)=="KO.MICE.EXPRESSION"] <- "Knockout Mice"


Mice.Human.expression.COUNT_ <- Mice.Human.expression.COUNT_ %>% 
  gather("Expression Type","Number of Genes",`Human Genes with Expression > 0 TPM`, `Human Genes with Expression >/= 0.1 TPM`,
        `Human Genes with Expression >/= 1 TPM`, `Knockout Mice`)%>%
  select(`Expression Type`,`Number of Genes`, HPO.MP, Expression)%>%
  distinct(`Expression Type`,`Number of Genes`, HPO.MP, Expression)


view(Mice.Human.expression.COUNT_)

#####################################################################################################################
################################################################################################################################

## PLOT

mouse.human.p1<-ggplot(Mice.Human.expression.COUNT_, aes(fill=Mice.Human.expression.COUNT_$Expression, 
                                                        y=Mice.Human.expression.COUNT_$`Number of Genes`,
                                                        x=Mice.Human.expression.COUNT_$HPO.MP))%>% 
  
  + labs(x= "HPO Toplevel (Human) + MP Toplevel (Mouse)", y="Number of Genes", 
         subtitle= "Number of Genes Expressed for each Human Phenotype Ontology (HPO) and Mammalian Phenotype (MP) Annotation")%>%
  + guides(fill=guide_legend(title="Gene Expression: ")) %>%
  +theme(axis.title=element_text(size=14,face="bold"),
        plot.caption=element_text(face = "italic", size=12, hjust = 0),
        text = element_text(size=10),
        legend.position = "bottom",legend.direction = "horizontal",
        legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
        axis.text.y= element_text(size=8),
        axis.text.x= element_text(size=10),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +geom_bar( stat="identity")%>%   
  +facet_wrap(~Mice.Human.expression.COUNT_$`Expression Type` , scales = "free_x")%>%
  + coord_flip()
mouse.human.p1



# save plot
library(cowplot)
save_plot("./Plots/Human+Mouse/mouse_human_p1.png",
         mouse.human.p1 ,base_height= 6 ,base_aspect_ratio = 3) 

 
 
 
 ## PLOT

mouse.human.p<-ggplot(Mice.Human.expression.COUNT_, aes(fill=Mice.Human.expression.COUNT_$`Expression Type`, 
                                                        y=Mice.Human.expression.COUNT_$`Number of Genes`,
                                                        x=Mice.Human.expression.COUNT_$HPO.MP))%>% 
 
 + labs(x= "HPO Toplevel (Human) + MP Toplevel (Mouse)", y="Number of Genes", 
        subtitle= "Number of Genes Expressed for each Human Phenotype Ontology (HPO) and Mammalian Phenotype (MP) Annotation")%>%
 + guides(fill=guide_legend(title="Gene Expression: ")) %>%
 +theme(axis.title=element_text(size=14,face="bold"),
        plot.caption=element_text(face = "italic", size=12, hjust = 0),
        text = element_text(size=10),
        legend.position = "bottom",legend.direction = "horizontal",
        legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
        axis.text.y= element_text(size=8),
        axis.text.x= element_text(size=10),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
 +geom_bar( stat="identity")%>%   
 +facet_wrap(~Mice.Human.expression.COUNT_$`Expression Type` , scales = "free_x")%>%
 + coord_flip()
mouse.human.p



# save plot
library(cowplot)
save_plot("./Plots/Human+Mouse/mouse_human_p.png",
         mouse.human.p ,base_height= 6 ,base_aspect_ratio = 3) 


 
#####################################################################################################################
################################################################################################################################
 
 
