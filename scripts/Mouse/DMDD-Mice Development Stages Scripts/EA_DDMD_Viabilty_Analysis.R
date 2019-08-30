
###############################################################################################################################
###############################################################################################################################

### Project: PhenoExpression ##################################################################################################
### Script: EA_DDMD_Viabilty_Analysis.R ###################################################################################################
### Purpose: To compare gene expression and viability of mice genes in different developmental stages######################################################
### Author: Friha Zafar ####################################################################################################
### Date: 17/07/2019 ##########################################################################################################

###############################################################################################################################
###############################################################################################################################

## load packages ##############################################################################################################

library(dplyr);library(tidyr);library(stringr);library(readr);library(tidyverse);library(ggpubr);library(ggplot2)

################################################################################################################################
################################################################################################################################

# Load Files

DMDD.Mouse.gene.expr<-read.csv("D:/MSC RESEARCH PROJECT/Output_Files/Mice/EA_DMDDs DATASET/sort_MEAN.Pheno.Mouse.Development.csv")

impc_primary_viability<-read.delim("D:/MSC RESEARCH PROJECT/GITHUB/data/mouse_viability_data_files/impc_primary_viability.txt")

impc_secondary_viability<-read.delim("D:/MSC RESEARCH PROJECT/GITHUB/data/mouse_viability_data_files/impc_secondary_viability.txt")

##################################################################################################################################################
################################################################################################################################################


######### ~~~~~~~~~~~Dividing the DMDD Gene Expression Dataset into three catagories:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#######
## TPM >/= 0.1 TPM (YES OR NO)  !!!!NOTE: THERE WERE NO GENES THAT HAD AN EXPRESSION <0.1 TPM!!! 
## TPM >/= 1TPM (YES OR NO)

## REMOVED X COLUMN

DMDD.Mouse.gene.expr<-DMDD.Mouse.gene.expr%>%
  select(MGI.Accession.ID, mp.description, 
         Expression,Development.Stage,mp.ancestors)%>%
  distinct()
View(DMDD.Mouse.gene.expr)


# TPM >/= 0.1 (yes, no)

tpm0.1.DMDD.Mouse.gene.expr<- DMDD.Mouse.gene.expr

tpm0.1.DMDD.Mouse.gene.expr$Expression[tpm0.1.DMDD.Mouse.gene.expr$Expression  >= 0.1 ] <- "Yes"
tpm0.1.DMDD.Mouse.gene.expr$Expression[tpm0.1.DMDD.Mouse.gene.expr$Expression  < 0.1 ] <- "No"

#write.csv(human_genes_TPM_0.1,'./Output_Files/Human/Expression_0.1/human_genes_TPM_0.1.csv')

names(tpm0.1.DMDD.Mouse.gene.expr)[names(tpm0.1.DMDD.Mouse.gene.expr)=="Expression"] <- "MOUSE TPM >/= 0.1 TPM"


View(tpm0.1.DMDD.Mouse.gene.expr)

# TPM>/= 1 (yes, no)

tpm1.DMDD.Mouse.gene.expr<- DMDD.Mouse.gene.expr

tpm1.DMDD.Mouse.gene.expr$Expression[tpm1.DMDD.Mouse.gene.expr$Expression >= 1 ] <- "Yes"

tpm1.DMDD.Mouse.gene.expr$Expression[tpm1.DMDD.Mouse.gene.expr$Expression < 1 ] <- "No"

names(tpm1.DMDD.Mouse.gene.expr)[names(tpm1.DMDD.Mouse.gene.expr)=="Expression"] <- "MOUSE TPM >/= 1 TPM"

#write.csv(human_genes_TPM_1,'./Output_Files/Human/Expression_1/human_genes_TPM_1.csv')

View(tpm1.DMDD.Mouse.gene.expr)


###################################################################################################

## MERGE THE THREHSOLDS INTO ONE DF


all_mouse_dmdd<-merge(tpm0.1.DMDD.Mouse.gene.expr, tpm1.DMDD.Mouse.gene.expr)

View(all_mouse_dmdd)

##TOTAL NUMBER OF GENES WITH A MP TERM LINKED

all_mouse_dmdd_total<-all_mouse_dmdd%>%
  select(mp.description, MGI.Accession.ID)%>%
  distinct(MGI.Accession.ID)

nrow(all_mouse_dmdd_total)
  
### THE TOTAL NUMBER OF GENES FOR EACH MP TERM

all_mouse_dmdd_num.genes_<-all_mouse_dmdd%>%
  select(mp.description,MGI.Accession.ID)%>%
  distinct()%>%
  group_by(mp.description )%>%
  summarise("Number of genes"=n()) %>%
  distinct()




all_mouse_dmdd_num.genes_$Percent<-(all_mouse_dmdd_num.genes_$`Number of genes`/8799)*100 
#View(all_mouse_dmdd_num.genes_)

################################################################################################################################
################################################################################################################################

## PLOT 

library(randomcoloR)
n <-70
palette <- distinctColorPalette(n)

all_mouse_dmdd_PLOT<-ggplot(all_mouse_dmdd_num.genes_, 
                            aes(x =reorder(all_mouse_dmdd_num.genes_$mp.description,
                                           -all_mouse_dmdd_num.genes_$`Number of genes`),
                                                y =all_mouse_dmdd_num.genes_$`Number of genes`,
                                                  fill=all_mouse_dmdd_num.genes_$`Number of genes`))%>%
  
  +labs(x= "Top-Level MP-Term",y="Number of Genes", 
        title= "The Number of Mice Genes for Each Top-Level MP-Term")%>%
  #+ scale_y_continuous(breaks = seq(0, 200, by = 10))%>%
  +scale_fill_gradient(low = "green", high = "red")%>%
  + theme(legend.title = element_blank(),
          legend.position = "none",
          axis.title=element_text(size=14,face="bold"),
          plot.caption=element_text(face = "italic", size=12, hjust = 0),
          text = element_text(size=10),
          axis.text.x= element_text(size=10),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +geom_bar(stat = "identity")%>%    # to create a stacked barchart
  +geom_text(show.legend = FALSE,aes(label=all_mouse_dmdd_num.genes_$`Number of genes`),size = 5, 
             position = position_stack(vjust = 0.5), colour=c("black"), fontface='bold') %>%
  
  + coord_flip()


# all_mouse_dmdd_PLOT

# !!!!UNHASH TO SAVE PLOT !!!!!!!!!

# library(cowplot)
# # 
# save_plot(paste("./Plots/Mouse/EA_Developmental_Stages/all_mouse_dmdd_PLOT.png"),
#           all_mouse_dmdd_PLOT ,base_height= 10 ,base_aspect_ratio = 2) 
# 
##################

#PERCENTAGE


all_mouse_dmdd_PLOT_percent<-ggplot(all_mouse_dmdd_num.genes_, 
                                    aes(x =reorder(all_mouse_dmdd_num.genes_$mp.description,
                                                   -all_mouse_dmdd_num.genes_$Percent),
                                        y =all_mouse_dmdd_num.genes_$Percent,
                                        fill=all_mouse_dmdd_num.genes_$Percent))%>%
  
  +labs(x= "Top-Level MP-Term",y="Percentage (%) of Genes with MP Annotation", 
        title= "Percentage (%) of Genes Associated to a Top-Level MP Term")%>%
  + scale_y_continuous(breaks = seq(0, 100, by = 10))%>%
  +scale_fill_gradient(low = "green", high = "red")%>%
  + theme(legend.title = element_blank(),
          legend.position = "none",
          axis.title=element_text(size=14,face="bold"),
          plot.caption=element_text(face = "italic", size=14, hjust = 0),
          text = element_text(size=14),
          axis.text.x= element_text(size=14),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +geom_bar(stat = "identity")%>%    # to create a stacked barchart
  +geom_text(show.legend = FALSE,aes(label=paste0(round(all_mouse_dmdd_num.genes_$Percent,0),"%")),size = 5, 
             position = position_stack(vjust = 0.5), colour=c("black"), fontface='bold') %>%
  
  + coord_flip()

# all_mouse_dmdd_PLOT_percent

#!!!!!!UNHASH TO SAVE PLOT!!!!
 
# library(cowplot)
# # 
# save_plot(paste("./Plots/Mouse/EA_Developmental_Stages/all_mouse_dmdd_PLOT_percent.png"),
#           all_mouse_dmdd_PLOT_percent ,base_height= 10 ,base_aspect_ratio = 2) 
# 



##################################################################################################
#################################################################################################################

########~~~~~~~~~~~~~~~~~~~~JOIN DATASET TO VIABILTY DATA~~~~~~~~~~~~~~~~~~~~~~~##################


View(all_mouse_dmdd)
str(impc_primary_viability$MGI_ID)

## TO ADD THE TPM VALUES TO THE DATAFRAME
all_mouse_dmdd_TPM<-right_join(DMDD.Mouse.gene.expr, all_mouse_dmdd)

View(all_mouse_dmdd_TPM)


## TO JOIN THE EXPRESSION DATA TO VIABILITY DATA
primary_viability_mouse_dmdd<-left_join(all_mouse_dmdd_TPM, 
                                        impc_primary_viability,
                                        by=c("MGI.Accession.ID"= "MGI_ID"))


#
View(primary_viability_mouse_dmdd)



## REMOVE GENES THAT HAD NO VIABILITY DATA

primary_viability_mouse_dmdd<-primary_viability_mouse_dmdd%>%
  select_all()%>%
  drop_na()%>%
  distinct()

View(primary_viability_mouse_dmdd)


Count_Primary<-primary_viability_mouse_dmdd%>%
  select(IMPC_Viability, MGI.Accession.ID )%>%
  distinct()

View(Count_Primary)

##UNHASH TO SAVE DF

#write.csv(primary_viability_mouse_dmdd,
#"D:/MSC RESEARCH PROJECT/Output_Files/Mice/EA_DMDDs DATASET/Viability/primary_viability_mouse_dmdd.csv")


########################################################################

## COLOURS
palette2 <- c("purple", "yellow", "green")

## SORT THE LEVELS OF THE FACTORS

primary_viability_mouse_dmdd$IMPC_Viability<-factor(primary_viability_mouse_dmdd$IMPC_Viability,levels=c("Lethal","Subviable","Viable"))


primary_viability_mouse_dmdd$Development.Stage <- factor(primary_viability_mouse_dmdd$Development.Stage,
                                                                 levels=c("4-somite stage", "5-somite stage", "6-somite stage",
                                                                                   "7-somite stage", "8-somite stage", "9-somite stage",
                                                                                   "10-somite stage", "11-somite stage","12-somite stage",
                                                                                   "13-somite stage", "14-somite stage","15-somite stage",
                                                                                   "16-somite stage", "17-somite stage","18-somite stage",
                                                                                   "19-somite stage", "20-somite stage","21-somite stage",
                                                                                   "22-somite stage", "23-somite stage","24-somite stage",
                                                                                   "25-somite stage", "26-somite stage","27-somite stage",
                                                                                   "28-somite stage", "29-somite stage","30-somite stage",
                                                                                   "31-somite stage","32-somite stage", "33-somite stage",
                                                                                   "34-somite stage","35-somite stage", "36-somite stage"))
# 
View(primary_viability_mouse_dmdd)



primary_viability_mouse_dmdd_PLOT1<-ggplot(primary_viability_mouse_dmdd,
                                           aes(x =primary_viability_mouse_dmdd$Development.Stage,
                                               y =primary_viability_mouse_dmdd$Expression,
                                               group=primary_viability_mouse_dmdd$MGI.Accession.ID ,
                                               colour=primary_viability_mouse_dmdd$IMPC_Viability))%>%
  +labs(x= "Developmental Stage", y="Gene Expression (TPM)", 
        title= "Mice Gene Expression (TPM) in Different Developmental Stages and their Associated Viabilities")%>%
  +scale_y_continuous(breaks = seq(0, 15000, by = 1000))%>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10, angle = 45, hjust = 1, vjust = 1)) %>%
  +geom_line()%>%
  +geom_point(size = 4)%>%
  + scale_color_manual(values =palette2, aesthetics = "colour")%>%
  +guides(color=guide_legend(title="Viabilty: "))%>%
  #+ xlim(rev(levels(primary_viability_mouse_dmdd$Development.Stage)))%>%
#+coord_flip()
  +facet_wrap( . ~ primary_viability_mouse_dmdd$IMPC_Viability)



#primary_viability_mouse_dmdd_PLOT1

# 
# # 
# library(cowplot)
# 
# save_plot(paste("./Plots/Mouse/EA_Developmental_Stages/Viability/primary_viability_mouse_dmdd_PLOT1_DevelopmentStage.png"),
#            primary_viability_mouse_dmdd_PLOT1 ,base_height= 10 ,base_aspect_ratio = 2) 


#################################################################################
## remove outliers

outlier<-boxplot.stats(primary_viability_mouse_dmdd$Expression)$out


primary_viability_mouse_dmdd_CLEANED <- primary_viability_mouse_dmdd[!primary_viability_mouse_dmdd$Expression %in% outlier, ]

View(primary_viability_mouse_dmdd_CLEANED)

# #UNHASH TO SAVE DF  

# write.csv(primary_viability_mouse_dmdd_CLEANED,
#           "D:/MSC RESEARCH PROJECT/Output_Files/Mice/EA_DMDDs DATASET/Viability/primary_viability_mouse_dmdd_CLEANED.csv")
# ###

######################################################################################
#####################################################################################################




primary_via_mouse_PLOT3_CLEANED<-ggplot(primary_viability_mouse_dmdd_CLEANED,
                                           aes(x =primary_viability_mouse_dmdd_CLEANED$Development.Stage,
                                               y =primary_viability_mouse_dmdd_CLEANED$Expression,
                                               group=primary_viability_mouse_dmdd_CLEANED$MGI.Accession.ID ,
                                               colour=primary_viability_mouse_dmdd_CLEANED$IMPC_Viability))%>%
  +labs(x= "Developmental Stage", y="Gene Expression (TPM)", 
        title= "Mice Gene Expression (TPM) in Different Developmental Stages and their Associated Viabilities")%>%
  +scale_y_continuous(breaks = seq(0, 15000, by = 10))%>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10, angle = 45, hjust = 1, vjust = 1)) %>%
  +geom_line()%>%
  +geom_point(size = 4)%>%
  + scale_color_manual(values =palette2, aesthetics = "colour")%>%
  +guides(color=guide_legend(title="Viabilty: "))%>%
  +facet_wrap( . ~ primary_viability_mouse_dmdd_CLEANED$IMPC_Viability)



#primary_viability_mouse_dmdd_PLOT1

#UNHASH TO SAVE PLOT  
#  
# library(cowplot)
#  
# save_plot(paste("./Plots/Mouse/EA_Developmental_Stages/Viability/primary_viability_mouse_dmdd_PLOT3_Develop_stage_cleaned.png"),
#           primary_via_mouse_PLOT3_CLEANED ,base_height= 10 ,base_aspect_ratio = 2) 


 



###############################################################################################################
################################################################################################################

## Calculating the mean, SE and SE  for each Lethal, Subviable and Viable



mean_primary_viab_mouse_CLEANED<-primary_viability_mouse_dmdd_CLEANED%>%
  select(IMPC_Viability, Expression, Development.Stage)%>%
  group_by(IMPC_Viability, Development.Stage)%>%
  mutate(Mean= mean(Expression), SD=sd(Expression),
         SE= sd(Expression)/sqrt(length(Expression)))%>%    
  distinct()

View(mean_primary_viab_mouse_CLEANED)


#########################################################################################################
###################################################################################################################


## PLOT MEAN GENE EXPRESSION +SD



library(reshape)
MEAN_primary_via_mouse_PLOT4_CLEANED<-ggplot(mean_primary_viab_mouse_CLEANED,
                                        aes(x =mean_primary_viab_mouse_CLEANED$Development.Stage,
                                            y =mean_primary_viab_mouse_CLEANED$Mean ,
                                            group=mean_primary_viab_mouse_CLEANED$IMPC_Viability ,
                                            colour=mean_primary_viab_mouse_CLEANED$IMPC_Viability))%>%
  +labs(x= "Developmental Stage", y="Mean Gene Expression (TPM)", 
        title= "Mean Mice Gene Expression (TPM) in Different Developmental Stages and their Associated Viabilities")%>%
  +scale_y_continuous(breaks = seq(0, 15000, by = 10))%>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10, angle = 45, hjust = 1, vjust = 1)) %>%
  +geom_line()%>%
  +geom_point(size = 4)%>%
  + scale_color_manual(values =palette2, aesthetics = "colour")%>%
  +geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=1.3)%>%
  +guides(color=guide_legend(title="Viabilty: "))%>%
  #+ xlim(rev(levels(primary_viability_mouse_dmdd$Development.Stage)))%>%
  #+coord_flip()
  +facet_wrap( . ~ mean_primary_viab_mouse_CLEANED$IMPC_Viability)



#primary_viability_mouse_dmdd_PLOT1



#### UNHASH TO SAVE PLOT!!!
# library(cowplot)
 
# save_plot(paste("./Plots/Mouse/EA_Developmental_Stages/Viability/MEAN_primary_viability_mouse_dmdd_PLOT4_Develop_stage_cleaned.png"),
#          MEAN_primary_via_mouse_PLOT4_CLEANED ,base_height= 10 ,base_aspect_ratio = 2) 
 
##########################################################################################################################


## PLOT MEAN GENE EXPRESSION +SE




MEAN_primary_via_mouse_PLOT4_CLEANED_SE<-ggplot(mean_primary_viab_mouse_CLEANED,
                                             aes(x =mean_primary_viab_mouse_CLEANED$Development.Stage,
                                                 y =mean_primary_viab_mouse_CLEANED$Mean ,
                                                 group=mean_primary_viab_mouse_CLEANED$IMPC_Viability ,
                                                 colour=mean_primary_viab_mouse_CLEANED$IMPC_Viability))%>%
  +labs(x= "Developmental Stage", y="Mean Gene Expression (TPM)", 
        title= "Mean Mice Gene Expression (TPM) in Different Developmental Stages and their Associated Viabilities")%>%
  +scale_y_continuous(breaks = seq(0, 50, by = 5))%>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10, angle = 45, hjust = 1, vjust = 1)) %>%
  +geom_line()%>%
  +geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=1, color="black")%>%
  +geom_point(size = 4)%>%
  + scale_color_manual(values =palette2, aesthetics = "colour")%>%
  +guides(color=guide_legend(title="Viabilty: "))%>%
  #+ xlim(rev(levels(primary_viability_mouse_dmdd$Development.Stage)))%>%
  #+coord_flip()
  +facet_wrap( . ~ mean_primary_viab_mouse_CLEANED$IMPC_Viability)




# 
# #### UNHASH TO SAVE PLOT!!!
# library(cowplot)
# 
# save_plot(paste("./Plots/Mouse/EA_Developmental_Stages/Viability/MEAN_primary_via_mouse_PLOT4_CLEANED_SE.png"),
#          MEAN_primary_via_mouse_PLOT4_CLEANED_SE ,base_height= 10 ,base_aspect_ratio = 2) 
# 

####################################################################################################
#NO FACETWRAP

MEAN_primary_via_mouse_PLOT5_CLEANED<-ggplot(mean_primary_viab_mouse_CLEANED,
                                             aes(x =mean_primary_viab_mouse_CLEANED$Development.Stage,
                                                 y =mean_primary_viab_mouse_CLEANED$Mean ,
                                                 group=mean_primary_viab_mouse_CLEANED$IMPC_Viability ,
                                                 colour=mean_primary_viab_mouse_CLEANED$IMPC_Viability))%>%
  +labs(x= "Developmental Stage", y="Mean Gene Expression (TPM)", 
        title= "Mean Mice Gene Expression (TPM) in Different Developmental Stages and their Associated Viabilities")%>%
  +scale_y_continuous(breaks = seq(0, 100, by = 5))%>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10, angle = 45, hjust = 1, vjust = 1)) %>%
  +geom_line()%>%
  +geom_point(size = 4)%>%
  + scale_color_manual(values =palette2, aesthetics = "colour")%>%
  +guides(color=guide_legend(title="Viabilty: "))



#primary_viability_mouse_dmdd_PLOT1

# 
# # 
# library(cowplot)
# # 
# save_plot(paste("./Plots/Mouse/EA_Developmental_Stages/Viability/MEAN_primary_viability_mouse_dmdd_PLOT5_Develop_stage_cleaned.png"),
#          MEAN_primary_via_mouse_PLOT5_CLEANED ,base_height= 10 ,base_aspect_ratio = 2) 
# # 





################################################################################################

## EACH MP.TERM





primary_viability_mouse_dmdd_PLOT2<-ggplot(primary_viability_mouse_dmdd,
                                           aes(x =primary_viability_mouse_dmdd$mp.description,
                                               y =primary_viability_mouse_dmdd$Expression,
                                               group=primary_viability_mouse_dmdd$MGI.Accession.ID ,
                                               colour=primary_viability_mouse_dmdd$IMPC_Viability))%>%
  +labs(x= "Top-Level MP-Term", y="Gene Expression (TPM)", 
        title= "Mice Gene Expression (TPM), their Linked Mammalian Phenotype (MP) and Associated Viabilities")%>%
  +scale_y_continuous(breaks = seq(0, 15000, by = 1000))%>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10, angle = 45, hjust = 1, vjust = 1)) %>%
  +geom_line()%>%
  +geom_point(size = 4)%>%
  + scale_color_manual(values =palette2, aesthetics = "colour")%>%
  +guides(color=guide_legend(title="Viabilty: "))%>%
  #+ xlim(rev(levels(primary_viability_mouse_dmdd$Development.Stage)))%>%
  #+coord_flip()
  +facet_wrap( . ~ primary_viability_mouse_dmdd$IMPC_Viability)



# primary_viability_mouse_dmdd_PLOT2

#UNHASH TO SAVE PLOT  

# # 
# library(cowplot)
# # 
# save_plot(paste("./Plots/Mouse/EA_Developmental_Stages/Viability/primary_viability_mouse_dmdd_MP_TERM.png"),
#            primary_viability_mouse_dmdd_PLOT2 ,base_height= 10 ,base_aspect_ratio = 2) 
# # 
# 
# 

















##############################################################################################################
######################################################################################################
## To merge the dataset with the secondary viability data



second_viability_mouse_dmdd<-left_join(primary_viability_mouse_dmdd, 
                                          impc_secondary_viability,
                                        by=c("MGI.Accession.ID"= "MGI_ID"))


#

second_viability_mouse_dmdd<-second_viability_mouse_dmdd%>%
  select_all()%>%
  select(-IMPC_Viability)%>%
  drop_na()%>%
  distinct()

View(second_viability_mouse_dmdd)


second_viability_mouse_dmdd$Window_Category <-factor(second_viability_mouse_dmdd$Window_Category,
                                                     levels=c("Early-gestation",
                                                              "Mid-gestation",
                                                              "Late-gestation"))


## save df

#write.csv(second_viability_mouse_dmdd,
#          "D:/MSC RESEARCH PROJECT/Output_Files/Mice/EA_DMDDs DATASET/Viability/second_viability_mouse_dmdd.csv")



## MEAN Expression

SE_mean_second_viability_mouse_dmdd<-second_viability_mouse_dmdd%>%
  select(Window_Category, Expression, Development.Stage)%>%
  group_by(Window_Category, Development.Stage)%>%
  mutate(Mean= mean(Expression), SD=sd(Expression),
         SE= sd(Expression)/sqrt(length(Expression)))%>% 
  distinct()


View(SE_mean_second_viability_mouse_dmdd)


# write.csv(SE_mean_second_viability_mouse_dmdd,
#           "D:/MSC RESEARCH PROJECT/Output_Files/Mice/EA_DMDDs DATASET/Viability/SE_mean_second_viability_mouse_dmdd.csv")




## COUNT THE NUMBER OF GENES FOR EACH CATEGORY
COUNT_SECONDARY<-second_viability_mouse_dmdd%>%
  select(Window_Category, MGI.Accession.ID)%>%
  distinct()

View(COUNT_SECONDARY)


###################################################################################################
##PLOT

######################################################################################################################

## Mean without SD 
palette3 <- c("red","orange","blue")

MEAN_second_viability_mouse_dmdd_plot<-ggplot(SE_mean_second_viability_mouse_dmdd,
                                               aes(x =SE_mean_second_viability_mouse_dmdd$Development.Stage,
                                                   y =SE_mean_second_viability_mouse_dmdd$Mean ,
                                                   group=SE_mean_second_viability_mouse_dmdd$Window_Category   ,
                                                   colour=SE_mean_second_viability_mouse_dmdd$Window_Category))%>%
  +labs(x= "Developmental Stage", y="Mean Gene Expression (TPM)", 
        title= "Mean Mice Gene Expression (TPM) in Different Developmental Stages and their Window of Lethality")%>%
  +scale_y_continuous(breaks = seq(-200, 400, by = 10))%>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10, angle = 45, hjust = 1, vjust = 1)) %>%
  +geom_line()%>%
  +geom_point(size = 4)%>%
  #+geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=1.3)%>%
  + scale_color_manual(values =palette3, aesthetics = "colour")%>%
  +guides(color=guide_legend(title="Window of Lethality: "))
  #+compare_means( Expression ~ Window_Category, data = SE_mean_second_viability_mouse_dmdd)
#+facet_wrap( . ~ SE_mean_second_viability_mouse_dmdd$Window_Category)




#primary_viability_mouse_dmdd_PLOT

 
# !!!!!!UNHASH TO SAVE PLOT  !!!!!!!!111

# library(cowplot)
 
# save_plot(paste("./Plots/Mouse/EA_Developmental_Stages/Viability/MEAN_second_viability_mouse_dmdd_plot.png"),
#             MEAN_second_viability_mouse_dmdd_plot ,base_height= 10 ,base_aspect_ratio = 2) 

################################################################################
## Mean +SD 
palette3 <- c("red","orange","blue")

MEAN_second_viability_mouse_dmdd_plot1<-ggplot(SE_mean_second_viability_mouse_dmdd,
                                             aes(x =SE_mean_second_viability_mouse_dmdd$Development.Stage,
                                                 y =SE_mean_second_viability_mouse_dmdd$Mean ,
                                                 group=SE_mean_second_viability_mouse_dmdd$Window_Category   ,
                                                 colour=SE_mean_second_viability_mouse_dmdd$Window_Category))%>%
  +labs(x= "Developmental Stage", y="Mean Gene Expression (TPM)", 
        title= "Mean Mice Gene Expression (TPM) in Different Developmental Stages and their Window of Lethality")%>%
  +scale_y_continuous(breaks = seq(-200, 400, by = 10))%>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10, angle = 45, hjust = 1, vjust = 1)) %>%
  +geom_line()%>%
  +geom_point(size = 4)%>%
  +geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=1.3)%>%
  + scale_color_manual(values =palette3, aesthetics = "colour")%>%
  +guides(color=guide_legend(title="Window of Lethality: "))%>%
  +compare_means( Expression ~ Window_Category, data = SE_mean_second_viability_mouse_dmdd)
  +facet_wrap( . ~ SE_mean_second_viability_mouse_dmdd$Window_Category)




#primary_viability_mouse_dmdd_PLOT1

# 
# # 
# library(cowplot)
# # 
# save_plot(paste("./Plots/Mouse/EA_Developmental_Stages/Viability/MEAN_second_viability_mouse_dmdd_plot1.png"),
#            MEAN_second_viability_mouse_dmdd_plot1 ,base_height= 10 ,base_aspect_ratio = 2) 

#########################################################################################################



## Mean +SE
palette3 <- c("red","orange","blue")

MEAN_second_viability_mouse_dmdd_plot1_SE<-ggplot(SE_mean_second_viability_mouse_dmdd,
                                               aes(x =SE_mean_second_viability_mouse_dmdd$Development.Stage,
                                                   y =SE_mean_second_viability_mouse_dmdd$Mean ,
                                                   group=SE_mean_second_viability_mouse_dmdd$Window_Category   ,
                                                   colour=SE_mean_second_viability_mouse_dmdd$Window_Category))%>%
  +labs(x= "Developmental Stage", y="Mean Gene Expression (TPM)", 
        title= "Mean Mice Gene Expression (TPM) in Different Developmental Stages and their Window of Lethality")%>%
  +scale_y_continuous(breaks = seq(-200, 400, by = 10))%>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10, angle = 45, hjust = 1, vjust = 1)) %>%
  +geom_line()%>%
  +geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=1, color="black")%>%
  +geom_point(size = 4)%>%
  + scale_color_manual(values =palette3, aesthetics = "colour")%>%
  +guides(color=guide_legend(title="Window of Lethality: "))%>%
  +facet_wrap( . ~ SE_mean_second_viability_mouse_dmdd$Window_Category)




#primary_viability_mouse_dmdd_PLOT1

# 
# # 
# library(cowplot)
# # 
# save_plot(paste("./Plots/Mouse/EA_Developmental_Stages/Viability/MEAN_second_viability_mouse_dmdd_plot1_SE.png"),
#           MEAN_second_viability_mouse_dmdd_plot1_SE ,base_height= 10 ,base_aspect_ratio = 2) 

############################################################################################################################
##########################################################################################################################

