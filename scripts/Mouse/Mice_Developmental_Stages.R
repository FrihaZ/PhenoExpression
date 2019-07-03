
###############################################################################################################################
###############################################################################################################################

### Project: PhenoExpression ##################################################################################################
### Script: Mice_Developmental_Stages.R ###################################################################################################
### Purpose: To compare gene expression of KO and WT mice in different developmental stages ########################################################
### Author: Friha Zafar ####################################################################################################
### Date: 26/06/2019 ##########################################################################################################

###############################################################################################################################
###############################################################################################################################

## load packages ##############################################################################################################

library(dplyr);library(tidyr);library(stringr);library(readr);library(tidyverse);library(ggpubr)

################################################################################################################################
################################################################################################################################

# Load Files
MGI_ALL_TS_knockout<-read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/MGIgeneExpressionQuery_20190626_131124_knockout.txt",delim = "\t")





################################################################################################################################
################################################################################################################################

###############~~~~~~~~~~~~~~~~~~~~~~~~CLEAN KO DATA AND MAP THEM TO MAMMALIAN PHENOTYPES~~~~~~~~~~~~#########

################################################################################################################################
################################################################################################################################




MGI_ALL_TS_knockout<-MGI_ALL_TS_knockout%>%
  select(`MGI Gene ID` ,Structure,Detected, `Theiler Stage`) %>% # select columns of interest 
  # we will try to always use MGI IDs for mouse genes / HGNC ids for human genes (these are stable identifiers)
  distinct(`MGI Gene ID`, Structure,Detected, `Theiler Stage`) %>% # retain unique rows
  arrange(Detected) %>% ## sort by this variable
  group_by(`MGI Gene ID`,Structure,  `Theiler Stage`) %>% # group data by gene and tissue %>%
  summarise(Detected_all = paste0(Detected,collapse="|")) %>% ## collpase all the expression data availabe for a given
  # gene and tissue  in one cell
  mutate(Detected_all = ifelse(Detected_all =="No","No",
                               ifelse(Detected_all == "Yes","Yes","Ambiguous_NotSpecified_Contradictory"))) %>%
  # recode variable based on values 
  rename(MGI.ID = `MGI Gene ID`, ALL.TS.mutant.expression.detected = Detected_all) %>% ## rename variables
  select(MGI.ID,Structure,ALL.TS.mutant.expression.detected, `Theiler Stage`) %>% # select columns of interest
  distinct(MGI.ID,Structure,ALL.TS.mutant.expression.detected, `Theiler Stage`) %>% # check again for duplicated rows
  mutate(Gene.Anatomy = paste0(MGI.ID,"-",Structure)) # create a new variable (Gene - Anatomy (Structure) term)
  

View(MGI_ALL_TS_knockout)



#write.csv(MGI_ALL_TS_knockout, "D:/MSC RESEARCH PROJECT/Output_Files/Mice/MGI_ALL_TS_knockout.csv")

################################################################################################################################
################################################################################################################################

############~~~~~~~~~~~~~~~SUMMARISE PER GENE, EMBRYONIC STAGE AND THE NUMBER OF YES/NO/AMBIGIOUS~~~~~~~~~~~~~~~~~~~~~~~~~######################## 


#per gene, embryonic stage and tissue the number of  Yes / No in order to classify each gene /embryonic stage / tissue as Yes/No/Ambiguous.
list_structures<- ungroup(MGI_ALL_TS_knockout)%>%
  select(Structure)%>%
  distinct()

View(list_structures)


MGI_ALL_TS_knockout_All<- MGI_ALL_TS_knockout%>%
  select(MGI.ID,Structure,ALL.TS.mutant.expression.detected, `Theiler Stage` )%>%
  distinct()


MGI_ALL_TS_knockout_All$MGI.ID<-as.factor(MGI_ALL_TS_knockout_All$MGI.ID )

MGI_ALL_TS_knockout_All$Structure<-as.factor(MGI_ALL_TS_knockout_All$Structure )

# Build the heat map from scratch

cols  <- c(
  "Yes" = "Green",
  "No" = "Red",
  "Ambiguous_NotSpecified_Contradictory" = "Black",
  na.value  = "Grey")

MGI_ALL_TS_knockout_All_p1<-ggplot(MGI_ALL_TS_knockout_All, aes(x =MGI_ALL_TS_knockout_All$MGI.ID, 
                                                            y =MGI_ALL_TS_knockout_All$Structure))%>%
  +labs(x= "MGI.ID", y="Structure", 
        subtitle= "KO Mice Developmental Gene Expression in Each Structure for Each Theiler Stage")%>%
  +theme(axis.title=element_text(size=10,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=1),
         axis.text.x= element_text(size=1),
         panel.border = element_rect(colour = "black", fill=NA, size=1),
         panel.background = element_rect(fill = 'lightgrey')) %>%
  + guides(fill=guide_legend(title="Gene Expression: ")) %>%
  +geom_tile(aes(fill = MGI_ALL_TS_knockout_All$ALL.TS.mutant.expression.detected)) %>%
  +scale_fill_manual(values = cols, na.value = "Grey") %>%# Adjust colors
  +facet_wrap( . ~ MGI_ALL_TS_knockout_All$`Theiler Stage`)%>% # Facet layer
  +coord_flip()



MGI_ALL_TS_knockout_All_p1




# UNHASH TO SAVE PLOT!!!!!!!!!!!!!!!!!
# 
# 
# library(cowplot)
#  save_plot("./Plots/Mouse/Developmental Stages/Heatmap/MGI_ALL_TS_knockout_All.pdf",
#            MGI_ALL_TS_knockout_All_p1 ,base_height= 6 ,base_aspect_ratio = 3) 
# # 







############~~~~~~~~~~~~~~~~~~CREATE FUNCTION FOR EACH STRUCTURE TO HAVE A HEATMAP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#############################




Heatmap_Structure_Gene<-function(strctr){
  
  
  ## create function
  MGI_ALL_TS_knockout_1<- MGI_ALL_TS_knockout%>%
    select(MGI.ID,Structure,ALL.TS.mutant.expression.detected, `Theiler Stage` )%>%
    filter(Structure  == strctr)%>%
    distinct()
  
  View(MGI_ALL_TS_knockout_1)
  
  
  MGI_ALL_TS_knockout_1$MGI.ID<-as.factor(MGI_ALL_TS_knockout_1$MGI.ID )
  
  MGI_ALL_TS_knockout_1$Structure<-as.factor(MGI_ALL_TS_knockout_1$Structure )
  
  # Build the heat map from scratch
  
  cols  <- c(
    "Yes" = "Green",
    "No" = "Red",
    "Ambiguous_NotSpecified_Contradictory" = "Black",
    na.value  = "Grey")
  
  MGI_ALL_TS_knockout_1_p2<-ggplot(MGI_ALL_TS_knockout_1, aes(x =MGI_ALL_TS_knockout_1$MGI.ID, 
                                    y =MGI_ALL_TS_knockout_1$Structure))%>%
    +labs(x= "MGI.ID", y="Structure", 
         subtitle= paste("KO Mice Developmental Gene Expression for Each Theiler Stage in:  \n ", strctr, "Structure"))%>%
    +theme(axis.title=element_text(size=10,face="bold"),
           plot.caption=element_text(face = "italic", size=12, hjust = 0),
           text = element_text(size=10),
           legend.position = "bottom",legend.direction = "horizontal",
           legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
           axis.text.y= element_text(size=1),
           axis.text.x= element_text(size=10),
           panel.border = element_rect(colour = "black", fill=NA, size=1),
           panel.background = element_rect(fill = 'lightgrey')) %>%
    + guides(fill=guide_legend(title="Gene Expression: ")) %>%
    +geom_tile(aes(fill = MGI_ALL_TS_knockout_1$ALL.TS.mutant.expression.detected)) %>%
    +scale_fill_manual(values = cols, na.value = "Grey") %>%# Adjust colors
    +facet_wrap( . ~ MGI_ALL_TS_knockout_1$`Theiler Stage`)%>% # Facet layer
    +coord_flip()
  
  
  
  #MGI_ALL_TS_knockout_1_p2
  
  return(MGI_ALL_TS_knockout_1_p2)
  
  
  }
# 
# 
# # !!!!!!! UNHASH TO USE THE FOR-LOOP !!!!!!!!!!
# 
# for(i in 1:nrow(list_structures)) {
#    ii <- list_structures[i,1]     #FIRST COLUMN
#    print(ii) # TO CHECK IF THE CORRECT ROW IS BEING TAKEN
#  
#    a= ii 
#    
#    
#    
#    MGI_ALL_TS_knockout_1_p2<-Heatmap_Structure_Gene(ii)
#     
#    
#     # REMOVE SPACES, /, -, SO THAT THE FILE CAN BE SAVED 
#    
#    a<-gsub(" ", "_", ii, fixed=TRUE)
#    a<-gsub("/", "_", a, fixed=TRUE)
#    a<-gsub("-", "_", a, fixed=TRUE)
#    
#    print(a)
# 
# # SAVE PLOT----CHANGE LOCATION!
#    library(cowplot)
#    
#    save_plot(paste("./Plots/Mouse/Developmental Stages/Heatmap/", a,".pdf"),
#               MGI_ALL_TS_knockout_1_p2 ,base_height= 5.5 ,base_aspect_ratio = 2) 
#    
#    
#  }
















#######################################################################################################
########################################################################################################

###~~~~~~~~~~~~~~~~~~~~~~~~~~Number of Genes Detected for Each Structure! ~~~~~~~~~#######################

Genes_Structure<-ungroup(MGI_ALL_TS_knockout)%>%
  select(Structure, MGI.ID, ALL.TS.mutant.expression.detected)%>%
  group_by(Structure, MGI.ID, ALL.TS.mutant.expression.detected )%>%
  summarise('Number.of.Genes'= n())%>%
  distinct(Structure, MGI.ID, Number.of.Genes, ALL.TS.mutant.expression.detected  )

Genes_Structure_sum<-aggregate(Genes_Structure$Number.of.Genes , 
                        by=list(Structure=Genes_Structure$Structure, 
                                ALL.TS.mutant.expression.detected=Genes_Structure$ALL.TS.mutant.expression.detected)
                        , FUN=sum)

View(Genes_Structure_sum)



Genes_Structure_sum_p1<-ggplot(Genes_Structure_sum, aes(fill=Genes_Structure_sum$ALL.TS.mutant.expression.detected, 
                                                       y=Genes_Structure_sum$x,
                                                       x=Genes_Structure_sum$Structure ))%>% 
  
  + labs(x= "Structure", y="Number of Genes", 
         subtitle="Number of Genes Detected for Each Structure for KO Mice")%>%
  + guides(fill=guide_legend(title="Gene Expression: ")) %>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=8),
         axis.text.x= element_text(size=10),
         panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  +geom_bar(stat="identity") %>%  
  +coord_flip()

Genes_Structure_sum_p1


# UNHASH TO SAVE PLOT!!!!!!!!!!!!!!!!!
# 
# 
# library(cowplot)
# save_plot(paste("./Plots/Mouse/Developmental Stages/Structure/Each_Structure.pdf"),
#           Genes_Structure_sum_p1 ,base_height= 5.5 ,base_aspect_ratio = 2) 





###~~~~~~~~~~~~~~~~~~~~~~~~~~Number of Genes Detected for Each TS! ~~~~~~~~~#######################

Genes_TS<-ungroup(MGI_ALL_TS_knockout)%>%
  select(`Theiler Stage`, MGI.ID, ALL.TS.mutant.expression.detected)%>%
  group_by(`Theiler Stage`, MGI.ID, ALL.TS.mutant.expression.detected )%>%
  summarise('Number.of.Genes'= n())%>%
  distinct(`Theiler Stage`, MGI.ID, Number.of.Genes, ALL.TS.mutant.expression.detected  )

Genes_TS_sum<-aggregate(Genes_TS$Number.of.Genes , 
                        by=list(Theiler.Stage=Genes_TS$`Theiler Stage`, 
                                ALL.TS.mutant.expression.detected=Genes_TS$ALL.TS.mutant.expression.detected)
                        , FUN=sum)

View(Genes_TS_sum)



Genes_TS_sum_p1<-ggplot(Genes_TS_sum, aes(fill=Genes_TS_sum$ALL.TS.mutant.expression.detected, 
                                          y=Genes_TS_sum$x,
                                          x=Genes_TS_sum$Theiler.Stage))%>% 
  
  + labs(x= "Theiler Stage", y="Number of Genes", 
         subtitle="Number of Genes Detected for Each Theiler Stage for KO Mice")%>%
  + guides(fill=guide_legend(title="Gene Expression: ")) %>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=8),
         axis.text.x= element_text(size=10),
         panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  + scale_x_continuous(breaks = seq(0, 28, by = 1))%>%
  +geom_bar(stat="identity") %>%  
  +geom_text(aes(label=paste0(Genes_TS_sum$x)),
             size = 3, position = position_stack(vjust = 0.5), 
             colour=c("black"), fontface='bold')%>%
  +coord_flip()

Genes_TS_sum_p1


# UNHASH TO SAVE PLOT!!!!!!!!!!!!!!!!!
# 
# 
# library(cowplot)
# save_plot(paste("./Plots/Mouse/Developmental Stages/Theiler_Stages/Each_TS.pdf"),
#           Genes_TS_sum_p1 ,base_height= 5.5 ,base_aspect_ratio = 2) 
# 

