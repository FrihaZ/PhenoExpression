
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

library(dplyr);library(tidyr);library(readr);library(ggplot2)

################################################################################################################################
################################################################################################################################

# Load Files
MGI_ALL_TS_knockout_<-read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/MGIgeneExpressionQuery_20190626_131124_knockout.txt",delim = "\t")

MGI_ALL_TS_WT<-read_delim("D:/MSC RESEARCH PROJECT/GITHUB/data/MGIgeneExpressionQuery_20190711_073852_WT.txt",delim = "\t")


################################################################################################################################
################################################################################################################################

###############~~~~~~~~~~~~~~~~~~~~~~~~CLEAN KO DATA AND MAP THEM TO MAMMALIAN PHENOTYPES~~~~~~~~~~~~#########

################################################################################################################################
################################################################################################################################

MGI_ALL_TS_knockout<-MGI_ALL_TS_knockout_%>%
  select(`MGI Gene ID` ,Structure,Detected, `Theiler Stage`) %>% # select columns of interest 
  # we will try to always use MGI IDs for mouse genes / HGNC ids for human genes (these are stable identifiers)
  distinct(`MGI Gene ID`, Structure,Detected, `Theiler Stage`) %>% # retain unique rows
  arrange(Detected) %>% ## sort by this variable
  group_by(`MGI Gene ID`,Structure,  `Theiler Stage`) %>% # group data by gene and tissue %>%
  summarise(Detected_all = paste0(Detected,collapse="|")) %>% ## collpase all the expression data availabe for a given
  # gene and tissue  in one cell
  mutate(Detected_all = ifelse(Detected_all =="No","No",
                               ifelse(Detected_all == "Yes","Yes","Ambiguous"))) %>%
  # recode variable based on values 
  rename(MGI.ID = `MGI Gene ID`,ALL.TS.expression.detected = Detected_all) %>% ## rename variables
  select(MGI.ID,Structure,ALL.TS.expression.detected, `Theiler Stage`) %>% # select columns of interest
  distinct(MGI.ID,Structure,ALL.TS.expression.detected, `Theiler Stage`) %>% # check again for duplicated rows
  mutate(Gene.Anatomy = paste0(MGI.ID,"-",Structure)) # create a new variable (Gene - Anatomy (Structure) term)

View(MGI_ALL_TS_knockout)

## count number of KO genes
count_MGI_ALL_TS_knockout<-ungroup(MGI_ALL_TS_knockout)%>%
  select(MGI.ID)%>%
  drop_na()%>%
  distinct()

nrow(count_MGI_ALL_TS_knockout)

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
  select(MGI.ID,Structure,ALL.TS.expression.detected, `Theiler Stage` )%>%
  distinct()


MGI_ALL_TS_knockout_All$MGI.ID<-as.factor(MGI_ALL_TS_knockout_All$MGI.ID )

MGI_ALL_TS_knockout_All$Structure<-as.factor(MGI_ALL_TS_knockout_All$Structure )

# Build the heat map from scratch

cols  <- c(
  "Yes" = "Green",
  "No" = "Red",
  "Ambiguous" = "Grey",
  na.value  = "Black")

MGI_ALL_TS_knockout_All_p1<-ggplot(MGI_ALL_TS_knockout_All, aes(x =MGI_ALL_TS_knockout_All$MGI.ID, 
                                                            y =MGI_ALL_TS_knockout_All$Structure))%>%
  +labs(x= "MGI.ID", y="Structure", 
        subtitle= "KO Mice Developmental Gene Expression in Each Structure for Each Theiler Stage")%>%
  +theme(axis.title=element_text(size=10,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=0.02),
         axis.text.x= element_text(size=0.02),
         panel.border = element_rect(colour = "black", fill=NA, size=1),
         panel.background = element_rect(fill = 'black')) %>%
  + guides(fill=guide_legend(title="Gene Expression: ")) %>%
  +geom_tile(aes(fill = MGI_ALL_TS_knockout_All$ALL.TS.expression.detected)) %>%
  +scale_fill_manual(values = cols, na.value = "black") %>%# Adjust colors
  +facet_wrap( . ~ MGI_ALL_TS_knockout_All$`Theiler Stage`)%>% # Facet layer
  +coord_flip()



MGI_ALL_TS_knockout_All_p1




# UNHASH TO SAVE PLOT!!!!!!!!!!!!!!!!!
# 
# # 
# library(cowplot)
# save_plot("./Plots/Mouse/MGI GXD Developmental Stages/KO/Heatmap/MGI_ALL_TS_knockout_All.pdf",
#           MGI_ALL_TS_knockout_All_p1 ,base_height= 6 ,base_aspect_ratio = 3) 
# # # 


##########################################################################################################
######################################################################################################################





############~~~~~~~~~~~~~~~~~~CREATE FUNCTION FOR EACH STRUCTURE TO HAVE A HEATMAP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#############################




Heatmap_Structure_Gene<-function(strctr, MGI_ALL_TS,GENO){
  
  
  ## create function
  MGI_ALL_TS_knockout_1<- MGI_ALL_TS%>%
    select(MGI.ID,Structure,ALL.TS.expression.detected, `Theiler Stage` )%>%
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
         subtitle= paste(GENO," Mice Developmental Gene Expression for Each Theiler Stage in:  \n ", strctr, "Structure"))%>%
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
    +geom_tile(aes(fill = MGI_ALL_TS_knockout_1$ALL.TS.expression.detected)) %>%
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
#   ii <- list_structures[i,1]     #FIRST COLUMN
#   print(ii) # TO CHECK IF THE CORRECT ROW IS BEING TAKEN
# 
#   a= ii 
#     
#     
#     
#   MGI_ALL_TS_knockout_1_p2<-Heatmap_Structure_Gene(strctr=ii, MGI_ALL_TS=MGI_ALL_TS_knockout, GENO="KO")
#      
#     
#    # REMOVE SPACES, /, -, SO THAT THE FILE CAN BE SAVED 
#     
#   a<-gsub(" ", "_", ii, fixed=TRUE)
#   a<-gsub("/", "_", a, fixed=TRUE)
#   a<-gsub("-", "_", a, fixed=TRUE)
# #    
#   print(a)
# # 
#       #SAVE PLOT----CHANGE LOCATION!
#   library(cowplot)
#     
#   save_plot(paste("./Plots/Mouse/MGI GXD Developmental Stages/KO/Heatmap/", a,".pdf"),
#             MGI_ALL_TS_knockout_1_p2 ,base_height= 5.5 ,base_aspect_ratio = 2) 
#    
#   
# }
















#######################################################################################################
########################################################################################################

###~~~~~~~~~~~~~~~~~~~~~~~~~~Number of Genes Detected for Each Structure! ~~~~~~~~~#######################

Genes_Structure<-ungroup(MGI_ALL_TS_knockout)%>%
  select(Structure, MGI.ID, ALL.TS.expression.detected)%>%
  distinct()%>%
  group_by(Structure, ALL.TS.expression.detected)%>%
  summarise(Number.of.Genes=n()) %>%
  mutate('Percentage' = Number.of.Genes/sum(Number.of.Genes)*100)%>%
  distinct(Structure, ALL.TS.expression.detected, Number.of.Genes,Percentage)

View(Genes_Structure)


### CREATE PLOT WITH HIGHEST AND SMALLEST NUMBER OF GENES

Genes_Structure_sum_sorted_high_low<-Genes_Structure%>%
  select_all()%>%
  filter(Structure =="embryo"|
           Structure =="heart"| 
           Structure=="lung"|
           Structure=="liver"|
           Structure=="brain" | 
           Structure =="metanephros"| 
           Structure== "neural tube" | 
           Structure=="cerebral cortex"|
           Structure=="telencephalon"|
           Structure== "forelimb bud"|
           Structure=="spinal cord")%>%
  distinct()




library(ggplot2)
Genes_Structure_sum_high_low<-ggplot(Genes_Structure_sum_sorted_high_low, 
                                     aes(fill=Genes_Structure_sum_sorted_high_low$ALL.TS.expression.detected, 
                                            y=Genes_Structure_sum_sorted_high_low$Percentage,
                                            x= reorder(Genes_Structure_sum_sorted_high_low$Structure,-Genes_Structure_sum_sorted_high_low$Percentage ) ))%>% 
  
  + labs(x= "Structure", y="Percentage of Genes (%)", 
         subtitle=paste("Percentage of Genes Detected for Each Structure for KO Mice"))%>%
  + guides(fill=guide_legend(title="Gene Expression: ")) %>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10, angle = 45, hjust = 1, vjust = 1),
         panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  + scale_y_continuous(breaks = seq(0, 100, by = 20))%>%
  +geom_bar(stat="identity") %>%
  +geom_text(aes(label=paste0(round(Genes_Structure_sum_sorted_high_low$Percentage,0),"%")),size = 4, position = position_stack(vjust = 0.5),
             colour=c("black"), fontface='bold') 

Genes_Structure_sum_high_low


# UNHASH TO SAVE PLOT!!!!!!!!!!!!!!!!!
# # 
# # # # # 
# library(cowplot)
# save_plot(paste("./Plots/Mouse/MGI GXD Developmental Stages/KO/Structure/Genes_Structure_sum_high_low.png"),
#         Genes_Structure_sum_high_low ,base_height= 6 ,base_aspect_ratio = 1) 
# # # # 
# 













###~~~~~~~~~~~~~~~~~~~~~~~~~~Number of Genes Detected for Each TS! ~~~~~~~~~#######################

Genes_TS<-ungroup(MGI_ALL_TS_knockout)%>%
  select(`Theiler Stage`, MGI.ID, ALL.TS.expression.detected)%>%
  group_by(`Theiler Stage`, MGI.ID, ALL.TS.expression.detected )%>%
  summarise('Number.of.Genes'= n())%>%
  distinct(`Theiler Stage`, MGI.ID, Number.of.Genes, ALL.TS.expression.detected  )

Genes_TS_sum<-aggregate(Genes_TS$Number.of.Genes , 
                        by=list(Theiler.Stage=Genes_TS$`Theiler Stage`, 
                                ALL.TS.expression.detected=Genes_TS$ALL.TS.expression.detected)
                        , FUN=sum)

View(Genes_TS_sum)

sum(Genes_TS_sum$x) # to get the total number of genes

Genes_TS_sum_p1<-ggplot(Genes_TS_sum, aes(fill=Genes_TS_sum$ALL.TS.expression.detected , 
                                          y=Genes_TS_sum$x,
                                          x=Genes_TS_sum$Theiler.Stage))%>% 
  
  + labs(x= "TS", y="Number of Genes", 
         subtitle="Number of Genes Detected for Each Theiler Stage (TS) in KO Mice")%>%
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
# save_plot(paste("./Plots/Mouse/MGI GXD Developmental Stages/KO/Theiler_Stages/Each_TS.png"),
#            Genes_TS_sum_p1 ,base_height= 5.5 ,base_aspect_ratio = 2) 
# 





































###############~~~~~~~~~~~~~~~~~~~~~~~~WILD TYPE MICE DEVELOPMENTAL SATGES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#########





################################################################################################################################
################################################################################################################################

###############~~~~~~~~~~~~~~~~~~~~~~~~CLEAN WILD TYPE DATA AND MAP THEM TO MAMMALIAN PHENOTYPES~~~~~~~~~~~~#########

################################################################################################################################
################################################################################################################################




MGI_ALL_TS_WT_<-MGI_ALL_TS_WT%>%
  select(`MGI Gene ID` ,Structure,Detected, `Theiler Stage`) %>% # select columns of interest 
  # we will try to always use MGI IDs for mouse genes / HGNC ids for human genes (these are stable identifiers)
  distinct(`MGI Gene ID`, Structure,Detected, `Theiler Stage`) %>% # retain unique rows
  arrange(Detected) %>% ## sort by this variable
  group_by(`MGI Gene ID`,Structure,  `Theiler Stage`) %>% # group data by gene and tissue %>%
  summarise(Detected_all = paste0(Detected,collapse="|")) %>% ## collpase all the expression data availabe for a given
  # gene and tissue  in one cell
  mutate(Detected_all = ifelse(Detected_all =="No","No",
                               ifelse(Detected_all == "Yes","Yes","Ambiguous"))) %>%
  # recode variable based on values 
  rename(MGI.ID = `MGI Gene ID`, ALL.TS.expression.detected = Detected_all) %>% ## rename variables
  select(MGI.ID,Structure,ALL.TS.expression.detected, `Theiler Stage`) %>% # select columns of interest
  distinct(MGI.ID,Structure,ALL.TS.expression.detected, `Theiler Stage`) %>% # check again for duplicated rows
  mutate(Gene.Anatomy = paste0(MGI.ID,"-",Structure)) # create a new variable (Gene - Anatomy (Structure) term)


View(MGI_ALL_TS_WT_)



count_MGI_ALL_TS_WT_<-ungroup(MGI_ALL_TS_WT_)%>%
  select(MGI.ID)%>%
  drop_na()%>%
  distinct()

nrow(count_MGI_ALL_TS_WT_)


## UNHASH TO SAVE CSV
#write.csv(MGI_ALL_TS_WT_, "D:/MSC RESEARCH PROJECT/Output_Files/Mice/MGI_ALL_TS_WT.csv")

################################################################################################################################
################################################################################################################################

############~~~~~~~~~~~~~~~SUMMARISE PER GENE, EMBRYONIC STAGE AND THE NUMBER OF YES/NO/AMBIGIOUS~~~~~~~~~~~~~~~~~~~~~~~~~######################## 


#per gene, embryonic stage and tissue the number of  Yes / No in order to classify each gene /embryonic stage / tissue as Yes/No/Ambiguous.
WT_list_structures<- ungroup(MGI_ALL_TS_WT_)%>%
  select(Structure)%>%
  distinct()

View(WT_list_structures)


MGI_ALL_TS_WT__All<- MGI_ALL_TS_WT_ %>%
  select(MGI.ID ,Structure,ALL.TS.expression.detected,
         `Theiler Stage` )%>%
  distinct()

View(MGI_ALL_TS_WT_)

MGI_ALL_TS_WT__All$MGI.ID<-as.factor(MGI_ALL_TS_WT__All$MGI.ID )

MGI_ALL_TS_WT__All$Structure<-as.factor(MGI_ALL_TS_WT__All$Structure )

View(MGI_ALL_TS_WT__All)

#########################################################################################################################

# Build the heat map from scratch


cols  <- c(
  "Yes" = "Green",
  "No" = "Red",
  "Ambiguous" = "Black",
  na.value  = "Grey")

MGI_ALL_TS_WT__All_P1<-ggplot(MGI_ALL_TS_WT__All, aes(x =MGI_ALL_TS_WT__All$MGI.ID, 
                                                                y =MGI_ALL_TS_WT__All$Structure ))%>%
  +labs(x= "MGI.ID", y="Structure", 
        subtitle= "WT Mice Developmental Gene Expression in Each Structure for Each Theiler Stage")%>%
  +theme(axis.title=element_text(size=10,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=0.02),
         axis.text.x= element_text(size=0.02),
         panel.border = element_rect(colour = "black", fill=NA, size=1),
         panel.background = element_rect(fill = 'lightgrey')) %>%
  + guides(fill=guide_legend(title="Gene Expression: ")) %>%
  +geom_tile(aes(fill = MGI_ALL_TS_WT__All$ALL.TS.expression.detected)) %>%
  +scale_fill_manual(values = cols, na.value = "Grey") %>%# Adjust colors
  +facet_wrap( . ~ MGI_ALL_TS_WT__All$`Theiler Stage`)%>% # Facet layer
  +coord_flip()



MGI_ALL_TS_WT__All_P1




# UNHASH TO SAVE PLOT!!!!!!!!!!!!!!!!!
# 
# # 
# # library(cowplot)
# save_plot("./Plots/Mouse/MGI GXD Developmental Stages/WT/Heatmap/MGI_ALL_TS_WT__All_P1.png",
#          MGI_ALL_TS_WT__All_P1 ,base_height= 6 ,base_aspect_ratio = 3) 
# # # # 







################################################################################################
###########################################################################################

########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~USE FUNCTION FOR WT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##########################


# # !!!!!!! UNHASH TO USE THE FOR-LOOP !!!!!!!!!!
# 
# for(t in 1:nrow(WT_list_structures)) {
#   tt <- WT_list_structures[t,1]     #FIRST COLUMN
#   print(tt) # TO CHECK IF THE CORRECT ROW IS BEING TAKEN
# #  
# #    b= tt
# #    
# #    
# #    
#   MGI_ALL_TS_WT__1_p2<-Heatmap_Structure_Gene(strctr=tt, MGI_ALL_TS=MGI_ALL_TS_WT_, GENO="KO")
# #     
# #    
# #     # REMOVE SPACES, /, -, SO THAT THE FILE CAN BE SAVED 
# #    
#   tb<-gsub(" ", "_", tt, fixed=TRUE)
#   tb<-gsub("/", "_", tb, fixed=TRUE)
#   tb<-gsub("-", "_", tb, fixed=TRUE)
#   
#   print(tb)
# # 
# # # SAVE PLOT----CHANGE LOCATION!
#   library(cowplot)
#  
#   save_plot(paste("./Plots/Mouse/MGI GXD Developmental Stages/WT/Heatmap/", tb,".pdf"),
#             MGI_ALL_TS_WT__1_p2 ,base_height= 5.5 ,base_aspect_ratio = 2) 
#     
#  }















#######################################################################################################
########################################################################################################

###~~~~~~~~~~~~~~~~~~~~~~~~~~Number of Genes Detected for Each Structure! ~~~~~~~~~#######################

WT_Genes_Structure<-MGI_ALL_TS_WT_%>%
  select(Structure, MGI.ID , ALL.TS.expression.detected)%>%
  distinct()%>%
  group_by(Structure, ALL.TS.expression.detected)%>%
  summarise(Number.of.Genes=n()) %>%
  mutate('Percentage' = Number.of.Genes/sum(Number.of.Genes)*100)%>%
  distinct(Structure, ALL.TS.expression.detected, Number.of.Genes,Percentage)

View(WT_Genes_Structure)
####################################################################################################################
##############################################################################################################

### CREATE PLOT WITH HIGHEST AND SMALLEST NUMBER OF GENES

WT_Genes_Structure_sum_high_low<-WT_Genes_Structure%>%
  select_all()%>%
  filter(Structure =="brain"|
           Structure =="spinal cord"| 
           Structure=="embryo"|
           Structure=="metanephros"|
           Structure=="cerebral cortex"|
           Structure=="lung" | 
           Structure =="midbrain"| 
           Structure== "testis" | 
           Structure=="liver"|
           Structure== "heart")%>%
  distinct()

View(WT_Genes_Structure_sum_high_low)




library(ggplot2)
WT_Genes_Structure_sum_high_low_P1<-ggplot(WT_Genes_Structure_sum_high_low, 
                                        aes(fill=WT_Genes_Structure_sum_high_low$ALL.TS.expression.detected, 
                                                                       y=WT_Genes_Structure_sum_high_low$Percentage ,
                                                                       x= reorder(WT_Genes_Structure_sum_high_low$Structure,-WT_Genes_Structure_sum_high_low$Percentage ) ))%>% 
  
  + labs(x= "Structure", y="Percentage of Genes (%)", 
         subtitle=paste("Percentage of Genes Detected for Each Structure for WT Mice"))%>%
  + guides(fill=guide_legend(title="Gene Expression: ")) %>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10, angle = 45, hjust = 1, vjust = 1),
         panel.border = element_rect(colour = "black", fill=NA, size=1)) %>%
  + scale_y_continuous(breaks = seq(0, 100, by = 20))%>%
  +geom_bar(stat="identity") %>%  
  +geom_text(aes(label=paste0(round(WT_Genes_Structure_sum_high_low$Percentage,0),"%")),size = 4, position = position_stack(vjust = 0.5),
             colour=c("black"), fontface='bold') 


#WT_Genes_Structure_sum_high_low_P1


# UNHASH TO SAVE PLOT!!!!!!!!!!!!!!!!!
# # 
# # # # 
# library(cowplot)
# save_plot(paste("./Plots/Mouse/MGI GXD Developmental Stages/WT/Structure/WT_Genes_Structure_sum_high_low.png"),
#          WT_Genes_Structure_sum_high_low_P1 ,base_height= 6 ,base_aspect_ratio = 1) 
# # # # 









###~~~~~~~~~~~~~~~~~~~~~~~~~~Number of Genes Detected for Each TS! ~~~~~~~~~#######################

View(MGI_ALL_TS_WT_)
MGI_ALL_TS_WT_<-ungroup(MGI_ALL_TS_WT_)%>%
  select(`Theiler Stage`, MGI.ID, ALL.TS.expression.detected)%>%
  distinct()%>%
  group_by(`Theiler Stage`, ALL.TS.expression.detected )%>%
  summarise('Number.of.Genes'= n())%>%
  distinct(`Theiler Stage`, MGI.ID, Number.of.Genes, ALL.TS.expression.detected  )

WT_Genes_TS_sum<-aggregate(MGI_ALL_TS_WT_$Number.of.Genes , 
                        by=list(Theiler.Stage=MGI_ALL_TS_WT_$`Theiler Stage`, 
                                ALL.TS.expression.detected=MGI_ALL_TS_WT_$ALL.TS.expression.detected)
                        , FUN=sum)

View(WT_Genes_TS_sum)

sum(WT_Genes_TS_sum$x)


WT_Genes_TS_sum_p1<-ggplot(WT_Genes_TS_sum, aes(fill=WT_Genes_TS_sum$ALL.TS.expression.detected, 
                                          y=WT_Genes_TS_sum$x,
                                          x=WT_Genes_TS_sum$Theiler.Stage ))%>% 
  
  + labs(x= "TS", y="Number of Genes", 
         subtitle="Number of Genes Detected for Each Theiler Stage (TS) for WT Mice")%>%
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
  +geom_text(aes(label=paste0(WT_Genes_TS_sum$x)),
             size = 3, position = position_stack(vjust = 0.5), 
             colour=c("black"), fontface='bold')%>%
  +coord_flip()

#WT_Genes_TS_sum_p1


# UNHASH TO SAVE PLOT!!!!!!!!!!!!!!!!!
# 
# # # 
# library(cowplot)
# save_plot(paste("./Plots/Mouse/MGI GXD Developmental Stages/WT/Theiler_Stages/Each_TS.png"),
#            WT_Genes_TS_sum_p1 ,base_height= 5.5 ,base_aspect_ratio = 2) 


