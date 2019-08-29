
###############################################################################################################################
###############################################################################################################################

### Project: PhenoExpression ##################################################################################################
### Script: PCA_DDMD_Viabilty.R ###################################################################################################
### Purpose: To conduct the Principle Component Analysis and t-SNE on the Viabilty dataset######################################################
### Author: Friha Zafar ####################################################################################################
### Date: 19/07/2019 ##########################################################################################################

###############################################################################################################################
###############################################################################################################################

## load packages ##############################################################################################################

library(dplyr);library(tidyr);library(stringr);library(readr);library(tidyverse);library(ggpubr);library(ggplot2)

################################################################################################################################
################################################################################################################################

# Load Files

prim_viab_mouse_dmdd<-read.csv("D:/MSC RESEARCH PROJECT/Output_Files/Mice/EA_DMDDs DATASET/Viability/primary_viability_mouse_dmdd.csv")

CLEANED_prim_viab_mouse_dmdd<-read.csv("D:/MSC RESEARCH PROJECT/Output_Files/Mice/EA_DMDDs DATASET/Viability/primary_viability_mouse_dmdd_CLEANED.csv")

second_viab_mouse_dmdd<-read.csv("D:/MSC RESEARCH PROJECT/Output_Files/Mice/EA_DMDDs DATASET/Viability/second_viability_mouse_dmdd.csv")

ST.ER.second_viab_mouse_dmdd<-read.csv("D:/MSC RESEARCH PROJECT/Output_Files/Mice/EA_DMDDs DATASET/Viability/SE_mean_second_viability_mouse_dmdd.csv")


impc_secondary_viability<-read.delim("D:/MSC RESEARCH PROJECT/GITHUB/data/mouse_viability_data_files/impc_secondary_viability.txt")

##################################################################################################################################################
################################################################################################################################################


CLEANED_prim_viab_mouse_dmdd<-CLEANED_prim_viab_mouse_dmdd%>%
  select_all()%>%
  select(-X)%>%
  distinct()

View(CLEANED_prim_viab_mouse_dmdd)


############~~~~~~~~~~~~~~~~~~PCA ON PRIMARY VIABILITY DATA~~~~~~~~~~~~~~~~~#################


pca_CLEANED_prim_viab<-CLEANED_prim_viab_mouse_dmdd%>%
  select_all()%>%
  distinct(MGI.Accession.ID, Expression, Development.Stage, IMPC_Viability)

View(pca_CLEANED_prim_viab)

pca_CLEANED_prim_viab$Expression = as.numeric(as.factor(pca_CLEANED_prim_viab$Expression ))

library(reshape2)

pca_2_CLEANED_prim_viab<-dcast(pca_CLEANED_prim_viab,
                               MGI.Accession.ID + IMPC_Viability  ~ Development.Stage, value.var="Expression")

pca_2_CLEANED_prim_viab<-column_to_rownames(pca_2_CLEANED_prim_viab, var = "MGI.Accession.ID")



##replace NAs with 0
pca_2_CLEANED_prim_viab[is.na(pca_2_CLEANED_prim_viab)] <- 0

str(pca_2_CLEANED_prim_viab)

######################################################

##CONDUCT PCA
pca_3_CLEANED_prim_viab<-pca_2_CLEANED_prim_viab

pca_3_CLEANED_prim_viab$IMPC_Viability<-as.numeric(pca_3_CLEANED_prim_viab$IMPC_Viability)

res.pca_3_CLEANED_prim_viab <- prcomp(pca_3_CLEANED_prim_viab,  scale = TRUE)

attributes(res.pca_3_CLEANED_prim_viab) #the attributes shows the information inside the object of class prcomp

summary_pca_prim_viab <- summary(res.pca_3_CLEANED_prim_viab)
str(summary_pca_prim_viab)


######################################################################
# EXPLAINED VARIANCE 
#The amount of varaiance that can be explained by the PRINCIPLE COMPONENTS
#PC1 explains 76.07% of the VARIATION IN THE data, PCA2 explains 6.591% of the data

expVar_pca_prim_viab  <- summary_pca_prim_viab$importance[2,] * 100
expVar_pca_prim_viab

##PLOT EXPLAINED VARIANCE

expVar_pca_prim_viab_barplot<-barplot(expVar_pca_prim_viab, xlab = "Principal Components", ylab = "Explained Variance (%)",
                        col ="steelblue", ylim=c(0,100),
                        main="The Explained Variance (%) of the Primary Viability Dataset")

## ADD LABELS
text(expVar_pca_prim_viab_barplot, expVar_pca_prim_viab, 
     labels=round(expVar_pca_prim_viab,1), xpd=TRUE, pos=3,) #this is the line of your interest in regards to the value labels


# Cumulative Variance
cumVar__pca_prim_viab <- summary_pca_prim_viab$importance[3,] * 100
cumVar__pca_prim_viab

plot(cumVar__pca_prim_viab, type="o" ,col="black", pch=21, bg="blue", cex=0.8, ylab="Cumulative Variance (%)",
     xlab="Principal Components", main="The Cumulative Variance (%) of the Primary Viability Dataset")


##################################################################################

#DEFINE COLOURS FOR ALL VIABILITIES

# # 
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GEOquery")   
# 
# library(GEOquery)
# pD_pca_2_CLEANED_prim_viab<-as(pca_2_CLEANED_prim_viab,"DataFrame")
# 
# #library(ballgown)

#pDat_pca_2_CLEANED_prim_viab <- pData(pca_2_CLEANED_prim_viab)

pClass_pca_2_CLEANED_prim_viab <- as.factor(pca_2_CLEANED_prim_viab$IMPC_Viability )
#names(pClass_pca_2_CLEANED_prim_viab) <- sampleNames(pca_2_CLEANED_prim_viab)
str(pClass_pca_2_CLEANED_prim_viab)


levels(pClass_pca_2_CLEANED_prim_viab)[levels(pClass_pca_2_CLEANED_prim_viab)=="Lethal"]
levels(pClass_pca_2_CLEANED_prim_viab)[levels(pClass_pca_2_CLEANED_prim_viab)=="Subviable"]
levels(pClass_pca_2_CLEANED_prim_viab)[levels(pClass_pca_2_CLEANED_prim_viab)=="Viable"]

table(pClass_pca_2_CLEANED_prim_viab)

LethalIndx <- which((pClass_pca_2_CLEANED_prim_viab == "Lethal") == TRUE)
SubviableIndx <- which((pClass_pca_2_CLEANED_prim_viab == "Subviable") == TRUE)
ViableIndx <- which((pClass_pca_2_CLEANED_prim_viab == "Viable") == TRUE)

pca_4_CLEANED_prim_viab<-pca_2_CLEANED_prim_viab%>%
  select(-IMPC_Viability)
colour <- NULL
colour[LethalIndx] <- "purple"
colour[SubviableIndx] <- "yellow"
colour[ViableIndx] <- "green"

boxplot(pca_4_CLEANED_prim_viab, col=colour, las=2)
legend("topright", legend=c("Lethal","Subviable", 
                              "Viable"), fill=c("purple", "yellow", "green"),
       lty=1:2, cex=0.8,box.lty=0)





#################################################################
## PLOT THE XSCORES OF PC 1 AND PC2


Xscores_res.pca_3_CLEANED_prim_viab <- res.pca_3_CLEANED_prim_viab$x

View(Xscores_res.pca_3_CLEANED_prim_viab)

plot(Xscores_res.pca_3_CLEANED_prim_viab, xlab="PC1", 
     ylab="PC2", pch=21,bg=colour, cex=1.5,cex.lab=0.7, cex.axis = 0.7,
     main="Principle Component 2 Score versus Principle Component 1 Score for the Primary Viability Dataset")

legend("bottomright", legen=c("Lethal","Subviable", 
                      "Viable"), fill=c("purple", "yellow", "green"),lty=1:2, cex=0.8,
       box.lty=0)

################################################################################

############~~~~~~~~~~~~~~~~~~ANOVA ON PRIMARY VIABILITY DATA~~~~~~~~~~~~~~~~~#################



anova_CLEANED_prim_viab_mouse_dmdd<-CLEANED_prim_viab_mouse_dmdd%>%
  select(IMPC_Viability, Expression, Development.Stage, 
         MGI.Accession.ID)%>%
  distinct()

View(anova_CLEANED_prim_viab_mouse_dmdd)



##### EXPRESSION ~ VIABILITY 
res.anova_CLEANED_prim_viab<- aov(Expression ~ IMPC_Viability, 
                                  data = anova_CLEANED_prim_viab_mouse_dmdd)


summary_anova_prim_viab<-summary(res.anova_CLEANED_prim_viab)

### Multiple pair-wise comparison 

options(scipen=999)

(tk.primary<-TukeyHSD(res.anova_CLEANED_prim_viab))


format(tk.primary$IMPC_Viability[,"p adj"], scientific = TRUE)
with(par(mai=c(1,2.5,1,1)),{plot(tk.primary, las=1,cex.axis=1)})

# 
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = Expression ~ IMPC_Viability, data = anova_CLEANED_prim_viab_mouse_dmdd)
# 
# $IMPC_Viability
#                    diff       lwr        upr      p adj
# Subviable-Lethal -10.07040 -10.77536  -9.365436     0
# Viable-Lethal    -21.35049 -21.81187 -20.889114     0
# Viable-Subviable -11.28010 -11.93782 -10.622372     0

















##### EXPRESSION ~ VIABILITY * DS
res.anova_CLEANED_prim_viab2<- aov(Expression ~ IMPC_Viability +Development.Stage , 
                                  data = anova_CLEANED_prim_viab_mouse_dmdd)


summary_anova_prim_viab<-summary(res.anova_CLEANED_prim_viab2)

### Multiple pair-wise comparison 

tuk<-TukeyHSD(res.anova_CLEANED_prim_viab2)
format(tuk$IMPC_Viability[,"p adj"], scientific = TRUE)



















# ############################################################################################################################################################################




anova_CLEANED_prim_viab_mouse_dmdd$Development.Stage <- factor(anova_CLEANED_prim_viab_mouse_dmdd$Development.Stage,
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

## DO ANOVA



## PLOT ANOVA RESULTS
library(ggplot2)

plot_anova_DS.VIAB <- ggplot(anova_CLEANED_prim_viab_mouse_dmdd, aes(x =anova_CLEANED_prim_viab_mouse_dmdd$Development.Stage,
                                                                 y = anova_CLEANED_prim_viab_mouse_dmdd$Expression,
                                                                 fill = anova_CLEANED_prim_viab_mouse_dmdd$IMPC_Viability ))%>%
  +labs(x= "Developmental Stage", y="Gene Expression (TPM)", 
        title= "Mice Gene Expression (TPM) in Different Developmental Stages and the Associated Viability ")%>%
  +scale_y_continuous(breaks = seq(0, 160, by = 20))%>%
  
  +scale_fill_manual(values=c("purple", "yellow", "green"))  %>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10, angle = 45, hjust = 1, vjust = 1)) %>%
  +guides(fill=guide_legend(title="Viability: "))%>%
  +geom_boxplot(position=position_dodge(1),horizontal=TRUE,
                axes=FALSE,outlier.shape=NA)%>%
  +coord_cartesian(ylim = c(0, 160))%>%
  +geom_hline(yintercept = mean(anova_CLEANED_prim_viab_mouse_dmdd$Expression ), 
                                linetype = 2, color="red")+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 140, label = "p.signif")+stat_compare_means(method = "anova", label.y = 160, mapping=aes(label=..p.adj..))



# 
# 
# library(cowplot)
#  
# save_plot(paste("./Plots/Mouse/EA_Developmental_Stages/Viability/plot_anova_DS_VIAB.png"),
#           plot_anova_DS.VIAB ,base_height= 10 ,base_aspect_ratio = 2) 


## TO GET THE P-VALUES FOR EACH ***
library(ggpubr)
library(ggsignif)
library(tidyverse)
anno_df <- compare_means(Expression ~ IMPC_Viability, method = "anova",group.by = "Development.Stage", data = anova_CLEANED_prim_viab_mouse_dmdd, 
                         p.adjust.method = "BH") %>%
  mutate(y_pos = 150, p.adj = format.pval(p.adj, digits = 2))

View(anno_df)
numb <- c(0.0000000000000002)
formatC(numb, format = "e", digits = 2)







#################################################################################
###############################################################################################






















####################~~~~~~~~~~~~~~~~~~~~~~~~~~~SECONDARY VIABILITY DATASET CONDUCTING PRINCIPLE COMPONENT ANALYSIS ~~~~~~~~~~~~~~~~~~~~~~~~##############

anova_second_viab_mouse_dmdd<-second_viab_mouse_dmdd%>%
  select(MGI.Accession.ID, Expression, Development.Stage, Window_Category)%>%
  distinct(Window_Category, Expression, Development.Stage, 
           MGI.Accession.ID)

View(anova_second_viab_mouse_dmdd)


anova_second_viab_mouse_dmdd$Development.Stage <- factor(anova_second_viab_mouse_dmdd$Development.Stage,
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

## DO ANOVA

res.anova_second_viab_mouse_dmdd <- aov(Expression ~ Window_Category, data = anova_second_viab_mouse_dmdd)


summary_<-summary(res.anova_second_viab_mouse_dmdd)
# 
# 
#                     Df    Sum Sq Mean Sq    F value Pr(>F)    
# Window_Category     2   4662145   2331072   118.1 <2e-16 ***
#   Residuals       10413 205515313   19736                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


### Multiple pair-wise comparison 
# 
(tuk.secondary<-TukeyHSD(res.anova_second_viab_mouse_dmdd))


format(tuk.secondary$Window_Category[,"p adj"], scientific = TRUE)

## PLOT THE TUKEY HSD RESULT TO FIND THE DIFFERENCES IN THE MEANS

with(par(mai=c(1,2.5,1,1)),{plot(tuk.secondary, las=1,cex.axis=1)})



# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = Expression ~ Window_Category, data = anova_second_viab_mouse_dmdd)
# 
# $Window_Category
#                            diff       lwr         upr     p adj
# Late-gestation-Early-gestation -39.39104 -46.36849 -32.4135759 0.0000000
# Mid-gestation-Early-gestation  -49.23404 -59.39400 -39.0740740 0.0000000
# Mid-gestation-Late-gestation    -9.84300 -20.34373   0.6577334 0.0716805 
# 
# 
# 

##There is a significant difference within the Window_catagory(p< 2e-16) , within the mp.descriptions (p< 2e-16)
# and within the Developmental Stages  (p< 2e-16)



###########################################################################################


###BOXPLOT TO SHOW 
library(ggpubr)

palette3 <- c("red","orange","blue")

anova_second_viab_mouse_dmdd$Window_Category <-factor(anova_second_viab_mouse_dmdd$Window_Category,
                                                      levels=c("Early-gestation",
                                                               "Mid-gestation",
                                                               "Late-gestation"))


## PLOT ANOVA RESULTS

plot_anova_DS.WINDOW <- ggplot(anova_second_viab_mouse_dmdd, aes(x =anova_second_viab_mouse_dmdd$Development.Stage,
            y = anova_second_viab_mouse_dmdd$Expression,
            fill = anova_second_viab_mouse_dmdd$Window_Category))%>%
  +labs(x= "Developmental Stage", y="Mean Gene Expression (TPM)", 
        title= "Mice Gene Expression (TPM) in Different Developmental Stages and the Window of Lethality ")%>%
  +scale_y_continuous(breaks = seq(0, 300, by = 20))%>%
  
  +scale_fill_manual(values=c("red", "orange", "blue"))  %>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10, angle = 45, hjust = 1, vjust = 1)) %>%
  +guides(fill=guide_legend(title="Window of Lethality: "))%>%
  +geom_boxplot(position=position_dodge(1),horizontal=TRUE,
                axes=FALSE,outlier.shape=NA)%>%
  +coord_cartesian(ylim = c(0, 300))%>%
  +geom_hline(yintercept = mean(anova_second_viab_mouse_dmdd$Expression ), 
             linetype = 2, color="green")+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 280, label = "p.signif")
#  stat_compare_means(method = "anova", label.y = 2500)


# 
# 
# library(cowplot)
# 
# save_plot(paste("./Plots/Mouse/EA_Developmental_Stages/Viability/plot_anova_DS.WINDOW.png"),
#           plot_anova_DS.WINDOW ,base_height= 10 ,base_aspect_ratio = 2) 

########################################################################################



## TO GET THE P-VALUES FOR EACH ***
library(ggpubr)
library(ggsignif)
library(tidyverse)
anno_df.sec <- compare_means(Expression ~ Window_Category ,
                         method = "anova",
                         group.by = "Development.Stage", 
                         data = anova_second_viab_mouse_dmdd, 
                         p.adjust.method = "BH") %>%
  mutate(y_pos = 400, p.adj = format.pval(p.adj, digits = 2))

View(anno_df.sec)
numb <- c(0.0000000000000002)
formatC(numb, format = "e", digits = 2)





#####################################################################################

## CLEAN DATA FOR PCA

pca_second_viab<-second_viab_mouse_dmdd%>%
  select_all()%>%
  distinct(MGI.Accession.ID , Expression, Development.Stage,
           Window_Category)


View(pca_second_viab)


library(reshape2)

pca_second_viab2<-dcast(pca_second_viab,
                        MGI.Accession.ID+ Window_Category ~ Development.Stage, value.var="Expression")

pca_second_viab2<-column_to_rownames(pca_second_viab2, var = "MGI.Accession.ID")


View(pca_second_viab2) #pca_try2





##replace NAs with 0
pca_second_viab2[is.na(pca_second_viab2)] <- 0


######################################################

##CONDUCT PCA
pca_3_secondary_viab<-pca_second_viab2

pca_3_secondary_viab$Window_Category  <-as.numeric(pca_3_secondary_viab$Window_Category)

res.pca_3_secondary_viab <- prcomp(pca_3_secondary_viab,  scale = TRUE)

attributes(res.pca_3_secondary_viab) #the attributes shows the information inside the object of class prcomp

summary_pca_3_secondary_viab <- summary(res.pca_3_secondary_viab)
str(summary_pca_3_secondary_viab)


######################################################################
# EXPLAINED VARIANCE 
#The amount of varaiance that can be explained by the PRINCIPLE COMPONENTS
#PC1 explains 95.376% of the VARIATION IN THE data, PCA2 explains 3.355% of the data

expVar_pca_second_viab  <- summary_pca_3_secondary_viab$importance[2,] * 100
expVar_pca_second_viab

##PLOT EXPLAINED VARIANCE

expVar_pca_second_viab_barplot<-barplot(expVar_pca_second_viab, xlab = "Principal Components", ylab = "Explained Variance (%)",
                                      col ="steelblue", ylim=c(0,100),
                                      main="The Explained Variance (%) of the Secondary Viability Dataset")

## ADD LABELS
text(expVar_pca_second_viab_barplot, expVar_pca_second_viab, 
     labels=round(expVar_pca_second_viab,1), xpd=TRUE, pos=3,) #this is the line of your interest in regards to the value labels


# Cumulative Variance
cumVar__pca_second_viab <- summary_pca_3_secondary_viab$importance[3,] * 100
cumVar__pca_second_viab

plot(cumVar__pca_second_viab, type="o" ,col="black", pch=21, bg="blue", cex=0.8, ylab="Cumulative Variance (%)",
     xlab="Principal Components", main="The Cumulative Variance (%) of the Secondary Viability Dataset")


##################################################################################

#DEFINE COLOURS FOR ALL VIABILITIES

# 
#pDat_pca_2_CLEANED_prim_viab <- pData(pca_2_CLEANED_prim_viab)

pClass_pca_second_viab2 <- as.factor(pca_second_viab2$Window_Category )
#names(pClass_pca_2_CLEANED_prim_viab) <- sampleNames(pca_2_CLEANED_prim_viab)
str(pClass_pca_second_viab2)


levels(pClass_pca_second_viab2)[levels(pClass_pca_second_viab2)=="Early-gestation"]
levels(pClass_pca_second_viab2)[levels(pClass_pca_second_viab2)=="Late-gestation"]
levels(pClass_pca_second_viab2)[levels(pClass_pca_second_viab2)=="Mid-gestation"]

table(pClass_pca_second_viab2)

EarlyIndx <- which((pClass_pca_second_viab2 == "Early-gestation") == TRUE)
LateIndx <- which((pClass_pca_second_viab2 == "Late-gestation") == TRUE)
MidIndx <- which((pClass_pca_second_viab2 == "Mid-gestation") == TRUE)

pca_4_second_viab2<-pca_second_viab2%>%
  select(-Window_Category)
colour2 <- NULL
colour2[EarlyIndx] <- "red"
colour2[LateIndx] <- "blue"
colour2[MidIndx] <- "orange"


### Plot

boxplot(pca_4_second_viab2, col=colour2, las=2)

legend("top", legend=c("Early-gestation","Late-gestation", 
                              "Mid-gestation"), 
       fill=c("red", "blue", "orange"),
       lty=1:2, cex=0.8,box.lty=0)




#################################################################
## PLOT THE XSCORES OF PC 1 AND PC2


Xscores_res.pca_3_secondary_viab <- res.pca_3_secondary_viab$x

View(Xscores_res.pca_3_secondary_viab)

plot(Xscores_res.pca_3_secondary_viab, xlab="PC1", 
     ylab="PC2", pch=21,bg=colour2, cex=1.5,cex.lab=0.7, cex.axis = 0.7,
     main="Principle Component 2 Score versus Principle Component 1 Score for the Secondary Viability Dataset")

legend(-48,2.0,legend=c("Early-gestation","Late-gestation", 
                             "Mid-gestation"), 
       fill=c("red", "blue", "orange"),
       lty=1:2, cex=0.8,
       box.lty=0, inset = 0.5)
################################################################################
# 
















###############################################################################



###################~~~~~~~~~~~~~~~~~~~~#PRIMARY VIABILTY: CONDUCT t-SNE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########



library(Rtsne)


TSNE_PRIMARY<-dcast(pca_CLEANED_prim_viab,
                               MGI.Accession.ID + IMPC_Viability  ~ Development.Stage, value.var="Expression")

View(TSNE_PRIMARY)


##replace NAs with 0
TSNE_PRIMARY[is.na(TSNE_PRIMARY)] <- 0


## MODIFIYING DATASET FOR BOTH TSNE AND PCA
LABELS<-TSNE_PRIMARY$IMPC_Viability

TSNE_PRIMARY$IMPC_Viability<-as.factor(TSNE_PRIMARY$IMPC_Viability)



TSNE_PRIMARY$IMPC_Viability <-factor(TSNE_PRIMARY$IMPC_Viability,
                                                      levels=c("Lethal",
                                                               "Subviable",
                                                               "Viable"))




## For Plotting

TSNE_COLOURS<-c("green", "purple", "yellow")

names(TSNE_COLOURS)= unique(TSNE_PRIMARY$IMPC_Viability)


View(TSNE_PRIMARY)

## PERFORM THE TEST
TSNE_TEST_PRIMARY <- Rtsne(TSNE_PRIMARY[,-2], dims = 2, 
              perplexity=30, verbose=TRUE, max_iter = 500, pca=TRUE)



exeTimeTsne<- system.time(Rtsne(TSNE_PRIMARY[,-2], dims = 2, 
                                perplexity=30, verbose=TRUE, max_iter = 500))



# ## PLOTTING
# # 
# plot(TSNE_TEST_PRIMARY$Y, t="n", main="t-SNE on Primary Viability Dataset", xlab="t-SNE 1", 
#       ylab="t-SNE 2", pch=18, cex=1.5,cex.lab=0.7, cex.axis = 0.7)
# text(TSNE_TEST_PRIMARY$Y, labels = TSNE_PRIMARY$IMPC_Viability ,
#       col=TSNE_COLOURS[TSNE_PRIMARY$IMPC_Viability ])
#      col=TSNE_COLOURS[TSNE_PRIMARY$IMPC_Viability])
# 
# legend("bottomright", legen=c("Lethal","Subviable", 
#                               "Viable"), fill=c("purple", "yellow", "green"),lty=1:2, cex=0.8,
#        box.lty=0)


library(ggplot2)
tsne_plot <- data.frame(x = TSNE_TEST_PRIMARY$Y[,1], y = TSNE_TEST_PRIMARY$Y[,2], 
                        col = TSNE_PRIMARY$IMPC_Viability)

tsne_prima.plot<-ggplot(tsne_plot,aes(x=x, y=y, 
                                     colour=TSNE_PRIMARY$IMPC_Viability))%>%
  + geom_point(aes(x=x, y=y, color=TSNE_PRIMARY$IMPC_Viability))%>%
  +labs(x= "t-SNE Dimension 1", y="t-SNE Dimension 2", 
        title= "t-SNE Analaysis for the Primary Viability Dataset")%>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10, angle = 45, hjust = 1, vjust = 1)) %>%
  +geom_point(size = 4)%>%
  + scale_color_manual(values =TSNE_COLOURS, aesthetics = "colour")%>%
  +guides(color=guide_legend(title="Viability: "))


# 
# # 
# library(cowplot)
# # # 
# save_plot(paste("./Plots/Mouse/EA_Developmental_Stages/Viability/TSNE RESULTS/tsne_prima.plot.png"),
#          tsne_prima.plot ,base_height= 10 ,base_aspect_ratio = 2) 
# # # 



###################~~~~~~~~~~~~~~~~~~~~#SECONDARY VIABILTY: CONDUCT t-SNE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########



library(Rtsne)


TSNE_SECONDARY<-dcast(pca_second_viab,
                          MGI.Accession.ID+ Window_Category ~ Development.Stage, value.var="Expression")
                    
View(TSNE_SECONDARY)

TSNE_SECONDARY$Window_Category <-factor(TSNE_SECONDARY$Window_Category,
                                        levels=c("Early-gestation",
                                                 "Mid-gestation",
                                                 "Late-gestation"))






##replace NAs with 0
TSNE_SECONDARY[is.na(TSNE_SECONDARY)] <- 0


## MODIFIYING DATASET FOR BOTH TSNE AND PCA
LABELS2<-TSNE_SECONDARY$Window_Category

TSNE_SECONDARY$Window_Category<-as.factor(TSNE_SECONDARY$Window_Category)

## For Plotting

TSNE2_COLOURS<-c("red", "blue","orange")

names(TSNE2_COLOURS)= unique(TSNE_SECONDARY$Window_Category)




TSNE_TEST_SECONDARY <- Rtsne(TSNE_SECONDARY[,-2], dims = 2, 
                           perplexity=30, verbose=TRUE, max_iter = 500, 
                           pca=TRUE)



exeTimeTsne_SECONDARY<- system.time(Rtsne(TSNE_SECONDARY[,-2], dims = 2, 
                                perplexity=30, verbose=TRUE, max_iter = 500))



## PLOTTING
# 
# plot(TSNE_TEST_SECONDARY$Y, t="n", main="t-SNE on Secondary Viability Dataset")
# text(TSNE_TEST_SECONDARY$Y, labels = TSNE_SECONDARY$Window_Category ,
#      col=TSNE2_COLOURS[TSNE_SECONDARY$Window_Category ])
# 
# 
# legend("topleft",legend=c("Early-gestation","Late-gestation", 
#                              "Mid-gestation"), 
#        fill=c("red", "green", "blue"),
#        lty=1:2, cex=0.8,box.lty=0, inset = 0.05)
# 


library(ggplot2)
tsne_sec.plot2 <- data.frame(x = TSNE_TEST_SECONDARY$Y[,1], y = TSNE_TEST_SECONDARY$Y[,2], 
                        col = TSNE_SECONDARY$Window_Category)

tsne_sec.plot<-ggplot(tsne_sec.plot2,aes(x=x, y=y, 
                                     colour=TSNE_SECONDARY$Window_Category))%>%
  + geom_point(aes(x=x, y=y, color=TSNE_SECONDARY$Window_Category))%>%
  +labs(x= "t-SNE Dimension 1", y="t-SNE Dimension 2", 
        title= "t-SNE Analaysis for the Secondary Viability Dataset")%>%
  +theme(axis.title=element_text(size=14,face="bold"),
         plot.caption=element_text(face = "italic", size=12, hjust = 0),
         text = element_text(size=10),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_rect(fill="lightyellow", size=0.5, linetype = "solid"),
         axis.text.y= element_text(size=10),
         axis.text.x= element_text(size=10, angle = 45, hjust = 1, vjust = 1)) %>%
  +geom_point(size = 4)%>%
  + scale_color_manual(values =TSNE2_COLOURS, aesthetics = "colour")%>%
  +guides(color=guide_legend(title="Window of Lethality: "))
  


# # 
# library(cowplot)
#  
save_plot(paste("./Plots/Mouse/EA_Developmental_Stages/Viability/TSNE RESULTS/tsne_sec.plot2.png"),
        tsne_sec.plot ,base_height= 10 ,base_aspect_ratio = 2) 
# # # 
