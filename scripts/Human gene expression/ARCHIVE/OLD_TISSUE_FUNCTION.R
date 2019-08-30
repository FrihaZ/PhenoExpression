


#################~~~~~~~~~~~~~~~~~~~~~~~~~~~~~THE TISSUE FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#################################



TISSUE_FISHER_TEST<-function(df,superclass, Tissue){  
  
  # For each threhsold the number of genes being expressed is counted 
  
  #change df to  long format using gather
  # keep gene id here
  pheno_cleaned <- df %>% 
    gather("GTEX.tissues","Expression", `Adipose - Subcutaneous`:`Whole Blood`)%>%
    select(GTEX.tissues,`HPO-Term-ID`,hpo.ancestors,hpo.description,GTEX.tissues,Expression)%>%
    distinct()
  #View(pheno_cleaned)
  
  
  
  pheno_cleaned_count<-pheno_cleaned%>%
    dplyr::group_by(GTEX.tissues,Expression) %>%  # group by gene here 
    dplyr::summarise(n=n()) 
  
  #View(pheno_cleaned_count)
  
  
  
  GTEX.Tissues.to.tissues<-GTEX.Tissues.to.HPO%>%
    select_all()%>%
    drop_na()%>%
    select(GTEX.tissues)
  
  ## There are some tissues that are found in different HPO.superclass.description tehrefore have been duplicated in the GTEX.Tissues.to.tissues
  
  # Merge the GTEX.Tissues.to.tissues and human_genes_TPM_greater0_pheno_count to get the duplicated tissues
  
  pheno_cleaned_Tissues<-merge(GTEX.Tissues.to.tissues, pheno_cleaned_count, 
                               by = c("GTEX.tissues", "GTEX.tissues"), all.x = TRUE)
  
  
  # The human_genes_TPM_greater0_pheno_count and GTEX.Tissues.to.HPO joind to get the HPO.superclass.description for each tissue
  
  pheno_cleaned_Tissues<-inner_join(GTEX.Tissues.to.HPO,pheno_cleaned_Tissues)
  pheno_cleaned_Tissues<-pheno_cleaned_Tissues%>%
    select_all()%>%
    drop_na()%>%    # the NAs are dropped
    distinct(GTEX.tissues,HPO.id, HPO.superclass.id,HPO.superclass.description,Expression,n ) # duplicates are removed that are not needed
  
  #View(pheno_cleaned_Tissues)
  
  
  # #View(human_genes_TPM_greater0_pheno_Tissues)
  
  pheno_cleaned_Tissues$Expression <- as.factor(pheno_cleaned_Tissues$Expression)
  pheno_cleaned_Tissues$GTEX.tissues <- as.factor(pheno_cleaned_Tissues$GTEX.tissues)
  pheno_cleaned_Tissues$HPO.superclass.description <- as.factor(pheno_cleaned_Tissues$HPO.superclass.description)
  
  pheno_cleaned_Tissues$HPO.id <- as.factor(pheno_cleaned_Tissues$HPO.id)
  pheno_cleaned_Tissues$HPO.superclass.id <- as.factor(pheno_cleaned_Tissues$HPO.superclass.id)
  
  
  ####~~~~~~~~~~~~~~~~~ To do a Fischer test for every Tissue vs Tissues not in the same System~~~~~~~
  
  
  
  
  greater0_pheno_fischer<-pheno_cleaned_Tissues%>%
    select(GTEX.tissues,HPO.superclass.description,Expression,n)%>%
    distinct()
  
  greater0_pheno_fischer_NO_TISSUE<-greater0_pheno_fischer%>%
    select(GTEX.tissues,HPO.superclass.description,Expression,n)%>%
    filter(HPO.superclass.description != superclass )%>%
    filter(GTEX.tissues != Tissue)%>%
    distinct(GTEX.tissues,Expression,n)
  
  
  
  #change them into factors
  greater0_pheno_fischer_NO_TISSUE$GTEX.tissues <- as.factor(greater0_pheno_fischer_NO_TISSUE$GTEX.tissues)
  greater0_pheno_fischer_NO_TISSUE$Expression <- as.factor(greater0_pheno_fischer_NO_TISSUE$Expression)
  
  
  
  
  greater0_pheno_fischer_tissue<-greater0_pheno_fischer%>%
    select(GTEX.tissues,HPO.superclass.description,Expression,n)%>%
    filter(HPO.superclass.description ==superclass)%>%
    filter(GTEX.tissues == Tissue )%>%
    distinct(GTEX.tissues,Expression,n)%>%
    spread(Expression,n)%>%
    distinct(GTEX.tissues,No, Yes)
  
  
  #All tissues without the Adipose_Subcutaneous
  
  All_tissues<-greater0_pheno_fischer_NO_TISSUE %>%
    spread(Expression,n)%>%
    distinct()
  
  ##View(all_tissues)
  
  
  # Sum all of the No and Yes rows to get total 
  
  All_tissues_<-All_tissues%>%
    select_all()%>%
    summarise_at(c("No","Yes"),sum) %>% #sums the columns
    distinct(No, Yes)%>%
    mutate(GTEX.tissues="Total")%>%
    distinct(GTEX.tissues,No,Yes)
  
  
  ##View(all_tissues_)
  
  
  # Bind the two datfarames to create a matrix
  Total_Tissues<-rbind(greater0_pheno_fischer_tissue,All_tissues_)
  
  
  #define the rownames 
  library(magrittr)
  Total_Tissues<-Total_Tissues %>%
    set_rownames(.$GTEX.tissues) %>% 
    select(-GTEX.tissues)
  
  
  #View(Total_Tissues)
  
  
  # perform Fisher test
  
  PVals<-fisher.test(Total_Tissues,alternative='two.sided')  # THIS WORKED
  
  PVALUES<-PVals$p.value
  
  return(PVALUES)
  
}


# #TEST FUNCTION
# TISSUE_FISHER_TEST(df=human_genes_TPM_greater0_pheno,
#                     superclass="Abnormality of the genitourinary system", 
#                     Tissue="Testis")

###########~~~~~~~~~~~~~~~~USE OF FOR-LOOP TO USE TISSUE_FISHER_TEST FUNCTION ~~~~~#################################

#THE FOR LOOP ITERATES THROUGH THE all_disease_genes_pheno_Tissues DATAFRAME WHICH HAS ALL TISSUES AND HPO.DESC(TOPLEVEL)

# CREATE A data frame with two columns and 55 rows, WHICH HAS THE TISSUE NAMES, HPO.SUPERCLASS AND AN EMPTY COLUMN


# NO GENE EXPRESSION 
PVALUE_Tissues_None <- data.frame(all_disease_genes_pheno_Tissues$GTEX.tissues,
                                  all_disease_genes_pheno_Tissues$HPO.superclass.description,
                                  P.VALUE=vector(length=55)) 


#View(human_genes_TPM_greater0_pheno)
for(i in 1:nrow(all_disease_genes_pheno_Tissues)) {
  ii <- all_disease_genes_pheno_Tissues[i,1]     #FIRST COLUMN
  jj <-all_disease_genes_pheno_Tissues[i,2]      #SECOND COLUMN
  
  print(ii)                                     # TO CHECK IF THE CORRECT ROW IS BEING TAKEN
  print(jj)
  
  
  ## THE OUTCOME OF THE FUNCTION IS SAVED AS PVAL WHICH IS THE FISHER TEST PVALUE 
  PVAL<-TISSUE_FISHER_TEST(df=human_genes_TPM_None_pheno, # THE FUNCTION NEEDS THE DF THAT HAS THE GENE GENE EXPRESSIONS AT SPECIFIC THRESHOLDS
                           
                           superclass= paste(jj), 
                           Tissue=paste(ii))
  
  print(PVAL)                  # THE PVALUE IS PRINTED
  
  PVALUE_Tissues_None$P.VALUE[i]<-paste(PVAL)              # THE PVALUE IS ADDED INTO THE DF CREATED ABOVE
  
}

# THE DATAFRAME IS #ViewED
#View(PVALUE_Tissues_None)



# > 0
PVALUE_Tissues_GREATER0 <- data.frame(all_disease_genes_pheno_Tissues$GTEX.tissues,
                                      all_disease_genes_pheno_Tissues$HPO.superclass.description,
                                      P.VALUE=vector(length=55)) 

for(i in 1:nrow(all_disease_genes_pheno_Tissues)) {
  ii <- all_disease_genes_pheno_Tissues[i,1]     #FIRST COLUMN
  jj <-all_disease_genes_pheno_Tissues[i,2]      #SECOND COLUMN
  
  print(ii)                                     # TO CHECK IF THE CORRECT ROW IS BEING TAKEN
  print(jj)
  
  
  
  ## THE OUTCOME OF THE FUNCTION IS SAVED AS PVAL WHICH IS THE FISHER TEST PVALUE 
  PVAL<-TISSUE_FISHER_TEST(df=human_genes_TPM_greater0_pheno, # THE FUNCTION NEEDS THE DF THAT HAS THE GENE GENE EXPRESSIONS AT SPECIFIC THRESHOLDS
                           superclass= paste(jj), 
                           Tissue=paste(ii))
  
  print(PVAL)                  # THE PVALUE IS PRINTED
  
  PVALUE_Tissues_GREATER0$P.VALUE[i]<-paste(PVAL)              # THE PVALUE IS ADDED INTO THE DF CREATED ABOVE
  
}

# THE DATAFRAME IS #ViewED
#View(PVALUE_Tissues_GREATER0)



# > 0.1
PVALUE_Tissues_0.1 <- data.frame(all_disease_genes_pheno_Tissues$GTEX.tissues,
                                 all_disease_genes_pheno_Tissues$HPO.superclass.description,
                                 P.VALUE=vector(length=55)) 

for(i in 1:nrow(all_disease_genes_pheno_Tissues)) {
  ii <- all_disease_genes_pheno_Tissues[i,1]     #FIRST COLUMN
  jj <-all_disease_genes_pheno_Tissues[i,2]      #SECOND COLUMN
  
  print(ii)                                     # TO CHECK IF THE CORRECT ROW IS BEING TAKEN
  print(jj)
  
  
  
  ## THE OUTCOME OF THE FUNCTION IS SAVED AS PVAL WHICH IS THE FISHER TEST PVALUE 
  PVAL<-TISSUE_FISHER_TEST(df=human_genes_TPM_0.1_pheno, # THE FUNCTION NEED THE DF THAT HAS TEH GENE GENE EXPRESSIONS AT SPECIFIC THRESHOLDS
                           superclass= paste(jj), 
                           Tissue=paste(ii))  # the default is two-sided fisher test if  you dont specify
  
  print(PVAL)                  # THE PVALUE IS PRINTED
  
  PVALUE_Tissues_0.1$P.VALUE[i]<-paste(PVAL)              # THE PVALUE IS ADDED INTO THE DF CREATED ABOVE
  
}

# THE DATAFRAME IS #ViewED
#View(PVALUE_Tissues_0.1)

# > 1
PVALUE_Tissues_1 <- data.frame(all_disease_genes_pheno_Tissues$GTEX.tissues,
                               all_disease_genes_pheno_Tissues$HPO.superclass.description,
                               P.VALUE=vector(length=55)) 

for(i in 1:nrow(all_disease_genes_pheno_Tissues)) {
  ii <- all_disease_genes_pheno_Tissues[i,1]     #FIRST COLUMN
  jj <-all_disease_genes_pheno_Tissues[i,2]      #SECOND COLUMN
  
  print(ii)                                     # TO CHECK IF THE CORRECT ROW IS BEING TAKEN
  print(jj)
  
  
  ## THE OUTCOME OF THE FUNCTION IS SAVED AS PVAL WHICH IS THE FISHER TEST PVALUE 
  PVAL<-TISSUE_FISHER_TEST(df=human_genes_TPM_1_pheno,# THE FUNCTION NEED THE DF THAT HAS TEH GENE GENE EXPRESSIONS AT SPECIFIC THRESHOLDS
                           superclass= paste(jj), 
                           Tissue=paste(ii))
  
  print(PVAL)                  # THE PVALUE IS PRINTED
  
  PVALUE_Tissues_1$P.VALUE[i]<-paste(PVAL)              # THE PVALUE IS ADDED INTO THE DF CREATED ABOVE
  
}

# THE DATAFRAME IS #ViewED
#View(PVALUE_Tissues_1)











