############################ Dataset preparation ##############################
### Libraries ####
library(readxl)
library(ggplot2)

### Functions####
merge_tables <- function(tables){
  data_table<-merge(x=tables[1],y=tables[2],by="ΚΩΔΙΚΟΣ.ΑΣΘΕΝΗ")
  for (i in 3:length(tables)){
    data_table<-merge(x=data_table,y=tables[i],by="ΚΩΔΙΚΟΣ.ΑΣΘΕΝΗ")
  }
  return(data_table)
}

columns_renaming<-function(old_columns,new_columns,dataset){
  for(i in 1:length(new_columns)){
    names(dataset)[names(dataset)==old_columns[i]] <- new_columns[i]
  }
  return(dataset)
}
### Parameters ###
path_data <- 'C:/Users/30697/Desktop/PhD/Data/GESS/PatientData_final.xlsx'
sheets <- readxl::excel_sheets(path_data)

### Main code here ####

### Patient, History, differential, entry, biochemical, syntax tables only
patient_table<-data.frame(read_excel(path_data,sheet=sheets[1]))
history_table<-data.frame(read_excel(path_data,sheet=sheets[2]))
differential_table<-data.frame(read_excel(path_data,sheet=sheets[3]))

entry_table<-data.frame(read_excel(path_data,sheet=sheets[4]))
biochemical_table<-data.frame(read_excel(path_data,sheet=sheets[5]))
blood_table<-data.frame(read_excel(path_data,sheet=sheets[6]))
syntax_table <-data.frame(read_excel(path_data,sheet="SYNTAX SCORE"))
subset_syntax<-subset(syntax_table,select=c("ΚΩΔΙΚΟΣ.ΑΣΘΕΝΗ","TOTAL.SS"))

##Merge table function
data_table <-merge_tables(list(patient_table,history_table,differential_table,entry_table,biochemical_table,
                               blood_table,subset_syntax))

##### Drop initial columns - uninformative columns
drop <-c("ΗΜΕΡΟΜΗΝΙΑ","...6","SYNTAX.SCORE.I","Οther","other.other","BMI...6","EXCLUDE.y","SYNTAX.SCORE.y","SYNTAX.SCORE","SYNTAX.SCORE.x")
data_init<-data_table[,!(names(data_table) %in% drop)]


### Rename columns
columns_to_rename<-c("X....FAMILY.HISTORY","BMI...5","TOTAL.SS","Death..1.Y..0.No.","ΚΑΡΔΙΑΚΗ.ΑΝΕΠΑΡΚΕΙΑ","Σοβαρού.βαθμού.στένωση.αορτικής.βαλβίδας","ΣΤΑΘΕΡΗ.ΣΤΗΘΑΓΧΗ","Unstable.Angina")
columns_new_names<-c("FAMILY.HISTORY","BMI","SYNTAX.SCORE","Death","HEART.FAILURE","SEVERE.AORTIC.STENOSIS","STABLE.ANGINA","UNSTABLE.ANGINA")
data_init <-columns_renaming(columns_to_rename,columns_new_names,data_init)

### New columns
new_columns <- c("AST","ALT")
data_init$AST <- data_init$SGOT
data_init$ALT <-data_init$SGPT

## Indexes
#### Index 1 ALT AST ratio####
## Formula: ALT/AST
data_init$index_ALT_AST_ratio <-data_init$ALT/data_init$AST





#### Index 2 HSI####
## Formula: 8*ALT/AST+BMI*(2*Diabetes Mellitus+ 2 if Female)    # Modification of formula in case of 1 Female and 0 Male
data_init$HSI <- 8 * (data_init$ALT / data_init$AST) + 
  data_init$BMI + (2 * data_init$DIABETES.MELLITUS + ifelse(data_init$Gender == 0, 2, 0))




### Index 3 TyG ####
## Formula: ln(TG*GLU/2)   ### Caution in case of fasting TG and fasting GLU is different from ours. 
data_init$TyG <- log(data_init$TG*data_init$GLU/2)




### Index 4 FSI ####
## Formula: 	FSI = -7.981+0.011*AGE-0.146*SEX[female (1) /male (0)] + 0.173*BMI +0.007*TG +0.593 *HYPERTENSION + 0.789 * DIABETES MELLITUS +1.1 * (SGPT/SGOT >= 1.33)
data_init$FSI<- -7.981 +   0.011 * data_init$AGE -  0.146 * ifelse(data_init$Gender==0,1,0) +  0.173 * data_init$BMI + 
  0.007 * data_init$TG + 
  0.593 * data_init$HYPERTENSION + 
  0.789 * data_init$DIABETES.MELLITUS + 
  1.1 * (data_init$ALT / data_init$AST >= 1.33)



### Index 5 acNASH ####
## Formula: acNASH = -0.031 × Age (years) – 2.352 × Ln(Creatinine (umol/L) # we have mg/dL # )  + 1.988 × Ln AST (U/L)
ratio<-88.4017
data_init$acNASH = -0.031*data_init$AGE -2.352*log(data_init$CREATININE * ratio)+ 1.988*log(data_init$AST)
data_init_acNash <-data_init[!is.na(data_init$acNASH),]



### Index 6 FIB-4 ####
## Formula: FIB-4=(Age (years) × AST (U/L)) / (Platelet Count (10^9) # (we have per 1000 or this is wrong) #  x √ALT (U/L))
data_init$FIB4=(data_init$AGE *data_init$AST) / ((data_init$PLT)*sqrt(data_init$ALT))


### Index 7 APRI ####
## Formula: APRI= [(AST (U/L)/upper limit of normal) / (Platelet Count (10^9/l)] × 100
AST_LIM<-40 # It is not provided in the paper ? I found a blog that proposes 40 https://www.webmd.com/hepatitis/what-is-apri-score
data_init$AST_LIM <- ifelse(data_init$Gender==1,40,32)
data_init$APRI=((data_init$AST/data_init$AST_LIM)/data_init$PLT)*100

data_init_APRI <-data_init[!is.na(data_init$APRI),]





################ data of interest ######################

cols_differential <- c('SYNTAX.SCORE','AGE','CHOL','TG','HYPERTENSION','ALT','AST','URIC.ACID','HDL','GLU','PLT','Gender','SMOKING',
                       'DIABETES.MELLITUS','AGE','BMI','GFR','LDL')

cols_indices <-c('index_ALT_AST_ratio','HSI','TyG','FSI','acNASH','FIB4','APRI')
cols_all<-cbind(c(cols_differential,cols_indices))
data_of_interest<-data_init[,(names(data_init) %in% cols_all)]

library(dplyr)


data_of_interest$SYNTAX.SCORE.BINARY<-as.factor(ifelse(data_of_interest$SYNTAX.SCORE>0,1,0))


save(data_of_interest,file='C:/Users/30697/Desktop/Project_Indices/Data/Data_NIT_with_selected_variables1.Rda')
