################## NEXT part of the code - Logistic regression analysis #################
#load('C:/Users/30697/Desktop/Project_Indices/280125/Clean_Data_NIT.Rda')
#libraries
library(stats)
#NULL model
load('C:/Users/30697/Desktop/Project_Indices/Data/Clean_Data_NIT.Rda')
data<-data_of_interest[!is.na(data_of_interest$APRI),]
data<-data[!is.na(data$FIB4),]
data<-data[!is.na(data$FSI),]
data<-data[!is.na(data$TyG),]
data<-data[!is.na(data$acNASH),]
data<-data[!is.na(data$GFR),]
data<-data[!is.na(data$LDL),]
#data<-data[!is.na(data$index_ALT_AST_ratio),]
#data<-data[!is.na(data$HSI),]
data<-data[!is.na(data$SYNTAX.SCORE.BINARY),]
data_of_interest<-data
data_of_interest$ratio<-data_of_interest$AST/data_of_interest$ALT
data_of_interest$target<-data_of_interest$SYNTAX.SCORE.BINARY
model_null <- glm(target~1,family=binomial,data=data_of_interest)

summary(model_null)

#### step 2. Model each NIT ####
nit_numerical<-c('APRI','FSI','TyG','acNASH','FIB4','HSI','ratio')
#nit_threshold<-c('APRI_threshold','FSI_threshold','TyG_threshold','acNASH_threshold','FIB4_threshold')
#numerical
important_nit<-c()
p_values<-c()
forms<-c('target ~ 1')
LRT_p<-c(NA)

for(nit in nit_numerical){
  data_nit<-data_of_interest[!is.na(data_of_interest[[nit]]),]
  model_null <- glm(target~1,family=binomial,data=data_nit)
  
  form <-as.formula(paste('target ~',nit))
  forms<-c(forms,paste('target ~',nit))
  model_nit <-glm(form,data=data_nit,family=binomial)
  
  res<- anova(model_null,model_nit,test="Chi")
  p_value<-res$`Pr(>Chi)`[[2]]
  LRT_p<-c(LRT_p,p_value)
  
    important_nit<-c(important_nit,nit)
    p_values<-c(p_values,p_value)
  
}


OR_values <- c()
CI_lower <- c()
CI_upper <- c()

for(nit in nit_numerical){
  data_nit <- data_of_interest[!is.na(data_of_interest[[nit]]),]
  model_null <- glm(target ~ 1, family = binomial, data = data_nit)
  
  form <- as.formula(paste('target ~', nit))
  forms <- c(forms, paste('target ~', nit))
  model_nit <- glm(form, data = data_nit, family = binomial)
  
  res <- anova(model_null, model_nit, test = "Chi")
  p_value <- res$`Pr(>Chi)`[[2]]
  LRT_p <- c(LRT_p, p_value)
  
  
    important_nit <- c(important_nit, nit)
    p_values <- c(p_values, p_value)
    
    # Extract OR and 95% CI
    OR <- exp(coef(model_nit)[2])  # Exponentiate coefficient to get OR
    CI <- exp(confint(model_nit))[2, ]  # Extract CI and exponentiate
    
    # Store results
    OR_values <- c(OR_values, OR)
    CI_lower <- c(CI_lower, CI[1])
    CI_upper <- c(CI_upper, CI[2])
  
}

# Create a summary dataframe
results <- data.frame(
  Variable = important_nit,
  P_Value = p_values,
  OR = OR_values,
  CI_Lower = CI_lower,
  CI_Upper = CI_upper
)

# Print results
print(results)
df_numerical<-data.frame(Model=forms,`LRT (p-value)`=LRT_p)
####### SAVE MODEL #################
save(df_numerical,file='C:/Users/30697/Desktop/Project_Indices/Working/Numerical_NIT_univariate.Rda')
########### Threshold #####################
#thresholds
important_nit_thresholds<-c()
p_values_thresholds<-c()
forms<-c('target ~ 1')
LRT_p<-c(NA)
for(nit in nit_threshold){
  data_nit<-data_of_interest[!is.na(data_of_interest[[nit]]),]
  model_null <- glm(target~1,family=binomial,data=data_nit)
  
  form <-as.formula(paste('target ~',nit))
  forms<-c(forms,paste('target ~',nit))
  
  
  model_nit <-glm(form,data=data_nit,family=binomial)
  
  res<- anova(model_null,model_nit,test="Chi")
  p_value<-res$`Pr(>Chi)`[[2]]
  LRT_p<-c(LRT_p,p_value)
  if(p_value < 0.05){
    important_nit_thresholds<-c(important_nit_thresholds,nit)
    p_values_thresholds<-c(p_values_thresholds,p_value)
    
  }
}

df_threshold<-data.frame(Model=forms,`LRT (p-value)`=LRT_p)

save(df_threshold,file='C:/Users/30697/Desktop/Project_Indices/Working/Threshold_NIT_univariate.Rda')


############ For each model --- add covariates ##############
#TyG equation involves TG and GLU
#APRI equation involves SGOT and PLT
#FSI equation involves AGE, Gender, BMI, TG, HYPERTENSION, DIABETES.MELLITUS, SGPT, SGOT
#FIB4 equation involves AGE, SGOT, PLT, SGPT
#acNASH equation involves AGE, CREATININE, SGOT
tyg_involved<-c('TG','GLU')
apri_involved<-c('SGOT','PLT')
fsi_involved<-c('AGE','Gender','BMI','TG','HYPERTENSION','DIABETES.MELLITUS','SGPT','SGOT')
fib4_involved<-c('AGE','SGOT','PLT','SGPT')
acNASH_involved<-c('AGE','CREATININE','SGOT')


covariates_all<-c('AGE','Gender','BMI','LDL','GFR','DIABETES.MELLITUS','SMOKING')

tyg_covariates<-setdiff(covariates_all,tyg_involved)
apri_covariates<-setdiff(covariates_all,apri_involved)
fsi_covariates<-setdiff(covariates_all,fsi_involved)
fib4_covariates<-setdiff(covariates_all,fib4_involved)
acNASH_covariates<-setdiff(covariates_all,acNASH_involved)


#important NITs 1. APRI, 2. FSI, 3. TyG, 4. acNASH, 5. FIB4
######### 1. APRI ###################
predictors <- apri_covariates
formulas<-c()
AIC_values<-c()
LRT_p<-c()
best_predictor<-NULL
final_predictors<-c()
dataset<-data_of_interest[!is.na(data_of_interest$APRI),]
predictors_numerical<-c('APRI')
while(length(predictors)>0){
  aic<-c()
  pred<-c()
  
  if(length(final_predictors)==0){
    string1<-paste(predictors_numerical,collapse=" + ")
    formulas<-c(formulas,paste('target ~ ',string1))
    init_model<-glm(as.formula(paste('target ~ ',string1)),data=dataset,family=binomial)
    AIC_values<-c(AIC_values,AIC(init_model))
    LRT_p<-c(LRT_p,NA)
    
  }else{
    string1<-paste(c(predictors_numerical,final_predictors),collapse=" + ")
    formulas<-c(formulas,paste('target ~',string1))
    
    init_model<-glm(as.formula(paste('target ~',string1)),data=dataset,family=binomial)
    
    AIC_values<-c(AIC_values,AIC(init_model))
    
    LRT_p<-c(LRT_p,NA)
  }
  for(predictor in predictors){
    
    data_pred<-dataset[!is.na(dataset[predictor]),]
    
    if(is.null(best_predictor)){
      string1<-paste(predictors_numerical,collapse=" + ")
      form<-as.formula(paste('target ~ ',string1))
      init_model<-glm(form,data=data_pred,family=binomial)
      
    }else{
      string1<-paste(c(predictors_numerical,final_predictors),collapse=" + ")
      
      form<-as.formula(paste('target ~',string1))
      
      init_model<-glm(form,data=data_pred,family=binomial)
    }
    if(is.null(best_predictor)){
      string1<-paste(c(predictors_numerical,predictor),collapse=" + ")
      form<-as.formula(paste('target ~',string1))
      
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      formulas<-c(formulas,paste('target ~',string1))
      
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
    }else{
      string1<-paste(c(predictors_numerical,final_predictors,predictor),collapse=" + ")
      
      form<-as.formula(paste('target ~',string1))
      
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      
      formulas<-c(formulas,paste('target ~',string1))
      
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
      
    }
    
    lrt_pred<-anova(init_model,model_pred,test="Chi")
    if(lrt_pred$`Pr(>Chi)`[[2]] < 0.05){
      pred<-c(pred,predictor)
      aic<-c(aic,AIC(model_pred))
    }
    
  }
  
  if(length(aic)!=0){
    best_predictor<-pred[which.min(aic)]
    
    final_predictors<-c(final_predictors,best_predictor)
    
    predictors<-setdiff(predictors,best_predictor)
    
    dataset<-dataset[!is.na(dataset[best_predictor]),]
  }else{
    predictors<-c()
  }
}
dataframe_numerical_with_covariates_apri<-data.frame(Model=formulas,`AIC values`=AIC_values,`LRT(p-values)`=LRT_p)

form<-as.formula('target ~ APRI + Gender + DIABETES.MELLITUS + GFR + SMOKING')
model_nit <-glm(form,data=data_pred,family=binomial)
summary(model_nit)

form <- as.formula('target ~ APRI + Gender + DIABETES.MELLITUS + GFR + SMOKING')
form <- as.formula('target ~ FSI + SMOKING + GFR')
form <- as.formula('target ~ TyG + Gender + GFR + SMOKING')

model_nit <- glm(form, data = data_pred, family = binomial)

# Extract Odds Ratios (OR), Confidence Intervals (95% CI), and p-values
OR <- exp(coef(model_nit))  # Odds Ratios
CI <- exp(confint(model_nit))  # 95% Confidence Intervals
p_values <- summary(model_nit)$coefficients[, 4]  # P-values

# Create a summary table
results <- data.frame(
  Variable = names(OR),
  OR = round(OR, 3),
  `95% CI` = paste0("(", round(CI[,1], 3), ", ", round(CI[,2], 3), ")"),
  `P-value` = round(p_values, 3)
)

# Print the table
print(results, row.names = FALSE)


save(dataframe_numerical_with_covariates_apri,file='C:/Users/30697/Desktop/Project_Indices/Supplementary/Numerical_NIT_LR_apri_with_covariates.Rda')


######### 2. FSI ###################
predictors <- fsi_covariates
formulas<-c()
AIC_values<-c()
LRT_p<-c()
best_predictor<-NULL
final_predictors<-c()
dataset<-data_of_interest[!is.na(data_of_interest$FSI),]
predictors_numerical<-c('FSI')
while(length(predictors)>0){
  aic<-c()
  pred<-c()
  
  if(length(final_predictors)==0){
    string1<-paste(predictors_numerical,collapse=" + ")
    formulas<-c(formulas,paste('target ~ ',string1))
    init_model<-glm(as.formula(paste('target ~ ',string1)),data=dataset,family=binomial)
    AIC_values<-c(AIC_values,AIC(init_model))
    LRT_p<-c(LRT_p,NA)
    
  }else{
    string1<-paste(c(predictors_numerical,final_predictors),collapse=" + ")
    formulas<-c(formulas,paste('target ~',string1))
    
    init_model<-glm(as.formula(paste('target ~',string1)),data=dataset,family=binomial)
    
    AIC_values<-c(AIC_values,AIC(init_model))
    
    LRT_p<-c(LRT_p,NA)
  }
  for(predictor in predictors){
    
    data_pred<-dataset[!is.na(dataset[predictor]),]
    
    if(is.null(best_predictor)){
      string1<-paste(predictors_numerical,collapse=" + ")
      form<-as.formula(paste('target ~ ',string1))
      init_model<-glm(form,data=data_pred,family=binomial)
      
    }else{
      string1<-paste(c(predictors_numerical,final_predictors),collapse=" + ")
      
      form<-as.formula(paste('target ~',string1))
      
      init_model<-glm(form,data=data_pred,family=binomial)
    }
    if(is.null(best_predictor)){
      string1<-paste(c(predictors_numerical,predictor),collapse=" + ")
      form<-as.formula(paste('target ~',string1))
      
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      formulas<-c(formulas,paste('target ~',string1))
      
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
    }else{
      string1<-paste(c(predictors_numerical,final_predictors,predictor),collapse=" + ")
      
      form<-as.formula(paste('target ~',string1))
      
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      
      formulas<-c(formulas,paste('target ~',string1))
      
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
      
    }
    
    lrt_pred<-anova(init_model,model_pred,test="Chi")
    if(lrt_pred$`Pr(>Chi)`[[2]] < 0.05){
      pred<-c(pred,predictor)
      aic<-c(aic,AIC(model_pred))
    }
    
  }
  
  if(length(aic)!=0){
    best_predictor<-pred[which.min(aic)]
    
    final_predictors<-c(final_predictors,best_predictor)
    
    predictors<-setdiff(predictors,best_predictor)
    
    dataset<-dataset[!is.na(dataset[best_predictor]),]
  }else{
    predictors<-c()
  }
}
dataframe_numerical_with_covariates_FSI<-data.frame(Model=formulas,`AIC values`=AIC_values,`LRT(p-values)`=LRT_p)

form <- as.formula('target ~ FSI + SMOKING + GFR')
model_nit <- glm(form, data = data_pred, family = binomial)

# Extract Odds Ratios (OR), Confidence Intervals (95% CI), and p-values
OR <- exp(coef(model_nit))  # Odds Ratios
CI <- exp(confint(model_nit))  # 95% Confidence Intervals
p_values <- summary(model_nit)$coefficients[, 4]  # P-values

# Create a summary table
results <- data.frame(
  Variable = names(OR),
  OR = round(OR, 2),
  `95% CI` = paste0("(", round(CI[,1], 2), ", ", round(CI[,2], 2), ")"),
  `P-value` = round(p_values, 3)
)

# Print the table
print(results, row.names = FALSE)



save(dataframe_numerical_with_covariates_FSI,file='C:/Users/30697/Desktop/Project_Indices/Supplementary/Numerical_NIT_LR_fsi_with_covariates.Rda')

######### 3. TyG ###################
predictors <- tyg_covariates
formulas<-c()
AIC_values<-c()
LRT_p<-c()
best_predictor<-NULL
final_predictors<-c()
dataset<-data_of_interest[!is.na(data_of_interest$TyG),]
predictors_numerical<-c('TyG')
while(length(predictors)>0){
  aic<-c()
  pred<-c()
  
  if(length(final_predictors)==0){
    string1<-paste(predictors_numerical,collapse=" + ")
    formulas<-c(formulas,paste('target ~ ',string1))
    init_model<-glm(as.formula(paste('target ~ ',string1)),data=dataset,family=binomial)
    AIC_values<-c(AIC_values,AIC(init_model))
    LRT_p<-c(LRT_p,NA)
    
  }else{
    string1<-paste(c(predictors_numerical,final_predictors),collapse=" + ")
    formulas<-c(formulas,paste('target ~',string1))
    
    init_model<-glm(as.formula(paste('target ~',string1)),data=dataset,family=binomial)
    
    AIC_values<-c(AIC_values,AIC(init_model))
    
    LRT_p<-c(LRT_p,NA)
  }
  for(predictor in predictors){
    
    data_pred<-dataset[!is.na(dataset[predictor]),]
    
    if(is.null(best_predictor)){
      string1<-paste(predictors_numerical,collapse=" + ")
      form<-as.formula(paste('target ~ ',string1))
      init_model<-glm(form,data=data_pred,family=binomial)
      
    }else{
      string1<-paste(c(predictors_numerical,final_predictors),collapse=" + ")
      
      form<-as.formula(paste('target ~',string1))
      
      init_model<-glm(form,data=data_pred,family=binomial)
    }
    if(is.null(best_predictor)){
      string1<-paste(c(predictors_numerical,predictor),collapse=" + ")
      form<-as.formula(paste('target ~',string1))
      
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      formulas<-c(formulas,paste('target ~',string1))
      
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
    }else{
      string1<-paste(c(predictors_numerical,final_predictors,predictor),collapse=" + ")
      
      form<-as.formula(paste('target ~',string1))
      
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      
      formulas<-c(formulas,paste('target ~',string1))
      
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
      
    }
    
    lrt_pred<-anova(init_model,model_pred,test="Chi")
    if(lrt_pred$`Pr(>Chi)`[[2]] < 0.05){
      pred<-c(pred,predictor)
      aic<-c(aic,AIC(model_pred))
    }
    
  }
  
  if(length(aic)!=0){
    best_predictor<-pred[which.min(aic)]
    
    final_predictors<-c(final_predictors,best_predictor)
    
    predictors<-setdiff(predictors,best_predictor)
    
    dataset<-dataset[!is.na(dataset[best_predictor]),]
  }else{
    predictors<-c()
  }
}
dataframe_numerical_with_covariates_TyG<-data.frame(Model=formulas,`AIC values`=AIC_values,`LRT(p-values)`=LRT_p)



form <- as.formula('target ~ TyG + Gender + GFR + SMOKING')
model_nit <- glm(form, data = data_pred, family = binomial)

# Extract Odds Ratios (OR), Confidence Intervals (95% CI), and p-values
OR <- exp(coef(model_nit))  # Odds Ratios
CI <- exp(confint(model_nit))  # 95% Confidence Intervals
p_values <- summary(model_nit)$coefficients[, 4]  # P-values

# Create a summary table
results <- data.frame(
  Variable = names(OR),
  OR = round(OR, 2),
  `95% CI` = paste0("(", round(CI[,1], 2), ", ", round(CI[,2], 2), ")"),
  `P-value` = round(p_values, 3)
)

# Print the table
print(results, row.names = FALSE)


save(dataframe_numerical_with_covariates_TyG,file='C:/Users/30697/Desktop/Project_Indices/Supplementary/Numerical_NIT_LR_tyg_with_covariates.Rda')




######### 4. acNASH ###################
predictors <- acNASH_covariates
formulas<-c()
AIC_values<-c()
LRT_p<-c()
best_predictor<-NULL
final_predictors<-c()
dataset<-data_of_interest[!is.na(data_of_interest$acNASH),]
predictors_numerical<-c('acNASH')
while(length(predictors)>0){
  aic<-c()
  pred<-c()
  
  if(length(final_predictors)==0){
    string1<-paste(predictors_numerical,collapse=" + ")
    formulas<-c(formulas,paste('target ~ ',string1))
    init_model<-glm(as.formula(paste('target ~ ',string1)),data=dataset,family=binomial)
    AIC_values<-c(AIC_values,AIC(init_model))
    LRT_p<-c(LRT_p,NA)
    
  }else{
    string1<-paste(c(predictors_numerical,final_predictors),collapse=" + ")
    formulas<-c(formulas,paste('target ~',string1))
    
    init_model<-glm(as.formula(paste('target ~',string1)),data=dataset,family=binomial)
    
    AIC_values<-c(AIC_values,AIC(init_model))
    
    LRT_p<-c(LRT_p,NA)
  }
  for(predictor in predictors){
    
    data_pred<-dataset[!is.na(dataset[predictor]),]
    
    if(is.null(best_predictor)){
      string1<-paste(predictors_numerical,collapse=" + ")
      form<-as.formula(paste('target ~ ',string1))
      init_model<-glm(form,data=data_pred,family=binomial)
      
    }else{
      string1<-paste(c(predictors_numerical,final_predictors),collapse=" + ")
      
      form<-as.formula(paste('target ~',string1))
      
      init_model<-glm(form,data=data_pred,family=binomial)
    }
    if(is.null(best_predictor)){
      string1<-paste(c(predictors_numerical,predictor),collapse=" + ")
      form<-as.formula(paste('target ~',string1))
      
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      formulas<-c(formulas,paste('target ~',string1))
      
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
    }else{
      string1<-paste(c(predictors_numerical,final_predictors,predictor),collapse=" + ")
      
      form<-as.formula(paste('target ~',string1))
      
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      
      formulas<-c(formulas,paste('target ~',string1))
      
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
      
    }
    
    lrt_pred<-anova(init_model,model_pred,test="Chi")
    if(lrt_pred$`Pr(>Chi)`[[2]] < 0.05){
      pred<-c(pred,predictor)
      aic<-c(aic,AIC(model_pred))
    }
    
  }
  
  if(length(aic)!=0){
    best_predictor<-pred[which.min(aic)]
    
    final_predictors<-c(final_predictors,best_predictor)
    
    predictors<-setdiff(predictors,best_predictor)
    
    dataset<-dataset[!is.na(dataset[best_predictor]),]
  }else{
    predictors<-c()
  }
}
dataframe_numerical_with_covariates_acNASH<-data.frame(Model=formulas,`AIC values`=AIC_values,`LRT(p-values)`=LRT_p)
save(dataframe_numerical_with_covariates_acNASH,file='C:/Users/30697/Desktop/Project_Indices/Working/Numerical_NIT_LR_acNASH_with_covariates.Rda')


######### 5. FIB4 ###################
predictors <- fib4_covariates
formulas<-c()
AIC_values<-c()
LRT_p<-c()
best_predictor<-NULL
final_predictors<-c()
dataset<-data_of_interest[!is.na(data_of_interest$FIB4),]
predictors_numerical<-c('FIB4')
while(length(predictors)>0){
  aic<-c()
  pred<-c()
  
  if(length(final_predictors)==0){
    string1<-paste(predictors_numerical,collapse=" + ")
    formulas<-c(formulas,paste('target ~ ',string1))
    init_model<-glm(as.formula(paste('target ~ ',string1)),data=dataset,family=binomial)
    AIC_values<-c(AIC_values,AIC(init_model))
    LRT_p<-c(LRT_p,NA)
    
  }else{
    string1<-paste(c(predictors_numerical,final_predictors),collapse=" + ")
    formulas<-c(formulas,paste('target ~',string1))
    
    init_model<-glm(as.formula(paste('target ~',string1)),data=dataset,family=binomial)
    
    AIC_values<-c(AIC_values,AIC(init_model))
    
    LRT_p<-c(LRT_p,NA)
  }
  for(predictor in predictors){
    
    data_pred<-dataset[!is.na(dataset[predictor]),]
    
    if(is.null(best_predictor)){
      string1<-paste(predictors_numerical,collapse=" + ")
      form<-as.formula(paste('target ~ ',string1))
      init_model<-glm(form,data=data_pred,family=binomial)
      
    }else{
      string1<-paste(c(predictors_numerical,final_predictors),collapse=" + ")
      
      form<-as.formula(paste('target ~',string1))
      
      init_model<-glm(form,data=data_pred,family=binomial)
    }
    if(is.null(best_predictor)){
      string1<-paste(c(predictors_numerical,predictor),collapse=" + ")
      form<-as.formula(paste('target ~',string1))
      
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      formulas<-c(formulas,paste('target ~',string1))
      
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
    }else{
      string1<-paste(c(predictors_numerical,final_predictors,predictor),collapse=" + ")
      
      form<-as.formula(paste('target ~',string1))
      
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      
      formulas<-c(formulas,paste('target ~',string1))
      
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
      
    }
    
    lrt_pred<-anova(init_model,model_pred,test="Chi")
    if(lrt_pred$`Pr(>Chi)`[[2]] < 0.05){
      pred<-c(pred,predictor)
      aic<-c(aic,AIC(model_pred))
    }
    
  }
  
  if(length(aic)!=0){
    best_predictor<-pred[which.min(aic)]
    
    final_predictors<-c(final_predictors,best_predictor)
    
    predictors<-setdiff(predictors,best_predictor)
    
    dataset<-dataset[!is.na(dataset[best_predictor]),]
  }else{
    predictors<-c()
  }
}
dataframe_numerical_with_covariates_FIB4<-data.frame(Model=formulas,`AIC values`=AIC_values,`LRT(p-values)`=LRT_p)
save(dataframe_numerical_with_covariates_FIB4,file='C:/Users/30697/Desktop/Project_Indices/Working/Numerical_NIT_LR_fib4_with_covariates.Rda')





























#important NIT threshold variables 1. FSI_threshold, 2. TyG threshold, 3. FIB4 threshold
########## 1. FSI threshold ############
predictors <- fsi_covariates
predictors_threshold<-c('FSI_threshold')
dataset<-data_of_interest[!is.na(data_of_interest$FSI_threshold),]
formulas<-c()
AIC_values<-c()
LRT_p<-c()
best_predictor<-NULL
final_predictors<-c()

while(length(predictors)>0){
  aic<-c()
  pred<-c()
  
  if(length(final_predictors)==0){
    string1<-paste(predictors_threshold,collapse=" + ")
    formulas<-c(formulas,paste('target ~ ',string1))
    init_model<-glm(as.formula(paste('target ~ ',string1)),data=dataset,family=binomial)
    AIC_values<-c(AIC_values,AIC(init_model))
    LRT_p<-c(LRT_p,NA)
    
  }else{
    string1<-paste(c(predictors_threshold,final_predictors),collapse=" + ")
    formulas<-c(formulas,paste('target ~',string1))
    
    init_model<-glm(as.formula(paste('target ~',string1)),data=dataset,family=binomial)
    
    AIC_values<-c(AIC_values,AIC(init_model))
    
    LRT_p<-c(LRT_p,NA)
  }
  for(predictor in predictors){
    
    data_pred<-dataset[!is.na(dataset[predictor]),]
    
    if(is.null(best_predictor)){
      string1<-paste(predictors_threshold,collapse=" + ")
      form<-as.formula(paste('target ~ ',string1))
      init_model<-glm(form,data=data_pred,family=binomial)
      
    }else{
      string1<-paste(c(predictors_threshold,final_predictors),collapse=" + ")
      
      form<-as.formula(paste('target ~',string1))
      
      init_model<-glm(form,data=data_pred,family=binomial)
    }
    if(is.null(best_predictor)){
      string1<-paste(c(predictors_threshold,predictor),collapse=" + ")
      form<-as.formula(paste('target ~',string1))
      
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      formulas<-c(formulas,paste('target ~',string1))
      
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
    }else{
      string1<-paste(c(predictors_threshold,final_predictors,predictor),collapse=" + ")
      
      form<-as.formula(paste('target ~',string1))
      
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      
      formulas<-c(formulas,paste('target ~',string1))
      
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
      
    }
    
    lrt_pred<-anova(init_model,model_pred,test="Chi")
    if(lrt_pred$`Pr(>Chi)`[[2]] < 0.05){
      pred<-c(pred,predictor)
      aic<-c(aic,AIC(model_pred))
    }
    
  }
  
  if(length(aic)!=0){
    best_predictor<-pred[which.min(aic)]
    
    final_predictors<-c(final_predictors,best_predictor)
    
    predictors<-setdiff(predictors,best_predictor)
    
    dataset<-dataset[!is.na(dataset[best_predictor]),]
  }else{
    predictors<-c()
  }
}

dataframe_threshold_with_covariates_fsi<-data.frame(Model=formulas,`AIC values`=AIC_values,`LRT(p-values)`=LRT_p)
save(dataframe_threshold_with_covariates_fsi,file='C:/Users/30697/Desktop/Project_Indices/Working/Threshold_NIT_LR_fsi_with_covariates.Rda')

######### 2. TyG Threshold ##############
predictors <- tyg_covariates
predictors_threshold<-c('TyG_threshold')
dataset<-data_of_interest[!is.na(data_of_interest$TyG_threshold),]
formulas<-c()
AIC_values<-c()
LRT_p<-c()
best_predictor<-NULL
final_predictors<-c()

while(length(predictors)>0){
  aic<-c()
  pred<-c()
  
  if(length(final_predictors)==0){
    string1<-paste(predictors_threshold,collapse=" + ")
    formulas<-c(formulas,paste('target ~ ',string1))
    init_model<-glm(as.formula(paste('target ~ ',string1)),data=dataset,family=binomial)
    AIC_values<-c(AIC_values,AIC(init_model))
    LRT_p<-c(LRT_p,NA)
    
  }else{
    string1<-paste(c(predictors_threshold,final_predictors),collapse=" + ")
    formulas<-c(formulas,paste('target ~',string1))
    
    init_model<-glm(as.formula(paste('target ~',string1)),data=dataset,family=binomial)
    
    AIC_values<-c(AIC_values,AIC(init_model))
    
    LRT_p<-c(LRT_p,NA)
  }
  for(predictor in predictors){
    
    data_pred<-dataset[!is.na(dataset[predictor]),]
    
    if(is.null(best_predictor)){
      string1<-paste(predictors_threshold,collapse=" + ")
      form<-as.formula(paste('target ~ ',string1))
      init_model<-glm(form,data=data_pred,family=binomial)
      
    }else{
      string1<-paste(c(predictors_threshold,final_predictors),collapse=" + ")
      
      form<-as.formula(paste('target ~',string1))
      
      init_model<-glm(form,data=data_pred,family=binomial)
    }
    if(is.null(best_predictor)){
      string1<-paste(c(predictors_threshold,predictor),collapse=" + ")
      form<-as.formula(paste('target ~',string1))
      
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      formulas<-c(formulas,paste('target ~',string1))
      
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
    }else{
      string1<-paste(c(predictors_threshold,final_predictors,predictor),collapse=" + ")
      
      form<-as.formula(paste('target ~',string1))
      
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      
      formulas<-c(formulas,paste('target ~',string1))
      
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
      
    }
    
    lrt_pred<-anova(init_model,model_pred,test="Chi")
    if(lrt_pred$`Pr(>Chi)`[[2]] < 0.05){
      pred<-c(pred,predictor)
      aic<-c(aic,AIC(model_pred))
    }
    
  }
  
  if(length(aic)!=0){
    best_predictor<-pred[which.min(aic)]
    
    final_predictors<-c(final_predictors,best_predictor)
    
    predictors<-setdiff(predictors,best_predictor)
    
    dataset<-dataset[!is.na(dataset[best_predictor]),]
  }else{
    predictors<-c()
  }
}

dataframe_threshold_with_covariates_tyg<-data.frame(Model=formulas,`AIC values`=AIC_values,`LRT(p-values)`=LRT_p)
save(dataframe_threshold_with_covariates_tyg,file='C:/Users/30697/Desktop/Project_Indices/Working/Threshold_NIT_LR_tyg_with_covariates.Rda')


######### 3. FIB4 threshold ###############
predictors <- fib4_covariates
predictors_threshold<-c('FIB4_threshold')
dataset<-data_of_interest[!is.na(data_of_interest$FIB4_threshold),]
formulas<-c()
AIC_values<-c()
LRT_p<-c()
best_predictor<-NULL
final_predictors<-c()

while(length(predictors)>0){
  aic<-c()
  pred<-c()
  
  if(length(final_predictors)==0){
    string1<-paste(predictors_threshold,collapse=" + ")
    formulas<-c(formulas,paste('target ~ ',string1))
    init_model<-glm(as.formula(paste('target ~ ',string1)),data=dataset,family=binomial)
    AIC_values<-c(AIC_values,AIC(init_model))
    LRT_p<-c(LRT_p,NA)
    
  }else{
    string1<-paste(c(predictors_threshold,final_predictors),collapse=" + ")
    formulas<-c(formulas,paste('target ~',string1))
    
    init_model<-glm(as.formula(paste('target ~',string1)),data=dataset,family=binomial)
    
    AIC_values<-c(AIC_values,AIC(init_model))
    
    LRT_p<-c(LRT_p,NA)
  }
  for(predictor in predictors){
    
    data_pred<-dataset[!is.na(dataset[predictor]),]
    
    if(is.null(best_predictor)){
      string1<-paste(predictors_threshold,collapse=" + ")
      form<-as.formula(paste('target ~ ',string1))
      init_model<-glm(form,data=data_pred,family=binomial)
      
    }else{
      string1<-paste(c(predictors_threshold,final_predictors),collapse=" + ")
      
      form<-as.formula(paste('target ~',string1))
      
      init_model<-glm(form,data=data_pred,family=binomial)
    }
    if(is.null(best_predictor)){
      string1<-paste(c(predictors_threshold,predictor),collapse=" + ")
      form<-as.formula(paste('target ~',string1))
      
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      formulas<-c(formulas,paste('target ~',string1))
      
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
    }else{
      string1<-paste(c(predictors_threshold,final_predictors,predictor),collapse=" + ")
      
      form<-as.formula(paste('target ~',string1))
      
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      
      formulas<-c(formulas,paste('target ~',string1))
      
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
      
    }
    
    lrt_pred<-anova(init_model,model_pred,test="Chi")
    if(lrt_pred$`Pr(>Chi)`[[2]] < 0.05){
      pred<-c(pred,predictor)
      aic<-c(aic,AIC(model_pred))
    }
    
  }
  
  if(length(aic)!=0){
    best_predictor<-pred[which.min(aic)]
    
    final_predictors<-c(final_predictors,best_predictor)
    
    predictors<-setdiff(predictors,best_predictor)
    
    dataset<-dataset[!is.na(dataset[best_predictor]),]
  }else{
    predictors<-c()
  }
}
dataframe_threshold_with_covariates_fib4<-data.frame(Model=formulas,`AIC values`=AIC_values,`LRT(p-values)`=LRT_p)
save(dataframe_threshold_with_covariates_fib4,file='C:/Users/30697/Desktop/Project_Indices/Working/Threshold_NIT_LR_fib4_with_covariates.Rda')


#### step 3. and 4. Sequentially Add Predictors #######

####### Numerical variables ############
predictors <- important_nit
selected_predictors <- c() # Track selected predictors
formulas<-c()
AIC_values<-c()
LRT_p<-c()
best_predictor<-NULL
dataset<-data_of_interest
final_predictors<-c()



while(length(predictors)>0){
  aic<-c()
  pred<-c()
  if(length(final_predictors)==0){
  formulas<-c(formulas,'target ~ 1')
  init_model<-glm(target ~ 1,data=dataset,family=binomial)
  AIC_values<-c(AIC_values,AIC(init_model))
  LRT_p<-c(LRT_p,NA)
  }else{
  string1<-paste(final_predictors,collapse=" + ")
  formulas<-c(formulas,paste('target ~',string1))
  init_model<-glm(as.formula(paste('target ~',string1)),data=dataset,family=binomial)
  AIC_values<-c(AIC_values,AIC(init_model))
  LRT_p<-c(LRT_p,NA)
  }
  for(predictor in predictors){
    data_pred<-dataset[!is.na(dataset[predictor]),]
    if(is.null(best_predictor)){
    init_model<-glm(target ~ 1 ,data=data_pred,family=binomial)
   
    }else{
      string1<-paste(final_predictors,collapse=" + ")
      form<-as.formula(paste('target ~',string1))
      init_model<-glm(form,data=data_pred,family=binomial)
    }
    if(is.null(best_predictor)){
    model_pred<-glm(as.formula(paste('target ~',predictor)),data=data_pred,family=binomial)
    
    formulas<-c(formulas,paste('target ~',predictor))
    AIC_values<-c(AIC_values,AIC(model_pred))
    lrt_pred<-anova(init_model,model_pred,test="Chi")
    LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
    }else{
      string1<-paste(c(final_predictors,predictor),collapse=" + ")
      form<-as.formula(paste('target ~',string1))
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      
      formulas<-c(formulas,paste('target ~',string1))
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
      
    }
    
    lrt_pred<-anova(init_model,model_pred,test="Chi")
    if(lrt_pred$`Pr(>Chi)`[[2]] < 0.05){
      pred<-c(pred,predictor)
      aic<-c(aic,AIC(model_pred))
    }
    
  }
  
  if(length(aic)!=0){
  best_predictor<-pred[which.min(aic)]
  
  final_predictors<-c(final_predictors,best_predictor)
  
  predictors<-setdiff(predictors,best_predictor)
  dataset<-dataset[!is.na(dataset[best_predictor]),]
  }else{
    predictors<-c()
  }
}


dataframe_numerical<-data.frame(Model=formulas,`AIC values`=AIC_values,`LRT(p-values)`=LRT_p)


### Final model 
dataset<-data_of_interest
dataset<-dataset[!is.na(dataset$APRI),]
dataset<-dataset[!is.na(dataset$TyG),]

final_model<- glm(target ~ TyG + APRI, data=dataset,family=binomial)
coeffs<-coef(final_model)
ORs<-exp(coeffs)
CI<-confint(final_model)
CI_OR<-exp(CI)
summary_table <- data.frame(
  Predictor = names(coeffs),
  Coefficient = round(coeffs, 3),
  OR = round(ORs, 3),
  `CI Lower` = round(CI_OR[, 1], 3),
  `CI Upper` = round(CI_OR[, 2], 3)
)

save(final_model,file='Numerical_NIT_mult_model_no_cov.Rda')
save(summary_table,file='Numerical_NIT_mult_model_summary_table.Rda')
#### Threshold variables ############
predictors <- important_nit_thresholds
selected_predictors <- c() # Track selected predictors
formulas<-c()
AIC_values<-c()
LRT_p<-c()
best_predictor<-NULL
dataset<-data_of_interest
final_predictors<-c()
names_coefs<-names(coeffs)[2:length(names(coeffs))]
while(length(predictors)>0){
  aic<-c()
  pred<-c()
  if(length(final_predictors)==0){
    formulas<-c(formulas,'target ~ 1')
    init_model<-glm(target ~ 1,data=dataset,family=binomial)
    AIC_values<-c(AIC_values,AIC(init_model))
    LRT_p<-c(LRT_p,NA)
  }else{
    string1<-paste(final_predictors,collapse=" + ")
    formulas<-c(formulas,paste('target ~',string1))
    init_model<-glm(as.formula(paste('target ~',string1)),data=dataset,family=binomial)
    AIC_values<-c(AIC_values,AIC(init_model))
    LRT_p<-c(LRT_p,NA)
  }
  for(predictor in predictors){
    data_pred<-dataset[!is.na(dataset[predictor]),]
    if(is.null(best_predictor)){
      init_model<-glm(target ~ 1 ,data=data_pred,family=binomial)
      
    }else{
      string1<-paste(final_predictors,collapse=" + ")
      form<-as.formula(paste('target ~',string1))
      init_model<-glm(form,data=data_pred,family=binomial)
    }
    if(is.null(best_predictor)){
      model_pred<-glm(as.formula(paste('target ~',predictor)),data=data_pred,family=binomial)
      
      formulas<-c(formulas,paste('target ~',predictor))
      AIC_values<-c(AIC_values,AIC(model_pred))
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
    }else{
      string1<-paste(c(final_predictors,predictor),collapse=" + ")
      form<-as.formula(paste('target ~',string1))
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      
      formulas<-c(formulas,paste('target ~',string1))
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
      
    }
    
    lrt_pred<-anova(init_model,model_pred,test="Chi")
    if(lrt_pred$`Pr(>Chi)`[[2]] < 0.05){
      pred<-c(pred,predictor)
      aic<-c(aic,AIC(model_pred))
    }
    
  }
  
  if(length(aic)!=0){
    best_predictor<-pred[which.min(aic)]
    
    final_predictors<-c(final_predictors,best_predictor)
    
    predictors<-setdiff(predictors,best_predictor)
    dataset<-dataset[!is.na(dataset[best_predictor]),]
  }else{
    predictors<-c()
  }
}
dataframe_thresholds<-data.frame(Model=formulas,`AIC values`=AIC_values,`LRT(p-values)`=LRT_p)

### Final model ###
dataset<-data_of_interest
dataset<-dataset[!is.na(dataset$FSI_threshold),]
dataset<-dataset[!is.na(dataset$FIB4_threshold),]
dataset<-dataset[!is.na(dataset$TyG_threshold),]

final_model<- glm(target ~ FSI_threshold + FIB4_threshold + TyG_threshold, data=dataset,family=binomial)
coeffs<-coef(final_model)
ORs<-exp(coeffs)
CI<-confint(final_model)
CI_OR<-exp(CI)
summary_table <- data.frame(
  Predictor = names(coeffs),
  Coefficient = round(coeffs, 3),
  OR = round(ORs, 3),
  `CI Lower` = round(CI_OR[, 1], 3),
  `CI Upper` = round(CI_OR[, 2], 3)
)
save(final_model,file='Threshold_NIT_mult_model_no_cov.Rda')
save(summary_table,file='Threshold_NIT_mult_model_summary_table.Rda')
save(dataframe_numerical,file='C:/Users/30697/Desktop/Project_Indices/Working/Numerical_NIT_mult_LR_without_covariates.Rda')
save(dataframe_thresholds,file='C:/Users/30697/Desktop/Project_Indices/Working/Threshold_NIT_LR_mult_without_covariates.Rda')










############### With covariates now ###################

load('C:/Users/30697/Desktop/Project_Indices/Working/Numerical_NIT_mult_LR_without_covariates.Rda')
load('C:/Users/30697/Desktop/Project_Indices/Working/Threshold_NIT_LR_mult_without_covariates.Rda')


predictors_numerical <- c('TyG','APRI')
predictors_threshold <-c('FSI_threshold','FIB4_threshold','TyG_threshold')






######## Numerical ############
load('C:/Users/30697/Desktop/Project_Indices/Working/Final_dataset.Rda')
dataset<-data_of_interest[!is.na(data_of_interest$TyG),]
dataset<-dataset[!is.na(dataset$APRI),]
covariates_numerical<-tyg_covariates[tyg_covariates %in% apri_covariates]
predictors <- covariates_numerical
formulas<-c()
AIC_values<-c()
LRT_p<-c()
best_predictor<-NULL
final_predictors<-c()

while(length(predictors)>0){
  aic<-c()
  pred<-c()
  
  if(length(final_predictors)==0){
    string1<-paste(predictors_numerical,collapse=" + ")
    formulas<-c(formulas,paste('target ~ ',string1))
    init_model<-glm(as.formula(paste('target ~ ',string1)),data=dataset,family=binomial)
    AIC_values<-c(AIC_values,AIC(init_model))
    LRT_p<-c(LRT_p,NA)
    
  }else{
    string1<-paste(c(predictors_numerical,final_predictors),collapse=" + ")
    formulas<-c(formulas,paste('target ~',string1))
    
    init_model<-glm(as.formula(paste('target ~',string1)),data=dataset,family=binomial)
    
    AIC_values<-c(AIC_values,AIC(init_model))
    
    LRT_p<-c(LRT_p,NA)
  }
  for(predictor in predictors){
    
    data_pred<-dataset[!is.na(dataset[predictor]),]
    
    if(is.null(best_predictor)){
      string1<-paste(predictors_numerical,collapse=" + ")
      form<-as.formula(paste('target ~ ',string1))
      init_model<-glm(form,data=data_pred,family=binomial)
      
    }else{
      string1<-paste(c(predictors_numerical,final_predictors),collapse=" + ")
      
      form<-as.formula(paste('target ~',string1))
      
      init_model<-glm(form,data=data_pred,family=binomial)
    }
    if(is.null(best_predictor)){
      string1<-paste(c(predictors_numerical,predictor),collapse=" + ")
      form<-as.formula(paste('target ~',string1))
      
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      formulas<-c(formulas,paste('target ~',string1))
      
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
    }else{
      string1<-paste(c(predictors_numerical,final_predictors,predictor),collapse=" + ")
      
      form<-as.formula(paste('target ~',string1))
      
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      
      formulas<-c(formulas,paste('target ~',string1))
      
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
      
    }
    
    lrt_pred<-anova(init_model,model_pred,test="Chi")
    if(lrt_pred$`Pr(>Chi)`[[2]] < 0.05){
      pred<-c(pred,predictor)
      aic<-c(aic,AIC(model_pred))
    }
    
  }
  
  if(length(aic)!=0){
    best_predictor<-pred[which.min(aic)]
    
    final_predictors<-c(final_predictors,best_predictor)
    
    predictors<-setdiff(predictors,best_predictor)
    
    dataset<-dataset[!is.na(dataset[best_predictor]),]
  }else{
    predictors<-c()
  }
}

dataframe_numerical_with_covariates<-data.frame(Model=formulas,`AIC values`=AIC_values,`LRT(p-values)`=LRT_p)

save(dataframe_numerical_with_covariates,file='C:/Users/30697/Desktop/Project_Indices/Working/Numerical_NIT_mult_LR_with_covariates.Rda')

### Final model 
dataset<-data_of_interest
dataset<-dataset[!is.na(dataset$TyG),]
dataset<-dataset[!is.na(dataset$APRI),]
dataset<-dataset[!is.na(dataset$ACS),]
dataset<-dataset[!is.na(dataset$Gender),]
dataset<-dataset[!is.na(dataset$GFR),]
dataset<-dataset[!is.na(dataset$DIABETES.MELLITUS),]


library(caret)
trainIndex<-createDataPartition(dataset$target,p=.67,list=FALSE,times=1)
train<-dataset[trainIndex,]
test<-dataset[-trainIndex,]

model_train<-glm(target ~ TyG + APRI + ACS + Gender + GFR + DIABETES.MELLITUS,data=train,family=binomial)

final_model<- glm(target ~ TyG + APRI + ACS + Gender + GFR + DIABETES.MELLITUS, data=dataset,family=binomial)
coeffs<-coef(final_model)
ORs<-exp(coeffs)
CI<-confint(final_model)
CI_OR<-exp(CI)
summary_table <- data.frame(
  Predictor = names(coeffs),
  Coefficient = round(coeffs, 3),
  OR = round(ORs, 3),
  `CI Lower` = round(CI_OR[, 1], 3),
  `CI Upper` = round(CI_OR[, 2], 3)
)

save(final_model,file='Numerical_NIT_model_with_cov.Rda')
save(summary_table,file='Numerical_NIT_model_with_cov_summary_table.Rda')





######## Threshold
load('C:/Users/30697/Desktop/Project_Indices/Working/Final_dataset.Rda')

dataset<-data_of_interest[!is.na(data_of_interest$FSI_threshold),]
dataset<-dataset[!is.na(dataset$FIB4_threshold),]
dataset<-dataset[!is.na(dataset$TyG_threshold),]

covs<-tyg_covariates[tyg_covariates %in% fib4_covariates]
covariates_thresholds<-covs[covs %in% fsi_covariates]

predictors <- covariates_thresholds
formulas<-c()
AIC_values<-c()
LRT_p<-c()
best_predictor<-NULL
final_predictors<-c()

while(length(predictors)>0){
  aic<-c()
  pred<-c()
  
  if(length(final_predictors)==0){
    string1<-paste(predictors_threshold,collapse=" + ")
    formulas<-c(formulas,paste('target ~ ',string1))
    init_model<-glm(as.formula(paste('target ~ ',string1)),data=dataset,family=binomial)
    AIC_values<-c(AIC_values,AIC(init_model))
    LRT_p<-c(LRT_p,NA)
    
  }else{
    string1<-paste(c(predictors_threshold,final_predictors),collapse=" + ")
    formulas<-c(formulas,paste('target ~',string1))
    
    init_model<-glm(as.formula(paste('target ~',string1)),data=dataset,family=binomial)
    
    AIC_values<-c(AIC_values,AIC(init_model))
    
    LRT_p<-c(LRT_p,NA)
  }
  for(predictor in predictors){
    
    data_pred<-dataset[!is.na(dataset[predictor]),]
    
    if(is.null(best_predictor)){
      string1<-paste(predictors_threshold,collapse=" + ")
      form<-as.formula(paste('target ~ ',string1))
      init_model<-glm(form,data=data_pred,family=binomial)
      
    }else{
      string1<-paste(c(predictors_threshold,final_predictors),collapse=" + ")
      
      form<-as.formula(paste('target ~',string1))
      
      init_model<-glm(form,data=data_pred,family=binomial)
    }
    if(is.null(best_predictor)){
      string1<-paste(c(predictors_threshold,predictor),collapse=" + ")
      form<-as.formula(paste('target ~',string1))
      
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      formulas<-c(formulas,paste('target ~',string1))
      
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
    }else{
      string1<-paste(c(predictors_threshold,final_predictors,predictor),collapse=" + ")
      
      form<-as.formula(paste('target ~',string1))
      
      model_pred<-glm(form,data=data_pred,family=binomial)
      
      
      formulas<-c(formulas,paste('target ~',string1))
      
      AIC_values<-c(AIC_values,AIC(model_pred))
      
      lrt_pred<-anova(init_model,model_pred,test="Chi")
      
      LRT_p<-c(LRT_p,lrt_pred$`Pr(>Chi)`[[2]])
      
    }
    
    lrt_pred<-anova(init_model,model_pred,test="Chi")
    if(lrt_pred$`Pr(>Chi)`[[2]] < 0.05){
      pred<-c(pred,predictor)
      aic<-c(aic,AIC(model_pred))
    }
    
  }
  
  if(length(aic)!=0){
    best_predictor<-pred[which.min(aic)]
    
    final_predictors<-c(final_predictors,best_predictor)
    
    predictors<-setdiff(predictors,best_predictor)
    
    dataset<-dataset[!is.na(dataset[best_predictor]),]
  }else{
    predictors<-c()
  }
}

dataframe_threshold_with_covariates<-data.frame(Model=formulas,`AIC values`=AIC_values,`LRT(p-values)`=LRT_p)

save(dataframe_threshold_with_covariates,file='C:/Users/30697/Desktop/Project_Indices/Working/Threshold_NIT_LR_mult_with_covariates.Rda')


################ MODEL #################
dataset<-data_of_interest
dataset<-dataset[!is.na(dataset$FSI_threshold),]
dataset<-dataset[!is.na(dataset$FIB4_threshold),]
dataset<-dataset[!is.na(dataset$TyG_threshold),]
dataset<-dataset[!is.na(dataset$ACS),]




final_model<- glm(target ~ FSI_threshold + FIB4_threshold + TyG_threshold + ACS, data=dataset,family=binomial)
coeffs<-coef(final_model)
ORs<-exp(coeffs)
CI<-confint(final_model)
CI_OR<-exp(CI)
summary_table <- data.frame(
  Predictor = names(coeffs),
  Coefficient = round(coeffs, 3),
  OR = round(ORs, 3),
  `CI Lower` = round(CI_OR[, 1], 3),
  `CI Upper` = round(CI_OR[, 2], 3)
)

save(final_model,file='Threshold_NIT_model_with_cov.Rda')
save(summary_table,file='Threshold_NIT_model_with_cov_summary_table.Rda')
