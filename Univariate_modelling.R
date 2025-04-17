############################ Logistic regression univariate ####################
library(knitr)
library(reshape2)
library(psych)
library(ggplot2)
library(ggExtra)
library(ggthemes)
library(dplyr)
library(arsenal)
library(kableExtra)


load('C:/Users/30697/Desktop/Project_Indices/Data/Data_NIT_with_selected_variables1.Rda')
cols_selected<-c('APRI','FSI','index_ALT_AST_ratio','HSI','TyG','acNASH','FIB4')
adjusted_iqr_outliers_df <- function(data,cols_selected, multiplier = 3) {
  cols_thresh<-c()
  for(col in cols_selected){
    if(col=="index_ALT_AST_ratio"){
      col_t<-"ALT_AST_ratio_threshold"
    }
    else{
      col_t<-paste0(col,'_threshold')}
    cols_thresh<-c(cols_thresh,col_t) 
    
  }
  clean_data <- data
  discard_data<-data
  for (num in 1:length(cols_selected)) {
    col<-cols_selected[num]
    colt<-cols_thresh[num]
    
    
    Q1 <- quantile(data[[col]], 0.25,na.rm =TRUE)
    Q3 <- quantile(data[[col]], 0.75,na.rm = TRUE)
    IQR_val <- IQR(data[[col]],na.rm=TRUE)
    
    lower_bound <- Q1 - (multiplier * IQR_val)
    upper_bound <- Q3 + (multiplier * IQR_val)
    
    # Filter out outliers for the column
    inds<-which(clean_data[[col]] <= lower_bound | clean_data[[col]] >= upper_bound)
    
    #Clean data - convert to NA
    clean_data[inds,c(col,colt)] <- NA
    discard_data[-inds,c(col,colt)]<-NA
  }
  return(list(clean_data,discard_data))
}
data1<-adjusted_iqr_outliers_df(data_of_interest,cols_selected,3)
data_of_interest <- data1[[1]]
save(data_of_interest,file='C:/Users/30697/Desktop/Project_Indices/Data/Clean_Data_NIT.Rda')

load('C:/Users/30697/Desktop/Project_Indices/Data/Clean_Data_NIT.Rda')

data<-data_of_interest[!is.na(data_of_interest$APRI),]
data<-data[!is.na(data$FIB4),]
data<-data[!is.na(data$FSI),]
data<-data[!is.na(data$TyG),]
data<-data[!is.na(data$acNASH),]
data<-data[!is.na(data$GFR),]
data<-data[!is.na(data$LDL),]
data<-data[!is.na(data$index_ALT_AST_ratio),]
data<-data[!is.na(data$HSI),]
data_of_interest<-data
data_of_interest$SYNTAX.SCORE.BINARY<-ifelse(data_of_interest$SYNTAX.SCORE>0,1,0)
data_of_interest$target<-data_of_interest$SYNTAX.SCORE.BINARY
data_of_interest<-data_of_interest[!is.na(data_of_interest$target),]

model_null <- glm(target~1,family=binomial,data=data_of_interest)

nit_numerical<-c('APRI','FSI','TyG','acNASH','FIB4','index_ALT_AST_ratio','HSI')



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



extract_logistic_results <- function(var) {
  formula <- as.formula(paste("SYNTAX.SCORE.BINARY ~", var))
  model <- glm(formula, data = data_of_interest, family = binomial)
  
  # Extract coefficients and confidence intervals
  coef_summary <- summary(model)$coefficients
  OR <- exp(coef_summary[2, 1])  # Odds ratio (exponent of the coefficient)
  p_value <- coef_summary[2, 4]  # P-value
  
  return(data.frame(Variable = var, OddsRatio = OR, PValue = p_value))
}

vars<-c('APRI','FSI','index_ALT_AST_ratio','HSI','TyG','acNASH','FIB4')
# Apply this extraction function to all variables
logistic_results <- do.call(rbind, lapply(vars, extract_logistic_results))

kable(logistic_results, caption = "Univariate Logistic Regression Results",
      col.names = c("Variable", "Odds Ratio", "P-value"), digits = 3)
