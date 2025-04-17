############## Modelling 10-02-2025 #################
library(caret)
library(stats)
library(PRROC)
library(precrec)
library(mltools)
library(precrec)
library(pROC)
library(PRROC)
library(ScottKnott)
library(arsenal)
library(kableExtra)
library(writexl)
library(xlsx)
set.seed(123)

load('C:/Users/30697/Desktop/Project_Indices/Data/Clean_Data_NIT.Rda')

data<-data_of_interest[!is.na(data_of_interest$APRI),]
data<-data[!is.na(data$FIB4),]
data<-data[!is.na(data$FSI),]
data<-data[!is.na(data$TyG),]
data<-data[!is.na(data$acNASH),]
data<-data[!is.na(data$GFR),]
data<-data[!is.na(data$LDL),]
data$SYNTAX.SCORE.BINARY<-ifelse(data$SYNTAX.SCORE>0,1,0)
data<-data[!is.na(data$SYNTAX.SCORE.BINARY),]
data$target<-data$SYNTAX.SCORE.BINARY
data_of_interest<-data
#####

df<-data_of_interest[c('AGE','CHOL','TG','HYPERTENSION','ALT','AST','URIC.ACID','HDL','GLU','PLT','Gender','SMOKING',
                       'DIABETES.MELLITUS','AGE','BMI','GFR','LDL','TyG','FSI','APRI','acNASH','FIB4','HSI','target')]
df$ratio<-data_of_interest$AST/data_of_interest$ALT
df$Gender<-as.factor(df$Gender)
df$DIABETES.MELLITUS<-as.factor(df$DIABETES.MELLITUS)
df$HYPERTENSION<-as.factor(df$HYPERTENSION)
df$SMOKING<-as.factor(df$SMOKING)
df$target<-as.factor(df$target)

#### normality results ########
numeric_cols <- sapply(df, is.numeric)  # Logical vector (TRUE for numeric columns)
df_numeric<-df[numeric_cols]
shapiro_results <- sapply(df_numeric, function(x) shapiro.test(x)$p.value)
shapiro_results
ks_results <- sapply(df_numeric, function(x) ks.test(x, "pnorm", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))$p.value)
ks_results



tab<-tableby(target ~., data=df,numeric.stats=c("mean","sd","median","q1q3"),numeric.test="wt",cat.test="chisq")
kable(summary(tab,text=TRUE, digits=2))


### Mann-Whitney U test
ps<-c()
for(col in colnames(df)[numeric_cols]){
  form<-as.formula(paste(col,'~','target'))
  res<-wilcox.test(form, data = df)
  ps<-c(ps,res$p.value)
}
####
results_formula<-c()
results_fold<-c()
results_recall<-c()
results_prec<-c()
results_acc<-c()
results_auc<-c()
results_roc<-c()
results_f1<-c()
results_balanced_accuracy<-c()
results_mcc<-c()
results_gmean<-c()
results_specificity<-c()
results_sensitivity<-c()
results_neg_pred_value<-c()
results_pos_pred_value<-c()

ls_for_roc<-list()
ls_for_pr<-list()
###### Create folds ###########
set.seed(123)  # For reproducibility
folds <- createFolds(data_of_interest$target, k = 10, list = TRUE)  # Creating indices for 10-fold CV
# Store train-test indices
cv_indices <- lapply(folds, function(test_idx) {
  train_idx <- setdiff(1:754, test_idx)
  list(train = train_idx, test = test_idx)
})

############ 1 predictor models without covariates ---numerical #####################
############ 1. APRI ###########################
formula<-'target ~ APRI'
form<-as.formula(formula)
fold<-1
dataset<-data_of_interest

# Storage for ROC/PR Curves (for Strategy 1)
roc_curves <- list()
pr_curves <- list()

# Storage for pooled predictions (for Strategy 2)
predicted_probabilities<-c()
test_cases<-factor(c(),levels=c(0,1))

probs_class_0<-c()
probs_class_1<-c()


for(cv in cv_indices){
  
  train<-dataset[cv$train,]
  test<-dataset[cv$test,]
  model_train<-glm(form,data=train,family=binomial)
  predictions <-data.frame('Probability'=stats::predict(model_train,test,type="response"))
  predictions$Class<-ifelse(predictions$Probability>0.5,1,0)
  actual<-test$target
  
  cm<- caret::confusionMatrix(factor(predictions$Class),factor(test$target),positive="1",mode=c("prec_recall"))
  predictions_performance<-data.frame(Results=c(cm$overall[1],
                                                cm$byClass[1:11]))
  
  
  
  
  predicted_probabilities<-c(predicted_probabilities,predictions$Probability)
  test_cases<-c(test_cases,actual)
  probs_class_0<-c(probs_class_0,predictions$Probability[actual == 1])
  probs_class_1<-c(probs_class_1,predictions$Probability[actual == 0])
  
  roc_obj <- roc(actual, predictions$Probability)
  # Compute PR-AUC using PRROC
  pr_curve <- pr.curve(scores.class0 = predictions$Probability[actual == 1],
                       scores.class1 = predictions$Probability[actual == 0], curve = TRUE)
  PR_auc <- pr_curve$auc.integral  # Extract PR-AUC
  roc_curves[[fold]] <- roc_obj
  pr_curves[[fold]] <- pr_curve
  ROC_auc <- as.numeric(auc(roc_obj))
  
  
  results_formula<-c(results_formula,formula)
  results_fold<-c(results_fold,fold)
  results_recall<-c(results_recall,predictions_performance['Recall',])
  results_prec<-c(results_prec,predictions_performance['Precision',])
  results_acc<-c(results_acc,predictions_performance['Accuracy',])
  results_auc<-c(results_auc,PR_auc)
  results_roc<-c(results_roc,ROC_auc)
  results_f1<-c(results_f1,predictions_performance['F1',])
  mcc_res<-mltools::mcc(factor(predictions$Class),factor(test$target))
  results_balanced_accuracy<-c(results_balanced_accuracy,predictions_performance['Balanced Accuracy',])
  results_mcc<-c(results_mcc,mcc_res)
  results_sensitivity<-c(results_sensitivity,predictions_performance['Sensitivity',])
  results_specificity<-c(results_specificity,predictions_performance['Specificity',])
  results_gmean<-c(results_gmean,sqrt(predictions_performance['Sensitivity',]*predictions_performance['Specificity',]))
  results_neg_pred_value<-c(results_neg_pred_value,ifelse(is.na(predictions_performance['Neg Pred Value',]),0,predictions_performance['Neg Pred Value',]))
  results_pos_pred_value<-c(results_pos_pred_value,predictions_performance['Pos Pred Value',])
  fold<-fold+1
  #all_preds <- c(all_preds, predictions$Probability)
  #all_labels <- c(all_labels, factor(actual))
  
}

ls_for_roc$`APRI_unadjusted`<-list(Probabilites=predicted_probabilities,Class=test_cases)
ls_for_pr$`APRI_unadjusted `<-list(Probabilities0=probs_class_0,Probabilites1=probs_class_1)


################# 2. TyG #####################
formula<-'target ~ TyG'
form<-as.formula(formula)
fold<-1
dataset<-data_of_interest

# Storage for ROC/PR Curves (for Strategy 1)
roc_curves <- list()
pr_curves <- list()

# Storage for pooled predictions (for Strategy 2)
predicted_probabilities<-c()
test_cases<-factor(c(),levels=c(0,1))

probs_class_0<-c()
probs_class_1<-c()


for(cv in cv_indices){
  
  train<-dataset[cv$train,]
  test<-dataset[cv$test,]
  model_train<-glm(form,data=train,family=binomial)
  predictions <-data.frame('Probability'=stats::predict(model_train,test,type="response"))
  predictions$Class<-ifelse(predictions$Probability>0.5,1,0)
  actual<-test$target
  
  cm<- caret::confusionMatrix(factor(predictions$Class),factor(test$target),positive="1",mode=c("prec_recall"))
  predictions_performance<-data.frame(Results=c(cm$overall[1],
                                                cm$byClass[1:11]))
  
  
  
  
  predicted_probabilities<-c(predicted_probabilities,predictions$Probability)
  test_cases<-c(test_cases,actual)
  probs_class_0<-c(probs_class_0,predictions$Probability[actual == 1])
  probs_class_1<-c(probs_class_1,predictions$Probability[actual == 0])
  
  roc_obj <- roc(actual, predictions$Probability)
  # Compute PR-AUC using PRROC
  pr_curve <- pr.curve(scores.class0 = predictions$Probability[actual == 1],
                       scores.class1 = predictions$Probability[actual == 0], curve = TRUE)
  PR_auc <- pr_curve$auc.integral  # Extract PR-AUC
  roc_curves[[fold]] <- roc_obj
  pr_curves[[fold]] <- pr_curve
  ROC_auc <- as.numeric(auc(roc_obj))
  
  
  results_formula<-c(results_formula,formula)
  results_fold<-c(results_fold,fold)
  results_recall<-c(results_recall,predictions_performance['Recall',])
  results_prec<-c(results_prec,predictions_performance['Precision',])
  results_acc<-c(results_acc,predictions_performance['Accuracy',])
  results_auc<-c(results_auc,PR_auc)
  results_roc<-c(results_roc,ROC_auc)
  results_f1<-c(results_f1,predictions_performance['F1',])
  mcc_res<-mltools::mcc(factor(predictions$Class),factor(test$target))
  results_balanced_accuracy<-c(results_balanced_accuracy,predictions_performance['Balanced Accuracy',])
  results_mcc<-c(results_mcc,mcc_res)
  results_sensitivity<-c(results_sensitivity,predictions_performance['Sensitivity',])
  results_specificity<-c(results_specificity,predictions_performance['Specificity',])
  results_gmean<-c(results_gmean,sqrt(predictions_performance['Sensitivity',]*predictions_performance['Specificity',]))
  results_neg_pred_value<-c(results_neg_pred_value,predictions_performance['Neg Pred Value',])
  results_pos_pred_value<-c(results_pos_pred_value,predictions_performance['Pos Pred Value',])
  fold<-fold+1
  #all_preds <- c(all_preds, predictions$Probability)
  #all_labels <- c(all_labels, factor(actual))
  
}

ls_for_roc$`TyG_unadjusted`<-list(Probabilites=predicted_probabilities,Class=test_cases)
ls_for_pr$`TyG_unadjusted `<-list(Probabilities0=probs_class_0,Probabilites1=probs_class_1)



################# 3. FSI #####################
formula<-'target ~ FSI'
form<-as.formula(formula)
fold<-1
dataset<-data_of_interest

# Storage for ROC/PR Curves (for Strategy 1)
roc_curves <- list()
pr_curves <- list()

# Storage for pooled predictions (for Strategy 2)
predicted_probabilities<-c()
test_cases<-factor(c(),levels=c(0,1))

probs_class_0<-c()
probs_class_1<-c()


for(cv in cv_indices){
  
  train<-dataset[cv$train,]
  test<-dataset[cv$test,]
  model_train<-glm(form,data=train,family=binomial)
  predictions <-data.frame('Probability'=stats::predict(model_train,test,type="response"))
  predictions$Class<-ifelse(predictions$Probability>0.5,1,0)
  actual<-test$target
  
  cm<- caret::confusionMatrix(factor(predictions$Class),factor(test$target),positive="1",mode=c("prec_recall"))
  predictions_performance<-data.frame(Results=c(cm$overall[1],
                                                cm$byClass[1:11]))
  
  
  
  
  predicted_probabilities<-c(predicted_probabilities,predictions$Probability)
  test_cases<-c(test_cases,actual)
  probs_class_0<-c(probs_class_0,predictions$Probability[actual == 1])
  probs_class_1<-c(probs_class_1,predictions$Probability[actual == 0])
  
  roc_obj <- roc(actual, predictions$Probability)
  # Compute PR-AUC using PRROC
  pr_curve <- pr.curve(scores.class0 = predictions$Probability[actual == 1],
                       scores.class1 = predictions$Probability[actual == 0], curve = TRUE)
  PR_auc <- pr_curve$auc.integral  # Extract PR-AUC
  roc_curves[[fold]] <- roc_obj
  pr_curves[[fold]] <- pr_curve
  ROC_auc <- as.numeric(auc(roc_obj))
  
  
  results_formula<-c(results_formula,formula)
  results_fold<-c(results_fold,fold)
  results_recall<-c(results_recall,predictions_performance['Recall',])
  results_prec<-c(results_prec,predictions_performance['Precision',])
  results_acc<-c(results_acc,predictions_performance['Accuracy',])
  results_auc<-c(results_auc,PR_auc)
  results_roc<-c(results_roc,ROC_auc)
  results_f1<-c(results_f1,predictions_performance['F1',])
  mcc_res<-mltools::mcc(factor(predictions$Class),factor(test$target))
  results_balanced_accuracy<-c(results_balanced_accuracy,predictions_performance['Balanced Accuracy',])
  results_mcc<-c(results_mcc,mcc_res)
  results_sensitivity<-c(results_sensitivity,predictions_performance['Sensitivity',])
  results_specificity<-c(results_specificity,predictions_performance['Specificity',])
  results_gmean<-c(results_gmean,sqrt(predictions_performance['Sensitivity',]*predictions_performance['Specificity',]))
  results_neg_pred_value<-c(results_neg_pred_value,ifelse(is.na(predictions_performance['Neg Pred Value',]),0,predictions_performance['Neg Pred Value',]))
  results_pos_pred_value<-c(results_pos_pred_value,predictions_performance['Pos Pred Value',])
  fold<-fold+1
  #all_preds <- c(all_preds, predictions$Probability)
  #all_labels <- c(all_labels, factor(actual))
  
}

ls_for_roc$`FSI_unadjusted`<-list(Probabilites=predicted_probabilities,Class=test_cases)
ls_for_pr$`FSI_unadjusted `<-list(Probabilities0=probs_class_0,Probabilites1=probs_class_1)


####### 1 predictor models with covariates -- numerical ###########
##### 1. APRI ###########
formula<-'target ~ APRI + Gender + DIABETES.MELLITUS + GFR + SMOKING '
form<-as.formula(formula)
fold<-1
dataset<-data_of_interest

# Storage for ROC/PR Curves (for Strategy 1)
roc_curves <- list()
pr_curves <- list()

# Storage for pooled predictions (for Strategy 2)
predicted_probabilities<-c()
test_cases<-factor(c(),levels=c(0,1))

probs_class_0<-c()
probs_class_1<-c()


for(cv in cv_indices){
  
  train<-dataset[cv$train,]
  test<-dataset[cv$test,]
  model_train<-glm(form,data=train,family=binomial)
  predictions <-data.frame('Probability'=stats::predict(model_train,test,type="response"))
  predictions$Class<-ifelse(predictions$Probability>0.5,1,0)
  actual<-test$target
  
  cm<- caret::confusionMatrix(factor(predictions$Class),factor(test$target),positive="1",mode=c("prec_recall"))
  predictions_performance<-data.frame(Results=c(cm$overall[1],
                                                cm$byClass[1:11]))
  
  
  

  predicted_probabilities<-c(predicted_probabilities,predictions$Probability)
  test_cases<-c(test_cases,actual)
  probs_class_0<-c(probs_class_0,predictions$Probability[actual == 1])
  probs_class_1<-c(probs_class_1,predictions$Probability[actual == 0])
  
  roc_obj <- roc(actual, predictions$Probability)
  # Compute PR-AUC using PRROC
  pr_curve <- pr.curve(scores.class0 = predictions$Probability[actual == 1],
                         scores.class1 = predictions$Probability[actual == 0], curve = TRUE)
  PR_auc <- pr_curve$auc.integral  # Extract PR-AUC
  roc_curves[[fold]] <- roc_obj
  pr_curves[[fold]] <- pr_curve
  ROC_auc <- as.numeric(auc(roc_obj))

  
  results_formula<-c(results_formula,formula)
  results_fold<-c(results_fold,fold)
  results_recall<-c(results_recall,predictions_performance['Recall',])
  results_prec<-c(results_prec,predictions_performance['Precision',])
  results_acc<-c(results_acc,predictions_performance['Accuracy',])
  results_auc<-c(results_auc,PR_auc)
  results_roc<-c(results_roc,ROC_auc)
  results_f1<-c(results_f1,predictions_performance['F1',])
  mcc_res<-mltools::mcc(factor(predictions$Class),factor(test$target))
  results_balanced_accuracy<-c(results_balanced_accuracy,predictions_performance['Balanced Accuracy',])
  results_mcc<-c(results_mcc,mcc_res)
  results_sensitivity<-c(results_sensitivity,predictions_performance['Sensitivity',])
  results_specificity<-c(results_specificity,predictions_performance['Specificity',])
  results_gmean<-c(results_gmean,sqrt(predictions_performance['Sensitivity',]*predictions_performance['Specificity',]))
  results_neg_pred_value<-c(results_neg_pred_value,predictions_performance['Neg Pred Value',])
  results_pos_pred_value<-c(results_pos_pred_value,predictions_performance['Pos Pred Value',])
  fold<-fold+1
  #all_preds <- c(all_preds, predictions$Probability)
  #all_labels <- c(all_labels, factor(actual))
  
}

ls_for_roc$`APRI`<-list(Probabilites=predicted_probabilities,Class=test_cases)
ls_for_pr$`APRI`<-list(Probabilities0=probs_class_0,Probabilites1=probs_class_1)
##### 3. FSI ###########
formula<-'target ~ FSI + SMOKING + GFR'
form<-as.formula(formula)
fold<-1
dataset<-data_of_interest

# Storage for ROC/PR Curves (for Strategy 1)
roc_curves <- list()
pr_curves <- list()

# Storage for pooled predictions (for Strategy 2)
predicted_probabilities<-c()
test_cases<-factor(c(),levels=c(0,1))

probs_class_0<-c()
probs_class_1<-c()


for(cv in cv_indices){
  
  train<-dataset[cv$train,]
  test<-dataset[cv$test,]
  model_train<-glm(form,data=train,family=binomial)
  predictions <-data.frame('Probability'=stats::predict(model_train,test,type="response"))
  predictions$Class<-ifelse(predictions$Probability>0.5,1,0)
  actual<-test$target
  
  cm<- caret::confusionMatrix(factor(predictions$Class),factor(test$target),positive="1",mode=c("prec_recall"))
  predictions_performance<-data.frame(Results=c(cm$overall[1],
                                                cm$byClass[1:11]))
  
  
  
  
  predicted_probabilities<-c(predicted_probabilities,predictions$Probability)
  test_cases<-c(test_cases,actual)
  probs_class_0<-c(probs_class_0,predictions$Probability[actual == 1])
  probs_class_1<-c(probs_class_1,predictions$Probability[actual == 0])
  
  roc_obj <- roc(actual, predictions$Probability)
  # Compute PR-AUC using PRROC
  pr_curve <- pr.curve(scores.class0 = predictions$Probability[actual == 1],
                       scores.class1 = predictions$Probability[actual == 0], curve = TRUE)
  PR_auc <- pr_curve$auc.integral  # Extract PR-AUC
  roc_curves[[fold]] <- roc_obj
  pr_curves[[fold]] <- pr_curve
  ROC_auc <- as.numeric(auc(roc_obj))
  
  
  results_formula<-c(results_formula,formula)
  results_fold<-c(results_fold,fold)
  results_recall<-c(results_recall,predictions_performance['Recall',])
  results_prec<-c(results_prec,predictions_performance['Precision',])
  results_acc<-c(results_acc,predictions_performance['Accuracy',])
  results_auc<-c(results_auc,PR_auc)
  results_roc<-c(results_roc,ROC_auc)
  results_f1<-c(results_f1,predictions_performance['F1',])
  mcc_res<-mltools::mcc(factor(predictions$Class),factor(test$target))
  results_balanced_accuracy<-c(results_balanced_accuracy,predictions_performance['Balanced Accuracy',])
  results_mcc<-c(results_mcc,mcc_res)
  results_sensitivity<-c(results_sensitivity,predictions_performance['Sensitivity',])
  results_specificity<-c(results_specificity,predictions_performance['Specificity',])
  results_gmean<-c(results_gmean,sqrt(predictions_performance['Sensitivity',]*predictions_performance['Specificity',]))
  results_neg_pred_value<-c(results_neg_pred_value,predictions_performance['Neg Pred Value',])
  results_pos_pred_value<-c(results_pos_pred_value,predictions_performance['Pos Pred Value',])
  fold<-fold+1
  #all_preds <- c(all_preds, predictions$Probability)
  #all_labels <- c(all_labels, factor(actual))
  
}

ls_for_roc$`FSI`<-list(Probabilites=predicted_probabilities,Class=test_cases)
ls_for_pr$`FSI`<-list(Probabilities0=probs_class_0,Probabilites1=probs_class_1)

##### 2. TyG ###########
formula<-'target ~ TyG + Gender + GFR + SMOKING '
form<-as.formula(formula)
fold<-1
dataset<-data_of_interest

# Storage for ROC/PR Curves (for Strategy 1)
roc_curves <- list()
pr_curves <- list()

# Storage for pooled predictions (for Strategy 2)
predicted_probabilities<-c()
test_cases<-factor(c(),levels=c(0,1))

probs_class_0<-c()
probs_class_1<-c()


for(cv in cv_indices){
  
  train<-dataset[cv$train,]
  test<-dataset[cv$test,]
  model_train<-glm(form,data=train,family=binomial)
  predictions <-data.frame('Probability'=stats::predict(model_train,test,type="response"))
  predictions$Class<-ifelse(predictions$Probability>0.5,1,0)
  actual<-test$target
  
  cm<- caret::confusionMatrix(factor(predictions$Class),factor(test$target),positive="1",mode=c("prec_recall"))
  predictions_performance<-data.frame(Results=c(cm$overall[1],
                                                cm$byClass[1:11]))
  
  
  
  
  predicted_probabilities<-c(predicted_probabilities,predictions$Probability)
  test_cases<-c(test_cases,actual)
  probs_class_0<-c(probs_class_0,predictions$Probability[actual == 1])
  probs_class_1<-c(probs_class_1,predictions$Probability[actual == 0])
  
  roc_obj <- roc(actual, predictions$Probability)
  # Compute PR-AUC using PRROC
  pr_curve <- pr.curve(scores.class0 = predictions$Probability[actual == 1],
                       scores.class1 = predictions$Probability[actual == 0], curve = TRUE)
  PR_auc <- pr_curve$auc.integral  # Extract PR-AUC
  roc_curves[[fold]] <- roc_obj
  pr_curves[[fold]] <- pr_curve
  ROC_auc <- as.numeric(auc(roc_obj))
  
  
  results_formula<-c(results_formula,formula)
  results_fold<-c(results_fold,fold)
  results_recall<-c(results_recall,predictions_performance['Recall',])
  results_prec<-c(results_prec,predictions_performance['Precision',])
  results_acc<-c(results_acc,predictions_performance['Accuracy',])
  results_auc<-c(results_auc,PR_auc)
  results_roc<-c(results_roc,ROC_auc)
  results_f1<-c(results_f1,predictions_performance['F1',])
  mcc_res<-mltools::mcc(factor(predictions$Class),factor(test$target))
  results_balanced_accuracy<-c(results_balanced_accuracy,predictions_performance['Balanced Accuracy',])
  results_mcc<-c(results_mcc,mcc_res)
  results_sensitivity<-c(results_sensitivity,predictions_performance['Sensitivity',])
  results_specificity<-c(results_specificity,predictions_performance['Specificity',])
  results_gmean<-c(results_gmean,sqrt(predictions_performance['Sensitivity',]*predictions_performance['Specificity',]))
  results_neg_pred_value<-c(results_neg_pred_value,predictions_performance['Neg Pred Value',])
  results_pos_pred_value<-c(results_pos_pred_value,predictions_performance['Pos Pred Value',])
  fold<-fold+1
  #all_preds <- c(all_preds, predictions$Probability)
  #all_labels <- c(all_labels, factor(actual))
  
}

ls_for_roc$`TyG`<-list(Probabilites=predicted_probabilities,Class=test_cases)
ls_for_pr$`TyG`<-list(Probabilities0=probs_class_0,Probabilites1=probs_class_1)




################### RESULTS #########################
res<-data.frame(Model=results_formula,Fold=results_fold,Recall=results_recall,Precision=results_prec,Accuracy=results_acc,F1=results_f1,AUCROC=results_roc,AUCPR=results_auc,
                Balanced_Accuarcy=results_balanced_accuracy,MCC=results_mcc,Gmean=results_gmean,Specificity=results_specificity,Sensitivity=results_sensitivity,
                Negative_Pred_Value=results_neg_pred_value,Positive_Pred_Value=results_pos_pred_value)
save(res,file='C:/Users/30697/Desktop/Project_Indices/Paper_results/Modelling_results_cv_full_metrics_all.Rda')











############################################################# NUMERICAL #####################################
################# Scott-Knott Results#################
metrics<-names(res)[3:15]
Comparison1 <- res
Comparison1$Model <- factor(Comparison1$Model)
Comparison1$Fold <-factor(Comparison1$Fold)
### Scott Knott  ####
scott_knott_ls<-list()
for(metric in metrics){
  formula<-paste(metric,'~ Model + Fold')
  form<-as.formula(formula)
  anv <- aov(form, data = Comparison1)
  anova(anv)
  
  # Perform the Scott-Knott clustering test
  sk_result <- SK(anv, which = "Model")
  
  # Print results
  print(sk_result)
  
  
  SKResultsStatistics.Algorithm <- data.frame(sk_result$info)
  
  Model <- SKResultsStatistics.Algorithm [,1]
  Measure <- SKResultsStatistics.Algorithm [,2]
  Lower <- SKResultsStatistics.Algorithm[,4]
  Upper <- SKResultsStatistics.Algorithm[,5]
  
  ClusterMatrix <- sk_result$out$Result
  ClusterMatrix <- ClusterMatrix[-1]
  ClusterMatrix <- cbind(ClusterMatrix, Group = apply(ClusterMatrix, 1, 
                                                      function(x) paste0(na.omit(x), collapse = "")))
  
  Groups.SK.Algorithm <- data.frame(Model,Measure,Lower,Upper,ClusterMatrix$Group)
  
  colnames(Groups.SK.Algorithm)[5] <- c("Cluster")
  scott_knott_ls[[metric]]<-Groups.SK.Algorithm
  }

metric<-metrics[1]
sk<-scott_knott_ls[[metric]]
sk_res<-sk[c('Model','Cluster')]
sk_res<-sk_res[order(sk_res$Model),]
colnames(sk_res)<-c('Model',metric)
for(metric in metrics[2:13]){
  sk<-scott_knott_ls[[metric]]
  cluster<-sk[order(sk$Model),]$Cluster
  sk_res[metric]<-cluster
}
save(sk_res,file='C:/Users/30697/Desktop/Project_Indices/Paper_results/SK_all_results_numerical.Rda')

################# Metrics Summary ####################
SummaryTables <- tableby(Model ~ ., data=Comparison1[,-2],
                         total=F,
                         control=tableby.control(numeric.stats=c(
                           'meansd',
                           "median",
                           "range")))

ResultsRQ1 <- summary(SummaryTables, text=T)
ResultsRQ1 <- data.frame(ResultsRQ1)
colnames(ResultsRQ1)[2:7] <- c(levels(Comparison1$Model))

save(ResultsRQ1,file='C:/Users/30697/Desktop/Project_Indices/Paper_results/Summary_Table_Scores_numerical.Rda')

################# Delong Tests ######################
names(ls_for_roc) <- c('APRI unadjusted','TyG unadjuted','FSI unadjusted','APRI adjusted','FSI adjusted','TyG adjusted')
ls_for_roc_numerical<-ls_for_roc

list_results_delong_numerical<-list()
names_all<-names(ls_for_roc_numerical)
for(name in names_all){
  roc1<-roc(ls_for_roc[[name]]$Class,ls_for_roc[[name]]$Probabilites)
  rest_names <-setdiff(names_all,name)
  for(second_name in rest_names){
    roc2<-roc(ls_for_roc[[second_name]]$Class,ls_for_roc[[second_name]]$Probabilites)
    res<-roc.test(roc1,roc2,method="delong",boot.n=1000)
    test<-paste(name,second_name,sep="~")
    list_results_delong_numerical[[test]]=res
  }
  names_all<-setdiff(names_all,name)
}

############# How can I find the threshold based on the predicted? Is there any way?


auc_1<-c()
auc_2<-c()
z<-c()
confint_1<-c()
confint_2<-c()
p_value<-c()
comparison_names<-c()
for(name in names(list_results_delong_numerical)){
  comp<-list_results_delong_numerical[[name]]
  comparison_names<-c(comparison_names,name)
  auc_1<-c(auc_1,as.numeric(comp$estimate[1]))
  auc_2<-c(auc_2,as.numeric(comp$estimate[2]))
  z<-c(z,as.numeric(comp$statistic))
  confint_1<-c(confint_1,comp$conf.int[1])
  confint_2<-c(confint_2,comp$conf.int[2])
  p_value<-c(p_value,comp$p.value)
}
names_all<-c('ROC Comparison','AUC (1)','AUC (2)','Statistic (Z)','P-Value','Lower CI','Upper CI')
res_delong_numerical<-data.frame(`ROC Comparison`=comparison_names,`AUC(1)`=auc_1,`AUC(2)`=auc_2,`Statistic (Z)`=z,`P-Value`=p_value,
                       `Lower CI`=confint_1,`Upper CI`=confint_2)
colnames(res_delong_numerical)<-names_all

res_delong_numerical$`P-Value`<-p.adjust(res_delong_numerical$`P-Value`,method='fdr')

res_delong_numerical[2:7]<-round(res_delong_numerical[2:7],3)

save(res_delong_numerical,file='C:/Users/30697/Desktop/Project_Indices/Paper_results/Delong_comparison_results_numerical.Rda')
write.xlsx(res_delong_numerical,file='C:/Users/30697/Desktop/Project_Indices/Paper_results/Delong_comparison_results_numerical.xlsx')

################# PR Tests #############
source('C:/Users/30697/Desktop/Project_Indices/Code/pr_test.R')
names(ls_for_roc) <- c('APRI unadjusted','TyG unadjuted','FSI unadjusted','APRI adjusted','FSI adjusted','TyG adjusted')
class<-ls_for_roc$`APRI unadjusted`$Class

list_results_pr_test<-list()
names_all<-names(ls_for_roc)


for(name in names_all){
  prob1<-ls_for_roc[[name]]$Probabilites
  rest_names <-setdiff(names_all,name)
  for(second_name in rest_names){
    prob2<-ls_for_roc[[second_name]]$Probabilites
    res<-pr.test(class,prob1,prob2,boot.n=1000)
    test<-paste(name,second_name,sep="~")
    list_results_pr_test[[test]]=res
  }
  names_all<-setdiff(names_all,name)
}

############################ Bootstrap #######################################


########## ROC bootstrap
names_all<-names(ls_for_roc)
library(pROC)
list_cis<-list()
for(name in names_all){
  prob<-ls_for_roc[[name]]$Probabilites
  class<-ls_for_roc[[name]]$Class
  roc<-pROC::roc(class,prob)
  list_cis[[name]]<-ci.auc(roc,boot.n=1000)
}

###### PR bootstrap
list_cis_pr<-list()
names_all<-names(ls_for_roc)
calculation_auc_ci<-function(pr_curve){
  x<-pr_curve$recall
  y1<-pr_curve$precision_low
  y2<-pr_curve$precision_high
  auc_low<-c()
  auc_high<-c()
  for(i in 1:length(x)-1){
    area_low_i <-(x[i+1]-x[i])*y1[i+1]
    area_high_i<-(x[i+1]-x[i])*y2[i+1]
    
    auc_low<-c(auc_low,area_low_i)
    auc_high<-c(auc_high,area_high_i)
  }
  auc_low<-round(sum(auc_low),3)
  auc_high<-round(sum(auc_high),3)
  
  res_ls<-paste(auc_low,auc_high,sep='-')
  }
for (name in names_all){
  prob<-ls_for_roc[[name]]$Probabilites
  class<-ls_for_roc[[name]]$Class
  pr_tbl = pr.boot(class, prob, boot.n = 1000, x_bins = 100,alpha=0.05) # default x_bin is 1000
  list_cis_pr[[name]]<-calculation_auc_ci(pr_tbl)
}

min_auc<-c()
auc_c<-c()
max_auc<-c()
names_all<-names(ls_for_roc)

# draw PR curve + add the bootstrap percentile confidence bands
library(ggplot2)
library(tibble)
for(name in names_all){
  prob<-ls_for_roc[[name]]$Probabilites
  class<-ls_for_roc[[name]]$Class
  pr_tbl = pr.boot(class, prob, boot.n = 100, x_bins = 100,alpha=0.05) # default x_bin is 1000
  pr_tbl |>
    ggplot(aes(x = recall, y = precision)) +
    geom_step() +
    ylim(c(0,1)) +
    geom_ribbon(aes(ymin = precision_low, ymax = precision_high), alpha = 0.2)
  pr_auc <- pr.curve(scores.class0 = pr_tbl$recall, weights.class0 = pr_tbl$precision, curve = TRUE)$auc.integral
  auc_c<-c(auc_c,pr_auc)
  pr_auc_min<-pr.curve(scores.class0=pr_tbl$recall,weights.class0 = pr_tbl$precision_low,curve=TRUE)$auc.integral
  min_auc<-c(min_auc,pr_auc_min)
  pr_auc_max<-pr.curve(scores.class0=pr_tbl$recall,weights.class0 = pr_tbl$precision_high,curve=TRUE)$auc.integral
  max_auc<-c(max_auc,pr_auc_max)
  
}

df<-data.frame(`Model`=names_all,`AUC`=auc_c,`Lower CI`=min_auc,`Upper CI`=max_auc)
df[2:4]<-round(df[2:4],3)

pr_auc <- pr.curve(scores.class0 = pr_tbl$recall, weights.class0 = pr_tbl$precision, curve = TRUE)

# Print the AUC value
print(pr_auc$auc.integral)
pr_auc_min<-pr.curve(scores.class0=pr_tbl$recall,weights.class0 = pr_tbl$precision_low,curve=TRUE)
print(pr_auc_min$auc.integral)
###################################################### ##################################################

auc_1<-c()
auc_2<-c()
z<-c()
confint_1<-c()
confint_2<-c()
p_value<-c()
comparison_names<-c()
for(name in names(list_results_pr_test)){
  comp<-list_results_pr_test[[name]]
  comparison_names<-c(comparison_names,name)
  auc_1<-c(auc_1,as.numeric(comp$auc1))
  auc_2<-c(auc_2,as.numeric(comp$auc2))
  #z<-c(z,as.numeric(comp$statistic))
  #confint_1<-c(confint_1,comp$conf.int[1])
  #confint_2<-c(confint_2,comp$conf.int[2])
  p_value<-c(p_value,comp$p.value)
}
names_all<-c('ROC Comparison','AUC (1)','AUC (2)','Statistic (Z)','P-Value','Lower CI','Upper CI')
names_all<-c('ROC Comparison','AUC (1)','AUC (2)','P-Value')

res_pr<-data.frame(`ROC Comparison`=comparison_names,`AUC(1)`=auc_1,`AUC(2)`=auc_2,`P-Value`=p_value)
                                 
colnames(res_pr)<-names_all
res_pr$`P-Value`<-p.adjust(res_pr$`P-Value`,method='fdr')
res_pr[2:4]<-round(res_pr[2:4],3)
save(res_pr,file='C:/Users/30697/Desktop/Project_Indices/Paper_results/PR_test_results.Rda')
write.xlsx(res_pr,file='C:/Users/30697/Desktop/Project_Indices/Paper_results/PR_test_results.xlsx')

################# ROC visualization ##################
library(dplyr)
roc_data <- lapply(names(ls_for_roc_numerical), function(name) {
  roc_obj <- roc(ls_for_roc[[name]]$Class, ls_for_roc[[name]]$Probabilites)
  auc_value <- round(auc(roc_obj), 3)  # Round AUC to 3 decimals
  data.frame(
    FPR = 1 - roc_obj$specificities,  # False Positive Rate
    TPR = roc_obj$sensitivities,      # True Positive Rate
    Model = paste0(name, " (AUC=", auc_value, ")")  # Append AUC to model name
  )
})

roc_df <- bind_rows(roc_data)
ggplot(roc_df, aes(x = FPR, y = TPR, group = Model, color = Model)) + 
  geom_line(linewidth = 1.5) +  
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), linewidth = 1, linetype = "dashed") + 
  theme_bw() +
  theme(legend.position = c(0.8, 0.2),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(face = "bold", color = "black", size = 12),
        axis.text.y = element_text(face = "bold", color = "black", size = 12),
        axis.title.x = element_text(face = "bold", color = "black", size = 12),
        axis.title.y = element_text(face = "bold", color = "black", size = 12),
        plot.title = element_text(face = "bold", color = "black", size = 16, hjust = 0.5)) +
  labs(y = "True Positive Rate (Sensitivity)",
       x = "False Positive Rate (1 - Specificity)",
       title = "ROC Curves")

############# PR Visualization ######################
pr_models <- list(
  "APRI" = ls_for_pr$APRI,
  "APRI_unadjusted" = ls_for_pr$`APRI_unadjusted `,
  "FSI" = ls_for_pr$FSI,
  "TyG" = ls_for_pr$TyG,
  "FSI_unadjusted" = ls_for_pr$`FSI_unadjusted `,
  "TyG_unadjusted"=ls_for_pr$`TyG_unadjusted `
)
names(pr_models) <- c('APRI adjusted','APRI unadjusted','FSI adjusted','TyG adjusted','FSI unadjusted','TyG unadjuted')

# Compute PR curves and extract Precision, Recall, and PR-AUC
pr_data <- lapply(names(pr_models), function(name) {
  model_data <- pr_models[[name]]
  
  # Compute PR curve
  pr_curve <- pr.curve(scores.class0 = model_data$Probabilities0,
                       scores.class1 = model_data$Probabilites1, 
                       curve = TRUE)
  
  pr_auc <- round(pr_curve$auc.integral, 3)  # Extract and round PR-AUC
  
  # Create a dataframe with Precision-Recall values
  data.frame(
    Recall = pr_curve$curve[, 1],    # Recall (x-axis)
    Precision = pr_curve$curve[, 2], # Precision (y-axis)
    Model = paste0(name, " (PR-AUC=", pr_auc, ")")  # Append PR-AUC to model name
  )
})

# Combine all PR curve data
pr_df <- bind_rows(pr_data)
ggplot(pr_df, aes(x = Recall, y = Precision, group = Model, color = Model)) + 
  geom_line(linewidth = 1.5) +  
  theme_bw() +
  theme(legend.position = c(0.8, 0.2),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(face = "bold", color = "black", size = 12),
        axis.text.y = element_text(face = "bold", color = "black", size = 12),
        axis.title.x = element_text(face = "bold", color = "black", size = 12),
        axis.title.y = element_text(face = "bold", color = "black", size = 12),
        plot.title = element_text(face = "bold", color = "black", size = 16, hjust = 0.5)) +
  labs(y = "Precision",
       x = "Recall",
       title = "Precision-Recall Curves")




ggplot(roc_df, aes(x = FPR, y = TPR, group = Model, color = Model)) + 
  geom_line(linewidth = 1.5) +  
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), linewidth = 1, linetype = "dashed") + 
  theme_bw() +
  theme(legend.position = c(0.8, 0.2),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(face = "bold", color = "black", size = 10),
        axis.text.y = element_text(face = "bold", color = "black", size = 9),
        axis.title.x = element_text(face = "bold", color = "black", size = 12),
        axis.title.y = element_text(face = "bold", color = "black", size = 12),
        plot.title = element_text(face = "bold", color = "black", size = 14, hjust = 0.5)) +
  labs(y = "True Positive Rate (Sensitivity)",
       x = "False Positive Rate (1 - Specificity)",
       title = "ROC Curves")



################################################### ROC-PR CI ######################################################

library(precrec)
library(ggplot2)




class <- ls_for_roc$`APRI adjusted`$Class
prob <- ls_for_roc$`APRI adjusted`$Probabilites

inds1<-1:75
inds2<-76:150
inds3<-152:226
inds4<-228:302
inds5<-303:377
inds6<-378:452
inds7<-450:524
inds8<-528:602
inds9<-602:676
inds10<-677:751
ind<-list(inds1,inds2,inds3,inds4,inds5,inds6,inds7,inds8,inds9,inds10)

scores1<-prob[inds1]
labels1<-class[inds1]

scores2<-prob[inds2]
labels2<-class[inds2]

scores3<-prob[inds3]
labels3<-class[inds3]


scores4<-prob[inds4]
labels4<-class[inds4]


scores5<-prob[inds5]
labels5<-class[inds5]

scores6<-prob[inds6]
labels6<-class[inds6]

scores7<-prob[inds7]
labels7<-class[inds7]

scores8<-prob[inds8]
labels8<-class[inds8]

scores9<-prob[inds9]
labels9<-class[inds9]

scores10<-prob[inds10]
labels10<-class[inds10]

scores<-join_scores(scores1,scores2,scores3,scores4,scores5,
                    scores6,scores7,scores8,scores9,scores10)
labels<-join_labels(labels1,labels2,labels3,labels4,labels5,
                    labels6,labels7,labels8,labels9,labels10)



# Create the mmdata object
mdat <- mmdata(scores, labels, dsids = c(1,2,3,4,5,6,7,8,9,10))
curves<-evalmod(mdat)
r<-auc_ci(curves,alpha=0.05)
r[3:6]<-round(r[3:6],3)

