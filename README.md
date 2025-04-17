# NIT_CAD_analysis
R scripts for the data pre-processing, feature selection and analysis of NIT indices aiming for the identification of patience with obstructive CAD.

# Dataset_preparation.R
Load the clinical dataset from the GEnetic SYNTAX Score clinical trial registered on ClinicalTrials.gov (NCT03150680). Required data-preprocessing techniques are employed, alongside the calculation and visualization of NITs indices:  i) HSI, ii) TyG,
iii) FSI, iv) acNASH, v) FIB-4 and vi) APRI.

# Univariate_modelling.R
Load the pre-processed dataset, applies IQR outlier removal and univariate logistic regression models are employed to identify the most informative NITs for the identification of obstructive CAD.

# LR_models_evaluation.R
Evaluates the important NITs with the potential covariates and selects the optimal multivariate model for each NIT.

# Modelling_predictive.R
Load the filtered dataset from the outlier instances and defines the optimal models identified from the LR_models_evaluation.R
10-fold stratified cross validation scheme is employed for the evaluation of multivariable and univariate NIT models. 
Scott-Knott algorithm, Delong Test for ROC curves and Bootstrap test for PR are applied to the results alongside with curve visualization.
