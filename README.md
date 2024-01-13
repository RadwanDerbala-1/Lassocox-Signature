# Prognostic Gene Signature for Late-Stage Gastric Cancer using Lasso Cox 
Project of the Machine Learning Course Offered in Fall 2023 @ Zewail City

This project aims to develop a prognostic model (gene signature) using [lassocox] based on the levels of certain genes as prognostic biomarkers. 


## Dataset:

The datsets utilized are the following:

- TCGA-Stomach Adenocarcinoma (STAD), while excluding patients with H.pylori infection, tpm normalized RNA-SEQ data.
- GSE62254, gastric cancer dataset that consists of RMA normalized microarray data 


## Data Preprocessing (R): 

- for the training/testing dataset (TCGA-STAD): tpm normalization, gene duplicates removal, and stage I / III stratification were applied
- Differential Gene expression was performed between STAGE I/III (using Raw count reads), and DEGs were adapted based on the parameters (LFC= 1 , p-value=0.05)
- For the GSE 62254, RMA normalized reads were adapted, with their corresponding clinical (survival) data.
## Data Preprocessing (Jupyter notebook):

- for both the TCGA-STAD / GSE62254, standardscaler function from sklearn was used
- No NA are present because of the pre-processing steps performed on R 
- The survival data events were transformed from 0-1 to true or false (boolean type)

## CoxnetSurvivalAnalysis Model and Evaluation:

- The conxnetsurvivalanalysis() function was used to predict the B (Beta) coefficients for each corresponding differentially expressed gene
- Gridsearch was applied to test different values of alphas (different penalties), and different values of l1-ratio (elastic net)
- The results were evaluated directly using a method called the concordance test; which indicates how concordant the hazard ratio was in respect to the survival time, calculated by c-index. 
- Gene signature formulation, this step is performed by multiplying each coefficient by the quantity of gene expression
- The gene signature accuracy is evaluated using the cumulative and dynamic AUC test

## Files:
- Source: contains the source code of jupyter
-     Machine_STAD.RMD: contains the TCGA data pre-processing
-     GSE_.RMD:contains the GSE data preprocessing
-     Project_STAD.ipynb: contains the actual machine learning pipeline
- Data : contains data files needed by R and/or by python
-         Up/Down DEGs: Differentially expressed genes resulted from STAGE I / III
-         TCGA_Clinical: from which data that doesn't have h.pylori was made
-         mrna_tpm: Transcriptome of the 100 patients from the TCGA to be used in python
-         GSE_63500.csv: Transcriptome of the 300 patients from GEO with their survival data
-         ACRG Survival.xlsx: Survival data of the GSE, for the R part

## Example (Gastric)

The following lines show the functions of the coxnetsurvivalanalysis provided by Sklearn Survival package and its main parameters
```python
estimator=CoxnetSurvivalAnalysis(alphas=[0.18], n_alphas=20, l1_ratio=1, alpha_min_ratio=0.5 ,verbose = 10)
estimator.fit(X,Y) #X will stand for the features matrix # Y will contain the corresponding survival data
```
- n_alphas refers to the number of steps that the software will take from the alpha you initialize to the alpha_min_ratio
- alphas refers to the value of alpha from which the software will start with
- l1_ratio is the parameter responsible for the elastic net function (1 --> Lasso only) / (0 --> Ridge only) / (0.1-0.9 --> mixed impact of the two penalities)

The previous function can be integrated into Gridsearch to choose the best parameters for the model

```python
#Steps in elastic net parameter
l1_ratios = np.arange(1, -1, -0.05)
#number of cross fold validations (testing will be 20% and will repeat 5 times to cover all the samples)
cv = KFold(n_splits=5, shuffle=True, random_state=0)
gcv = GridSearchCV((CoxnetSurvivalAnalysis(alpha_min_ratio=1)),
    param_grid={"alphas": [[v] for v in alphas], "l1_ratio" : l1_ratios},
    cv=cv, 
    error_score=0.5,
    n_jobs=-1, #to work on the GPU
).fit(X, Y)
cv_results = pd.DataFrame(gcv.cv_results_)

# Display the relevant information
print("Best Parameters:", gcv.best_params_)
print("Best Score:", gcv.best_score_)
print("Grid Search Results:")
print(cv_results[['params', 'mean_test_score', 'std_test_score']])

```

## Dependencies

- scikit-survival           0.22.1 
- scikit-learn-intelex      2023.1.1
- scikit-learn              1.3.0 
- lifelines                 0.27.8
- matplotlib                3.7.2      
- matplotlib-base           3.7.2 
- matplotlib-inline         0.1.6  
- pandas                    2.0.3
- numpy                     1.24.3
