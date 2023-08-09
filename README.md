![image](https://github.com/ybryan95/UCSC_Xena_Dimensionality_Reduction_Demo/assets/123009743/a5a09993-63e0-4331-845e-35ab276cb7bb)

# TCGA Liver Cancer Data Analysis

This repository contains R code for analyzing gene expression and phenotype data from the GDC TCGA Liver Cancer (LIHC) dataset. The code performs data preprocessing, exploratory data analysis, and applies machine learning techniques such as Decision Trees and Random Forests for predictive modeling.

## Table of Contents
- [Data Description](#data-description)
- [Data Preprocessing](#data-preprocessing)
- [Exploratory Data Analysis](#exploratory-data-analysis)
- [Machine Learning Models](#machine-learning-models)
- [Dependencies](#dependencies)

## Data Description <a name = "data-description"></a>

The data used in this project are two files from the GDC TCGA Liver Cancer (LIHC) dataset:

1. `TCGA-LIHC.htseq_fpkm.tsv`: This file contains gene expression data obtained through RNA sequencing. The data are in the form of Fragments Per Kilobase of transcript per Million mapped reads (FPKM), which is a common method for estimating gene expression levels.

2. `TCGA-LIHC.GDC_phenotype.tsv`: This file contains phenotype data for the liver cancer samples, including information such as fibrosis stage, sample type, and disease stage.

## Data Preprocessing <a name = "data-preprocessing"></a>

The code first loads the data and performs some basic exploratory data analysis, such as checking the dimensions of the data and the names of the columns. It then preprocesses the data by converting it into a numeric matrix, calculating the standard deviation for each row (gene), and selecting the top 1000 genes with the highest standard deviation. The phenotype data is also preprocessed to align with the gene expression data, and labels are defined based on fibrosis and disease stage.

## Exploratory Data Analysis <a name = "exploratory-data-analysis"></a>

The code performs exploratory data analysis by examining the distribution of fibrosis stage, sample type, and disease stage in the phenotype data. This helps to understand the characteristics of the data and inform the subsequent machine learning modeling.

## Machine Learning Models <a name = "machine-learning-models"></a>

The code applies two machine learning models to the preprocessed data:

1. **Decision Tree**: A decision tree model is trained on a subset of the data and used to predict the labels of the test data. The performance of the model is evaluated using a confusion matrix.

2. **Random Forest**: A random forest model is trained on the data, with the number of trees ranging from 1 to 50. The error rate of each model is calculated and plotted against the number of trees to visualize the performance of the models. The model is also evaluated using a confusion matrix.

In addition, the code performs 3-fold cross-validation to assess the robustness of the models. The data is divided into three folds, and the models are trained on two folds and tested on the third fold. The performance of the models is evaluated using precision, recall, and accuracy.

## Dependencies <a name = "dependencies"></a>

- R version 3.6.0 or above
- `rpart` package for decision tree modeling
- `rpart.plot` package for visualizing decision trees
- `randomForest` package for random forest modeling
- `caret` package for cross-validation and model evaluation

## Usage

1. Set the directory path to where the data files are located.
2. Load the data files into R.
3. Perform data preprocessing and exploratory data analysis.
4. Train the decision tree and random forest models, and evaluate their performance.
5. Perform 3-fold cross-validation to assess the robustness of the models.

## Note

This code is provided as is without any guarantees. It is for educational purposes and is free to be used and modified. The data used in this project is publicly available from the GDC TCGA Liver Cancer (LIHC) dataset.
```
