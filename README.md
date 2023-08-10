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

<a href="https://xenabrowser.net/datapages/?host=https%3A%2F%2Fgdc.xenahubs.net&dataset=TCGA-LIHC.htseq_fpkm.tsv&allSamples=true&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443"><img src="https://img.shields.io/badge/RNASeq-E7E1E1?style=flat-square&logo=RNASeq&logoColor=white"/>

'https://xenabrowser.net/datapages/?host=https%3A%2F%2Fgdc.xenahubs.net&dataset=TCGA-LIHC.htseq_fpkm.tsv&allSamples=true&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443'

2. `TCGA-LIHC.GDC_phenotype.tsv`: This file contains phenotype data for the liver cancer samples, including information such as fibrosis stage, sample type, and disease stage.

<a href="https://xenabrowser.net/datapages/?dataset=TCGA-LIHC.GDC_phenotype.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443"><img src="https://img.shields.io/badge/phenotype-FFCA28?style=flat-square&logo=phenotype&logoColor=white"/>
*https://xenabrowser.net/datapages/?dataset=TCGA-LIHC.GDC_phenotype.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443*

## Data Preprocessing <a name = "data-preprocessing"></a>

The code first loads the data and performs some basic exploratory data analysis, such as checking the dimensions of the data and the names of the columns. It then preprocesses the data by converting it into a numeric matrix, calculating the standard deviation for each row (gene), and selecting the top 1000 genes with the highest standard deviation. The phenotype data is also preprocessed to align with the gene expression data, and labels are defined based on fibrosis and disease stage.

## Exploratory Data Analysis <a name = "exploratory-data-analysis"></a>

The code performs exploratory data analysis by examining the distribution of fibrosis stage, sample type, and disease stage in the phenotype data. This helps to understand the characteristics of the data and inform the subsequent machine learning modeling.

## Machine Learning Models <a name = "machine-learning-models"></a>
Both Decision Tree and Random Forest models are used for the same task: predicting the fibrosis condition of liver cancer samples based on gene expression data.
The code applies two machine learning models to the preprocessed data: 

1. **Decision Tree**: A decision tree model is trained on a subset of the data and used to predict the labels of the test data. The performance of the model is evaluated using a confusion matrix.
```bash
#Decision Tree
library(rpart)
library(rpart.plot)
n_sample <- dim(data2)[1]
idx_train <- sample(c(1:n_sample), round(n_sample * 0.67))
idx_test <- setdiff(c(1:n_sample), idx_train)
data_train <- data[idx_train, ]
data_test <- data[idx_test,]
# training
fit1 <- rpart(label3 ~., data=data_train, method = 'class')
rpart.plot(fit1)

# prediction
prediction1 <-predict(fit1, data_test, type = "class")

table_mat <- table(data_test$label3, prediction1)
table_mat
prediction2 <-predict(fit1, data_train, type = 'class')
table_mat <- table(data_train$label3, prediction2)
table_mat

```
2. **Random Forest**: A random forest model is trained on the data, with the number of trees ranging from 1 to 50. The error rate of each model is calculated and plotted against the number of trees to visualize the performance of the models. The model is also evaluated using a confusion matrix.
```bash
#Random forest
#install.packages('randomForest')
library('randomForest')
rf <- randomForest(label3~., data=data_train, ntree=10)
print(rf)
pre1 <- predict(rf, data_test)
table_mat_rf <- table(data_test$label3, t(pre1))
print(table_mat_rf)
pre1 <- predict(rf, data_train)
table_mat_rf <- table(data_train$label3, t(pre1))
print(table_mat_rf)

# Create a trees vs error plot
error_rate <- rep(0, 50) # record error rate for each model
for (i in 1:50) {
  rf <- randomForest(label3~., data=data_train, ntree=i)
  pre1 <- predict(rf, data_test)
  error_rate[i] <- 1 - sum(diag(table(data_test$label3, pre1))) / sum(table(data_test$label3, pre1))
}

plot(1:50, error_rate, type='b', xlab='Number of Trees', ylab='Error Rate')

```
In addition, the code performs 3-fold cross-validation to assess the robustness of the models. The data is divided into three folds, and the models are trained on two folds and tested on the third fold. The performance of the models is evaluated using precision, recall, and accuracy.
```bash
# 3 fold validation
#install.packages('caret')
library(caret)
# define 3 folds for cross-validation
#divide dataset into 3 folds randomly
set.seed(123) # set a seed for reproducibility
n_sample <- dim(data2)[1]
folds <- sample(1:3, n_sample, replace = TRUE) # randomly divide the data into 3 folds

#train the model using the first and second fold
train_idx <- which(folds != 3) # select the indices of the first and second fold for training
data_train <- data2[train_idx, ]
label_train <- label3[train_idx]

#Test the model on the 3rd fold
test_idx <- which(folds == 3) # select the indices of the third fold for testing
data_test <- data2[test_idx, ]
label_test <- label3[test_idx]

#evaluate models on the test datasets with confusion matrix, precision, recall and accuracy
library(class)
predicted_labels <- knn(data_train, data_test, label_train, k = 5) # make predictions on the test data using the trained model
conf_mat <- table(predicted_labels, label_test) # compute confusion matrix
conf_mat
precision <- diag(conf_mat)/colSums(conf_mat) # compute precision
precision
recall <- diag(conf_mat)/rowSums(conf_mat) # compute recall
recall
accuracy <- sum(diag(conf_mat))/sum(conf_mat) # compute accuracy
accuracy
```
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
