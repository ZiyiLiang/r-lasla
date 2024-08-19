# Set working directory to the LASLA folder
setwd("C:/Users/liang/Documents/GitHub/r-lasla")
source('./methods/lasla_funcs.R')
source('./methods/utils.R')

library(adaptMT)

# Data source
# https://www.kaggle.com/competitions/loan-default-prediction/data

# Helper functions to preprocess the data
#-------------------------------------------------------------------------------
# Function to classify variable type
classify_variable <- function(x) {
  if (is.numeric(x)) {
    return("Continuous")
  } else if (is.factor(x) || is.character(x)) {
    return("Categorical")
  } else {
    return("Other")
  }
}


# Function to identify and factorize categorical-like numeric columns
factorize_categoricals <- function(df, threshold) {
  for (col in names(df)) {
    if (is.numeric(df[[col]]) && length(unique(df[[col]])) < threshold) {
      df[[col]] <- as.factor(df[[col]])
    }
  }
  return(df)
}


# Function to identify variables with fewer than two levels
find_single_level_factors <- function(df) {
  single_level_vars <- sapply(df, function(x) {
    if (is.factor(x) || is.character(x)) {
      return(length(unique(x)) < 2)
    }
    return(FALSE)
  })
  return(names(single_level_vars[single_level_vars]))
}
#-------------------------------------------------------------------------------
load = FALSE

if (load){
  zipfile <- "./data/loan_default.zip"
  filename_in_zip <- "train_v2.csv"
  
  loan_data_raw <- read.csv(unz(zipfile, filename_in_zip))
  
  # Remove all rows with any NA values
  loan_data <- na.omit(loan_data_raw)
  loan_data <- loan_data[, -which(names(loan_data) == "id")]
}
# # Convert all non-zero values in the loss column to 1
# loan_data$loss <- ifelse(loan_data$loss != 0, 1, 0)

# Set a threshold for numeric columns to be converted to factors
threshold <- 10
loan_data <- factorize_categoricals(loan_data, threshold)

# Remove all categorical variables 
numeric_vars <- sapply(loan_data, is.numeric)
loan_data_clean <- loan_data[, numeric_vars]

# # Find single level factors
# single_level_vars <- find_single_level_factors(loan_data)
# print(single_level_vars)
# 
# # Remove factors with fewer than two levels
# loan_data <- loan_data[, !names(loan_data) %in% single_level_vars]
# 
# # Count the occurrences of each classification
# variable_types <- sapply(loan_data, classify_variable)
# classification_counts <- table(variable_types)
# print(classification_counts)


# Determine the fold size
num_folds <- 5
fold_size <- nrow(loan_data) %/% num_folds
remainder <- nrow(loan_data) %% num_folds

# Create the folds dynamically
for (i in 1:num_folds) {
  if (i < num_folds) {
    start_idx <- (i - 1) * fold_size + 1
    end_idx <- i * fold_size
  } else {
    start_idx <- (i - 1) * fold_size + 1
    end_idx <- i * fold_size + remainder
  }
  fold_name <- paste0("fold", i)
  assign(fold_name, loan_data_clean[start_idx:end_idx, ])
}


# Use data-splitting to remove the highly correlated variables
y_fold1 = fold1["loss"]
x_fold1 = fold1[, names(fold1) != "loss"]

# Calculate the correlation matrix
cor_matrix <- cor(x_fold1, use = "pairwise.complete.obs")

# Set a correlation threshold
cor_threshold <- 0.8

# Find highly correlated pairs
high_cor_pairs <- which(abs(cor_matrix) > cor_threshold & abs(cor_matrix) < 1, arr.ind = TRUE)

# Print the highly correlated pairs
high_cor_pairs <- data.frame(
  Var1 = rownames(cor_matrix)[high_cor_pairs[, 1]],
  Var2 = colnames(cor_matrix)[high_cor_pairs[, 2]],
  Correlation = cor_matrix[high_cor_pairs]
)


# Create a vector to hold features to remove
features_to_remove <- c()

# Loop through the pairs and add one feature from each pair to the removal list
for (i in 1:nrow(high_cor_pairs)) {
  if (!high_cor_pairs$Var1[i] %in% features_to_remove & !high_cor_pairs$Var2[i] %in% features_to_remove){
    features_to_remove <- c(features_to_remove, high_cor_pairs$Var2[i])
  }
}

reduced_dim <- dim(x_fold1)[2]-length(features_to_remove)
t_list <- matrix(rep(0, (num_folds-1)*reduced_dim), num_folds-1, reduced_dim)
pv_list <- matrix(rep(0, (num_folds-1)*reduced_dim), num_folds-1, reduced_dim)

for (i in 2:num_folds){
  fold_name = paste0("fold",i)
  fold_data <- get(fold_name)
  
  y = fold_data["loss"]
  x = fold_data[, names(fold_data) != "loss"]
  x_reduced <- x[, !names(x) %in% features_to_remove]
  
  model <- lm(loss ~., data=cbind(x_reduced,y))
  summary <- summary(model)
  pv_list[i-1,] <- summary$coefficients[, "Pr(>|t|)"][2:(reduced_dim+1)]
  t_list[i-1,] <- summary$coefficients[, "t value"][2:(reduced_dim+1)]
}


m<-reduced_dim
mu0<-rep(0,m)
sd0<-rep(1,m)
fdr_level<-c(0.1,0.2)
nrep<-length(fdr_level)

bon.nr<-rep(0,nrep)
bh.nr<-rep(0,nrep)
lasla.dd.nr<-rep(0,nrep)
adapt.nr<-rep(0,nrep)

bon.de<-matrix(rep(0, m*nrep), nrep, m)
bh.de<-matrix(rep(0, m*nrep), nrep, m)
lasla.dd.de<-matrix(rep(0, m*nrep), nrep, m)
adapt.de<-matrix(rep(0, m*nrep), nrep, m)

# Primiary p-values 
pv <- pv_list[1,]
t <- t_list[1,]

# Auxiliary sequences
S <- t(t_list[2:dim(t_list)[1],])

# Create distance matrix from the auxiliary sequences
R<-cov(S)
d_lasla<-matrix(rep(0,m^2),m,m)
for (k in 1:m) {
  d_lasla[k,]=mahalanobis(S,S[k,],R)
}
d_lasla[lower.tri(d_lasla)] <- t(d_lasla)[lower.tri(d_lasla)]
d_lasla <- normalize_distance(t, d_lasla)

for (i in 1:nrep) {
  cat("Running with FDR level:", fdr_level[i], "\n")
  q<-fdr_level[i]
  
  bon.de[i,which(pv<=q/m)]<-1
  bon.nr[i]<-length(which(pv<=q/m))
  
  bh.res<-bh.func(pv, q)
  bh.de[i,]<-bh.res$de
  bh.nr[i]<-bh.res$nr
  
  bw=density(t,bw="SJ-ste")$bw
  pis_lasla <- lasla_pis(t, d_lasla, pv, tau=bh.func(pv,0.8)$th, h=bw)
  
  weight <- lasla_weights(t, d_lasla, pis_lasla, mu0, sd0, q, h=bw)
  lasla.dd.res <- lasla_thres(pv, pis_lasla, weight, q)
  lasla.dd.de[i,]<-lasla.dd.res$de
  lasla.dd.nr[i]<-lasla.dd.res$nr
  

  adapt.res <-  adapt_xgboost(S ,pv,
                              verbose = list(print = FALSE,
                                             fit = FALSE,
                                             ms = FALSE),
                              piargs = list("nrounds" = 50,
                                            "max_depth" = 1,
                                            "nthread" = 1,
                                            "verbose" = 0),
                              muargs = list("nrounds" = 50,
                                            "max_depth" = 1,
                                            "nthread" = 1,
                                            "verbose" = 0),
                              alphas = c(q),
                              nfits = 5)
  adapt.de[i, which(adapt.res$qvals <= q)]<-1
  adapt.nr[i]<-length(which(adapt.res$qvals <= q))
}

