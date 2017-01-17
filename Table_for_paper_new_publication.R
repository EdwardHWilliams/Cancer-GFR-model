source("rmse_median_IQR_function.R")
source("CKD-EPI_model_functions.R")
source("Cockcroft_MDRD_formula.R")

GFR_fitted_values <- function(data){
  MDRD <- LM_revised_equation(data$Age, data$Creat, data$Sex) 
  MDRD_adj <- LM_revised_adj_equation(data$Age, data$Creat, data$Sex, data$SufA) 
  Cockcroft <- Cockcroft_equation_2(data$Age, data$Wt, data$Sex, data$Creat)
  Jelliffe <- Jelliffe_equation(data$Age, data$Creat, data$Sex)
  Jelliffe_adj <- Jelliffe_adj_equation(data$Age, data$Creat, data$Sex, data$SufA)
  Wright <- Wright_equation(data$Age, data$Creat, data$Sex, data$SufA)
  CKD <- Original_CKD_model(data$Sex, data$Creat, data$Age)
  CKD_adj <- Original_CKD_model_adjusted(data$Sex, data$Creat, data$Age, BSA = data$SufA)
  Mayo <- Mayo_equation(data$Age, data$Creat, data$Sex)
  Mayo_adj <- Mayo_adj_equation(data$Age, data$Creat, data$Sex, data$SufA)
  Martin <- Martin_equation(data$Age, data$Creat, data$Sex, data$Wt)
  fitted_values <- data.frame(CKD, CKD_adj, MDRD, MDRD_adj, Cockcroft, Jelliffe, Jelliffe_adj,
                              Wright, Mayo, Mayo_adj, Martin)
  return(fitted_values)
}

apply_func_to_fitted_vals <- function(func, actual, fitted, Dose=F){
  thefunc <- function(x){
    func(actual, x, Dose=Dose)
  }
  apply(fitted, 2, thefunc)
}
apply_func_to_fitted_vals_under60 <- function(func, actual, fitted, Dose=F){
  thefunc <- function(x){
    res = x<60
    func(actual, x, restriction=res, Dose=Dose)
  }
  apply(fitted, 2, thefunc)
}
apply_func_to_fitted_vals_60more <- function(func, actual, fitted, Dose=F){
  thefunc <- function(x){
    res = x>=60
    func(actual, x, restriction=res, Dose=Dose)
  }
  apply(fitted, 2, thefunc)
}



