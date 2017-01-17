library(dplyr)
library(foreach)

load("GFR_data_original.Rdata")
thedata[2428,c("Weight..kg.", "Height..cm.")] <- 
  thedata[2428,c("Height..cm.", "Weight..kg.")]

DeBois <- function(Ht, Wt){
  0.007184 * Ht^0.725 * Wt^0.425
}
thedata$Surface.area..m2..DuBois <- DeBois(thedata$Height..cm., thedata$Weight..kg.)

thedata <- thedata[, c("EDTA.GFR", "Creat..mg.dL.", "Age", "Sex", "Weight..kg.", 
                       "Height..cm.", "Surface.area..m2..DuBois")]
names(thedata) <- c("GFR", "Creat", "Age", "Sex", "Wt", "Ht", "SufA")

thedata <- thedata[thedata$Age >= 18,]
thedata$log_Creat <- log(thedata$Creat)

transdata <- thedata %>% mutate(log_Creat = log(Creat)) %>%
  mutate(GFR=sqrt(GFR)) %>% dplyr::select(-Creat)

rmse <- function(x,y){
  sqrt(mean((x-y)^2))
}

################################################################################
################################################################################
################################################################################

n <- 100
split_func <- function(data=transdata, val.frac=1/5){
  n <- dim(data)[1]
  val_index <- sample(1:n, n*val.frac)
  data_val <- data[val_index,]
  data_dev <- data[-val_index,]
  return(list(data_dev, data_val))
}
set.seed(12345)
split_list <- lapply(1:n, function(i) split_func() )

################################################################################

models_BIC <- list()

foreach(j=1:n) %do% {
  full <- glm(GFR ~ .*. + I(log_Creat^2) + I(log_Creat^3), data = split_list[[j]][[1]], 
              family = gaussian(link = 'identity'), x=T, y=T)
  null <- glm(GFR ~ 1, data = split_list[[j]][[1]], 
              family = gaussian(link = 'identity'), x=T, y=T)
  simple <- glm(GFR ~Age + log_Creat + SufA + Sex, data = split_list[[j]][[1]])
  
  BIC_model <- step(null, scope = list(lower=null, upper=full), trace=F, 
                    k=log(dim(split_list[[j]][[1]])[1]))
  print(paste0("DONE Split ", j, " BIC"))
  models_BIC[[j]] <- BIC_model
}

# save(models_BIC, file= "Different_splits_BIC_models.RData")

################################################################################

source('Functions_for_step_cross_selection_gaussian_leave-one-out.R')

models_leave_one_out <- list()

foreach(j=1:n) %do% {
  full <- glm(GFR ~ .*. + I(log_Creat^2) + I(log_Creat^3), data = split_list[[j]][[1]], 
              family = gaussian(link = 'identity'), x=T, y=T)
  null <- glm(GFR ~ 1, data = split_list[[j]][[1]], 
              family = gaussian(link = 'identity'), x=T, y=T)
  simple <- glm(GFR ~Age + log_Creat + SufA + Sex, data = split_list[[j]][[1]])
  
  model <- RMSE_step(simple, direction = "both", 
                             scope=list(lower=null, upper=full), 
                             seed=seed, 
                             frac_improv = 0.999, trace=T)
  print(paste0("DONE Split ", j, " for leave-out-one"))
  models_leave_one_out[[j]] <- model
}

# save(models_leave_one_out, file= "Different_splits_leaveoneout_models.RData")


################################################################################




source('Functions_for_step_cross_selection_gaussian.R')
source("Table_list_glm_function.R")

models_5_fold <- list()


foreach(j=1:n) %do% {
  full <- glm(GFR ~ .*. + I(log_Creat^2) + I(log_Creat^3), data = split_list[[j]][[1]], 
              family = gaussian(link = 'identity'), x=T, y=T)
  null <- glm(GFR ~ 1, data = split_list[[j]][[1]], 
              family = gaussian(link = 'identity'), x=T, y=T)
  simple <- glm(GFR ~Age + log_Creat + SufA + Sex, data = split_list[[j]][[1]])
  
  n <- 10
  models <- list()
  for(i in 1:n){
    seed = sample(1:1000000, 1)
    models[[i]] <- RMSE_step(simple, direction = "both", 
                             scope=list(lower=null, upper=full), 
                             seed=seed, 
                             frac_improv = 0.999, trace=1)
  }
  model <- Table_list_glm(models)[[1]]$model
  print(paste0("DONE Split ", j, " for 5 fold"))
  models_5_fold[[j]] <- model
}
# save(models_5_fold, file= "Different_splits_5fold_models.RData")


################################################################################
################################################################################
################################################################################

freq_table_func <- function(models){
  coeficients <- lapply(models, function(x) names(x$coef) )
  
  vvv <- unique(unlist(coeficients))
  coef_mat <- matrix(NA, nrow=length(coeficients), ncol=length(vvv))
  for(i in 1:length(coeficients)){
    coef_mat[i,] <- vvv %in% coeficients[[i]]
  }
  colnames(coef_mat) <- vvv
  
  coef_mat <- coef_mat[,-c(1:7)]
  
  mm <- mgcv::uniquecombs(coef_mat)
  class(mm)
  ind <- attr(mm,"index")
  table(ind)
  
  df <- data.frame(mm)
  colnames(df) <- vvv[-c(1:7)]
  df$freq <- table(ind)
  df[order(df$freq, decreasing = T),]
}
# freq_table_func(models_BIC)
# freq_table_func(models_5_fold)
# freq_table_func(models_leave_one_out)
# 
# 
# 
# coeficients_BIC <- lapply(models_BIC, function(x) sort(names(x$coef)) )
# coeficients_5fold <- lapply(models_5_fold, function(x) sort(names(x$coef)) )
# coeficients_leave <- lapply(models_leave_one_out, function(x) sort(names(x$coef)) )
# final_coef <- coeficients_BIC[[96]]
# final_coef <- coeficients_5fold[[90]]
# 
# func <- function(y){
#   final_BIC <- unlist(lapply(coeficients_BIC, function(x) identical(x, y)))
#   final_5fold <- unlist(lapply(coeficients_5fold, function(x) identical(x, y)))
#   final_leave <- unlist(lapply(coeficients_leave, function(x) identical(x, y)))
#   final_df <- data.frame(final_BIC, final_5fold, final_leave)
#   sum(rowSums(final_df)==3)
# }
# sort(unlist(lapply(coeficients_5fold, func)))
# 
# final_BIC <- unlist(lapply(coeficients_BIC, function(x) identical(x, final_coef)))
# final_5fold <- unlist(lapply(coeficients_5fold, function(x) identical(x, final_coef)))
# final_leave <- unlist(lapply(coeficients_leave, function(x) identical(x, final_coef)))
# final_df <- data.frame(final_BIC, final_5fold, final_leave)
# sum(rowSums(final_df)==3)
# 
# 
# vvv <- unique(unlist(coeficients))
# coef_mat <- matrix(NA, nrow=length(coeficients), ncol=length(vvv))
# for(i in 1:length(coeficients)){
#   coef_mat[i,] <- vvv %in% coeficients[[i]]
# }
# colnames(coef_mat) <- vvv
# 
# coef_mat <- coef_mat[,-c(1:7)]
# 
# mm <- mgcv::uniquecombs(coef_mat)
# class(mm)
# ind <- attr(mm,"index")
# table(ind)
# 
# df <- data.frame(mm)
# colnames(df) <- vvv[-c(1:7)]
# df$freq <- table(ind)
# df[order(df$freq, decreasing = T),]
# 

