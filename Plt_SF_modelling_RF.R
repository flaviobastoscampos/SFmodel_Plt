# Evapotranspiration dynamics and partitioning in a grassed vineyard: 
# ecophysiological and computational modelling approaches
# 
# Flávio Bastos Campos
# 
# code for SFmodel
# model Tv in Plantaditsch vineyard using Tv(obs) and environmental variables. 


# Libraries & functions -------------------------------

library(tidyverse)
library(egg)          # to use ggarrange function
library(scales)       # to use the function date_format() in ggplot
library(caret)        # to use the train function in the modelling process


# computes the RMSE between two vectors of same length
RMSE_fc <- function(Actual, Predicted) {
  diff_x_y <- (Actual - Predicted)^2
  RMSE <- (sum(diff_x_y)/length(Actual))^0.5
  return(RMSE)
}



# modelling Tv ----------------------

# Using Random forest
# L.O.M.O. C.V. (Leave-one-month-out cross validation)

R2 <- vector()
R2_test <- vector()
RMSE_SF_fit <- vector()
M_index <- c(6,7,8)

V_Imp_list <- list()


# start LOOP  ----------------------------------

for (i in c(1:3)) {
  
  # load SF_ud (Tv_obs)
  setwd("C:/Users/bastosca/Desktop/data_SFmodel")
  SF_vines_h <-  readRDS("SF_vines_h.rds")

  # load Plt 2021
  data_2021 <- readRDS("Plt_data_21_h_gf_EBC1m_Fp_ET.rds")
  data_2021 <- data_2021 %>% filter(date>=as.POSIXct("2021-05-01 00:00:00", tz="GMT") & 
                                      date<=as.POSIXct("2021-12-31 23:59:00", tz="GMT"))
  
  # load Plt 2022
  data_2022 <- readRDS("Plt_data_22_h_gf_EBC1m_Fp_ET.rds")
  data_2022 <- data_2022 %>% filter(date>=as.POSIXct("2022-01-01 00:00:00", tz="GMT") & 
                                      date<=as.POSIXct("2022-08-31 23:59:00", tz="GMT"))  
  
  # merge 2022
  data_2022 <- data_2022 %>% 
    left_join(SF_vines_h[,c("date", "SF_B", "SF_C", "SF_F", "SF_ud_mean")], by="date") %>% 
    mutate(
      hour= hour(date),
      DoY= yday(date))
  
  # select columns and adjust names
  data_2022 <- data_2022[,c("date", "Year", "month", "day", "hour", "DoY",
                            "LE_f", "H_f",
                            "Rn_Wm2", "G_gf", "usoil_10cm", "Precip",
                            "Tair", "rH", "VPD_KPa", "Ustar", 
                            "ET_EC_mmh_f",
                            "SF_B", "SF_C", "SF_F", "SF_ud_mean")]
  
  colnames(data_2022) <- c("date", "Year", "month", "day", "hour", "DoY",
                           "LE", "H",
                           "Rn", "G", "usoil_10cm", "Precip",
                           "Tair", "rH", "VPD", "Ustar", 
                           "ET",
                           "SF_B", "SF_C", "SF_F", "SF_ud_mean")
  
  
  data <- data_2022 %>% filter(month>=6 & month<=8)
  
  # fit using data_2022
  data_2022_short <- data_2022[,c("date", 
                                  "LE", "H",
                                  "Rn", "G", "usoil_10cm", "Precip",
                                  "Tair", "rH", "VPD", "Ustar", 
                                  "ET")]
  data_2022_short <- data_2022_short[!apply(data_2022_short, 1, anyNA),]
  


  dat <- data[!apply(data[,c("ET", 
                             "LE", "H", "Rn", "G", "usoil_10cm", "Precip",
                             "Tair", "rH", "VPD", "Ustar",
                             "SF_B", "SF_C", "SF_F", "SF_ud_mean")], 1, anyNA),] 
  nrow(dat)
  levels(as_factor(dat$hour))
  
  # set dir
  dir_path <- file.path(getwd(), "RF")
  
  if (!dir.exists(dir_path)) {     # If the directory doesn't exist, create it
    dir.create(dir_path)
  }
  
  setwd(dir_path)  # Set the working directory to "RF"
  
  
  # join & merge
  data_2022_mean <- dat %>% select(-SF_B, - SF_C, -SF_F)
  names(data_2022_mean)[names(data_2022_mean)=="SF_ud_mean"] <- "SF_ud"
  dat <- data_2022_mean
  
  
  # RF ----------------------
  
  # step 1: get data_train ----------------------------------
  index= M_index[i]
  data_train <- dat %>% filter(month!=index) %>% as.data.frame(.)
  data_test <- dat %>% filter(month==index) %>% as.data.frame(.)
  
  
  # step 2: Train the model ----------------------------------
  model <- train(
    x= data_train[,c("LE", "H", "Rn", "G",
                     "usoil_10cm", "Precip", 
                     "Tair", "VPD", "Ustar")],
    y= data_train[,c("SF_ud")],
    method= "rf",
    ntree= 500,
    nodesize= 30,
    metric= "RMSE",
    trControl= trainControl(method= "boot", number= 15))
  
  df <- model$results %>% as.data.frame()
  R2[i] <- max(df[,3])
  
  SF_fit_test <- predict(model, data_test) %>% as.vector(.)
  reg <- lm(SF_fit_test ~ data_test$SF_ud)
  R2_test[i] <- round(summary(reg)$r.squared, digits= 7)
  RMSE_SF_fit[i]= round(RMSE_fc(SF_fit_test, data_test$SF_ud), 4)

  
  # step 3: Predict using the model ----------------------------------
  
  # fit using data_test
  SF_fit_test <- predict(model, data_test) %>% as.vector(.)
  data_test <- data_test %>% mutate(SF_fit_test= SF_fit_test) %>% as_tibble(.)
  
  # fit using dat (dat has months 6-8)
  SF_fit_dat <- predict(model, dat) %>% as.vector(.)
  dat <- dat %>% mutate(SF_fit_dat= SF_fit_dat) %>% as_tibble(.)
  
  SF_fit_data_2022 <- predict(model, data_2022_short) %>% as.vector(.)
  data_2022_short <- data_2022_short %>% mutate(SF_fit_data_2022= SF_fit_data_2022) %>% as_tibble(.)
  
  data_2022 <- data_2022 %>% left_join(data_2022_short[,c("date", "SF_fit_data_2022")], by="date")
  

  # set dir
  setwd("C:/Users/bastosca/Desktop/data_SFmodel")
  data_2021 <- readRDS("Plt_data_21_h_gf_EBC1m_Fp_ET.rds") %>% filter(Year==2021)
  
  data_2021_short <- data_2021[,c("date", "month", 
                                  "LE_f", "H_f", "Rn_Wm2", "G_gf",
                                  "usoil_10cm", "Precip",
                                  "Tair", "rH", "VPD_KPa", "Ustar", 
                                  "ET_EC_mmh_f")]
  
  colnames(data_2021_short) <- c("date", "month",
                                 "LE", "H", "Rn", "G", 
                                 "usoil_10cm", "Precip",
                                 "Tair", "rH", "VPD", "Ustar", 
                                 "ET")
  
  data_2021_short <- data_2021_short[!apply(data_2021_short, 1, anyNA),]
  
  SF_fit_data_2021 <- predict(model, data_2021_short) %>% as.vector(.)
  data_2021_short <- data_2021_short %>% mutate(SF_fit_data_2021= SF_fit_data_2021) %>% as_tibble(.)
  
  data_2021 <- data_2021 %>% left_join(data_2021_short[,c("date", "SF_fit_data_2021")], by="date")
  
  
  # set dir
  dir_path <- file.path(getwd(), "RF")
  
  if (!dir.exists(dir_path)) {     # If the directory doesn't exist, create it
    dir.create(dir_path)
  }
  
  setwd(dir_path)  # Set the working directory to "RF"
  
  saveRDS(model, file= paste0("model_", index, ".rds"))
  saveRDS(dat, file=paste0("dat_NAfree_", index, ".rds"))
  saveRDS(data_train, file=paste0("data_train_", index, ".rds"))
  saveRDS(data_test, file=paste0("data_test_", index, ".rds"))
  
  saveRDS(dat, file= paste0("dat_", index, ".rds"))
  saveRDS(data, file= paste0("data_", index, ".rds"))
  saveRDS(data_2021_short, file= paste0("data_2021_short_", index, ".rds"))
  
  saveRDS(data_2022, file=paste0("Plt_SF_fit_2022_", index, ".rds"))
  saveRDS(data_2021, file=paste0("Plt_SF_fit_2021_", index, ".rds"))
  
  
  # step 4: V_Imp of model ----------------------------------

  # All measures of importance are scaled to have a maximum value of 100, - Liaw and Wiener (2002):
  # Liaw A, Wiener M (2002). “Classification and Regression by randomForest.” R News, 2(3),
  # 18–22. URL http://CRAN.R-project.org/doc/Rnews/.
  
  V_Imp <- varImp(model)$importance
  V_Imp <- tibble::rownames_to_column(V_Imp)
  colnames(V_Imp) <- c("var", "Imp")
  V_Imp$Imp <- round(V_Imp$Imp, 0)
  # store in list
  V_Imp_list[[i]] <- V_Imp
}


# step 5: Test the model ----------------------------------

# R2 of the model
R2
mean(R2)

# R2 of the data_test
R2_test
mean(R2_test)
sd(R2_test)

# RMSE of the data_test
RMSE_SF_fit
mean(RMSE_SF_fit)
sd(RMSE_SF_fit)


# test model ----------------------

# scatter plots
par(mfrow=c(1,3), mai = c(0.35, 0.35, 0.35, 0.35))
M_index <- c(6,7,8)

for (i in c(1:3)) {
  index= M_index[i]

  # set dir
  dir_path <- file.path(getwd(), "RF")
  
  if (!dir.exists(dir_path)) {     # If the directory doesn't exist, create it
    dir.create(dir_path)
  }
  
  setwd(dir_path)  # Set the working directory to "RF"
  
  data_test <- readRDS(paste0("data_test_", index, ".rds"))
  
  reg <- lm(SF_ud ~ SF_fit_test, data_test)   # y~x
  R2_vec <- round(summary(reg)$r.squared, digits= 3)
  p_value_vec <- round(summary(reg)$coefficients[,4][[2]], 8)
  a <- round(unname(coef(reg)[1]), digits= 3)
  b <- round(unname(coef(reg)[2]), digits= 3)
  
  # plot SF_fit_test x SF_ud (test the model using data_test)
  Lscale_0 <- min(quantile(data_test$SF_ud, 0), quantile(data_test$SF_fit_test, 0)) 
  Lscale_1 <- max(quantile(data_test$SF_ud, 1), quantile(data_test$SF_fit_test, 1)) 
  
  plot(data_test$SF_fit_test, data_test$SF_ud,
       xlim=c(Lscale_0, Lscale_1), ylim=c(Lscale_0, Lscale_1),  
       ylab= "SF_ud", xlab="SF_fitted",
       mtext(paste0("SF_ud= ", a, "+ ", b, "* SF_fitted", "\n R²= ", R2_vec), side=3))
  abline(a=0, b=1, col="black")
  abline(a=a, b=b, col="blue")
}


# metrics

# gives GCV, GRSq, nrpune, n_predictors and selected terms +eq.
# Summarizes the results and gives us the metrics about the model
res <- model$results %>% 
  filter(mtry== model$bestTune$mtry)

mtry <- res$mtry
model_RMSE <- res$RMSE
model_RMSE_sd <- res$RMSESD
model_R2 <- res$Rsquared
model_R2_sd <- res$RsquaredSD

tib <- tibble(
  R2= R2_vec,
  p_value= p_value_vec,
  mtry= mtry,
  model_RMSE= model_RMSE,
  model_RMSE_sd= model_RMSE_sd,
  model_R2= model_R2,
  model_R2_sd= model_R2_sd)

write.table(tib, file= "tib_Plt.txt", sep= "\t", row.names=F)


# done.
 


