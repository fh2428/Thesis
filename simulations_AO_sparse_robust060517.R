setwd("~/College/FINISHLINE/Thesis/Thesis/THESIS/R and Simulations")
setwd("H:/R/THESIS")


# rm(list=ls())

BICglmnet<-function(x,y,beta){
  n<-length(y)
  df<-length(which(beta!=0))
  sigma<-sqrt(mean((y-x%*%beta)^2))
  BIC<-log(sigma)+df*log(n)/n
}

library(glmnet)
library(robustbase)
library(robustHD)
library(MASS)
library(data.table)
library(forecast)

####################################################################################################################################
# Simulate true AR(p), p = 3  model. robust ------   ###############################################################################
####################################################################################################################################

# arguments
p <- 3
nsim <- 10 # nr. of simulation
burnin <- 200 # Burn in period
n <- burnin + 50 + p 
constant <- 0.5
truebeta <- c(0.5, 0.3, 0.2, 0.1)


# SIMULATION
# stores results ----
MAEE.OLS <- rep(NA, nsim) 
MAEE.LASSO <- rep(NA, nsim)
MAEE.LTS <- rep(NA, nsim)
MAEE.LTS.LASSO <- rep(NA, nsim)

MASEE.OLS <- rep(NA, nsim) 
MASEE.LASSO <- rep(NA, nsim) 
MASEE.LTS <- rep(NA, nsim) 
MASEE.LTS.LASSO <- rep(NA, nsim) 

VAR.OLS <- rep(NA, nsim) 
VAR.LASSO <- rep(NA, nsim) 
VAR.LTS <- rep(NA, nsim) 
VAR.LTS.LASSO <- rep(NA, nsim) 

mean.ols2 <- c()
mean.lasso2 <- c()
mean.lts2 <- c()
mean.lts.lasso2 <- c()

mase.ols2 <- c()
mase.lasso2 <- c()
mase.lts2 <- c()
mase.lts.lasso2 <- c()



# contamination -----
epsilon <- seq(0, 0.5, by = 0.01)
#epsilon = 0; 0.01; 0.05; 0.1; 0.05

for (j in 1:51){

 for (isim in 1:nsim) { # outer loop nsim times  
  cat("Simulation run", isim, "\n") # how long code runs once executed
  
  # STEP 1: Generate data from an AR(3) model
  Y <- matrix(NA, nrow = n, ncol = 1)   # matrix Y with 1 column, and n rows: column values -> response AR(3) model [before: NA]
  
  Y[1:p,] <- constant + rnorm(p, mean = 0, sd = 0.9)    # first p rows of Y; previous lags not available -> use constant + error term

    for(t in (p+1):n){
    Y[t,] = truebeta[1] + truebeta[2]*Y[t-1,] + truebeta[3]*Y[t-2,] + truebeta[4]*Y[t-3,] + rnorm(1, mean = 0, sd = 0.1) # Y is now a matrix, not a column extract the elements from truebeta
  } 

  Y <- Y[- (1:burnin),] # remove the simulations of the burn in period

  #### Add AO
  outliers <- sample((p+1):(n-burnin), epsilon[j]*(n-burnin-p)) # randomly draw observations numbers to put the outliers (don't use first p=3 observations because those will be lost)
  Y[outliers] <- rnorm(length(outliers), mean = 15, sd = 0)

  # Data matrix
  DATA <- embed(Y, dimension = p+1) # set p=10 for order selection simulation
    
  # Response: first column of DATA
  Ymatrix <- DATA[,1]
  Xmatrix <- DATA[,-1]
 
   
  # STEP 2: Estimation -----
  # OLS
  fit.ols <- auto.arima(Ymatrix, xreg = Xmatrix, 
                        lambda = NULL, 
                        ic = "bic",   # order selection
                        max.q = 0, 
                        trace = FALSE,
                        max.p = 3, max.order = 3)
  
  
  
  # LTS
  fit.lts <- ltsReg(x = Xmatrix, 
                    y = Ymatrix, 
                    alpha = 0.75)
  
  
  # LASSO
  fit.glmnet <- glmnet(x = Xmatrix,      # Apply lasso for grid of lambda values (grid is internally chosen)
                       y = Ymatrix, 
                       family = "gaussian", 
                       intercept = TRUE,
                       alpha = 1) # alpha = 1 indicates Lasso penalty 
  BICvalues <- apply(fit.glmnet$beta, 2, BICglmnet, x = Xmatrix, y = Ymatrix) # BIC values for different lamdba values
  lambdanbr.opt <- which.min(BICvalues) #number of optimal lambda value, the one that minimizes BIC criterion
  beta.lasso <- c(fit.glmnet$a0[lambdanbr.opt], fit.glmnet$beta[ ,lambdanbr.opt]) #betahat lasso: intercept and predictor coefficients for optimal lambda value
  
  
  # LTS LASSO
  l0 <- lambda0(Xmatrix, Ymatrix) # smaller steps could be considered 
  l0_grid <- seq(0, l0, by = 0.1*l0)
  fit.lts.lasso <- sparseLTS(x = Xmatrix,
                             y = Ymatrix, 
                             lambda = l0_grid, 
                             mode = "fraction", 
                             alpha = 0.75, 
                             crit = "BIC",
                             initial = "sparse", pmax = 3)
  
  # STEP 3: Estimation accuracy ------     
  MAEE.OLS[isim] <- mean(abs(truebeta - fit.ols$coef))
  MAEE.LASSO[isim] <- mean(abs(truebeta - (beta.lasso)))
  MAEE.LTS[isim] <- mean(abs(truebeta - summary(fit.lts)$coefficients[,1]))
  MAEE.LTS.LASSO[isim] <- mean(abs(truebeta - fit.lts.lasso$coefficients[,1]))
  
  MASEE.OLS[isim] <- mean((truebeta - fit.ols$coef)^2)
  MASEE.LASSO[isim] <- mean((truebeta - (beta.lasso)^2))
  MASEE.LTS[isim] <- mean((truebeta - summary(fit.lts)$coefficients[,1])^2)
  MASEE.LTS.LASSO[isim] <- mean((truebeta - fit.lts.lasso$coefficients[,1])^2)
  
  # STEP 4: Diff for 1st parameter: true - coeff 
  VAR.OLS[isim] <- (truebeta - fit.ols$coef)
  VAR.LASSO[isim] <- (truebeta - (beta.lasso))
  VAR.LTS[isim] <- (truebeta - summary(fit.lts)$coefficients[2,1])
  VAR.LTS.LASSO[isim] <- (truebeta - fit.lts.lasso$coefficients[2,1])
  
  }

  # Mean of MAEE
  mean.ols <- mean(MAEE.OLS)
  mean.lasso <- mean(MAEE.LASSO)
  mean.lts <- mean(MAEE.LTS)
  mean.lts.lasso <- mean(MAEE.LTS.LASSO)
  
  mean.ols2 <- cbind(mean.ols2, mean.ols)
  mean.lasso2 <- cbind(mean.lasso2, mean.lasso)
  mean.lts2 <- cbind(mean.lts2, mean.lts)
  mean.lts.lasso2 <- cbind(mean.lts.lasso2, mean.lts.lasso)
  
  # Mean of MASE
  mase.ols <- mean(MASEE.OLS)
  mase.lasso <- mean(MASEE.LASSO)
  mase.lts <- mean(MASEE.LTS)
  mase.lts.lasso <- mean(MASEE.LTS.LASSO)
  
  mase.ols2 <- cbind(mase.ols2, mase.ols)
  mase.lasso2 <- cbind(mase.lasso2, mase.lasso)
  mase.lts2 <- cbind(mase.lts2, mase.lts)
  mase.lts.lasso2 <- cbind(mase.lts.lasso2, mase.lts.lasso)
  
  # Variance truebeta - 1st coef over isim
  var.ols <- var(VAR.OLS)
  var.lasso <- var(VAR.LASSO)
  var.lts <- var(VAR.LTS)
  var.lts.lasso <- var(VAR.LTS.LASSO)
  
  
}

mean.ols2t <- t(mean.ols2)
mean.lasso2t <- t(mean.lasso2)
mean.lts2t <- t(mean.lts2)
mean.lts.lasso2t <- t(mean.lts.lasso2)

mase.ols2t <- t(mase.ols2)
mase.lasso2t <- t(mase.lasso2)
mase.lts2t <- t(mase.lts2)
mase.lts.lasso2t <- t(mase.lts.lasso2)

cont <- seq(from = 0, to = 0.5, by = 0.01)
contamination <- rep(cont, times = 4)
maee <- as.data.frame(as.table(rbind(mean.ols2t, mean.lasso2t, mean.lts2t, mean.lts.lasso2t)))
maee_data <- cbind(maee, cont)
maee_data$Var2 <- NULL
setnames(maee_data, old = c('Var1','Freq'), new = c('Estimator','MAEE'))

mase <- as.data.frame(as.table(rbind(mase.ols2t, mase.lasso2t, mase.lts2t, mase.lts.lasso2t)))
mase$Var2 <- NULL
setnames(mase, old = c('Var1','Freq'), new = c('Estimator2','MASE'))

data <- cbind(maee_data, mase)
write.csv(data, file = "PE_data_AO_AR3.csv")


# Efficiency ----
eff.lts <- (var.ols/var.lts)
eff.lasso <- (var.ols/var.lasso)
eff.lts.lasso <- (var.ols/var.lts.lasso)


# save.image("AO_simulations_AR3.rda")  

# Boxplot BIAS  -----(only for single contamination levels)
boxplot(data.frame(MAEE.OLS, MAEE.LASSO, MAEE.LTS, MAEE.LTS.LASSO),
        names = c("OLS-estimator","Lasso-estimator", "LTS-estimator", "LTS-LASSO-estimator"), 
        main = "Estimation accuracy MAEE", 
        ylab = "Bias (MAEE)", 
        ylim = c(0, 5))

boxplot(data.frame(MASEE.OLS, MASEE.LASSO, MASEE.LTS, MASEE.LTS.LASSO),
        names = c("OLS-estimator","Lasso-estimator", "LTS-estimator", "LTS-LASSO-estimator"), 
        main = "Estimation accuracy MASE", 
        ylab = "Bias (MASE)", 
        ylim = c(0, 5))


# Performance versus contamination ------


data_maee_ao_truear <- read.csv("PE_data_AO_AR3.csv")

# plot IO AR3
plot(data_maee_ao_truear$cont[data_maee_ao_truear$Estimator == "mean.ols"], data_maee_ao_truear$MAEE[data_maee_ao_truear$Estimator == "mean.ols"],
     type = "l",
     xlab = "contamination (in %)",
     ylab = "MAEE",
     ylim = c(0, 5),
     xlim = c(0, 0.5),
     main = "Prediction performance under additive outliers", frame.plot = TRUE, 
     col = "black", lwd = 2)
#axis(1, at = c(0.0, 0.01, 0.05, 0.1), labels = c("0.0", "0.01", "0.05", "0.1"))
lines(data_maee_ao_truear$cont[data_maee_ao_truear$Estimator == "mean.lasso"], data_maee_ao_truear$MAEE[data_maee_ao_truear$Estimator == "mean.lasso"], col = "orange", lwd = 2)
lines(data_maee_ao_truear$cont[data_maee_ao_truear$Estimator == "mean.lts"], data_maee_ao_truear$MAEE[data_maee_ao_truear$Estimator == "mean.lts"], col = "blue", lwd = 2)
lines(data_maee_ao_truear$cont[data_maee_ao_truear$Estimator == "mean.lts.lasso"], data_maee_ao_truear$MAEE[data_maee_ao_truear$Estimator == "mean.lts.lasso"], col = "red", lwd = 2)
#legend("top", c("OLS", "Lasso", "LTS", "LTS.Lasso"), lty=c("solid","dotdash","twodash","longdash"), cex = 0.75, horiz = TRUE)
legend("top", c("OLS", "Lasso", "LTS", "LTS.Lasso"), col = c("black","orange","blue","red"), xjust=0, yjust=0, 
       xpd = TRUE, horiz = TRUE, inset = c(0, 0),  text.width=c(0.01, 0.01, 0.01, 0.01), lwd=3, cex = 0.75, bty ="n")


plot(data_maee_ao_truear$cont[data_maee_ao_truear$Estimator == "mean.ols"], data_maee_ao_truear$MASE[data_maee_ao_truear$Estimator == "mean.ols"],
     type = "l",
     xlab = "contamination (in %)",
     ylab = "MSEE",
     ylim = c(0, 1),
     xlim = c(0, 0.5),
     main = "Prediction Performance", frame.plot = TRUE,
     col = "black", lwd = 2)
#axis(1, at = c(0.0, 0.01, 0.05, 0.1), labels = c("0.0", "0.01", "0.05", "0.1"))
lines(data_maee_ao_truear$cont[data_maee_ao_truear$Estimator == "mean.lasso"], data_maee_ao_truear$MASE[data_maee_ao_truear$Estimator == "mean.lasso"], col = "orange", lwd = 2)
lines(data_maee_ao_truear$cont[data_maee_ao_truear$Estimator == "mean.lts"], data_maee_ao_truear$MASE[data_maee_ao_truear$Estimator == "mean.lts"], col = "blue", lwd = 2)
lines(data_maee_ao_truear$cont[data_maee_ao_truear$Estimator == "mean.lts.lasso"], data_maee_ao_truear$MASE[data_maee_ao_truear$Estimator == "mean.lts.lasso"], col = "red", lwd = 2)
#legend("top", c("OLS", "Lasso", "LTS", "LTS.Lasso"), lty=c("solid","dotdash","twodash","longdash"), cex = 0.75, horiz = TRUE)
legend("top", c("OLS", "Lasso", "LTS", "LTS.Lasso"), col = c("black","orange","blue","red"), xjust=0, yjust=0, 
       xpd = TRUE, horiz = TRUE, inset = c(0, 0),  text.width=c(0.01, 0.01, 0.01, 0.01), lwd=3, cex = 0.75, bty ="n")


###################################################################################################################################################
# Simulate AR(p), p = 10       sparserobust ----- #################################################################################################
###################################################################################################################################################

rm(list=ls())

BICglmnet<-function(x,y,beta){
  n<-length(y)
  df<-length(which(beta!=0))
  sigma<-sqrt(mean((y-x%*%beta)^2))
  BIC<-log(sigma)+df*log(n)/n
}

library(glmnet)
library(robustbase)
library(robustHD)
library(MASS)
library(data.table)
library(forecast)

# arguments -----
p <- 10
nsim <- 10  # nr. of simulation
burnin <- 200 # Burn in period
n <- burnin + 100 + p 
constant <- 0.5
truebeta <- c(0.5, 0.3, 0.2, 0.1)


# SIMULATION
# stores results -----
MAEE.OLS <- rep(NA, nsim) 
MAEE.LASSO <- rep(NA, nsim)
MAEE.LTS <- rep(NA, nsim)
MAEE.LTS.LASSO <- rep(NA, nsim)

MASEE.OLS <- rep(NA, nsim) 
MASEE.LASSO <- rep(NA, nsim) 
MASEE.LTS <- rep(NA, nsim) 
MASEE.LTS.LASSO <- rep(NA, nsim) 

VAR.OLS <- rep(NA, nsim) 
VAR.LASSO <- rep(NA, nsim) 
VAR.LTS <- rep(NA, nsim) 
VAR.LTS.LASSO <- rep(NA, nsim) 

TP.OLS <- rep(NA, nsim) 
TP.LASSO <- rep(NA, nsim)
TP.LTS <- rep(NA, nsim)
TP.LTS.LASSO <- rep(NA, nsim)

FP.OLS <- rep(NA, nsim) 
FP.LASSO <- rep(NA, nsim)
FP.LTS <- rep(NA, nsim)
FP.LTS.LASSO <- rep(NA, nsim)

mean.ols2 <- c()
mean.lasso2 <- c()
mean.lts2 <- c()
mean.lts.lasso2 <- c()

mase.ols2 <- c()
mase.lasso2 <- c()
mase.lts2 <- c()
mase.lts.lasso2 <- c()

tp.ols2 <- c()
tp.lasso2 <- c()
tp.lts2 <- c()
tp.lts.lasso2 <- c()

fp.ols2 <- c()
fp.lasso2 <- c()
fp.lts2 <- c()
fp.lts.lasso2 <- c()

# contamination levels -----
epsilon <- seq(0, 0.5, by = 0.01)
#epsilon = 0; 0.01; 0.05; 0.1; 0.05



for (j in 1:5){
  
  for (isim in 1:nsim) { # outer loop nsim times  
    cat("Simulation run", isim, "\n") # how long code runs once executed
    
    # STEP 1: Generate data from an AR(3) model
    Y <- matrix(NA, nrow = n, ncol = 1)   # matrix Y with 1 column, and n rows: column values -> response AR(3) model [before: NA]
    
    Y[1:p,] <- constant + rnorm(p, mean = 0, sd = 0.9)    # first p rows of Y; previous lags not available -> use constant + error term
    
    for(t in (p+1):n){
      Y[t,] = truebeta[1] + truebeta[2]*Y[t-1,] + truebeta[3]*Y[t-2,] + truebeta[4]*Y[t-3,] + rnorm(1, mean = 0, sd = 0.1) # Y is now a matrix, not a column extract the elements from truebeta
    } 
    
    Y <- Y[- (1:burnin),] # remove the simulations of the burn in period
    
    #### Add AO (AR10) -----
    outliers <- sample((p+1):(n-burnin), epsilon[j]*(n-burnin-p)) # randomly draw observations numbers to put the outliers (don't use first p=3 observations because those will be lost)
    Y[outliers] <- rnorm(length(outliers), mean = 15, sd = 0)
    
    # Data matrix
    DATA <- embed(Y, dimension = p+1) # set p=10 for order selection simulation
    
    # Response: first column of DATA
    Ymatrix <- DATA[,1]
    Xmatrix <- DATA[,-1]
    
    
    # STEP 2: Estimation ------
    # OLS
    fit.ols <- auto.arima(Ymatrix,
                          d = 0,
                          lambda = NULL, 
                          ic = "bic",   # order selection
                          max.q = 0, 
                          trace = FALSE,
                          max.p = 3, max.order = 3)
    
    
    
    # LTS
    fit.lts <- ltsReg(x = Xmatrix, 
                      y = Ymatrix, 
                      alpha = 0.75)
    
    
    # LASSO
    fit.glmnet <- glmnet(x = Xmatrix,      # Apply lasso for grid of lambda values (grid is internally chosen)
                         y = Ymatrix, 
                         family = "gaussian", 
                         intercept = TRUE,
                         alpha = 1) # alpha = 1 indicates Lasso penalty 
    BICvalues <- apply(fit.glmnet$beta, 2, BICglmnet, x = Xmatrix, y = Ymatrix) # BIC values for different lamdba values
    lambdanbr.opt <- which.min(BICvalues) #number of optimal lambda value, the one that minimizes BIC criterion
    beta.lasso <- c(fit.glmnet$a0[lambdanbr.opt], fit.glmnet$beta[ ,lambdanbr.opt]) #betahat lasso: intercept and predictor coefficients for optimal lambda value
    
    
    # LTS LASSO
    l0 <- lambda0(Xmatrix, Ymatrix) # smaller steps could be considered 
    l0_grid <- seq(0, l0, by = 0.1*l0)
    fit.lts.lasso <- sparseLTS(x = Xmatrix,
                               y = Ymatrix, 
                               lambda = l0_grid, 
                               mode = "fraction", 
                               alpha = 0.75, 
                               crit = "BIC",
                               initial = "sparse", pmax = 3)
    
    # lasso
    beta.lasso.sparse <- names(beta.lasso)[which(beta.lasso!=0)]
    beta.lasso.sparse <- beta.lasso.sparse[!grepl("s", beta.lasso.sparse)]
    # True positives rate
    TP_lasso <- (sum(beta.lasso.sparse %in% c("V1", "V2", "V3")) 
                + sum(grepl("s", names(beta.lasso))))/length(truebeta)
    # False positive rate
    FP_lasso <- sum(! beta.lasso.sparse %in% c("V1", "V2", "V3"))/(p-length(truebeta)+1)
    
    #lts lasso
    lts.lasso = names(coef(fit.lts.lasso)[which(coef(fit.lts.lasso)!=0)])
    lts.lasso.sparse <- lts.lasso[-grepl("Intercept", lts.lasso)]
    TP_lts.lasso <- (sum(lts.lasso.sparse %in% c("x1", "x2", "x3")) 
                     + sum(grepl("Intercept", lts.lasso)))/length(truebeta)
    FP_lts.lasso <- sum(! lts.lasso.sparse %in% c("x1", "x2", "x3"))/(p-length(truebeta)+1)
    
    
    # auto arima
    arima.coef <- sub("ar", "",names(fit.ols$coef))
    arima.coef.nointerc <- as.integer(arima.coef[!grepl("intercept", arima.coef)])
    TP_ols <- (sum(arima.coef.nointerc %in% 1:3) + sum(grepl("intercept", arima.coef)))/length(truebeta)
    FP_ols <- sum(!arima.coef.nointerc %in% 1:3)/(p - length(truebeta) + 1)
    
    # get TP and fP for the sim run
    TP.OLS[isim] <- TP_ols
    FP.OLS[isim] <- FP_ols
    
    TP.LASSO[isim] <- TP_lasso
    FP.LASSO[isim] <- FP_lasso
    
    TP.LTS.LASSO[isim] <- TP_lts.lasso
    FP.LTS.LASSO[isim] <- FP_lts.lasso
    
    # STEP 3: Estimation accuracy and store results ------      
    MAEE.OLS[isim] <- mean(abs(truebeta - fit.ols$coef))
    MAEE.LASSO[isim] <- mean(abs(truebeta - (beta.lasso)))
    MAEE.LTS[isim] <- mean(abs(truebeta - summary(fit.lts)$coefficients[,1]))
    MAEE.LTS.LASSO[isim] <- mean(abs(truebeta - fit.lts.lasso$coefficients[,1]))
    
    MASEE.OLS[isim] <- mean((truebeta - fit.ols$coef)^2)
    MASEE.LASSO[isim] <- mean((truebeta - (beta.lasso)^2))
    MASEE.LTS[isim] <- mean((truebeta - summary(fit.lts)$coefficients[,1])^2)
    MASEE.LTS.LASSO[isim] <- mean((truebeta - fit.lts.lasso$coefficients[,1])^2)
    
    # STEP 4: Diff for 1st parameter: true - coeff ------
    VAR.OLS[isim] <- (truebeta - fit.ols$coef)
    VAR.LASSO[isim] <- (truebeta - (beta.lasso))
    VAR.LTS[isim] <- (truebeta - summary(fit.lts)$coefficients[,1])
    VAR.LTS.LASSO[isim] <- (truebeta - fit.lts.lasso$coefficients[,1])
    
  }
  
  # mean TP and FP rates
  tp.ols2 <- cbind(tp.ols2, mean(TP.OLS))
  tp.lasso2 <- cbind(tp.lasso2, mean(TP.LASSO))
  tp.lts2 <- cbind(tp.lts2, mean(TP.LTS))
  tp.lts.lasso2 <- cbind(tp.lts.lasso2, mean(TP.LTS.LASSO))
  
  fp.ols2 <- cbind(fp.ols2, mean(FP.OLS))
  fp.lasso2 <- cbind(fp.ols2, mean(FP.LASSO))
  fp.lts2 <- cbind(fp.lts2, mean(FP.LTS))
  fp.lts.lasso2 <- cbind(fp.lts.lasso2, mean(FP.LTS.LASSO))
  
  # Mean of MAEE
  mean.ols <- mean(MAEE.OLS)
  mean.lasso <- mean(MAEE.LASSO)
  mean.lts <- mean(MAEE.LTS)
  mean.lts.lasso <- mean(MAEE.LTS.LASSO)
  
  mean.ols2 <- cbind(mean.ols2, mean.ols)
  mean.lasso2 <- cbind(mean.lasso2, mean.lasso)
  mean.lts2 <- cbind(mean.lts2, mean.lts)
  mean.lts.lasso2 <- cbind(mean.lts.lasso2, mean.lts.lasso)
  
  # Mean of MASE
  mase.ols <- mean(MASEE.OLS)
  mase.lasso <- mean(MASEE.LASSO)
  mase.lts <- mean(MASEE.LTS)
  mase.lts.lasso <- mean(MASEE.LTS.LASSO)
  
  mase.ols2 <- cbind(mase.ols2, mase.ols)
  mase.lasso2 <- cbind(mase.lasso2, mase.lasso)
  mase.lts2 <- cbind(mase.lts2, mase.lts)
  mase.lts.lasso2 <- cbind(mase.lts.lasso2, mase.lts.lasso)
  
  # Variance truebeta - 1st coef over isim
  var.ols <- var(VAR.OLS)
  var.lasso <- var(VAR.LASSO)
  var.lts <- var(VAR.LTS)
  var.lts.lasso <- var(VAR.LTS.LASSO)
  
  
}

tp.ols2 
tp.lasso2
tp.lts2
tp.lts.lasso2

fp.ols2 
fp.lasso2
fp.lts2
fp.lts.lasso2


# store data for plot -----
mean.ols2t <- t(mean.ols2)
mean.lasso2t <- t(mean.lasso2)
mean.lts2t <- t(mean.lts2)
mean.lts.lasso2t <- t(mean.lts.lasso2)

mase.ols2t <- t(mase.ols2)
mase.lasso2t <- t(mase.lasso2)
mase.lts2t <- t(mase.lts2)
mase.lts.lasso2t <- t(mase.lts.lasso2)

cont <- seq(from = 0, to = 0.5, by = 0.01)
contamination <- rep(cont, times = 4)
maee <- as.data.frame(as.table(rbind(mean.ols2t, mean.lasso2t, mean.lts2t, mean.lts.lasso2t)))
maee_data <- cbind(maee, cont)
maee_data$Var2 <- NULL
setnames(maee_data, old = c('Var1','Freq'), new = c('Estimator','MAEE'))

mase <- as.data.frame(as.table(rbind(mase.ols2t, mase.lasso2t, mase.lts2t, mase.lts.lasso2t)))
mase$Var2 <- NULL
setnames(mase, old = c('Var1','Freq'), new = c('Estimator2','MASE'))

data <- cbind(maee_data, mase)
write.csv(data, file = "PE_data_AO_AR3.csv")


# Efficiency -----
eff.lts <- (var.ols/var.lts)
eff.lasso <- (var.ols/var.lasso)
eff.lts.lasso <- (var.ols/var.lts.lasso)


# save.image("AO_simulations_AR10.rda")  

# Boxplot BIAS -----
boxplot(data.frame(MAEE.OLS, MAEE.LASSO, MAEE.LTS, MAEE.LTS.LASSO),
        names = c("OLS-estimator","Lasso-estimator", "LTS-estimator", "LTS-LASSO-estimator"), 
        main = "Estimation accuracy MAEE", 
        ylab = "Bias (MAEE)", 
        ylim = c(0, 2))

boxplot(data.frame(MASEE.OLS, MASEE.LASSO, MASEE.LTS, MASEE.LTS.LASSO),
        names = c("OLS-estimator","Lasso-estimator", "LTS-estimator", "LTS-LASSO-estimator"), 
        main = "Estimation accuracy MASE", 
        ylab = "Bias (MASE)", 
        ylim = c(0, 2))


data_maee_ao_ar10 <- read.csv("PE_data_AO_AR10.csv")

# plot IO AR10 ------
plot(data_maee_ao_ar10$cont[data_maee_ao_ar10$Estimator == "mean.ols"], data_maee_ao_ar10$MAEE[data_maee_ao_ar10$Estimator == "mean.ols"],
     type = "l",
     xlab = "contamination (in %)",
     ylab = "MAEE",
     ylim = c(0, 5),
     xlim = c(0, 0.5),
     main = "Prediction performance under additive outliers", frame.plot = TRUE, 
     col = "black", lwd = 2)
#axis(1, at = c(0.0, 0.01, 0.05, 0.1), labels = c("0.0", "0.01", "0.05", "0.1"))
lines(data_maee_ao_ar10$cont[data_maee_ao_ar10$Estimator == "mean.lasso"], data_maee_ao_ar10$MAEE[data_maee_ao_ar10$Estimator == "mean.lasso"], col = "orange", lwd = 2)
lines(data_maee_ao_ar10$cont[data_maee_ao_ar10$Estimator == "mean.lts"], data_maee_ao_ar10$MAEE[data_maee_ao_ar10$Estimator == "mean.lts"], col = "blue", lwd = 2)
lines(data_maee_ao_ar10$cont[data_maee_ao_ar10$Estimator == "mean.lts.lasso"], data_maee_ao_ar10$MAEE[data_maee_ao_ar10$Estimator == "mean.lts.lasso"], col = "red", lwd = 2)
#legend("top", c("OLS", "Lasso", "LTS", "LTS.Lasso"), lty=c("solid","dotdash","twodash","longdash"), cex = 0.75, horiz = TRUE)
legend("top", c("OLS", "Lasso", "LTS", "LTS.Lasso"), col = c("black","orange","blue","red"), xjust=0, yjust=0, 
       xpd = TRUE, horiz = TRUE, inset = c(0, 0),  text.width=c(0.01, 0.01, 0.01, 0.01), lwd=3, cex = 0.75, bty ="n")


plot(data_maee_ao_ar10$cont[data_maee_ao_ar10$Estimator == "mean.ols"], data_maee_ao_ar10$MASE[data_maee_ao_ar10$Estimator == "mean.ols"],
     type = "l",
     xlab = "contamination (in %)",
     ylab = "MSEE",
     ylim = c(0, 1),
     xlim = c(0, 0.5),
     main = "Prediction Performance", frame.plot = TRUE,
     col = "black", lwd = 2)
#axis(1, at = c(0.0, 0.01, 0.05, 0.1), labels = c("0.0", "0.01", "0.05", "0.1"))
lines(data_maee_ao_ar10$cont[data_maee_ao_ar10$Estimator == "mean.lasso"], data_maee_ao_ar10$MASE[data_maee_ao_ar10$Estimator == "mean.lasso"], col = "orange", lwd = 2)
lines(data_maee_ao_ar10$cont[data_maee_ao_ar10$Estimator == "mean.lts"], data_maee_ao_ar10$MASE[data_maee_ao_ar10$Estimator == "mean.lts"], col = "blue", lwd = 2)
lines(data_maee_ao_ar10$cont[data_maee_ao_ar10$Estimator == "mean.lts.lasso"], data_maee_ao_ar10$MASE[data_maee_ao_ar10$Estimator == "mean.lts.lasso"], col = "red", lwd = 2)
#legend("top", c("OLS", "Lasso", "LTS", "LTS.Lasso"), lty=c("solid","dotdash","twodash","longdash"), cex = 0.75, horiz = TRUE)
legend("top", c("OLS", "Lasso", "LTS", "LTS.Lasso"), col = c("black","orange","blue","red"), xjust=0, yjust=0, 
       xpd = TRUE, horiz = TRUE, inset = c(0, 0),  text.width=c(0.01, 0.01, 0.01, 0.01), lwd=3, cex = 0.75, bty ="n")

