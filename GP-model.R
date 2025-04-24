###------------------------------------------------###
###      Code for prior sensitivity analysis       ###
###      of Gaussian process regression model      ###
###------------------------------------------------###
library(MASS)
library(rstan)

## settings 
n <- 50
set.seed(1)

## generation of synthetic data 
x <- sort( runif(n, 0, 3) )
mu <- sin(pi*x) + x
sig <- 0.5
y <- mu + sig*rnorm(n)

## MCMC 
model <- stan_model(file="GP-model.stan")
data <- list(n=n, x=x, y=y)
fit <- sampling(model, data=data, iter=2000, chains=1, seed=123)
pos <- extract(fit)

## posterior samples
sig_pos <- pos$sig
tau_pos <- pos$tau_mu
psi_pos <- pos$psi_mu
eta_pos <- pos$mu

## Index for simulating conditional posterior 
MC <- length(tau_pos)
num_local <- 50
Ind <- matrix(NA, MC, num_local)
for(k in 1:MC){
  eta <- eta_pos[k,]
  ord <- rank( apply((t(eta_pos)-eta)^2, 2, mean) )
  Ind[k,] <- which(ord<=num_local)
}

## set of alternative hyperparameters 
L <- 50
para_set <- seq(0, 10, length=L+1)[-1]


## Sensitivity measures 
Score1_H <- Score2_H <- matrix(NA, L, L)
Score1_KL <- Score2_KL <- matrix(NA, L, L)
for(l1 in 1:L){
  for(l2 in 1:L){
    log_ratio <- function(mm1, mm2){
      val1 <- dgamma(mm1, para_set[l1], para_set[l1], log=T) - dgamma(mm1, 1, 1, log=T)
      val2 <- dgamma(mm2, para_set[l2], para_set[l2], log=T) - dgamma(mm2, 1, 1, log=T)
      return(val1 + val2)
    }
    # joint posterior (H-sensitivity) 
    numer <- mean( sqrt(exp(log_ratio(tau_pos, psi_pos))) )
    denom <- sqrt( mean(exp(log_ratio(tau_pos, psi_pos))) )
    Score1_H[l1, l2] <- 1 - ( numer/denom ) 
    
    # joint posterior (KL-sensitivity) 
    term1 <- log(mean( exp(log_ratio(tau_pos, psi_pos)) ))
    term2 <- mean( log_ratio(tau_pos, psi_pos) )
    Score1_KL[l1, l2] <- term1 - term2 
    
    # marginal posterior (H- & KL-sensitivity)
    val <- c()
    for(k in 1:MC){
      val[k] <- mean( (exp(log_ratio(tau_pos[Ind[k,]], psi_pos[Ind[k,]]))) )
    }
    numer <- mean( sqrt(val) )
    Score2_H[l1, l2] <- 1 - ( numer/denom ) 
    term2 <- mean( log(val) )
    Score2_KL[l1, l2] <- term1 - term2
  }
}

# adjustment 
Score2_H[Score2_H<0] <- 0   
Score2_KL[Score2_KL<0] <- 0   


## Visualization by ggplot 
library(scales)
library(akima)
library(ggplot2)
library(gridExtra)

data1 <- data.frame(x=rep(para_set, times=L), y=rep(para_set, each=L), value=as.vector(Score1_H))
data2 <- data.frame(x=rep(para_set, times=L), y=rep(para_set, each=L), value=as.vector(Score2_H))
data3 <- data.frame(x=rep(para_set, times=L), y=rep(para_set, each=L), value=as.vector(Score1_KL))
data4 <- data.frame(x=rep(para_set, times=L), y=rep(para_set, each=L), value=as.vector(Score2_KL))

interp_result1 <- with(data1, interp(x=x, y=y, z=value, xo=seq(min(x), max(x), length=100), 
                                     yo=seq(min(y), max(y), length=100), linear=T))
interp_result2 <- with(data2, interp(x=x, y=y, z=value, xo=seq(min(x), max(x), length=100), 
                                     yo=seq(min(y), max(y), length=100), linear=T))
interp_result3 <- with(data3, interp(x=x, y=y, z=value, xo=seq(min(x), max(x), length=100), 
                                     yo=seq(min(y), max(y), length=100), linear=T))
interp_result4 <- with(data4, interp(x=x, y=y, z=value, xo=seq(min(x), max(x), length=100), 
                                     yo=seq(min(y), max(y), length=100), linear=T))

interp_data1 <- as.data.frame(expand.grid(x=interp_result1$x, y=interp_result1$y))
interp_data1$z <- as.vector(interp_result1$z)
interp_data2 <- as.data.frame(expand.grid(x=interp_result2$x, y=interp_result2$y))
interp_data2$z <- as.vector(interp_result2$z)
interp_data3 <- as.data.frame(expand.grid(x=interp_result3$x, y=interp_result3$y))
interp_data3$z <- as.vector(interp_result3$z)
interp_data4 <- as.data.frame(expand.grid(x=interp_result4$x, y=interp_result4$y))
interp_data4$z <- as.vector(interp_result4$z)


pdf("GP-result.pdf", height=14, width=14, pointsize=15)
plots <- list()
plots[[1]] <- ggplot(interp_data1, aes(x=x, y=y, fill=z)) +
  geom_tile() +
  scale_fill_gradient(low="blue", high="red", na.value="white", name="", limits=c(0, 0.25)) +
  theme_minimal() +
  labs(title="H-sensitivity (Joint posterior)", x="Hyperparameter for variance", y="Hyperparameter for range") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.text.y = element_text(), 
        plot.title = element_text(hjust = 0.5, size = 25),  
        axis.title = element_text(size = 14),             
        legend.title = element_text(size = 13),           
        legend.text = element_text(size = 14)) +
  geom_point(aes(x=1, y=1), color="black", size=5, pch=4) 

plots[[2]] <- ggplot(interp_data2, aes(x=x, y=y, fill=z)) +
  geom_tile() +
  scale_fill_gradient(low="blue", high="red", na.value="white", name="", limits=c(0, 0.25)) +
  theme_minimal() +
  labs(title="H-sensitivity (Marginal posterior)", x="Hyperparameter for variance", y = "Hyperparameter for range") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.text.y = element_text(), 
        plot.title = element_text(hjust = 0.5, size = 25),  
        axis.title = element_text(size = 14),             
        legend.title = element_text(size = 13),           
        legend.text = element_text(size = 14)) +
  geom_point(aes(x=1, y=1), color="black", size=5, pch=4) 

plots[[3]] <- ggplot(interp_data3, aes(x=x, y=y, fill=z)) +
  geom_tile() +
  scale_fill_gradient(low="blue", high="red", na.value="white", name="", limits=c(0, 2.05)) +
  theme_minimal() +
  labs(title="KL-sensitivity (Joint posterior)", x="Hyperparameter for variance", y = "Hyperparameter for range") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.text.y = element_text(), 
        plot.title = element_text(hjust = 0.5, size = 25),  
        axis.title = element_text(size = 14),             
        legend.title = element_text(size = 13),           
        legend.text = element_text(size = 14)) +
  geom_point(aes(x=1, y=1), color="black", size=5, pch=4) 

plots[[4]] <- ggplot(interp_data4, aes(x=x, y=y, fill=z)) +
  geom_tile() +
  scale_fill_gradient(low="blue", high="red", na.value="white", name="", limits=c(0, 2.05)) +
  theme_minimal() +
  labs(title="KL-sensitivity (Marginal posterior)", x="Hyperparameter for variance", y = "Hyperparameter for range") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.text.y = element_text(), 
        plot.title = element_text(hjust = 0.5, size = 25),  
        axis.title = element_text(size = 14),             
        legend.title = element_text(size = 13),           
        legend.text = element_text(size = 14)) +
  geom_point(aes(x=1, y=1), color="black", size=5, pch=4) 

grid.arrange(grobs=plots, ncol=2, nrow=2)
dev.off()


