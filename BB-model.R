###------------------------------------------------###
###      Code for prior sensitivity analysis       ###
###             of binomial-beta model             ###
###------------------------------------------------###
library(rstan)

## Data  
y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,
       2,1,5,2,5,3,2,7,7,3,3,2,9,10,4,4,4,4,4,4,4,10,4,4,4,5,11,12,
       5,5,6,5,6,6,6,6,16,15,15,9,4)
N <- c(20,20,20,20,20,20,20,19,19,19,19,18,18,17,20,20,20,20,19,19,18,18,25,24,
       23,20,20,20,20,20,20,10,49,19,46,27,17,49,47,20,20,13,48,50,20,20,20,20,
       20,20,20,48,19,19,19,22,46,49,20,20,23,19,22,20,20,20,52,46,47,24,14)

data <- list()
data$J <- length(N)
data$N <- N
data$y <- y


## Fitting binomial-beta model with two parameterizations
model1 <- stan_model(file="BB-model-para1.stan")
fit1 <- sampling(model1, data=data)
pos1 <- extract(fit1)

model2 <- stan_model(file="BB-model-para2.stan")
fit2 <- sampling(model2, data=data)
pos2 <- extract(fit2)


# posterior of transformed parameters
m_pos1 <- pos1$trans1
v_pos1 <- pos1$trans2

m_pos2 <- pos2$alpha
v_pos2 <- pos2$beta


## set of alternative hyperparameters 
L <- 50
Score1_H <- Score2_H <- matrix(NA, L, L)
Score1_KL <- Score2_KL <- matrix(NA, L, L)
para_set <- seq(0, 10, length=L+1)[-1]


## Sensitivity measures 
for(l1 in 1:L){
  for(l2 in 1:L){
    log_ratio <- function(mm1, mm2){
      val1 <- dgamma(mm1, para_set[l1], para_set[l1], log=T) - dgamma(mm1, 1, 1, log=T)
      val2 <- dgamma(mm2, para_set[l2], para_set[l2], log=T) - dgamma(mm2, 1, 1, log=T)
      return(val1 + val2)
    }
    
    # H-sensitivity
    numer1 <- mean( sqrt(exp(log_ratio(m_pos1, v_pos1))) )
    denom1 <- sqrt( mean(exp(log_ratio(m_pos1, v_pos1))) )
    Score1_H[l1, l2] <- 1 - ( numer1/denom1 ) 
    
    numer2 <- mean( sqrt(exp(log_ratio(m_pos2, v_pos2))) )
    denom2 <- sqrt( mean(exp(log_ratio(m_pos2, v_pos2))) )
    Score2_H[l1, l2] <- 1 - ( numer2/denom2 ) 
    
    # KL-sensitivity
    term1 <- log(mean( exp(log_ratio(m_pos1, v_pos1)) ))
    term2 <- mean( log_ratio(m_pos1, v_pos1) )
    Score1_KL[l1, l2] <- term1 - term2 
    
    term1 <- log(mean( exp(log_ratio(m_pos2, v_pos2)) ))
    term2 <- mean( log_ratio(m_pos2, v_pos2) )
    Score2_KL[l1, l2] <- term1 - term2 
  }
}

mKL <- round( max(Score1_KL, Score2_KL) )


## Visualization by ggplot 
library(scales)
library(akima)
library(ggplot2)
library(gridExtra)

pdf("BB-result.pdf", height=14, width=14, pointsize=13)
plots <- list()

# H-sensitivity (parameterization 1)
data1 <- data.frame(x=rep(para_set, times=L), y=rep(para_set, each=L), value=as.vector(Score1_H))
interp_result1 <- with(data1, interp(x=x, y=y, z=value, xo=seq(min(x), max(x), length=100), 
                                     yo=seq(min(y), max(y), length=100), linear=T))
interp_data1 <- as.data.frame(expand.grid(x=interp_result1$x, y=interp_result1$y))
interp_data1$z <- as.vector(interp_result1$z)
plots[[1]] <- ggplot(interp_data1, aes(x=x, y=y, fill=z)) +
  geom_tile() +
  scale_fill_gradient(low="blue", high="red", na.value="white", name="", limits=c(0, 1)) +
  theme_minimal() +
  labs(title="H-Sensitivity (Parameterization 1)", x="Hyperparameter for delta", y = "Hyperparameter for gamma") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.text.y = element_text(), 
        plot.title = element_text(hjust = 0.5, size = 25),  
        axis.title = element_text(size = 14),             
        legend.title = element_text(size = 13),           
        legend.text = element_text(size = 14)) +
  geom_point(aes(x=1, y=1), color="black", size=5, pch=4) 

# H-sensitivity (parameterization 2)
data2 <- data.frame(x=rep(para_set, times=L), y=rep(para_set, each=L), value=as.vector(Score2_H))
interp_result2 <- with(data2, interp(x=x, y=y, z=value, xo=seq(min(x), max(x), length=100), 
                                     yo=seq(min(y), max(y), length=100), linear=T))
interp_data2 <- as.data.frame(expand.grid(x=interp_result2$x, y=interp_result2$y))
interp_data2$z <- as.vector(interp_result2$z)
plots[[2]] <- ggplot(interp_data2, aes(x=x, y=y, fill=z)) +
  geom_tile() +
  scale_fill_gradient(low="blue", high="red", na.value="white", name="", limits=c(0, 1)) +
  theme_minimal() +
  labs(title="H-Sensitivity (Parameterization 2)", x="Hyperparameter for alpha", y = "Hyperparameter for beta") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.text.y = element_text(), 
        plot.title = element_text(hjust = 0.5, size = 25),  
        axis.title = element_text(size = 14),             
        legend.title = element_text(size = 13),           
        legend.text = element_text(size = 14)) +
  geom_point(aes(x=1, y=1), color="black", size=5, pch=4) 

# KL-sensitivity (parameterization 1)
data3 <- data.frame(x=rep(para_set, times=L), y=rep(para_set, each=L), value=as.vector(Score1_KL))
interp_result3 <- with(data3, interp(x=x, y=y, z=value, xo=seq(min(x), max(x), length=100), 
                                     yo=seq(min(y), max(y), length=100), linear=T))
interp_data3 <- as.data.frame(expand.grid(x=interp_result3$x, y=interp_result3$y))
interp_data3$z <- as.vector(interp_result3$z)
plots[[3]] <- ggplot(interp_data3, aes(x=x, y=y, fill=z)) +
  geom_tile() +
  scale_fill_gradient(low="blue", high="red", na.value="white", name="", limits=c(0, mKL)) +
  theme_minimal() +
  labs(title="KL-Sensitivity (Parameterization 1)", x="Hyperparameter for delta", y = "Hyperparameter for gamma") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.text.y = element_text(), 
        plot.title = element_text(hjust = 0.5, size = 25),  
        axis.title = element_text(size = 14),             
        legend.title = element_text(size = 13),           
        legend.text = element_text(size = 14)) +
  geom_point(aes(x=1, y=1), color="black", size=5, pch=4) 

# KL-sensitivity (parameterization 2)
data4 <- data.frame(x=rep(para_set, times=L), y=rep(para_set, each=L), value=as.vector(Score2_KL))
interp_result4 <- with(data4, interp(x=x, y=y, z=value, xo=seq(min(x), max(x), length=100), 
                                     yo=seq(min(y), max(y), length=100), linear=T))
interp_data4 <- as.data.frame(expand.grid(x=interp_result4$x, y=interp_result4$y))
interp_data4$z <- as.vector(interp_result4$z)
plots[[4]] <- ggplot(interp_data4, aes(x=x, y=y, fill=z)) +
  geom_tile() +
  scale_fill_gradient(low="blue", high="red", na.value="white", name="", limits=c(0, mKL)) +
  theme_minimal() +
  labs(title="KL-Sensitivity (Parameterization 2)", x="Hyperparameter for alpha", y = "Hyperparameter for beta") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.text.y = element_text(), 
        plot.title = element_text(hjust = 0.5, size = 25),  
        axis.title = element_text(size = 14),             
        legend.title = element_text(size = 13),           
        legend.text = element_text(size = 14)) +
  geom_point(aes(x=1, y=1), color="black", size=5, pch=4) 

grid.arrange(grobs=plots, ncol=2, nrow=2)
dev.off()



