###------------------------------------------------###
###      Code for prior sensitivity analysis       ###
###               of normal model                  ###
###------------------------------------------------###

## settings 
Y <- c(-2, -1, -0.5, 0, 0.5, 1, 2)   # observed data
n <- length(Y)
sig <- 1       # variance parameter (fixed)

## prior 
prior_mean <- 0
prior_prec <- 0.001

## posterior
MC <- 1000     # number of Monte Carlo samples
A <- 1/(n/sig^2+prior_prec)
B <- sum(Y)/sig^2 + prior_prec*prior_mean 
mu_pos <- rnorm(MC, A*B, sqrt(A))


## Computation of sensitivity measures
L <- 50     # number of alternative values for hyperparameters 
prior_mean_set <- seq(prior_mean, 5, length=L+1)[-1] 
prior_prec_set <- seq(prior_prec, 10, length=L+1)[-1] 

PM <- matrix(NA, L, L)          # Posterior mean
Score_H <- matrix(NA, L, L)     # H-sensitivity
Score_KL <- matrix(NA, L, L)    # KL-sensitivity
for(l1 in 1:L){
  for(l2 in 1:L){
    # Posterior mean 
    A <- 1/(n/sig^2+prior_prec_set[l2])
    B <- sum(Y)/sig^2 + prior_prec_set[l2]*prior_mean_set[l1] 
    PM[l1, l2] <- A*B
    
    # function for computing log-ratio of alternative and base priors
    log_ratio <- function(mm){
      dnorm(mm, prior_mean_set[l1], 1/sqrt(prior_prec_set[l2]), log=T) - dnorm(mm, prior_mean, 1/sqrt(prior_prec), log=T)
    }
    
    # H-sensitivity
    numer <- mean( sqrt(exp(log_ratio(mu_pos))) )
    denom <- sqrt( mean(exp(log_ratio(mu_pos))) )
    Score_H[l1, l2] <- 1 - ( numer/denom ) 
    
    # KL-sensitivity
    term1 <- log(mean( exp(log_ratio(mu_pos)) ))
    term2 <- mean( log_ratio(mu_pos) )
    Score_KL[l1, l2] <- term1 - term2 
  }
}


## Visualization by ggplot 
library(akima)
library(ggplot2)
library(gridExtra)

pdf("Normal-result.pdf", height=7, width=21, pointsize=15)
plots <- list()

# Posterior mean
data <- data.frame(x=rep(prior_mean_set, each=L), y=rep(prior_prec_set, times=L), 
                   value=as.vector(t(PM)))
interp_result <- with(data, interp(x = x, y = y, z = value, 
                                   xo = seq(min(x), max(x), length = 100), 
                                   yo = seq(min(y), max(y), length = 100), linear=T))
interp_data <- as.data.frame(expand.grid(x = interp_result$x, y = interp_result$y))
interp_data$z <- as.vector(interp_result$z)
plots[[1]] <- ggplot(interp_data, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_gradient(low="blue", high="red", na.value="white", name="") +
  theme_minimal() +
  labs(title = "Difference of posterior mean", x = "Prior mean", y = "Prior precision") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.text.y = element_text(), 
        plot.title = element_text(hjust = 0.5, size = 25),  
        axis.title = element_text(size = 14),             
        legend.title = element_text(size = 13),           
        legend.text = element_text(size = 14))

# H-sensitivity
data <- data.frame(x=rep(prior_mean_set, each=L), y=rep(prior_prec_set, times=L), 
                   value=as.vector(t(Score_H)))
interp_result <- with(data, interp(x = x, y = y, z = value, 
                                   xo = seq(min(x), max(x), length = 100), 
                                   yo = seq(min(y), max(y), length = 100), linear=T))
interp_data <- as.data.frame(expand.grid(x = interp_result$x, y = interp_result$y))
interp_data$z <- as.vector(interp_result$z)
data <- data.frame(x=rep(prior_mean_set, each=L), y=rep(prior_prec_set, times=L), 
                   value=as.vector(Score_H))
interp_result <- with(data, interp(x = x, y = y, z = value, 
                                   xo = seq(min(x), max(x), length = 100), 
                                   yo = seq(min(y), max(y), length = 100), linear=T))
interp_data <- as.data.frame(expand.grid(x = interp_result$x, y = interp_result$y))
interp_data$z <- as.vector(interp_result$z)
plots[[2]] <- ggplot(interp_data, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_gradient(low="blue", high="red", na.value="white", name="") +
  theme_minimal() +
  labs(title = "Sensitivity (H)", x = "Prior mean", y = "Prior precision") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.text.y = element_text(), 
        plot.title = element_text(hjust = 0.5, size = 25),  
        axis.title = element_text(size = 14),             
        legend.title = element_text(size = 13),           
        legend.text = element_text(size = 14))

# KL-sensitivity
data <- data.frame(x=rep(prior_mean_set, each=L), y=rep(prior_prec_set, times=L), 
                   value=as.vector(t(Score_KL)))
interp_result <- with(data, interp(x = x, y = y, z = value, 
                                   xo = seq(min(x), max(x), length = 100), 
                                   yo = seq(min(y), max(y), length = 100), linear=T))
interp_data <- as.data.frame(expand.grid(x = interp_result$x, y = interp_result$y))
interp_data$z <- as.vector(interp_result$z)
data <- data.frame(x=rep(prior_mean_set, each=L), y=rep(prior_prec_set, times=L), 
                   value=as.vector(Score_KL))
interp_result <- with(data, interp(x = x, y = y, z = value, 
                                   xo = seq(min(x), max(x), length = 100), 
                                   yo = seq(min(y), max(y), length = 100), linear=T))
interp_data <- as.data.frame(expand.grid(x = interp_result$x, y = interp_result$y))
interp_data$z <- as.vector(interp_result$z)
plots[[3]] <- ggplot(interp_data, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_gradientn(colors=c("blue", "red"), na.value = "white", name="") +
  theme_minimal() +
  labs(title = "Sensitivity (KL)", x = "Prior mean", y = "Prior precision") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.text.y = element_text(),  
        plot.title = element_text(hjust = 0.5, size = 25),  
        axis.title = element_text(size = 14),             
        legend.title = element_text(size = 13),           
        legend.text = element_text(size = 14))

grid.arrange(grobs=plots, ncol=3, nrow=1)
dev.off()

