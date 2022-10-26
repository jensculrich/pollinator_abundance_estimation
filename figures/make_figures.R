# library(rstan)
library(tidyverse)
# library(bayesplot)

#-------------------------------------------------
### Model Output
out_real_M5 <- readRDS("./model_outputs/M5_out_real.rds")
list_of_draws <- as.data.frame(out_real_M5)
print(out_real_M5, digits = 3)

out_real_M6_ymax <- readRDS("./model_outputs/M6_out_real_ymax.rds")
list_of_draws_ymax <- as.data.frame(out_real_M6_ymax)
print(out_real_M6_ymax, digits = 3)

## Some required functions 
ilogit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-x))  

n_species <- 7 

#-------------------------------------------------
### New data

## M5
X <- as.numeric(c(0, 1))
Y <- as.numeric(c(-3, 3))
df <- as.data.frame(cbind(X, Y))


#-------------------------------------------------
### Effect of habitat on detection probability 

(q <- ggplot(df) +
    geom_jitter(aes(x=X, y=Y),
                size = 0) +
    theme_bw() +
    # scale_color_viridis(discrete=TRUE) +
    scale_x_continuous(name="", breaks = c(0, 1),
                       labels=c("mowed", "meadow"),
                       limits = c(0,1.05)) +
    scale_y_continuous(name="Detection Probability",
                       limits = c(0,1)) +
    guides(color = guide_legend(title = "")) +
    theme(legend.text=element_text(size=10),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 2),
          axis.title.y = element_text(size = 18),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
)

# add species specific intercepts and slopes
for(i in 1:n_species){
  q <- q + geom_abline(intercept = ilogit(median(list_of_draws[,26]) + # row 26 is estimate for beta0
                                            median(list_of_draws[,(26+i)])), 
                       # when X=1, the detection is det[1] = ilogit(median(list_of_draws[,26+i]) + median(list_of_draws[,35+1+i]))
                       # when x=0, the detection is det[0] = ilogit(median(list_of_draws[,26]))
                       # so slope is det[1] - det[0] / (1 - 0)
                       slope = ilogit(median(list_of_draws[,26]) + median(list_of_draws[,(26+i)]) + median(list_of_draws[,35+1+i])) -
                         ilogit(median(list_of_draws[,26]) + median(list_of_draws[,26+i])), # row 35 is estimate for mu_beta1
                       col = "grey",
                       lwd = 1, alpha = 0.75)
}
q

# add community mean
(q <- q +
    geom_abline(intercept = ilogit(median(list_of_draws[,26])), 
                # slope is the change in detection when the equation includes beta1*1.
                # when X=1, the detection is det[1] = ilogit(median(list_of_draws[,26]) + median(list_of_draws[,35]))
                # when x=0, the detection is det[0] = ilogit(median(list_of_draws[,26]))
                # so slope is det[1] - det[0] / (1 - 0)
                slope = ilogit(median(list_of_draws[,26]) + median(list_of_draws[,35])) -
                            ilogit(median(list_of_draws[,26])), # 
                lwd = 2, alpha = 0.5,
                col = "black") 
)

# 95% confidence ints
(q <- q +
    geom_abline(intercept = ilogit(median(list_of_draws[,26]) +
                                     2 * sd(list_of_draws[,26])), 
                slope = ilogit(median(list_of_draws[,26]) + 
                                 median(list_of_draws[,35]) -
                                 ilogit(median(list_of_draws[,26])) +
                                 2 * sd(list_of_draws[,35])), 
                lwd = 1, alpha = 0.5,
                lty = "dashed",
                col = "blue") +
    geom_abline(intercept = ilogit(median(list_of_draws[,26]) -
                                     2 * sd(list_of_draws[,26])), 
                slope = ilogit(median(list_of_draws[,26]) + 
                                 median(list_of_draws[,35]) -
                                 ilogit(median(list_of_draws[,26])) -
                                 2 * sd(list_of_draws[,35])), 
                lwd = 1, alpha = 0.5,
                lty = "dashed",
                col = "blue") 
)

#-------------------------------------------------
### Effect of habitat on abundance 

(p <- ggplot(df) +
    geom_jitter(aes(x=X, y=Y),
                size = 0) +
    theme_bw() +
    # scale_color_viridis(discrete=TRUE) +
    scale_x_continuous(name="", breaks = c(0, 1),
                       labels=c("mowed", "meadow"),
                       limits = c(0,1.05)) +
    scale_y_continuous(name="Mean Abundance",
                       limits = c(0, 60)) +
    guides(color = guide_legend(title = "")) +
    theme(legend.text=element_text(size=10),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 2),
          axis.title.y = element_text(size = 18),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
)

# add species specific intercepts and slopes
for(i in 1:n_species){
  p <- p + geom_abline(intercept = exp(median(list_of_draws[,8]) + # row 8 is estimate for alpha0
                                            median(list_of_draws[,(8+4)])), 
                       slope = exp(median(list_of_draws[,17+1+4])), # row 17 is estimate for mu_alpha1
                       col = "grey",
                       lwd = 1, alpha = 0.75)
}
p

# add median community mean slope
(p <- p +
    geom_abline(intercept = exp(median(list_of_draws[,8])), 
                slope = exp(median(list_of_draws[,17])), 
                lwd = 2, alpha = 0.5,
                col = "black") 
)

# 95% confidence ints
(p <- p +
    geom_abline(intercept = exp(median(list_of_draws[,8]) +
                                     2 * sd(list_of_draws[,8])), 
                slope = exp(median(list_of_draws[,17]) + 
                                 2 * sd(list_of_draws[,17])), 
                lwd = 1, alpha = 0.5,
                lty = "dashed",
                col = "blue") +
    geom_abline(intercept = exp(median(list_of_draws[,8]) -
                                     2 * sd(list_of_draws[,8])), 
                slope = exp(median(list_of_draws[,17]) -
                                 2 * sd(list_of_draws[,17])), 
                lwd = 1, alpha = 0.5,
                lty = "dashed",
                col = "blue") 
)

#-------------------------------------------------
### Effect of habitat on abundance WITHOUT accounting for detection

(r <- ggplot(df) +
   geom_jitter(aes(x=X, y=Y),
               size = 0) +
   theme_bw() +
   # scale_color_viridis(discrete=TRUE) +
   scale_x_continuous(name="", breaks = c(0, 1),
                      labels=c("mowed", "meadow"),
                      limits = c(0,1.05)) +
   scale_y_continuous(name="Mean Abundance",
                      limits = c(0, 60)) +
   guides(color = guide_legend(title = "")) +
   theme(legend.text=element_text(size=10),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18),
         axis.title.x = element_text(size = 2),
         axis.title.y = element_text(size = 18),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"))
)

# add species specific intercepts and slopes
for(i in 1:n_species){
  r <- r + geom_abline(intercept = exp(median(list_of_draws_ymax[,8]) + # row 8 is estimate for alpha0
                                         median(list_of_draws_ymax[,(8+i)])), 
                       slope = exp(median(list_of_draws_ymax[,17+1+i])), # row 17 is estimate for mu_alpha1
                       col = "grey",
                       lwd = 1, alpha = 0.75)
}
r

# add median community mean slope
(r <- r +
    geom_abline(intercept = exp(median(list_of_draws_ymax[,8])), 
                slope = exp(median(list_of_draws_ymax[,17])), 
                lwd = 2, alpha = 0.5,
                col = "black") 
)

# 95% confidence ints
(r <- r +
    geom_abline(intercept = exp(median(list_of_draws_ymax[,8]) +
                                  2 * sd(list_of_draws_ymax[,8])), 
                slope = exp(median(list_of_draws_ymax[,17]) + 
                              2 * sd(list_of_draws_ymax[,17])), 
                lwd = 1, alpha = 0.5,
                lty = "dashed",
                col = "blue") +
    geom_abline(intercept = exp(median(list_of_draws_ymax[,8]) -
                                  2 * sd(list_of_draws_ymax[,8])), 
                slope = exp(median(list_of_draws_ymax[,17]) -
                              2 * sd(list_of_draws_ymax[,17])), 
                lwd = 1, alpha = 0.5,
                lty = "dashed",
                col = "blue") 
)

#-------------------------------------------------
### Posterior model estimates for effect of habitat on detection and occupancy
X <- c(0.25, 0.75) 
Y <- c(mean(list_of_draws[,35]), mean(list_of_draws[,17])) # median detection mean and mean abundance mean
# lower is mean - 1.96*sd
lower_95 <- c(mean(list_of_draws[,35]) - 1.96*sd(list_of_draws[,35]), 
              mean(list_of_draws[,17]) - 1.96*sd(list_of_draws[,17])) 
upper_95 <- c(mean(list_of_draws[,35]) + 1.96*sd(list_of_draws[,35]), 
              mean(list_of_draws[,17]) + 1.96*sd(list_of_draws[,17])) 

df_estimates <- as.data.frame(cbind(X, Y, lower_95, upper_95))

species_estimates <- data.frame()

for(i in 1: n_species){
  
  species_estimates[1,i] <- mean(list_of_draws[,35+1+i])
  species_estimates[2,i] <- mean(list_of_draws[,17+1+i])
  
}

df_estimates <- cbind(df_estimates, species_estimates)

(s <- ggplot(df_estimates) +
   geom_point(aes(x=X, y=Y),
               size = 3) +
  geom_errorbar(aes(x=X, ymin=lower_95, ymax =upper_95),
                color="black",width=0.05,size=1,alpha=0.5) +
   theme_bw() +
   # scale_color_viridis(discrete=TRUE) +
   scale_x_continuous(name="Community Mean Effect of Meadow Habitat", breaks = c(0.25, 0.75),
                      labels=c("detection", "abundance"),
                      limits = c(0,1)) +
   scale_y_continuous(name="Posterior Model Estimate",
                      limits = c(-3, 3)) +
   guides(color = guide_legend(title = "")) +
   geom_hline(yintercept = 0, lty = "dashed") +
   theme(legend.text=element_text(size=10),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18),
         axis.title.x = element_text(size = 18),
         axis.title.y = element_text(size = 18),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"))
)

for(i in 1: n_species){
  
  test <- as.data.frame(cbind(X, df_estimates[1:2,4+i]))
  colnames(test) <- c("X", "Y")
  
  s <- s + geom_point(data = test, aes(x=X, y=Y), 
             col = "skyblue", size = 12, shape = "-", alpha = 0.7)
  
}

s

#-------------------------------------------------
### Posterior model estimates for naive regression
X <- c(0.5) 
Y <- c(mean(list_of_draws_ymax[,17])) # median detection mean and mean abundance mean
# lower is mean - 1.96*sd
lower_95 <- c(mean(list_of_draws_ymax[,17]) - 1.96*sd(list_of_draws_ymax[,17])) 
upper_95 <- c(mean(list_of_draws_ymax[,17]) + 1.96*sd(list_of_draws_ymax[,17])) 

df_estimates <- as.data.frame(cbind(X, Y, lower_95, upper_95))

species_estimates <- data.frame()

for(i in 1: n_species){
  
  species_estimates[1,i] <- mean(list_of_draws_ymax[,17+1+i])
  
}

df_estimates <- cbind(df_estimates, species_estimates)

(t <- ggplot(df_estimates) +
    geom_point(aes(x=X, y=Y),
               size = 3) +
    geom_errorbar(aes(x=X, ymin=lower_95, ymax =upper_95),
                  color="black",width=0.05,size=1,alpha=0.5) +
    theme_bw() +
    # scale_color_viridis(discrete=TRUE) +
    scale_x_continuous(name="Community Mean Effect of Meadow Habitat", breaks = c(0.5),
                       labels=c("abundance"),
                       limits = c(0,1)) +
    scale_y_continuous(name="Posterior Model Estimate",
                       limits = c(-3, 3)) +
    guides(color = guide_legend(title = "")) +
    geom_hline(yintercept = 0, lty = "dashed") +
    theme(legend.text=element_text(size=10),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
)

for(i in 1: n_species){
  
  test <- as.data.frame(cbind(X, df_estimates[1,4+i]))
  colnames(test) <- c("X", "Y")
  
  t <- t + geom_point(data = test, aes(x=X, y=Y), 
                      col = "skyblue", size = 12, shape = "-", alpha = 0.7)
  
}

t
