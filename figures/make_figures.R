library(rstan)
library(tidyverse)
library(bayesplot)

#-------------------------------------------------
### Model Output
out_real_M5 <- readRDS("./model_outputs/M5_out_real.rds")
list_of_draws <- as.data.frame(out_real_M5)
print(out_real_M5, digits = 3)

out_real_M6_ymax <- readRDS("./model_outputs/M6_out_real_ymax.rds")
list_of_draws_ymax <- as.data.frame(out_real_M6_ymax)

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

for(i in 1:n_species){
  q <- q + geom_abline(intercept = ilogit(median(list_of_draws[,26]) + # row 26 is estimate for beta0
                                            median(list_of_draws[,(26+i)])), 
                       slope = ilogit(median(list_of_draws[,26]) +
                                        median(list_of_draws[,35+i])), # row 35 is estimate for mu_beta1
                       col = "grey",
                       lwd = 1, alpha = 0.75)
}

# add species specific intercepts and slopes
(q <- q +
    geom_abline(intercept = ilogit(median(list_of_draws[,26])), 
                slope = ilogit(median(list_of_draws[,26]) + 
                                 median(list_of_draws[,35])), 
                lwd = 2, alpha = 0.5,
                col = "black") 
)

# 95% confidence ints
(q <- q +
    geom_abline(intercept = ilogit(median(list_of_draws[,26]) +
                                     2 * sd(list_of_draws[,26])), 
                slope = ilogit(median(list_of_draws[,26]) + 
                                 median(list_of_draws[,35]) + 
                                 2 * sd(list_of_draws[,35])), 
                lwd = 1, alpha = 0.5,
                lty = "dashed",
                col = "blue") +
    geom_abline(intercept = ilogit(median(list_of_draws[,26]) -
                                     2 * sd(list_of_draws[,26])), 
                slope = ilogit(median(list_of_draws[,26]) + 
                                 median(list_of_draws[,35])-
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
                       limits = c(0, 100)) +
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
  p <- p + geom_abline(intercept = exp(median(list_of_draws[,8]) + # row 8 is estimate for beta0
                                            median(list_of_draws[,(8+i)])), 
                       slope = exp(median(list_of_draws[,8]) +
                                        median(list_of_draws[,17+i])), # row 17 is estimate for mu_alpha1
                       col = "grey",
                       lwd = 1, alpha = 0.75)
}
p

# add median community mean slope
(p <- p +
    geom_abline(intercept = exp(median(list_of_draws[,8])), 
                slope = exp(median(list_of_draws[,8]) + 
                                 median(list_of_draws[,17])), 
                lwd = 2, alpha = 0.5,
                col = "black") 
)

# 95% confidence ints
(p <- p +
    geom_abline(intercept = exp(median(list_of_draws[,8]) +
                                     2 * sd(list_of_draws[,8])), 
                slope = exp(median(list_of_draws[,8]) + 
                                 median(list_of_draws[,17]) + 
                                 2 * sd(list_of_draws[,17])), 
                lwd = 1, alpha = 0.5,
                lty = "dashed",
                col = "blue") +
    geom_abline(intercept = exp(median(list_of_draws[,8]) -
                                     2 * sd(list_of_draws[,8])), 
                slope = exp(median(list_of_draws[,8]) + 
                                 median(list_of_draws[,17])-
                                 2 * sd(list_of_draws[,17])), 
                lwd = 1, alpha = 0.5,
                lty = "dashed",
                col = "blue") 
)

