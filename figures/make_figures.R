library(rstan)
library(tidyverse)
library(bayesplot)

#-------------------------------------------------
### Model Output
out_real_M5 <- readRDS("out_real_M5.rds")
list_of_draws <- as.data.frame(out_real_M5)

out_real_M6_ymax <- readRDS("out_real_M6_ymax.rds")
list_of_draws_ymax <- as.data.frame(out_real_M6_ymax)

## Some required functions 
ilogit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-x))  



#-------------------------------------------------
### Effect of habitat on 

newdata <- seq(0, 1, .01)
preds <- predict(out_real_M5, newdata = data.frame(x=newdata))
X <- as.numeric(c(0, 1))
Y <- as.numeric(c(-3, 3))
df <- as.data.frame(cbind(X, Y))

(q <- ggplot(df) +
    geom_jitter(aes(x=X, y=Y),
                size = 0) +
    theme_bw() +
    # scale_color_viridis(discrete=TRUE) +
    scale_x_continuous(name="Site Category", breaks = c(0, 1),
                       labels=c("mowed", "no-mow"),
                       limits = c(0,1)) +
    scale_y_continuous(name="Detection Probability",
                       limits = c(0,1)) +
    guides(color = guide_legend(title = "")) +
    theme(legend.text=element_text(size=10),
          axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
)

J <- 7
for(i in 1:J){
  q <- q + geom_abline(intercept = ilogit(mean(list_of_draws[,26]) +
                                            mean(list_of_draws[,(26+i)])), 
                       slope = ilogit(mean(list_of_draws[,26]) +
                                        mean(list_of_draws[,35+i])), 
                       col = "grey",
                       lwd = 1, alpha = 0.75)
}

# add species specific intercepts and slopes
(q <- q +
    geom_abline(intercept = ilogit(mean(list_of_draws[,26])), 
                slope = ilogit(mean(list_of_draws[,26]) + 
                                 mean(list_of_draws[,35])), 
                lwd = 2, alpha = 0.5,
                col = "black") 
)

# 95% confidence ints
(q <- q +
    geom_abline(intercept = ilogit(mean(list_of_draws[,26]) +
                                     2 * sd(list_of_draws[,26])), 
                slope = ilogit(mean(list_of_draws[,26]) + 
                                 mean(list_of_draws[,35]) + 
                                 2 * sd(list_of_draws[,35])), 
                lwd = 1, alpha = 0.5,
                lty = "dashed",
                col = "blue") +
    geom_abline(intercept = ilogit(mean(list_of_draws[,26]) -
                                     2 * sd(list_of_draws[,26])), 
                slope = ilogit(mean(list_of_draws[,26]) + 
                                 mean(list_of_draws[,35])-
                                 2 * sd(list_of_draws[,35])), 
                lwd = 1, alpha = 0.5,
                lty = "dashed",
                col = "blue") 
)

