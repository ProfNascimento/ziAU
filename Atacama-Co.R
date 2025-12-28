## Libraries
library(tidyverse)
library(sp)
library(GWmodel)
library(caret)
library(factoextra)
library(ggplot2)
library(scales)

set.seed(123)

## DATA IMPORT
DB=read.csv("https://raw.githubusercontent.com/ProfNascimento/ziAU/refs/heads/main/Cobalt.csv",sep=";")
View(DB)
str(DB)

# GENERAL DESCRIPTIVE 
par(mfrow = c(2, 1), mar=c(0,4,1,2))
hist(DB$Co/10, probability = T, xlim=c(0,0.3), main = "")
boxplot(DB$Co/10, horizontal = TRUE, ylim=c(0,0.3),
        frame.plot=F )

summary(DB$Co);sd(DB$Co)

# CORRELATION BETWEEN THE 40 ELEMENTS (COR PLOT)
graficar_heatmap <- function(data) {
  ggcorrplot::ggcorrplot(round(cor(data[sapply(data, is.numeric)], use = "complete.obs"), 1),
             method = "square", type = "lower", lab = TRUE, lab_size = 2,
             colors = c("blue", "white", "red"), tl.cex = 10,
             legend.title = "Corr", outline.color = "gray")
}
graficar_heatmap(DB[,-c(1:4)])

# CORRELATION COBALTO & SOME ELEMENTS (GGPAIRS)
elementos <- DB[, c("Co", "Bi", "Ni", "Fe", "Rb")]
GGally::ggpairs(elementos, 
                title = "Pairwise Minerals (Co, Bi, Ni, Fe, Rb)",
                upper = list(continuous = "cor"),  
                lower = list(continuous = "points"),
                diag = list(continuous = "barDiag")) 

# SPATIAL SAMPLE VISUALIZATION
DB$ZONE <- as.factor(DB$ZONE)
ggplot(DB, aes(LON, LAT, color = ZONE)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("indianred", "olivedrab", "steelblue", "mediumpurple")) +
  labs(title = "Spatial Distribution per Zone", x = "Longitude", y = "Latitude", color = "Zone") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = NA))


## DATA PRE-PROCESSING (PCA)
fit.pca <- prcomp(DB[,-c(1:4,18)], scale. = T, center = T)

fviz_screeplot(fit.pca,ncp=40, addlabels = TRUE) + coord_flip()
fviz_pca_var(fit.pca,axes = c(1, 2), col.var = "red")
fviz_pca_var(fit.pca,axes = c(1, 3), col.var = "red")
fviz_pca_var(fit.pca,axes = c(2, 3), col.var = "red")

# HOLD-OUT DATA SET
DB0=cbind(fit.pca$x,                        # PCs
          scale(DB[,-c(1:4,18)]),   # Std Geochemical Elements
          scale(DB[,1:2]),          # Std Latitude & Longitude
          Zone=DB[,4],                   # Zone
          Co=DB$Co/10)  %>%          # Y TRANSFORM
    as.data.frame()
colnames(DB0[,1:40])=paste0("PC", 1:40)

#str(DB0)
indices <- caret::createDataPartition(DB0$Y, p = 0.7, list = FALSE)
DB_train <- DB0[indices, ]
DB_test <- DB0[-indices, ]

#############################################################################
source("GWR-ziAU.R", local = TRUE)
#-----------------------------------------------------------------------------#
# CREACION SPATIALPOINTS DATAFRAME
spCob <- SpatialPointsDataFrame(
  coords      = data.frame(LON = DB_train$LON, LAT = DB_train$LAT),
  data        = data.frame(
    y     = DB_train$Co,
    LON  = DB_train$LON,
    LAT = DB_train$LAT,
    DB_train[, paste0("PC", 1:40)]
  ),
  proj4string = CRS("")
)

###############################################################################
# FORMULAS
form_coord <- ~ LON + LAT
form_FULL_PCA <- as.formula(paste("~", paste(paste0("PC", 1:40), collapse = " + ")))

################################################################################
## MODEL (ESPECIFICACION)
# GWR-ZIAU FULL 4 (NU & Alpha = FULL PCA); Enlace: Logit(nu) & log(alpha)
set.seed(123)
bw_cob1 <- bw.ggwr.ziau(
  formula    = update(form_FULL_PCA,y ~ .),
  data       = spCob,
  approach   = "CV",
  kernel     = "gaussian", 
  adaptive   = FALSE,
  link.alpha = "log",
  link.nu    = "logit",
  formula.nu = form_FULL_PCA,
  c_inflate  = 0
)

# FIT GWR-ziAU
fit_cob5 <- ggwr.ziau.basic(
  formula     = update(form_FULL_PCA,y ~ .),
  data        = spCob,
  bw          = bw_cob1,
  kernel      = "gaussian",
  adaptive    = FALSE,
  cv          = TRUE,
  link.alpha  = "log",
  link.nu     = "logit", 
  formula.nu  = form_FULL_PCA,
  nu.cut      = 0.7,
  c_inflate   = 0
)
print(fit_cob5)

# MODELO FULL 4 (TRAIN SET)
boxplot(fit_cob5$SDF$nu ~ ifelse(DB_train$Co==0,0,1), 
        main="GWR-ziAU Full #5", ylab="P(Ŷ=0 | PCs)", xlab="Y_{OBS}")

roc.theta = pROC::roc(as.factor(ifelse(DB_train$Co==0,0,1)), 
                      fit_cob5$SDF$nu)
plot(roc.theta)
BEST.p = pROC::coords(roc.theta, "best", ret="threshold", transpose = FALSE)
roc.theta$auc

yhat <- ifelse(fit_cob5$SDF$nu > as.numeric(BEST.p), fit_cob5$SDF$yhat, 0)

hist(yhat,xlim=c(0,0.20),ylim=c(0,250),breaks = 18)
hist(fit_cob5$SDF$yhat,xlim=c(0,0.20),ylim=c(0,250),xlab="Cobalto %wt (per 10)")

hist(fit_cob5$SDF$y,xlim=c(0,0.20),probability = T, xlab="Co %wt (per 10)",main="")

# Density fitted.01
mean_nu_full1 <- mean(fit_cob$SDF$nu) # promedio nu 
mean_alpha_full1 <- mean(fit_cob$SDF$alpha)
curve((1-mean_nu_full1)*(2/(mean_alpha_full1*x))*(log(x)/mean_alpha_full1)^2*(dnorm(log(x)/mean_alpha_full1)),0,1,add=T,col="blue",lty = 2, lwd = 1)

# Density fitted.02
mean_nu_full3 <- mean(fit_cob3$SDF$nu) # promedio nu 
mean_alpha_full3 <- mean(fit_cob3$SDF$alpha)
curve((1-mean_nu_full3)*(2/(mean_alpha_full3*x))*(log(x)/mean_alpha_full3)^2*(dnorm(log(x)/mean_alpha_full3)),0,1,add=T,col="red",lty = 2, lwd = 1)

# Density fitted.03
mean_nu_full4 <- mean(fit_cob4$SDF$nu) # promedio nu 
mean_alpha_full4 <- mean(fit_cob4$SDF$alpha)
curve((1-mean_nu_full4)*(2/(mean_alpha_full4*x))*(log(x)/mean_alpha_full4)^2*(dnorm(log(x)/mean_alpha_full4)),0,1,add=T,col="green",lty = 2, lwd = 2)

# Density fitted.04
mean_nu_full5 <- mean(fit_cob5$SDF$nu) # promedio nu 
mean_alpha_full5 <- mean(fit_cob5$SDF$alpha)
curve((1-mean_nu_full5)*(2/(mean_alpha_full5*x))*(log(x)/mean_alpha_full5)^2*(dnorm(log(x)/mean_alpha_full5)),0,1,add=T,col="purple",lty = 2, lwd = 1)

legend("top",c("GWR-ziAU.#1","GWR-ziAU.#2","GWR-ziAU.#3","GWR-ziAU.#4"),col=c("blue","red","green","purple"),lwd=rep(1,4))

#------------------------------------------------------------------------------#
# MODELS' PERFORMANCE (TEST SET)
# GWR-ZIAU FULL 4
## Ancla fórmulas y TRAIN dentro del fit
attr(fit_cob, "form_alpha") <- form_alpha_FULL_PCA
attr(fit_cob, "form_nu")    <- form_nu_FULL_PCA
attr(fit_cob, "sp_train")   <- spCob

pred_out <- predict(
  fit_cob5,
  newdata      = spTe,        
  type         = "mean",      
  c_inflate    = 0,           
  return_coefs = FALSE,       
  metrics      = TRUE
)

pred_out$yhat <- ifelse(pred_out$nu_hat >= 0.4676828, 0, pred_out$pred) # TRESHOLD 0.5 
resid_test <- pred_out$y_obs - pred_out$yhat

## Performance Metrics
# RMSE
round(sqrt(mean(resid_test^2)),4)
# MAE
round(mean(abs(resid_test)),4)
# NULL PRED
sum(ifelse(pred_out$y_obs == pred_out$yhat,1,0))/129 *100

# MSE solo para casos donde y == 0
zero_indices_m1 <- DB_test$Co == 0 # zero indices
mse_null1 <- mean((DB_test$Co[zero_indices_m1] - pred_out$yhat[zero_indices_m1])^2)


### INFORVIS PER ZONES
# PREDICTED NULLs (Co %wt)
tapply(fit_cob5$SDF$nu, as.factor(DB_train$Zone), mean)

fit_cob5$SDF %>% as.data.frame() %>% mutate(Zone=DB_train$Zone) %>% filter(Zone == 1) %>%  
  ggplot(aes(x = LON, y = LAT, color = nu, label = nu)) +
  geom_point(size = 2, alpha = 0.8) + xlab("WEST") + ylab("SOUTH")+
  #geom_text_repel() +
  scale_color_gradient(
    low = "blue", high = "red",
    name   = "P(Y=0|PCs)",
    limits = c(0, 1),                          
    breaks = c(0, 0.25, 0.50, 0.75, 1)
  ) +
  coord_equal()

fit_cob5$SDF %>% as.data.frame() %>% mutate(Zone=DB_train$Zone) %>% filter(Zone == 2) %>%  
  ggplot(aes(x = LON, y = LAT, color = nu, label = nu)) +
  geom_point(size = 2, alpha = 0.8) + xlab("WEST") + ylab("SOUTH")+
  #geom_text_repel() +
  scale_color_gradient(
    low = "blue", high = "red",
    name   = "P(Y=0|PCs)",
    limits = c(0, 1),                          
    breaks = c(0, 0.25, 0.50, 0.75, 1)
  ) +
  coord_equal()

fit_cob5$SDF %>% as.data.frame() %>% mutate(Zone=DB_train$Zone) %>% filter(Zone == 3) %>%  
  ggplot(aes(x = LON, y = LAT, color = nu, label = nu)) +
  geom_point(size = 2, alpha = 0.8) + xlab("WEST") + ylab("SOUTH")+
  #geom_text_repel() +
  scale_color_gradient(
    low = "blue", high = "red",
    name   = "P(Y=0|PCs)",
    limits = c(0, 1),                          
    breaks = c(0, 0.25, 0.50, 0.75, 1)
  ) +
  coord_equal()

fit_cob5$SDF %>% as.data.frame() %>% mutate(Zone=DB_train$Zone) %>% filter(Zone == 4) %>%  
  ggplot(aes(x = LON, y = LAT, color = nu, label = nu)) +
  geom_point(size = 2, alpha = 0.8) + xlab("WEST") + ylab("SOUTH")+
  #geom_text_repel() +
  scale_color_gradient(
    low = "blue", high = "red",
    name   = "P(Y=0|PCs)",
    limits = c(0, 1),                          
    breaks = c(0, 0.25, 0.50, 0.75, 1)
  ) +
  coord_equal()

## PCs EXPLICABILITY (DB => RAW DATA  &  DB0 => SCALED DATA)
# Quantiles -- Colbalt’s geochemical affinity elements 
vars <- c("Ni", "Cu", "Fe", "As", "Bi", "MnO", "Pb", "Zn")

DB %>%
  select(all_of(vars)) %>%
  reframe(
    across(
      everything(),
      ~ qv <- quantile(.x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
    )
  )

##########################
vars  <- c("Ni","Cu","Fe","As","Bi","MnO","Pb","Zn")
probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)

Q_raw <- DB0 %>% #filter(V84==3)  %>%     # FILTER PER ZONE 
  select(all_of(vars)) %>%
  summarise(
    across(everything(),
           ~ quantile(.x, probs = probs, na.rm = TRUE))
  ) %>%
  as.matrix()

L <- fit.pca$rotation[vars, ]
PCA_quantiles <- Q_raw %*% L

nu_hat <- 1 / (1 + exp(-(PCA_quantiles %*% tail(fit_cob5$glms$beta.nu, -1) + fit_cob5$glms$beta.nu[1])))
a=exp(PCA_quantiles %*% colMeans(fit_cob5$glms$beta.alpha[, -1]) + mean(fit_cob5$glms$beta.alpha[,1]))

#cbind(seq(0.1, 0.9,by=0.01),round(nu_hat,4))

# Conditional Mean (ziAU w/ Co correction - x10)
ifelse(nu_hat > 0.4676828, 0,  #based on the threshold
       (1-nu_hat)*2*exp(a^2/2)*( a*dnorm(a)+(1-pnorm(a))-2*a*dnorm(a)+a^2*(1-pnorm(a)) )*10)

##########################
vars2  <- c("Hg","Ag","Se","S","Tl","Cd")
probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)

Q_raw2 <- DB0 %>%
  select(all_of(vars2)) %>%
  summarise(
    across(everything(),
           ~ quantile(.x, probs = probs, na.rm = TRUE))
  ) %>%
  as.matrix()

L2 <- fit.pca$rotation[vars2, ]
PCA_quantiles2 <- Q_raw2 %*% L2

nu_hat2 <- 1 / (1 + exp(-(PCA_quantiles2 %*% tail(fit_cob5$glms$beta.nu, -1) + fit_cob5$glms$beta.nu[1])))

##########################
vars3  <- c("Au","P2O5","Sn","Mo")
probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)

Q_raw3 <- DB0 %>%
  select(all_of(vars)) %>%
  summarise(
    across(everything(),
           ~ quantile(.x, probs = probs, na.rm = TRUE))
  ) %>%
  as.matrix()

L3 <- fit.pca$rotation[vars, ]
PCA_quantiles3 <- Q_raw3 %*% L3

nu_hat3 <- 1 / (1 + exp(-(PCA_quantiles3 %*% tail(fit_cob5$glms$beta.nu, -1) + fit_cob5$glms$beta.nu[1])))

##########################
vars4  <- c("Cl","V","TiO2","Nb","Th","U")
Q_raw4 <- DB0 %>%
  select(all_of(vars)) %>%
  summarise(
    across(everything(),
           ~ quantile(.x, probs = probs, na.rm = TRUE))
  ) %>%
  as.matrix()
L4 <- fit.pca$rotation[vars, ]
PCA_quantiles4 <- Q_raw4 %*% L4
nu_hat4 <- 1 / (1 + exp(-(PCA_quantiles4 %*% tail(fit_cob5$glms$beta.nu, -1) + fit_cob5$glms$beta.nu[1])))

vars5  <- c("La","Y")
Q_raw5 <- DB0 %>%
  select(all_of(vars)) %>%
  summarise(
    across(everything(),
           ~ quantile(.x, probs = probs, na.rm = TRUE))
  ) %>%
  as.matrix()
L5 <- fit.pca$rotation[vars, ]
PCA_quantiles5 <- Q_raw5 %*% L5
nu_hat5 <- 1 / (1 + exp(-(PCA_quantiles5 %*% tail(fit_cob5$glms$beta.nu, -1) + fit_cob5$glms$beta.nu[1])))

vars6  <- c("K2O","Rb")
Q_raw6 <- DB0 %>%
  select(all_of(vars)) %>%
  summarise(
    across(everything(),
           ~ quantile(.x, probs = probs, na.rm = TRUE))
  ) %>%
  as.matrix()
L6 <- fit.pca$rotation[vars, ]
PCA_quantiles6 <- Q_raw6 %*% L6
nu_hat6 <- 1 / (1 + exp(-(PCA_quantiles6 %*% tail(fit_cob5$glms$beta.nu, -1) + fit_cob5$glms$beta.nu[1])))

vars7  <- c("CaO","MgO","Sr","Ba")
Q_raw7 <- DB0 %>%
  select(all_of(vars)) %>%
  summarise(
    across(everything(),
           ~ quantile(.x, probs = probs, na.rm = TRUE))
  ) %>%
  as.matrix()
L7 <- fit.pca$rotation[vars, ]
PCA_quantiles7 <- Q_raw7 %*% L7
nu_hat7 <- 1 / (1 + exp(-(PCA_quantiles7 %*% tail(fit_cob5$glms$beta.nu, -1) + fit_cob5$glms$beta.nu[1])))

##########################
vars8  <- c("Al2O3","SiO2","Cr","Ta","Hf","Zr")
probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)

Q_raw8 <- DB0 %>%
  select(all_of(vars)) %>%
  summarise(
    across(everything(),
           ~ quantile(.x, probs = probs, na.rm = TRUE))
  ) %>%
  as.matrix()

L8 <- fit.pca$rotation[vars, ]
PCA_quantiles4 <- Q_raw8 %*% L8

nu_hat8 <- 1 / (1 + exp(-(PCA_quantiles8 %*% tail(fit_cob5$glms$beta.nu, -1) + fit_cob5$glms$beta.nu[1])))
