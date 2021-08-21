#########################################################################
#### Xue Sun, Zexu Long, Jingbo Jia
#### A multi-scale Maxent approach to model habitat suitability for the giant pandas in the Qionglai mountain, China
#### Global Ecology and Conservation
#########################################################################

##### Load packages (install from CRAN if necessary)
library(dismo)
library(corrplot)
library(ecospat)
library(usdm)
library(raster)
library(ENMeval)
library(dplyr)
library(ggplot2)
library(ENMTools)

### set work dataspace
setwd("path/multiscale variables")
# read filtered panda occurrence
pandaocc <- read.csv(file = "path/pandaocc_1000m_filter.csv",
                 header = TRUE) # panda occurrence file was not available due to panda conservation reason
head(pandaocc)
pandaocc <- pandaocc[,-1]
# all multiscale variables
env.path <- list.files(path = "path/multiscale variables",
                       pattern = "asc$", 
                       full.names = TRUE) # path of variables
env.path
env <- stack(env.path)  # 读入变量

env.name <- list.files(path = "path/multiscale variables",
                       pattern = "asc$", 
                       full.names = FALSE) # name of variables
env.name

# univariate maxent modeling，calculate mean testAUC
AUC <- c()  # create a null vector to store AUC value
bg <- randomPoints(env[[1]],1000)
fold <- kfold(pandaocc, k=10)
for ( i in 1:length(env.name)){
  predictor <- stack(env.path[i])
  testAUC <- c()
  for (j in 1:10){
  occtest <- pandaocc[fold == j, ]
  occtrain <- pandaocc[fold != j, ]
  me <- maxent(predictor, occtrain, args = c("-P"))
  e <- evaluate(me, p = occtest, a = bg, x = predictor)
  testAUC[j] <- e@auc}
  AUC[i] <- mean(testAUC) 
}
AUC
write.csv(cbind(env.name, AUC), file = "/UnivariateAUC.csv")  # AUC for univariate Maxent modeling


#  plot variables at their best scale (Figure 3 in MS)
table <- read.csv(file = "path/optimized scale for variables.csv",
                  header = T)  # 
table
names(table) <- c("Variable", "Scale", "mean.testAUC")
tiff(file = "path/scatter point plot of variables at their best scale.tiff",
     width = 6, height = 4, units = "in", res = 300)
ggplot(data=table,aes(x=factor(Variable), y=factor(Scale)))+ 
  geom_point(size = 3, color = "blue")+
  xlab("Variable")+
  ylab("Scale(m)")+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  geom_text(label = factor(round(table$mean.testAUC, digits = 2)), 
    nudge_x = 0.25, nudge_y = 0.25, 
    check_overlap = T)
dev.off()

# optimized scale variables 
optvar <- raster::stack(list.files(path = "path/variable_at_optimized_scale",
                           pattern = "asc",
                           full.names = TRUE))
optvar
names(optvar) <- c("AI", "Asp", "Bam", "CLUM_cf", "CLUM_of", "Disfarm", "Disroad", "Disvil",
                   "Slp", "ED_cf", "ED_of", "Ele", "LPI_cf", "LPI_of", "NPP", "SHDI")
randombg <- randomPoints(optvar[[1]], n = 20000)  
data <- raster::extract(optvar, randombg)
tiff(file = "path/correlation plot between paired variables.tiff",
     width = 6, height = 6, units = "in", res = 300)
corrplot(corr = cor(data, method = "spearman"), method = "square", type = "lower",
         diag = FALSE, addCoef.col = "black", number.cex = 0.7, tl.cex = 0.7)  
dev.off()
optvar <- optvar[[c(-7,-9,-14,-16)]]   # delete 4 variables: disroad, LPI_of, SHDI, Slp
####################################################################

# multi-scale model calibrated using ENMeval
tune.args <- list(fc = c("L", "Q" ,"LQP"), rm = 1:4)
partitions <- "block"
partition.settings <- list(orientation = c("lat_lon"))
pandaENM <- ENMevaluate(occ = pandaocc, 
                        envs = optvar, 
                        tune.args = list(fc = c("L", "Q" ,"LQP"), rm = 1:4),
                        partitions = "randomkfold",
                        partition.settings = list(kfolds = 5),
                        algorithm = "maxent.jar", 
                        doClamp = T,
                        bin.output = FALSE,
                        parallel = TRUE,
                        taxon.name = "panda")
pandaENM@results
tiff(filename = "path/AICc vs rm and fc.tiff",
     width = 6, height = 4, unit = "in", res = 300)    # plot AICc varied according to different combinations of regularized multiplier and feature types
evalplot.stats(pandaENM, stats = "AICc", x.var = "rm", color.var = "fc", error.bars = T)
dev.off()
optmodel <- pandaENM@models[[which(pandaENM@results$delta.AICc == 0)]]    # optimized model selected by the smallest AICc
var.imp <- ecospat.maxentvarimport(optmodel, dfvar = optmodel@absence, nperm = 10) # variable importance
var.imp
pred <- predict(pandaENM@models[[which(pandaENM@results$delta.AICc == 0)]],
                   optvar,  args = c("outputformat=logistic"))     # model prediction
plot(pred)
AUCmulti <- pandaENM@results[which(pandaENM@results$delta.AICc == 0), 8]  # mean testAUC for the optimized Maxent model
CBImulti <- pandaENM@results[which(pandaENM@results$delta.AICc == 0), 10] # mean testCBI for the optimized Maxent model

e <- dismo::evaluate(p = pandaocc, a = randomPoints(optvar[[1]],10000),
                     model = optmodel, x = optvar)
e
t <- cbind(e@t, e@TNR, e@TPR, e@TPR+e@TNR-1)
max(e@TPR+e@TNR-1)        # threshold and TSS for the convert from continous layer to binary layer

###### single scale models
# 250m
env250 <- stack(list.files(path = "path/scale250",
                              pattern = "asc", full.names = T))
panda250 <- ENMevaluate(occ = pandaocc, 
                        envs = env250, 
                        tune.args = list(fc = c("LQP"), rm = 1),
                        partitions = "randomkfold",
                        partition.settings = list(kfolds = 5),
                        doClamp = T,
                        algorithm = "maxent.jar", 
                        bin.output = FALSE,
                        parallel = TRUE)
panda250@results
optmodel250 <- panda250@models[[which(panda250@results$delta.AICc == 0)]]
pred250 <- predict(panda250@models[[which(panda250@results$delta.AICc == 0)]],
                env250,  args = c("outputformat=logistic"))
plot(pred250)
AUC250 <- panda250@results[which(panda250@results$delta.AICc == 0), 8]
CBI250 <- panda250@results[which(panda250@results$delta.AICc == 0), 10]

# 500m
env500 <- stack(list.files(path = "path/scale500",
                           pattern = "asc", full.names = T))
panda500 <- ENMevaluate(occ = pandaocc, 
                        envs = env500, 
                        tune.args = list(fc = c("LQP"), rm = 1),
                        partitions = "randomkfold",
                        partition.settings = list(kfolds = 5),
                        doClamp = T,
                        algorithm = "maxent.jar", 
                        bin.output = FALSE,
                        parallel = TRUE)
panda500@results
optmodel500 <- panda500@models[[which(panda500@results$delta.AICc == 0)]]
pred500 <- predict(panda500@models[[which(panda500@results$delta.AICc == 0)]],
                   env500,  args = c("outputformat=logistic"))
plot(pred500)
AUC500 <- panda500@results[which(panda500@results$delta.AICc == 0), 8]
CBI500 <- panda500@results[which(panda500@results$delta.AICc == 0), 10]

# 1000m
env1000 <- stack(list.files(path = "path/scale1000",
                           pattern = "asc", full.names = T))
panda1000 <- ENMevaluate(occ = pandaocc, 
                        envs = env1000, 
                        tune.args = list(fc = c("LQP"), rm = 1),
                        partitions = "randomkfold",
                        partition.settings = list(kfolds = 5),
                        doClamp = T,
                        algorithm = "maxent.jar", 
                        bin.output = FALSE,
                        parallel = TRUE)
panda1000@results
optmodel1000 <- panda1000@models[[which(panda1000@results$delta.AICc == 0)]]
pred1000 <- predict(panda1000@models[[which(panda1000@results$delta.AICc == 0)]],
                   env1000,  args = c("outputformat=logistic"))
plot(pred1000)
AUC1000 <- panda1000@results[which(panda1000@results$delta.AICc == 0), 8]
CBI1000 <- panda1000@results[which(panda1000@results$delta.AICc == 0), 10]

# 2000m
env2000 <- stack(list.files(path = "path/scale2000",
                           pattern = "asc", full.names = T))
panda2000 <- ENMevaluate(occ = pandaocc, 
                        envs = env2000, 
                        tune.args = list(fc = c("LQP"), rm = 1),
                        partitions = "randomkfold",
                        partition.settings = list(kfolds = 5),
                        doClamp = T,
                        algorithm = "maxent.jar", 
                        bin.output = FALSE,
                        parallel = TRUE)
panda2000@results
optmodel2000 <- panda2000@models[[which(panda2000@results$delta.AICc == 0)]]
pred2000 <- predict(panda2000@models[[which(panda2000@results$delta.AICc == 0)]],
                   env2000,  args = c("outputformat=logistic"))
plot(pred2000)
AUC2000 <- panda2000@results[which(panda2000@results$delta.AICc == 0), 8]
CBI2000 <- panda2000@results[which(panda2000@results$delta.AICc == 0), 10]

# 3000m
env3000 <- stack(list.files(path = "path/scale3000",
                           pattern = "asc", full.names = T))
panda3000 <- ENMevaluate(occ = pandaocc, 
                        envs = env3000, 
                        tune.args = list(fc = c("LQP"), rm = 1),
                        partitions = "randomkfold",
                        partition.settings = list(kfolds = 5),
                        doClamp = T,
                        algorithm = "maxent.jar", 
                        bin.output = FALSE,
                        parallel = TRUE)
panda3000@results
optmodel3000 <- panda3000@models[[which(panda3000@results$delta.AICc == 0)]]
pred3000 <- predict(panda3000@models[[which(panda3000@results$delta.AICc == 0)]],
                   env3000,  args = c("outputformat=logistic"))
plot(pred3000)
AUC3000 <- panda3000@results[which(panda3000@results$delta.AICc == 0), 8]
CBI3000 <- panda3000@results[which(panda3000@results$delta.AICc == 0), 10]

# 4000m
env4000 <- stack(list.files(path = "path/scale4000",
                           pattern = "asc", full.names = T))
panda4000 <- ENMevaluate(occ = pandaocc, 
                        envs = env4000, 
                        tune.args = list(fc = c("LQP"), rm = 1),
                        partitions = "randomkfold",
                        partition.settings = list(kfolds = 5),
                        doClamp = T,
                        algorithm = "maxent.jar", 
                        bin.output = FALSE,
                        parallel = TRUE)
panda4000@results
optmodel4000 <- panda4000@models[[which(panda4000@results$delta.AICc == 0)]]
pred4000 <- predict(panda4000@models[[which(panda4000@results$delta.AICc == 0)]],
                   env4000,  args = c("outputformat=logistic"))
plot(pred4000)
AUC4000 <- panda4000@results[which(panda4000@results$delta.AICc == 0), 8]
CBI4000 <- panda4000@results[which(panda4000@results$delta.AICc == 0), 10]

# 5000m
env5000 <- stack(list.files(path = "path/scale5000",
                           pattern = "asc", full.names = T))
panda5000 <- ENMevaluate(occ = pandaocc, 
                        envs = env5000, 
                        tune.args = list(fc = c("LQP"), rm = 1),
                        partitions = "randomkfold",
                        partition.settings = list(kfolds = 5),
                        doClamp = T,
                        algorithm = "maxent.jar", 
                        bin.output = FALSE,
                        parallel = TRUE)
panda5000@results
optmodel5000 <- panda5000@models[[which(panda5000@results$delta.AICc == 0)]]
pred5000 <- predict(panda5000@models[[which(panda5000@results$delta.AICc == 0)]],
                   env5000,  args = c("outputformat=logistic"))
plot(pred5000)
AUC5000 <- panda5000@results[which(panda5000@results$delta.AICc == 0), 8]
CBI5000 <- panda5000@results[which(panda5000@results$delta.AICc == 0), 10]

# 6000m
env6000 <- stack(list.files(path = "path/scale6000",
                           pattern = "asc", full.names = T))
panda6000 <- ENMevaluate(occ = pandaocc, 
                        envs = env6000, 
                        tune.args = list(fc = c("LQP"), rm = 1),
                        partitions = "randomkfold",
                        partition.settings = list(kfolds = 5),
                        doClamp = T,
                        algorithm = "maxent.jar", 
                        bin.output = FALSE,
                        parallel = TRUE)
panda6000@results
optmodel6000 <- panda6000@models[[which(panda6000@results$delta.AICc == 0)]]
pred6000 <- predict(panda6000@models[[which(panda6000@results$delta.AICc == 0)]],
                   env6000,  args = c("outputformat=logistic"))
plot(pred6000)
AUC6000 <- panda6000@results[which(panda6000@results$delta.AICc == 0), 8]
CBI6000 <- panda6000@results[which(panda6000@results$delta.AICc == 0), 10]

## comparing multi scale model and single scale models
auc <- c(AUCmulti,AUC250, AUC500, AUC1000, AUC2000, AUC3000, AUC4000, AUC5000, AUC6000) 
CBI <- c(CBImulti, CBI250, CBI500, CBI1000, CBI2000, CBI3000, CBI4000, CBI5000, CBI6000)
AICcvalue <- c(
  min(pandaENM@results$AICc),
  min(panda250@results$AICc),
  min(panda500@results$AICc),
  min(panda1000@results$AICc),
  min(panda2000@results$AICc),
  min(panda3000@results$AICc),
  min(panda4000@results$AICc),
  min(panda5000@results$AICc),
  min(panda6000@results$AICc)) 
overlap <- c(
raster.overlap(pred, pred)$D,
raster.overlap(pred, pred250)$D,
raster.overlap(pred, pred500)$D,
raster.overlap(pred, pred1000)$D,
raster.overlap(pred, pred2000)$D,
raster.overlap(pred, pred3000)$D,
raster.overlap(pred, pred4000)$D,
raster.overlap(pred, pred5000)$D,
raster.overlap(pred, pred6000)$D)   # overlap between multiscale prediction and each single scale model prediction

# variable importance between multiscale model and single scale models
var.imp # multiscale model 
varimp250 <- ecospat.maxentvarimport(optmodel250, dfvar = optmodel250@absence, nperm = 5)
varimp500 <- ecospat.maxentvarimport(optmodel500, dfvar = optmodel500@absence, nperm = 5)
varimp1000 <- ecospat.maxentvarimport(optmodel1000, dfvar = optmodel1000@absence, nperm = 5)
varimp2000 <- ecospat.maxentvarimport(optmodel2000, dfvar = optmodel2000@absence, nperm = 5)
varimp3000 <- ecospat.maxentvarimport(optmodel3000, dfvar = optmodel3000@absence, nperm = 5)
varimp4000 <- ecospat.maxentvarimport(optmodel4000, dfvar = optmodel4000@absence, nperm = 5)
varimp5000 <- ecospat.maxentvarimport(optmodel5000, dfvar = optmodel5000@absence, nperm = 5)
varimp6000 <- ecospat.maxentvarimport(optmodel6000, dfvar = optmodel6000@absence, nperm = 5)
varimp.df <- cbind(var.imp, varimp250, varimp500, varimp1000, varimp2000, varimp3000, varimp4000, varimp5000, varimp6000)
write.csv(varimp.df, file = "variableimportance_for_single_and_multi_scale_model.csv")
#####################################################################################################