## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup, include = FALSE---------------------------------------------------
library(ClustImpute)
required_packages <- Hmisc::Cs(psych,ggplot2,tidyr,Hmisc,tictoc,ClusterR,copula,dplyr,corrplot)
lapply(required_packages, require, character.only = TRUE)

## -----------------------------------------------------------------------------
### Random Dataset
set.seed(739)
n <- 7500 # numer of points
nr_other_vars <- 4
mat <- matrix(rnorm(nr_other_vars*n),n,nr_other_vars)
me<-4 # mean
x <- c(rnorm(n/3,me/2,1),rnorm(2*n/3,-me/2,1)) 
y <- c(rnorm(n/3,0,1),rnorm(n/3,me,1),rnorm(n/3,-me,1))
true_clust <- c(rep(1,n/3),rep(2,n/3),rep(3,n/3)) # true clusters
dat <- cbind(mat,x,y)
dat<- as.data.frame(scale(dat)) # scaling
summary(dat)

## -----------------------------------------------------------------------------
library(ggExtra)
dat4plot <- dat
dat4plot$true_clust_fct <- factor(true_clust)
p_base <- ggplot(dat4plot,aes(x=x,y=y,color=true_clust_fct)) + geom_point()
ggMarginal(p_base, groupColour = TRUE, groupFill = TRUE)

## -----------------------------------------------------------------------------
dat_with_miss <- miss_sim(dat,p=.2,seed_nr=120)
summary(dat_with_miss)
mis_ind <- is.na(dat_with_miss) # missing indicator

## -----------------------------------------------------------------------------
corrplot(cor(mis_ind),method="number")

## -----------------------------------------------------------------------------
dat_median_imp <- dat_with_miss
for (j in 1:dim(dat)[2]) {
  dat_median_imp[,j] <- Hmisc::impute(dat_median_imp[,j],fun=median)
}
imp <- factor(pmax(mis_ind[,5],mis_ind[,6]),labels=c("Original","Imputed")) # point is imputed if x or y is imputed
p_median_imp <- ggplot(dat_median_imp) + geom_point(aes(x=x,y=y,color=imp))
ggMarginal(p_median_imp,groupColour = TRUE, groupFill = TRUE)

## -----------------------------------------------------------------------------
dat_random_imp <- dat_with_miss
for (j in 1:dim(dat)[2]) {
  dat_random_imp[,j] <- impute(dat_random_imp[,j],fun="random")
}
imp <- factor(pmax(mis_ind[,5],mis_ind[,6]),labels=c("Original","Imputed")) # point is imputed if x or y is imputed
p_random_imp <- ggplot(dat_random_imp) + geom_point(aes(x=x,y=y,color=imp))
ggMarginal(p_random_imp,groupColour = TRUE, groupFill = TRUE)

## -----------------------------------------------------------------------------
tic("Clustering based on random imputation")
cl_compare <- KMeans_arma(data=dat_random_imp,clusters=3,n_iter=100,seed=751)
toc()
dat_random_imp$pred <- predict_KMeans(dat_random_imp,cl_compare)
p_random_imp <- ggplot(dat_random_imp) + geom_point(aes(x=x,y=y,color=factor(pred)))
ggMarginal(p_random_imp,groupColour = TRUE, groupFill = TRUE)

## -----------------------------------------------------------------------------
nr_iter <- 10 # iterations of procedure
n_end <- 10 # step until convergence of weight function to 1
nr_cluster <- 3 # number of clusters
c_steps <- 50 # numer of cluster steps per iteration
tic("Run ClustImpute")
res <- ClustImpute(dat_with_miss,nr_cluster=nr_cluster, nr_iter=nr_iter, c_steps=c_steps, n_end=n_end) 
toc()

## ----eval=FALSE---------------------------------------------------------------
#  res
#  summary(res)
#  attributes(res)

## -----------------------------------------------------------------------------
p_clustimpute <- ggplot(res$complete_data,aes(x,y,color=factor(res$clusters))) + geom_point()
ggMarginal(p_clustimpute,groupColour = TRUE, groupFill = TRUE)

## ----fig.width=10, fig.height=7-----------------------------------------------
plot(res)+xlim(-2.5,2.5)

## ----fig.width=10, fig.height=7-----------------------------------------------
plot(res, type="box")

## -----------------------------------------------------------------------------
res2 <- ClustImpute(dat_with_miss,nr_cluster=nr_cluster, nr_iter=nr_iter, c_steps=c_steps, n_end=n_end,seed_nr = 2)
res3 <- ClustImpute(dat_with_miss,nr_cluster=nr_cluster, nr_iter=nr_iter, c_steps=c_steps, n_end=n_end,seed_nr = 3)
mean_all <- rbind(res$imp_values_mean,res2$imp_values_mean,res3$imp_values_mean)
sd_all <- rbind(res$imp_values_sd,res2$imp_values_sd,res3$imp_values_sd)

## -----------------------------------------------------------------------------
mean_all <- cbind(mean_all,seed=rep(c(150519,2,3),each=11))
sd_all <- cbind(sd_all,seed=rep(c(150519,2,3),each=11))

## -----------------------------------------------------------------------------
ggplot(as.data.frame(mean_all)) + geom_line(aes(x=iter,y=V1,color=factor(seed))) + ggtitle("Mean")
ggplot(as.data.frame(sd_all)) + geom_line(aes(x=iter,y=V1,color=factor(seed))) + ggtitle("Std. dev.")

## -----------------------------------------------------------------------------
external_validation(true_clust, res$clusters)

## -----------------------------------------------------------------------------
class(dat_random_imp$pred) <- "numeric"
external_validation(true_clust, dat_random_imp$pred)

## -----------------------------------------------------------------------------
## complete cases
idx <- which(complete.cases(dat_with_miss)==TRUE)
sprintf("Number of complete cases is %s",length(idx))
sprintf("Rand index for this case %s", external_validation(true_clust[idx], res$clusters[idx]))

## -----------------------------------------------------------------------------
external_validation(true_clust, res$clusters,summary_stats = TRUE)

## -----------------------------------------------------------------------------
res_var <- var_reduction(res)
res_var$Variance_reduction
res_var$Variance_by_cluster

## -----------------------------------------------------------------------------
res <- ClustImpute(dat_with_miss,nr_cluster=10, nr_iter=nr_iter, c_steps=c_steps, n_end=n_end)
res_var <- var_reduction(res)
res_var$Variance_reduction
res_var$Variance_by_cluster

## -----------------------------------------------------------------------------
ClustImpute2 <- function(dataFrame,nr_cluster, nr_iter=10, c_steps=1, wf=default_wf, n_end=10, seed_nr=150519) {
  return(ClustImpute(dataFrame,nr_cluster, nr_iter, c_steps, wf, n_end, seed_nr))
}
res_list <- lapply(X=1:10,FUN=ClustImpute2,dataFrame=dat_with_miss, nr_iter=nr_iter, c_steps=c_steps, n_end=n_end)

## -----------------------------------------------------------------------------
tmp <- var_reduction(res_list[[1]])
var_by_clust <- tmp$Variance_by_cluster
for (k in 2:10) {
  tmp <- var_reduction(res_list[[k]])
  var_by_clust <- rbind(var_by_clust,tmp$Variance_by_cluster)
}
var_by_clust$nr_clusters <- 1:10

## -----------------------------------------------------------------------------
data2plot <- tidyr::gather(var_by_clust,key = "variable", value = "variance", -dplyr::one_of("nr_clusters"))
ggplot(data2plot,aes(x=nr_clusters,y=variance,color=variable)) + geom_line() + scale_x_continuous(breaks=1:10)

