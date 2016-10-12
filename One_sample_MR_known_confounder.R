#Simulation example for one sample MR with known confounder

#Load used libraries
library(MASS)
library(ggplot2)
library(GA)
library(metafor)

##Simulation parameters
n_sample=10000
met_corr=0.35

#Generate x1* (metabolite 1) and x2* (Metabolite 2)
x1StarAndx2Star <- mvrnorm(n_sample, c(20, 15), matrix(c(1, met_corr, met_corr, 1), 2, 2))

#Simulate 20 independent SNPs with MAF 0.1-0.3 and 10,000 individuals
SNPs1<-t(matrix(rbinom(n_sample*4,2,0.10),ncol=n_sample))
SNPs2<-t(matrix(rbinom(n_sample*4,2,0.15),ncol=n_sample))
SNPs3<-t(matrix(rbinom(n_sample*4,2,0.20),ncol=n_sample))
SNPs4<-t(matrix(rbinom(n_sample*4,2,0.25),ncol=n_sample))
SNPs5<-t(matrix(rbinom(n_sample*4,2,0.30),ncol=n_sample))
#Merge all simualted SNPs matrices
SNPs<-as.data.frame(cbind(SNPs1,SNPs2,SNPs3,SNPs4,SNPs5))

#Generate x1 and x2
#x1 affected by all SNPs
x1 <- x1StarAndx2Star[,1]+rowSums(SNPs)+rnorm(nrow(x1StarAndx2Star),0,0.1)
#x2 affected by half of the SNPs (odd)
x2 <- x1StarAndx2Star[,2]+rowSums(SNPs[c(1,3,5,7,9,11,13,15,17,19)])+rnorm(nrow(x1StarAndx2Star),0,0.1)

#Generate the risk factor
y<-1+x1+x2+rnorm(n_sample, 0, 0.5)

#Run IV two-stage IV for all SNPs
#First step IV regress all SNPs to x1
IV1 <- lm(x1 ~ . , data = SNPs)
summary(IV1)
#Extract fitted values from IV1 model
x1Hat <- IV1$fitted.values
#Second step IV regress fitted values to y
mrx1<-lm(y ~ x1Hat)
summary(mrx1)

#####Genetic algorithm#####

#####Both target and confounder are observed

#Model 1 to be used for fitness estimation
mod1<-lm(x1~. , data=SNPs)
summary(mod1)
#Save x and y for model1
mod1_x <- model.matrix(mod1)[, -1]
mod1_y <- model.response(model.frame(mod1))

#Model 2 to be used for fitness estimation
mod2<-lm(x2~. , data=SNPs)
summary(mod2)
#Save x and y for model2
mod2_x <- model.matrix(mod2)[, -1]
mod2_y <- model.response(model.frame(mod2))


####Adjusted R-square difference####

#Fitness function R-squared difference
fitness_rqsdiff <- function(string) {
  #SNPs included in teh chromosomes
  inc <- which(string == 1)
  
  #Subset of SNPs to be used
  X1 <- cbind(1, mod1_x[,inc])
  #Rerun model for specific combination
  model1 <- lm(mod1_y ~ X1)
  #Extract Rsq
  rsq1<-summary(model1)$adj.r.squared
  
  #Subset of SNPs to be used
  X2 <- cbind(1, mod2_x[,inc])
  #Rerun model for specific combination
  model2 <- lm(mod2_y ~ X2)
  #Extract Rsq
  rsq2<-summary(model2)$adj.r.squared
  
  #Adjusted R square selection. Suitable after a subset of SNPs has
  #been selected for the target measure
  #return(rsq2)
  #Difference between adjusted R squares. Suitable to slect SNPs maximising
  #target measure but minimising confounder information
  return(rsq1-rsq2)
}

#Running the genetic algorithm
GA_rqsdiff <- ga("binary",fitness=fitness_rqsdiff,nBits=ncol(SNPs),names=colnames(SNPs),
                 popSize = 40,maxiter=100,monitor=plot)
#Summary of GA
summary(GA_rqsdiff)

#Run IV two-stage IV with selected SNPs
mod_fin1<-lm(x1~.,data=data.frame(mod1_x[,GA_rqsdiff@solution==1]))
x1Hat.GA_rqsdiff<-mod_fin1$fitted.values
mrx1.GA_rqsdiff<-lm(y ~ x1Hat.GA_rqsdiff)
summary(mrx1.GA_rqsdiff)


####F-statistic difference####

###With individual SNPs

#Fitness function F ratio
fitness_fratio <- function(string) {
  #SNPs included in the chromosomes
  inc <- which(string == 1)
  
  #Subset of SNPs to be used
  X1 <- cbind(1, mod1_x[,inc])
  #Rerun model for specific combination
  model1 <- lm(mod1_y ~ X1)
  #Extract F1
  f1<-summary(model1)$fstatistic[1]
  #If null set to 0
  f1<-ifelse(is.null(f1),NA,f1)
  
  #Subset of SNPs to be used
  X2 <- cbind(1, mod2_x[,inc])
  #Rerun model for specific combination
  model2 <- lm(mod2_y ~ X2)
  #Extract F1
  f2<-summary(model2)$fstatistic[1]
  #If null set to 0
  f2<-ifelse(is.null(f2),NA,f2)
  
  #F of overall models 1 and 2
  f1_all<-summary(mod1)$fstatistic[1]
  f2_all<-summary(mod2)$fstatistic[1]
  
  #Ratio of F-statistic between the two models
  #return(f1/f2)
  #Ratio of the F-statistic of the models scaled for F of full model
  return((f1/f1_all)/(f2/f2_all))
  
}

#Running the genetic algorithm
GA_fratio <- ga("binary",fitness=fitness_fratio,nBits=ncol(SNPs),names=colnames(SNPs),
                popSize = 40,maxiter=100,monitor=plot)
#Summary of GA
summary(GA_fratio)

#Run IV two-stage IV with selected SNPs
#Fitting selected model
mod_fin1<-lm(x1~.,data=data.frame(mod1_x[,GA_fratio@solution==1]))
x1Hat.GA_fratio<-mod_fin1$fitted.values
mrx1.GA_fratio<-lm(y ~ x1Hat.GA_fratio)
summary(mrx1.GA_fratio)


###With unweighted score

#Fitness function for F difference
fitness_fdiff <- function(string) {
  #SNPs included in the chromosomes
  inc <- which(string == 1)
  
  #Subset of SNPs to be used
  X1 <- cbind(1, mod1_x[,inc])
  #Generate unweighted gene score
  X1_score<-rowSums(X1)
  #Rerun model for specific combination
  model1 <- lm(mod1_y ~ X1_score)
  #Extract F1
  f1<-summary(model1)$fstatistic[1]
  #If null set to 0
  f1<-ifelse(is.null(f1),NA,f1)
  
  #Subset of SNPs to be used
  X2 <- cbind(1, mod2_x[,inc])
  #Generate unweighted gene score
  X2_score<-rowSums(X2)
  #Rerun model for specific combination
  model2 <- lm(mod2_y ~ X2_score)
  #Extract F1
  f2<-summary(model2)$fstatistic[1]
  #If null set to 0
  f2<-ifelse(is.null(f2),NA,f2)
  
  #F of overall models 1 and 2
  f1_all<-summary(mod1)$fstatistic[1]
  f2_all<-summary(mod2)$fstatistic[1]
  
  #Selecting SNPs maximising the ratio of F-statistics for the target
  #and confounding measures 
  return(f1/f2)
  
}

#Running the genetic algorithm
GA_fdiff <- ga("binary",fitness=fitness_fdiff,nBits=ncol(SNPs),names=colnames(SNPs),
               popSize = 40,maxiter=100,monitor=plot)
#Summary of GA
summary(GA_fdiff)

#Run IV two-stage IV with selected SNPs
mod_fin1<-lm(x1~.,data=data.frame(mod1_x[,GA_fdiff@solution==1]))
x1Hat.GA_fdiff<-mod_fin1$fitted.values
mrx1.GA_fdiff<-lm(y ~ x1Hat.GA_fdiff)
summary(mrx1)
summary(mrx1.GA_fdiff)
