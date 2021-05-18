library(Hmisc)
library(ggplot2)
library(lme4)
library(gtools)
library(glmnet)

##rm(list=ls())
load("processed.rda")
## The processed data contains the following objects:
## NSB.v11, a 993x106 dimensional matrix of selected nasal
## gene expressions;
## CD19.v11, a 662x66 dimensional matrix of selected gene
## expressions obtained from purified CD19 cells;
## CD4.v11, a 454x92 dimensional matrix of selected gene
## expressions obtained from purified CD4 cells;
## CD8.v12, a 333x51 dimensional matrix of selected gene
## expressions obtained from purified CD8 cells;
## OTU.v11, a 15x104 dimensional matrix of relative abundance of
## selected OTUs from the nasal microbiome data
## ss, a 139x3 dimensional matrix that contains participant ID, visit
## age, and Global respiratory Severity Scores (GRSS) of all subjects.

## In this integrative analysis (CD4 + CD8 + CD19), all three sets of
## data were summarized by top principal components in PCA. The number
## of PCs were determined by a prespecified threshold of the
## proportion of variance explained by these PCs. Three choices of
## such thresholds were considered: 70%, 80%, and 90%. The third
## choice (90% of explained variance) achieved the best prediction
## accuracy therefore was selected and reported in the manuscript.
varprop.thresholds <- c(0.7, 0.8, 0.9)
for (s in 1:3)
  { 
    threshold<-varprop.thresholds[s]
    fit.CD8 <- prcomp(t(CD8.v12), scale=T) 
    cv<-round(cumsum(fit.CD8$sdev^2/sum(fit.CD8$sdev^2)), 3)
    n.cut<-length(cv[cv<threshold])+1  
    CD8.rtt<-fit.CD8$rotation
    if (s==3){CD8.pc.2<-fit.CD8$x[, 1:n.cut] }
    CD8.pc<-fit.CD8$x[, 1:n.cut] 
    
    fit.CD4 <- prcomp(t(CD4.v11), scale=T) 
    cv<-round(cumsum(fit.CD4$sdev^2/sum(fit.CD4$sdev^2)), 3)
    n.cut<-length(cv[cv<threshold])+1
    CD4.rtt<-fit.CD4$rotation
    
    if (s==3){CD4.pc.2<-fit.CD4$x[, 1:n.cut]}
    CD4.pc<-fit.CD4$x[, 1:n.cut]  
    
    fit.CD19<- prcomp(t(CD19.v11), scale=T)
    cv<-round(cumsum(fit.CD19$sdev^2/sum(fit.CD19$sdev^2)), 3)
    n.cut<-length(cv[cv<threshold])+1
    CD19.rtt<-fit.CD19$rotation
    
    if (s==3){CD19.pc.2<-fit.CD19$x[, 1:n.cut]}
    CD19.pc<-fit.CD19$x[, 1:n.cut] 
    
    comm.id<-intersect(intersect(rownames(CD4.pc), rownames(CD8.pc)), rownames(CD19.pc))
    CD4.pc.t<-t(CD4.pc[comm.id, ]); rownames(CD4.pc.t)<-paste("CD4", rownames(CD4.pc.t), sep=".")
    CD8.pc.t<-t(CD8.pc[comm.id, ]); rownames(CD8.pc.t)<-paste("CD8", rownames(CD8.pc.t), sep=".")
    CD19.pc.t<-t(CD19.pc[comm.id, ]); rownames(CD19.pc.t)<-paste("CD19", rownames(CD19.pc.t), sep=".")
    all.pc<-rbind(CD4.pc.t, CD19.pc.t, CD8.pc.t)
    all.pc.t<-t(all.pc)
    
    ss.35<-ss[which(ss$ParticipantId %in% rownames(all.pc.t)), ]
    GRSS<-ss.35$GRSS
    
    ## we fixed the random seed to ensure that the results obtained
    ## from cv.glmnet() are reproducible
    set.seed(42)
    MODEL<-cv.glmnet(all.pc.t, GRSS, alpha=0.9)
    best.lambda.1se<-MODEL$lambda.1se 
    best.lambda.min<-MODEL$lambda.min
    
    bestcf.1se <- coef(MODEL, s=best.lambda.1se)
    bestcf.min <- coef(MODEL, s=best.lambda.min)
    
    sel.PC.1se<-rownames(bestcf.1se)[as.numeric(bestcf.1se) !=0 ][-1] 
    sel.PC.min<-rownames(bestcf.min)[as.numeric(bestcf.min) !=0 ][-1] 
    
    mod.name.1se <- formula(paste("GRSS ~", paste(sel.PC.1se, collapse=" + ")))
    fit.2 <-lm(mod.name.1se, data=as.data.frame(all.pc.t))
    summary(fit.2)
      
    mod.name.min <- formula(paste("GRSS ~", paste(sel.PC.min, collapse=" + ")))
    fit.3 <-lm(mod.name.min, data=as.data.frame(all.pc.t))
    summary(fit.3)
    
    if (s==3) 
    {
      selected.PC.min<-sel.PC.min
      ## save the glmnet predictive model for later use
      mod3.glmnet <- fit.3
      beta<-summary(fit.3)$coefficients[, "Estimate"]
      beta.cd4<-beta[grepl("CD4", names(beta)) ]
      beta.cd19<-beta[grepl("CD19", names(beta)) ]
      beta.cd8<-beta[grepl("CD8", names(beta)) ]
    }
    
    N <- nrow(all.pc.t)
  
    nmod2 <- formula(paste("grss ~", paste(sel.PC.1se, collapse=" + ")))
    nmod3 <- formula(paste("grss ~", paste(sel.PC.min, collapse=" + ")))
 
    yhat2 <- yhat3 <- rep(0, N)
    for (j in 1:N){
      grss <- GRSS[-j]
      dat.training <- cbind(as.data.frame(all.pc.t[-j, ]), grss)
  
      fit2.train <- lm(nmod2, data=dat.training)
      fit3.train <- lm(nmod3, data=dat.training)
      cc2 <- coef(fit2.train)
      cc3 <- coef(fit3.train)
      
      yhat2[j] <- sum(cc2 * c(1, unlist(all.pc.t[j, sel.PC.1se])))
      yhat3[j] <- sum(cc3 * c(1, unlist(all.pc.t[j, sel.PC.min])))
      
    }
    
    sumtab <- matrix(0, 2, 7)
    colnames(sumtab) <- c("number of selected", "Naive RSS", "Naive Cor.", "Naive misclass.",
                          "CV RSS", "CV Cor.", "CV misclass.")
    rownames(sumtab) <- c("Model 2 (1se)", "Model 3 (min)")
    n<-c(length(sel.PC.1se), length(sel.PC.min))
    GRSS.bin <- ifelse(GRSS<=3.5, "Mild", "Severe")
    for (k in 2:3){
      GRSS.pred.naive <- predict(get(paste0("fit.", k)))
      GRSS.pred.cv <- get(paste0("yhat", k))
      sumtab[k-1, "number of selected"] <- n[k-1] 
      sumtab[k-1, "Naive RSS"] <- sum((GRSS - GRSS.pred.naive)^2)/(N-1)
      sumtab[k-1, "CV RSS"] <- sum((GRSS - GRSS.pred.cv)^2)/(N-1)
      sumtab[k-1, "Naive Cor."] <- cor(GRSS, GRSS.pred.naive)
      sumtab[k-1, "CV Cor."] <- cor(GRSS, GRSS.pred.cv)
      GRSS.pred.bin.naive <- ifelse(GRSS.pred.naive<=3.5, "Mild", "Severe")
      GRSS.pred.bin.cv <- ifelse(GRSS.pred.cv<=3.5, "Mild", "Severe")
      sumtab[k-1, "Naive misclass."] <- N-sum(diag(table(GRSS.bin, GRSS.pred.bin.naive)))
      sumtab[k-1, "CV misclass."] <- N-sum(diag(table(GRSS.bin, GRSS.pred.bin.cv)))
      
      if (s==3 & k==3){
        dat.result<-data.frame(PID=names(GRSS.pred.naive), GRSS=as.numeric(GRSS), 
                               GRSS.pred.naive.integrated=as.numeric(GRSS.pred.naive), 
                               GRSS.pred.cv.integrated=as.numeric(GRSS.pred.cv),
                               GRSS.pred.naive.NT=as.numeric(rep(-99, length(GRSS))), 
                               GRSS.pred.cv.NT=as.numeric(rep(-99, length(GRSS))), 
                               GRSS.pred.naive.CD4=as.numeric(rep(-99, length(GRSS))), 
                               GRSS.pred.cv.CD4=as.numeric(rep(-99, length(GRSS))),
                               GRSS.pred.naive.OTU=as.numeric(rep(-99, length(GRSS))), 
                               GRSS.pred.cv.OTU=as.numeric(rep(-99, length(GRSS))))
      }
    }
    Name <- paste("tab", s, sep=".")
    assign(Name, sumtab)
}
## The results of the above for loop are stored in the following three
## tables: tab.1 (the initial PCA explains about 70% of total
## variance), tab.2 (80% of variance; this choice has the best
## cross-validated prediction accuracy), and tab.3 (90% of variance).
  

######################################################################
## Predictive models using only an individual type of data. To be
## consistent with the integrative prective model, we selected the
## number of PCs based on 90% of explained variance for all predictive
## models based on individual type of data.
######################################################################

###########################
## individual data type CD4
##########################
CD4.pc <- CD4.pc.2[comm.id , ]
GRSS<-ss.35$GRSS
set.seed(42)
MODEL<-cv.glmnet(CD4.pc, GRSS, alpha=0.9)
best.lambda.min<-MODEL$lambda.min
bestcf.min <- coef(MODEL, s=best.lambda.min)
sel.PC.min<-rownames(bestcf.min)[as.numeric(bestcf.min) !=0 ][-1] 
mod.name.min <- formula(paste("GRSS ~", paste(sel.PC.min, collapse=" + ")))
fit.3 <-lm(mod.name.min, data=as.data.frame(CD4.pc))

N <- nrow(CD4.pc)
nmod3 <- formula(paste("grss ~", paste(sel.PC.min, collapse=" + ")))
yhat2 <-  rep(0, N)
for (j in 1:N){
  grss <- GRSS[-j]
  dat.training <- cbind(as.data.frame(CD4.pc[-j, ]), grss)
  fit3.train <- lm(nmod3, data=dat.training)
  cc2 <- coef(fit3.train)
  yhat2[j] <- sum(cc2 * c(1, unlist(CD4.pc[j, sel.PC.min])))
}

GRSS.pred.naive <- predict(fit.3)
GRSS.pred.cv <- yhat2
dat.result[, "GRSS.pred.naive.CD4"]<-GRSS.pred.naive
dat.result[, "GRSS.pred.cv.CD4"]<-GRSS.pred.cv


###########################
## individual data type CD8
##########################
CD8.pc <- CD8.pc.2[comm.id , ]
GRSS<-ss.35$GRSS
set.seed(42)
MODEL<-cv.glmnet(CD8.pc, GRSS, alpha=0.9)
best.lambda.min<-MODEL$lambda.min
bestcf.min <- coef(MODEL, s=best.lambda.min)
sel.PC.min<-rownames(bestcf.min)[as.numeric(bestcf.min) !=0 ][-1] #16 LPC
mod.name.min <- formula(paste("GRSS ~", paste(sel.PC.min, collapse=" + ")))
fit.3 <-lm(mod.name.min, data=as.data.frame(CD8.pc))

N <- nrow(CD8.pc)
nmod3 <- formula(paste("grss ~", paste(sel.PC.min, collapse=" + ")))

yhat2 <-  rep(0, N)
for (j in 1:N){
  grss <- GRSS[-j]
  dat.training <- cbind(as.data.frame(CD8.pc[-j, ]), grss)
  fit3.train <- lm(nmod3, data=dat.training)
  cc2 <- coef(fit3.train)
  yhat2[j] <- sum(cc2 * c(1, unlist(CD8.pc[j, sel.PC.min])))
}

GRSS.pred.naive <- predict(fit.3)
GRSS.pred.cv <- yhat2
dat.result[, "GRSS.pred.naive.CD8"]<-GRSS.pred.naive
dat.result[, "GRSS.pred.cv.CD8"]<-GRSS.pred.cv


###########################
## individual data type CD19
##########################
CD19.pc <- CD19.pc.2[comm.id , ]
GRSS<-ss.35$GRSS
set.seed(42)
MODEL<-cv.glmnet(CD19.pc, GRSS, alpha=0.9)
best.lambda.min<-MODEL$lambda.min
bestcf.min <- coef(MODEL, s=best.lambda.min)
sel.PC.min<-rownames(bestcf.min)[as.numeric(bestcf.min) !=0 ][-1] #16 LPC

mod.name.min <- formula(paste("GRSS ~", paste(sel.PC.min, collapse=" + ")))
fit.3 <-lm(mod.name.min, data=as.data.frame(CD19.pc))

N <- nrow(CD19.pc)
nmod3 <- formula(paste("grss ~", paste(sel.PC.min, collapse=" + ")))

yhat2 <-  rep(0, N)
for (j in 1:N){
  grss <- GRSS[-j]
  dat.training <- cbind(as.data.frame(CD19.pc[-j, ]), grss)
  fit3.train <- lm(nmod3, data=dat.training)
  cc2 <- coef(fit3.train)
  yhat2[j] <- sum(cc2 * c(1, unlist(CD19.pc[j, sel.PC.min])))
}

GRSS.pred.naive <- predict(fit.3)
GRSS.pred.cv <- yhat2
dat.result[, "GRSS.pred.naive.CD19"]<-GRSS.pred.naive
dat.result[, "GRSS.pred.cv.CD19"]<-GRSS.pred.cv


dat.result.1<-dat.result
dat.result.1[, "resi.cv.integrated"]<-dat.result.1[, "GRSS.pred.cv.integrated"]-dat.result.1[, "GRSS"]
dat.result.1[, "resi.cv.CD4"]<-dat.result.1[, "GRSS.pred.cv.CD4"]-dat.result.1[, "GRSS"]
dat.result.1[, "resi.cv.CD8"]<-dat.result.1[, "GRSS.pred.cv.CD8"]-dat.result.1[, "GRSS"]
dat.result.1[, "resi.cv.CD19"]<-dat.result.1[, "GRSS.pred.cv.CD19"]-dat.result.1[, "GRSS"]


dat.resi.1<-data.frame(residual=c(dat.result.1$resi.cv.integrated, dat.result.1$resi.cv.CD4, 
                                  dat.result.1$resi.cv.CD19, dat.result.1$resi.cv.CD8) , 
                       type=c(rep("Integrated", nrow(dat.result.1)), rep("CD4", nrow(dat.result.1)), 
                              rep("CD19", nrow(dat.result.1)), rep("CD8", nrow(dat.result.1))   ),
                       o=c(rep(1, nrow(dat.result.1)), rep(2, nrow(dat.result.1)), rep(3, nrow(dat.result.1)),
                           rep(4, nrow(dat.result.1))))

ggplot(dat.resi.1, aes(x=factor(o), y=residual, fill=type)) + 
  scale_y_continuous(breaks=c(-5.5, -3.5, -1.5, 0, 1.5, 3.5, 5.5))+
  scale_x_discrete(breaks=c("1","2","3", "4"),
                   labels=c("Integrated", "CD4", "CD19", "CD8"))+
  geom_boxplot(alpha=0.75, outlier.size = 0, width=0.7) + 
  scale_fill_manual(values=c("#DA9694", "#B8CCE4", "#D8E4BC", "white"))+
  theme_bw()+
  geom_point(position=position_jitter(width=0.1))+
  theme(axis.text.x=element_text(size=13))+
  theme(axis.text.y=element_text(size=13))+
  xlab(" ")+
  ylab(paste("Predicted \n - Observed GRSS"))+
  theme(axis.title.y = element_text(size = rel(1.5), vjust=1.8))+
  theme(axis.text.x = element_text(hjust = 0.5, size=rel(1.5)))


######################################################################
## Computing the weights for genes based on a backpropagation
## algorithm. The absolute values of See Supplemental Text, section
## "Feature Weight Calculation" for more details.
######################################################################

## note that the fitted linear coefficients related to NT and CD4 are
## stored in "beta.nt" and "beta.cd4", respectively.
pc.cd4 <- sapply(names(beta.cd4), function(s) strsplit(s, "\\.")[[1]][2])
pc.cd8 <- sapply(names(beta.cd8), function(s) strsplit(s, "\\.")[[1]][2])
pc.cd19 <- sapply(names(beta.cd19), function(s) strsplit(s, "\\.")[[1]][2])

## 1. CD4 weights
fit.CD4 <- prcomp(t(CD4.v11), scale=T)
means <- rowMeans(CD4.v11)
sds <- apply(CD4.v11, 1, sd)
CD4.param <- diag(1/sds) %*% fit.CD4$rotation
CD4.weights <- drop(CD4.param[, pc.cd4] %*% beta.cd4)
dat.CD4.weights <- data.frame(gene.id=rownames(fit.CD4$rotation), means=means, weights=CD4.weights)

## 2. CD8 weights
fit.CD8 <- prcomp(t(CD8.v12), scale=T)
means <- rowMeans(CD8.v12)
sds <- apply(CD8.v12, 1, sd)
CD8.param <- diag(1/sds) %*% fit.CD8$rotation
CD8.weights <- drop(CD8.param[, pc.cd8] %*% beta.cd8)
dat.CD8.weights <- data.frame(gene.id=rownames(fit.CD8$rotation), means=means, weights=CD8.weights)

## 3. CD19 weights
fit.CD19 <- prcomp(t(CD19.v11), scale=T)
means <- rowMeans(CD19.v11)
sds <- apply(CD19.v11, 1, sd)
CD19.param <- diag(1/sds) %*% fit.CD19$rotation
CD19.weights <- drop(CD19.param[, pc.cd19] %*% beta.cd19)
dat.CD19.weights <- data.frame(gene.id=rownames(fit.CD19$rotation), means=means, weights=CD19.weights)


######################################################################
## check the equivalence between the PCs and these weights. Note that we
## also need to compensate for the centering effects.
######################################################################

## 1. GRSS predicted by the glmnet model (based on PCs)
beta0 <- coef(mod3.glmnet)["(Intercept)"]
GRSS.pred1 <- predict(mod3.glmnet)

## Now let us compute the predicted GRSS from the equivalent
## weights. Note that we need to center NSB and CD4 first; then select
## those 61 subjects that have NSB, CD4, and OTU data
CD4c <- sweep(CD4.v11, 1, rowMeans(CD4.v11))[, ss.35$ParticipantId]
CD8c <- sweep(CD8.v12, 1, rowMeans(CD8.v12))[, ss.35$ParticipantId]
CD19c <- sweep(CD19.v11, 1, rowMeans(CD19.v11))[, ss.35$ParticipantId]
GRSS.pred2 <- beta0 +drop(t(CD4c)%*%CD4.weights +t(CD8c)%*%CD8.weights +t(CD19c)%*%CD19.weights)

## check if the two approaches are equivalent:
all.equal(GRSS.pred2,GRSS.pred1) #TRUE

