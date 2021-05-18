library(Hmisc)
library(ggplot2)
library(lme4)
library(gtools)
library(glmnet)

rm(list=ls())
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

## In this integrative analysis (NT + CD4 + OTU), both gene expression
## data (NT and CD4) were summarized by top principal components in
## PCA. The number of PCs were determined by a prespecified threshold
## of the proportion of variance explained by these PCs. Three choices
## of such thresholds were considered: 70%, 80%, and 90%. The second
## choice (80% of explained variance) achieved the best prediction
## accuracy therefore was selected and reported in the manuscript.
varprop.thresholds <- c(0.7, 0.8, 0.9)
for (s in 1:3)
  {
    threshold<-varprop.thresholds[s]
    fit.NT <- prcomp(t(NSB.v11), scale=T) 
    cv<-round(cumsum(fit.NT$sdev^2/sum(fit.NT$sdev^2)), 3)
    n.cut<-length(cv[cv<threshold])+1
    NT.rtt<-fit.NT$rotation
    if (s==2) {NT.pc.2<-fit.NT$x[, 1:n.cut]}
    NT.pc<-fit.NT$x[, 1:n.cut] 
    
    fit.CD4 <- prcomp(t(CD4.v11), scale=T)
    cv<-round(cumsum(fit.CD4$sdev^2/sum(fit.CD4$sdev^2)), 3)
    n.cut<-length(cv[cv<threshold])+1
    CD4.rtt<-fit.CD4$rotation
    
    if (s==2){CD4.pc.2<-fit.CD4$x[, 1:n.cut]}
    CD4.pc<-fit.CD4$x[, 1:n.cut]  
    
    comm.id<-intersect(intersect(rownames(NT.pc), rownames(CD4.pc)), colnames(OTU.v11)) 
    NT.pc.t<-t(NT.pc[comm.id, ]); rownames(NT.pc.t)<-paste("NT", rownames(NT.pc.t), sep=".")
    CD4.pc.t<-t(CD4.pc[comm.id, ]); rownames(CD4.pc.t)<-paste("CD4", rownames(CD4.pc.t), sep=".")
    
    map <- match(colnames((NT.pc.t)), colnames(CD4.pc.t))
    NT.pc.t <- NT.pc.t[, map ]
    OTU.v11<-OTU.v11[, colnames(NT.pc.t)]
    
    all.pc<-rbind(NT.pc.t, CD4.pc.t, OTU.v11) 
    all.pc.t<-t(all.pc)
    
    ss.61<-ss[which(ss$ParticipantId %in% rownames(all.pc.t)), ]
    all.pc.t<-all.pc.t[ss.61$ParticipantId,]
    GRSS.61<-ss.61$GRSS
   
    ## we fixed the random seed to ensure that the results obtained
    ## from cv.glmnet() are reproducible
    set.seed(42)
    MODEL<-cv.glmnet(all.pc.t, GRSS.61, alpha=0.9)
    best.lambda.1se<-MODEL$lambda.1se  
    best.lambda.min<-MODEL$lambda.min
    bestcf.1se <- coef(MODEL, s=best.lambda.1se)
    bestcf.min <- coef(MODEL, s=best.lambda.min)
    
    sel.PC.1se<-rownames(bestcf.1se)[as.numeric(bestcf.1se) !=0 ][-1]
    sel.PC.min<-rownames(bestcf.min)[as.numeric(bestcf.min) !=0 ][-1] 
    
    set.seed(42)
    mod.name.1se <- formula(paste("GRSS.61 ~", paste(sel.PC.1se, collapse=" + ")))
    fit.2 <-lm(mod.name.1se, data=as.data.frame(all.pc.t))
    summary(fit.2)
    
    if (s==2) 
      {
       selected.PC.1se<-sel.PC.1se
       ## save the glmnet predictive model for later use
       mod2.glmnet <- fit.2
       beta<-summary(fit.2)$coefficients[, "Estimate"]
       ## we divide the fitted linear coefficients according to the
       ## sources of the data. "beta.nt", "beta.cd4", and "beta.otu"
       ## will be used for computing weights for the original features
       ## (e.g., genes and OTUs).
       beta.nt<-beta[grepl("NT", names(beta)) ]
       beta.cd4<-beta[grepl("CD4", names(beta)) ]
       beta.otu<-beta[grepl("g", names(beta)) ]
       }
    
    mod.name.min <- formula(paste("GRSS.61 ~", paste(sel.PC.min, collapse=" + ")))
    fit.3 <-lm(mod.name.min, data=as.data.frame(all.pc.t))
    summary(fit.3)
    
    N <- nrow(all.pc.t)
   
    nmod2 <- formula(paste("grss ~", paste(sel.PC.1se, collapse=" + ")))
    nmod3 <- formula(paste("grss ~", paste(sel.PC.min, collapse=" + ")))
    
    yhat2 <- yhat3 <- rep(0, N)
    for (j in 1:N){
      grss <- GRSS.61[-j]
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
    GRSS.bin <- ifelse(GRSS.61<=3.5, "Mild", "Severe")
    for (k in 2:3){
      GRSS.pred.naive <- predict(get(paste0("fit.", k)))
      GRSS.pred.cv <- get(paste0("yhat", k))
      sumtab[k-1, "number of selected"] <- n[k-1] 
      sumtab[k-1, "Naive RSS"] <- sum((GRSS.61 - GRSS.pred.naive)^2)/(N-1)
      sumtab[k-1, "CV RSS"] <- sum((GRSS.61 - GRSS.pred.cv)^2)/(N-1)
      sumtab[k-1, "Naive Cor."] <- cor(GRSS.61, GRSS.pred.naive)
      sumtab[k-1, "CV Cor."] <- cor(GRSS.61, GRSS.pred.cv)
      GRSS.pred.bin.naive <- ifelse(GRSS.pred.naive<=3.5, "Mild", "Severe")
      GRSS.pred.bin.cv <- ifelse(GRSS.pred.cv<=3.5, "Mild", "Severe")
      sumtab[k-1, "Naive misclass."] <- N-sum(diag(table(GRSS.bin, GRSS.pred.bin.naive)))
      sumtab[k-1, "CV misclass."] <- N-sum(diag(table(GRSS.bin, GRSS.pred.bin.cv)))
      
      if (s==2 & k==2){
        dat.result<-data.frame(PID=names(GRSS.pred.naive), GRSS=as.numeric(GRSS.61), 
                               GRSS.pred.naive.integrated=as.numeric(GRSS.pred.naive), 
                               GRSS.pred.cv.integrated=as.numeric(GRSS.pred.cv),
                               GRSS.pred.naive.NT=as.numeric(rep(-99, length(GRSS.61))), 
                               GRSS.pred.cv.NT=as.numeric(rep(-99, length(GRSS.61))), 
                               GRSS.pred.naive.CD4=as.numeric(rep(-99, length(GRSS.61))), 
                               GRSS.pred.cv.CD4=as.numeric(rep(-99, length(GRSS.61))),
                               GRSS.pred.naive.OTU=as.numeric(rep(-99, length(GRSS.61))), 
                               GRSS.pred.cv.OTU=as.numeric(rep(-99, length(GRSS.61))))
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
## number of PCs based on 80% of explained variance for all predictive
## models based on individual type of data.
######################################################################

###########################
## individual data: NSB
###########################

NT.pc<-NT.pc.2[rownames(NT.pc.2) %in% comm.id , ] 
ss<-ss[which(ss$ParticipantId %in% rownames(NT.pc)), ]
NT.pc<-NT.pc[ss$ParticipantId, ]
GRSS<-ss$GRSS
set.seed(42)
MODEL<-cv.glmnet(NT.pc, GRSS, alpha=0.9)
best.lambda.1se<-MODEL$lambda.1se  
bestcf.1se <- coef(MODEL, s=best.lambda.1se)
sel.PC.1se<-rownames(bestcf.1se)[as.numeric(bestcf.1se) !=0 ][-1] 
mod.name.1se <- formula(paste("GRSS ~", paste(sel.PC.1se, collapse=" + ")))
fit.2 <-lm(mod.name.1se, data=as.data.frame(NT.pc))

N <- nrow(NT.pc)
nmod2 <- formula(paste("grss ~", paste(sel.PC.1se, collapse=" + ")))
yhat2 <-  rep(0, N)
for (j in 1:N){
  grss <- GRSS[-j]
  dat.training <- cbind(as.data.frame(NT.pc[-j, ]), grss)
  fit2.train <- lm(nmod2, data=dat.training)
  cc2 <- coef(fit2.train)
  yhat2[j] <- sum(cc2 * c(1, unlist(NT.pc[j, sel.PC.1se])))
}

GRSS.pred.naive <- predict(fit.2)
GRSS.pred.cv <- yhat2
dat.result[, "GRSS.pred.naive.NT"]<-GRSS.pred.naive
dat.result[, "GRSS.pred.cv.NT"]<-GRSS.pred.cv


###########################
## individual data type CD4
##########################

CD4.pc<-CD4.pc.2[rownames(CD4.pc) %in% comm.id , ] 
ss<-ss[which(ss$ParticipantId %in% rownames(CD4.pc)), ]
CD4.pc<-CD4.pc[ss$ParticipantId, ]
GRSS<-ss$GRSS
set.seed(42)
MODEL<-cv.glmnet(CD4.pc, GRSS, alpha=0.9)
best.lambda.1se<-MODEL$lambda.1se  
bestcf.1se <- coef(MODEL, s=best.lambda.1se)
sel.PC.1se<-rownames(bestcf.1se)[as.numeric(bestcf.1se) !=0 ][-1] 
mod.name.1se <- formula(paste("GRSS ~", paste(sel.PC.1se, collapse=" + ")))
fit.2 <-lm(mod.name.1se, data=as.data.frame(CD4.pc))

N <- nrow(CD4.pc)
nmod2 <- formula(paste("grss ~", paste(sel.PC.1se, collapse=" + ")))
yhat2 <-  rep(0, N)
for (j in 1:N){
  grss <- GRSS[-j]
  dat.training <- cbind(as.data.frame(CD4.pc[-j, ]), grss)
  fit2.train <- lm(nmod2, data=dat.training)
  cc2 <- coef(fit2.train)
  
  yhat2[j] <- sum(cc2 * c(1, unlist(CD4.pc[j, sel.PC.1se])))
}

GRSS.pred.naive <- predict(fit.2)
GRSS.pred.cv <- yhat2
dat.result[, "GRSS.pred.naive.CD4"]<-GRSS.pred.naive
dat.result[, "GRSS.pred.cv.CD4"]<-GRSS.pred.cv


############################
## individual data type OTU
############################
OTU.v11<-t(OTU.v11)
OTU.v11<-OTU.v11[rownames(OTU.v11) %in% comm.id , ] 
ss<-ss[which(ss$ParticipantId %in% rownames(OTU.v11)), ]
OTU.v11<-OTU.v11[ss$ParticipantId, ]
GRSS<-ss$GRSS
set.seed(42)
MODEL<-cv.glmnet(OTU.v11, GRSS, alpha=0.9)
best.lambda.1se<-MODEL$lambda.1se  
bestcf.1se <- coef(MODEL, s=best.lambda.1se)
sel.PC.1se<-rownames(bestcf.1se)[as.numeric(bestcf.1se) !=0 ][-1] 
mod.name.1se <- formula(paste("GRSS ~", paste(sel.PC.1se, collapse=" + ")))
fit.2 <-lm(mod.name.1se, data=as.data.frame(OTU.v11))

N <- nrow(OTU.v11)
nmod2 <- formula(paste("grss ~", paste(sel.PC.1se, collapse=" + ")))
yhat2 <-  rep(0, N)
for (j in 1:N){
  grss <- GRSS[-j]
  dat.training <- cbind(as.data.frame(OTU.v11[-j, ]), grss)
  fit2.train <- lm(nmod2, data=dat.training)
  cc2 <- coef(fit2.train)
  yhat2[j] <- sum(cc2 * c(1, unlist(OTU.v11[j, sel.PC.1se])))
}

GRSS.pred.naive <- predict(fit.2)
GRSS.pred.cv <- yhat2
dat.result[, "GRSS.pred.naive.OTU"]<-GRSS.pred.naive
dat.result[, "GRSS.pred.cv.OTU"]<-GRSS.pred.cv

## Visualization
dat.result.1<-dat.result
dat.result.1[, "resi.cv.integrated"]<-dat.result.1[, "GRSS.pred.cv.integrated"]-dat.result.1[, "GRSS"]
dat.result.1[, "resi.cv.NT"]<-dat.result.1[, "GRSS.pred.cv.NT"]-dat.result.1[, "GRSS"]
dat.result.1[, "resi.cv.CD4"]<-dat.result.1[, "GRSS.pred.cv.CD4"]-dat.result.1[, "GRSS"]
dat.result.1[, "resi.cv.OTU"]<-dat.result.1[, "GRSS.pred.cv.OTU"]-dat.result.1[, "GRSS"]

dat.resi.1<-data.frame(residual=c(dat.result.1$resi.cv.integrated, dat.result.1$resi.cv.NT,
                                  dat.result.1$resi.cv.CD4, dat.result.1$resi.cv.OTU), 
                       type=c(rep("Integrated", nrow(dat.result.1)),  rep("NT", nrow(dat.result.1)),
                              rep("CD4", nrow(dat.result.1)), rep("OTU", nrow(dat.result.1))   ),
                       o=c(rep(1, nrow(dat.result.1)), rep(2, nrow(dat.result.1)), rep(3, nrow(dat.result.1)),
                           rep(4, nrow(dat.result.1))))

ggplot(dat.resi.1, aes(x=factor(o), y=residual, fill=type)) + 
  scale_y_continuous(breaks=c(-7.5, -5, -2.5, 0, 2.5, 5))+
  scale_x_discrete(breaks=c("1","2", "3", "4"),
                   labels=c("Integrated", "NT", "CD4", "OTU"))+
  geom_boxplot(alpha=0.75, outlier.size = 0, width=0.7) + 
  scale_fill_manual(values=c("#B8CCE4", "white", "#B1A0C7", "#BF8E5C"  ))+
  theme_bw()+
  geom_point(position=position_jitter(width=0.1))+
  theme(axis.text.x=element_text(size=13))+
  theme(axis.text.y=element_text(size=13))+
  xlab(" ")+
  ylab(paste("Predicted \n - Observed GRSS"))+
  theme(axis.title.y = element_text(size = rel(1.5), vjust=1.8))+
  theme(axis.text.x = element_text(hjust = 0.5, size=rel(1.5)))

 
######################################################################
## Computing the weights for genes/OTUs based on a backpropagation
## algorithm. The absolute values of See Supplemental Text, section
## "Feature Weight Calculation" for more details.
######################################################################

## note that the fitted linear coefficients related to NT and CD4 are
## stored in "beta.nt" and "beta.cd4", respectively.
pc.nt <- sapply(names(beta.nt), function(s) strsplit(s, "\\.")[[1]][2])
pc.cd4 <- sapply(names(beta.cd4), function(s) strsplit(s, "\\.")[[1]][2])

## 1. Nasal transcriptome weights.
## Note that our PCA is applied to centered and SCALED data. So we
## need to compensate the scaling effects in computing weights.
fit.NT <- prcomp(t(NSB.v11), scale=T)
means <- rowMeans(NSB.v11)
sds <- apply(NSB.v11, 1, sd)
NT.param <- diag(1/sds) %*% fit.NT$rotation
NT.weights <- drop(NT.param[, pc.nt] %*% beta.nt)
dat.NT.weights <- data.frame(gene.id=rownames(fit.NT$rotation), means=means, weights=NT.weights)

## 2. CD4 weights
fit.CD4 <- prcomp(t(CD4.v11), scale=T)
means <- rowMeans(CD4.v11)
sds <- apply(CD4.v11, 1, sd)
CD4.param <- diag(1/sds) %*% fit.CD4$rotation
CD4.weights <- drop(CD4.param[, pc.cd4] %*% beta.cd4)
dat.CD4.weights <- data.frame(gene.id=rownames(fit.CD4$rotation), means=means, weights=CD4.weights)

## 3. OTU weights. No centering was applied to the OTU data
OTU.weights <- beta.otu
dat.OTU.weights <- data.frame(OTU=names(beta.otu), means=0, weights=OTU.weights)

######################################################################
## check the equivalence between the PCs and these weights. Note that we
## also need to compensate for the centering effects.
######################################################################

## 1. GRSS predicted by the glmnet model (based on PCs)
beta0 <- coef(mod2.glmnet)["(Intercept)"]
GRSS.pred1 <- predict(mod2.glmnet)

## Now let us compute the predicted GRSS from the equivalent
## weights. Note that we need to center NSB and CD4 first; then select
## those 61 subjects that have NSB, CD4, and OTU data
NSBc <- sweep(NSB.v11, 1, rowMeans(NSB.v11))[, ss.61$ParticipantId]
CD4c <- sweep(CD4.v11, 1, rowMeans(CD4.v11))[, ss.61$ParticipantId]
## select and transpose the OTU data
OTU <- t(OTU.v11)[names(OTU.weights), ss.61$ParticipantId]
GRSS.pred2 <- beta0 +drop(t(NSBc)%*%NT.weights +t(CD4c)%*%CD4.weights +t(OTU)%*%OTU.weights)

## check if the two approaches are equivalent:
all.equal(GRSS.pred2,GRSS.pred1) #TRUE

