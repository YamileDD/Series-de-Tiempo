
install.packages(c("ggplot2", "ggpubr", "tidyverse", "broom", "AICcmodavg"))
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)

Detergentes<-matrix(c(77,72,76,81,58,85,71,74,82,76,66,80,80,70,77), nrow=5, byrow=T)
colnames(Detergentes)<-c("A","B","C")
rownames(Detergentes)<-c("P1","P2","P3","P4","P5")
Detergentes
apply(Detergentes,2, FUN=mean)


blancura<-c(77,81,71,76,80,72,58,74,66,70,76,85,82,80,77)
detergentes<-as.factor(c(rep(c("A","B","C"), each=5)))
boxplot(blancura~detergentes, col=c("red", "blue","yellow"), ylab= "% de blancura")
tapply(blancura, detergentes, mean)

ANOVA<-aov(lm(blancura~detergentes))
summary(ANOVA)
Fcrit<-qf(0.01, 3-1, 15-3, lower.tail = F) #alfa, k-1, N-k
#Fcrit<-qf(0.01, 3-1, 15-3) 
Fcrit

