#Ex1. Seminario de Investigacion 25 Enero 2022

# Lectura de la base de datos desde Excel
combustible<-read.csv("C:\\Users\\Areliuuu\\Dropbox\\PC\\Desktop\\Seminario de Investigacion I\\Ejemplo1_Consumo_de_Combustible.csv")

#Grafico de dispersion de Temperatura vs Combustible
plot(combustible$Temperatura, combustible$Combustible)
plot(combustible)

#Regresion lm=ajusta modelos lineales
regresion<-lm(formula=Combustible~Temperatura, data=combustible)#y(Comb)=x(Temp); dataframe=archivo
regresion

summary(regresion) #Valor min, quartil, valorp, t, etc/ ver asteriscos
#R-squares: Coeficiente de determinacion 
#R= raiz cuadrada de Rsquare

plot(combustible$Temperatura, combustible$Combustible)
abline(regresion)

cor(combustible) #correlacion
pairs(combustible)
anova(regresion)

#Verificación de supuestos
residuos<-regresion$residuals
length(residuos)

#Normalidad de los residuos
plot(residuos)
library(nortest) #lillie & Anderson-Da
library(moments) #D'Agostino
library(ggplots2)

shapiro.test(residuos) #Shapiro-Wilk normality test
ad.test(residuos) #Anderson-Darling normality test
lillie.test(residuos) #Lilliefors (Kolmogorov-Smirnov) normality test
agostino.test(residuos) #D'Agostino

Desv=sqrt(var(residuos))
Desv
LI2D=-2*Desv
LS2D=2*Desv
LI2D
LS2D
plot(residuos)
abline(h=LS2D, col="red")
abline(h=LI2D, col="red")

#Residuos con media cero
TestMedias<-t.test(residuos, alternative = 'two.sided', conf.level = 0.95, mu=0) #One Sample t-test
TestMedias
Media=mean(residuos)
Media
Desv=sqrt(var(residuos))
Desv
N=length(combustible)
cociente=(sqrt(N))*(Media/Desv)
cociente

#Independiencia de errores
library(forecast)
library(lmtest)

checkresiduals(regresion)  #Breusch-Godfrey test
dwtest(regresion) #Durbin-Watson test

pairs(combustible)
anova(regresion)
Box.test(residuos, lag=3, type="Ljung") #Box-Ljung test

#help(plot) //ejemplo de ayuda
#length(x) tamaño 


