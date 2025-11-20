install.packages("fpp2")
install.packages("tseries")
install.packages("lmtest")
install.packages("car")
install.packages("readr")
library(fpp2)
library(ggplot2)
library(forecast)
library(fma)
library(expsmooth)
library(tseries)
library(nortest)
library(lmtest)
library(car)
library(readr)

show_col_types=FALSE
dataBase<-read_csv("serie BIE.csv")
head(dataBase)

### Convertimos en serie de tiempo
# Primero del 2025
tsP1<-ts(dataBase, start=c(2005,01),end=c(2024,04), frequency = 4)
tsP1

var(tsP1) #84.87063

#BoxCox
lambdaP<-BoxCox.lambda(tsP1) #-0.9999242

BoxCoxTs1<-BoxCox(tsP1,lambdaP)
var(BoxCoxTs1) #1.037444e-06 disminuye considerablemente la varianza


#Gráficas
autoplot(tsP1)
autoplot(BoxCoxTs1)
#Por lo que vemos en las gráficas no es normal.

#Dickey-Fuller
adf.test(BoxCoxTs1)
#0.01344 no es estacionario
# Casi es estacionaro pero levemente el p-value es mayor a 0.01 por lo que agregamos una diferencia 
shapiro.test(BoxCoxTs1) #0.0783


#Debemos aplicar diferencias
#1RA Diferencia
BoxCoxdiff1<-diff(BoxCoxTs1,differences=1)
adf.test(BoxCoxdiff1)
# p-value ya es menor a 0.01.
#Con una diferencia ya es estacionario
shapiro.test(BoxCoxdiff1) #1.778e-07

var(BoxCoxdiff1) #1.787235e-07 menor varianza con una diferencia
par(mfrow=c(1,1))

modelo1<-auto.arima(tsP1) #ARIMA(0,1,1),theta1=-0.5764 & el drift de 0.4028 
modelo1 #ARIMA(0,1,1) con drift 0.4028

#FAC y FACP 1ra dif.
autoplot(acf(BoxCoxdiff1,plot=FALSE))+labs(title="Diff1 Autocorrelaciones Simples")
autoplot(pacf(BoxCoxdiff1,plot=FALSE))+labs(title="Diff1 Autocorrelaciones Parciales")
#De acuerdo a la gráfica de autocorrelaciones consideramos el modelo ARIMA(0,1,1) 

#2DA Diferencia
BoxCoxdiff2<-diff(BoxCoxTs1,differences=2)

autoplot(acf(BoxCoxdiff2,plot=FALSE))+labs(title="Diff2 Autocorrelaciones Simples")
autoplot(pacf(BoxCoxdiff2,plot=FALSE))+labs(title="Diff2 Autocorrelaciones Parciales")
#De acuerdo a las gráficas pareciera que tenemos un ARIMA(0,2,1)

#3RA Diferencia
BoxCoxdiff3<-diff(BoxCoxTs1,differences=3)

autoplot(acf(BoxCoxdiff3,plot=FALSE))+labs(title="Diff3 Autocorrelaciones Simples")
autoplot(pacf(BoxCoxdiff3,plot=FALSE))+labs(title="Diff3 Autocorrelaciones Parciales")
#De acuerdo a la gráficas pareciera que tenemos un ARIMA(0,3,1)

c(var_diff1 = var(BoxCoxdiff1),
  var_diff2 = var(BoxCoxdiff2),
  var_diff3 = var(BoxCoxdiff3))
#tienes que sacar la varianza de cada uno y te quedes con el de la varianza más chica

#agregar el drift 0.4028 y como es constante sale ###
#aquí le cambié el drift que ya esta aplicado al usar el modelo ya transformado
#Probar supuestos
modelo1<-Arima(BoxCoxdiff1,c(0,0,1), include.mean = TRUE)

#Supuesto
resModelo1<-modelo1$residuals

nResM1<-length(resModelo1)

nTotal<-length(tsP1)

#Parámetros

mediaResModelo1<-mean(resModelo1) #-2.930871e-06
desvResModelo1<-sqrt(sum(resModelo1^2)/(nResM1-1)) #0.0003963612

modelo1
coef(modelo1)
modelo1[["var.coef"]]

theta1M1 <- -1 * as.numeric(coef(modelo1)["ma1"])
muM1    <- as.numeric(coef(modelo1)["intercept"])

resModelo1 <- residuals(modelo1)
nResM1 <- length(resModelo1)
nTotal <- length(tsP1)

#### Supuesto 1:  media ~ 0
sumaMediaResM1 <- sum(resModelo1)       
mediaResModelo1 <- sumaMediaResM1 / nResM1
desvResModelo1 <- sqrt(sum(resModelo1^2) / (nResM1 - 1))
abs( sqrt(nResM1) * mediaResModelo1 / desvResModelo1 )

t.test(resModelo1)$p.value

#### Supuesto 2: Residuos varianza constante
autoplot(resModelo1) + ggtitle("Residuos del modelo")  #en general existe varianza monótona

##### Supuesto 3: Independencia, residuos no correlacionados
#Pruebas Box-Pierce, Ljung-Box o DW
#En este caso utilizamos Box-Pierce o Ljung-Box ya que no tenemos AR(1)

#H0: pk = 0 para todo k
#HA: existe k, pk=!0

Box.test(resModelo1,type=c("Box-Pierce"), lag = 4) #0.47
Box.test(resModelo1,type=c("Ljung-Box"), lag = 4) #0.4373

#No rechazamos H0, por lo que los residuos son independientes.

##### Supuesto 4: Normalidad

lillie.test(resModelo1) #0.0003541
shapiro.test(resModelo1) #4.835e-08
ad.test(resModelo1) #8.723e-06
#Los residuos no son normales

#Regla empírica ±2σ
supModelo1 <- mediaResModelo1 + 2 * desvResModelo1
infModelo1 <- mediaResModelo1 - 2 * desvResModelo1
mean(resModelo1 > infModelo1 & resModelo1 < supModelo1) #0.9620253
#Se cumple que los residuos se encuentran dentro de +-2 desv std

##### Supuesto 5: Modelo Parsimonioso
#Construir intervalos alrededor de los parámetros
#Sabemos que nuestros parámetros(coeficientes) son theta1=0.4137609
vc <- modelo1[["var.coef"]]
varTheta1 <- as.numeric(vc["ma1","ma1"])
IC_theta1 <- c(theta1M1 - 2*sqrt(varTheta1), theta1M1 + 2*sqrt(varTheta1))
IC_theta1
#(0.1586426 , 0.6688793)
#No contiene al cero

#### Supuesto 6: Modelo admisible
#Verificar estacionariedad e invertibilidad
#Para MA(1) invertible: |ma1(param de R)| < 1
abs(as.numeric(coef(modelo1)["ma1"])) < 1
#Al ser un MA(1) es estacionario y comprobamos que es invertible.

#### Supuesto 7: Estabilidad (correlaciones entre parámetros)
# Solo hay 1 parámetro MA1 en este ajuste sin drift;
#si incluimos al drift podemos revisar correlaciones y cambia


#### Supuesto 8: Valores atípicos ±3σ
sup_m1 <- mediaResModelo1 + 3 * desvResModelo1 #0.001194728
inf_m1 <- mediaResModelo1 - 3 * desvResModelo1 #-0.001199485
sum(!(resModelo1 > inf_m1 & resModelo1 < sup_m1))
# # Solo 1 residuo queda fuera de ±3σ; no hay evidencia de outliers importantes.

################
#Modelo 2
#ARIMA(0,2,1)
#Probar supuestos
###Supuesto 1: Residuos tienen media 0

modelo2<-Arima(BoxCoxdiff1,c(0,1,1))
modelo2 #theta=-1

resModelo2<-modelo2$residuals

nResM2<-length(resModelo2) #79

#Párametros

mediaResModelo2<-mean(resModelo2) #-1.385077e-05
desvResModelo2<-sqrt(sum(resModelo2^2)/(nResM2-1)) #0.000422757

modelo2
coef(modelo2)
modelo2[["var.coef"]]

theta1M2 <- -1 * as.numeric(coef(modelo2)["ma1"])

resModelo2 <- residuals(modelo2)
nResM2 <- length(resModelo2)
nTotal <- length(tsP1)

#Prueba 
abs(sqrt(nTotal-1-1)*mediaResModelo2/desvResModelo2) #0.2893545
#Ésto es menor a 2, los residuos tiene media cero
t.test(resModelo2)$p.value

####Supuesto 2: residuos varianza constante
autoplot(resModelo2)

plot(resModelo2, col="blue",main="Residuos Modelo Arima(0,2,1)",xlab="Tiempo",ylab="Residuos")
#en general existe varianza monótona

###Supuesto 3: Independencia, residuos no correlacionados
#Pruebas Box-Pierce, Ljung-Box o DW

#H0: pk = 0 para todo k
#HA: existe k, pk=!0
###lag = 4
Box.test(resModelo2,type=c("Box-Pierce"), lag = 4) #pvalue = 0.2732
Box.test(resModelo2,type=c("Ljung-Box"), lag = 4) #pvalue= 0.2518

#no se rechaza H0, por lo que los residuos son independientes

### Supuesto 4:Normalidad
lillie.test(resModelo2) #4.414e-05
shapiro.test(resModelo2) #3.542e-08
ad.test(resModelo2) #1.412e-08
#Dado que los p-values son muy pequeños, se rechaza H0 de normalidad,
#por lo que los residuos NO siguen una distribución normal.

#Regla empírica ±2σ
supModelo2 <- mediaResModelo2 + 2 * desvResModelo2
infModelo2 <- mediaResModelo2 - 2 * desvResModelo2
mean(resModelo2 > infModelo2 & resModelo2 < supModelo2) #0.9493671
#Se cumple que los residuos se encuentran dentro de +-2 desv std

##### Supuesto 5: Modelo Parsimonioso
#Construir intervalos alrededor de los parámetros
#Sabemos que nuestros parámetros(coeficientes) son theta1=0.4137609
vc <- modelo2[["var.coef"]]
varTheta2 <- as.numeric(vc["ma1","ma1"])
IC_theta2 <- c(theta1M2 - 2*sqrt(varTheta2), theta1M2 + 2*sqrt(varTheta2))
IC_theta2
#(0.9297714 , 1.0702286)
#No contiene al cero

#### Supuesto 6: Modelo admisible
#Verificar estacionariedad e invertibilidad
#Para MA(1) invertible: |ma1(param de R)| < 1
abs(as.numeric(coef(modelo2)["ma1"])) < 1
#Al ser un MA(1) es estacionario y comprobamos que es invertible.

#### Supuesto 7: Estabilidad (correlaciones entre parámetros)
# Solo hay 1 parámetro MA1 en este ajuste sin drift;
#si incluimos al drift podemos revisar correlaciones y cambia

#### Supuesto 8: Valores atípicos ±3σ
sup_m2 <- mediaResModelo2 + 3 * desvResModelo2 #0.00125442
inf_m2 <- mediaResModelo2 - 3 * desvResModelo2 #-0.001282122
sum(!(resModelo2 > inf_m2 & resModelo2 < sup_m2))
#Solo 2 residuos quedan fuera de ±3σ; pocos outliers y 
#sin afectar el modelo de forma relevante.


#################
#Modelo 3

modelo3<-Arima(BoxCoxTs1,c(0,3,1))
modelo3 #theta=-0.9999

resModelo3<-modelo3$residuals
nResM3<-length(resModelo3)

#Parámetros
mediaResModelo3<-mean(resModelo3) #-5.55924e-05
desvResModelo3<-sqrt(sum(resModelo3^2)/(nResM3-1)) #0.0006903497

modelo3
coef(modelo3)
modelo3[["var.coef"]]

#### Supuesto 1: media cero

theta1M3 <- -1 * as.numeric(coef(modelo3)["ma1"])

resModelo3 <- residuals(modelo3)
nResM3 <- length(resModelo3)
nTotal <- length(tsP1)

t.test(resModelo3)$p.value

#Prueba 
abs(sqrt(nTotal-1-1)*mediaResModelo3/desvResModelo3) #0.711203
#esto es menor a 2, los residuos tiene media cero

t.test(resModelo3)$p.value

#### Supuesto 2: residuos varianza constante
autoplot(resModelo3)
plot(resModelo3, col="blue",main="Residuos Modelo Arima(0,3,1)",xlab="Tiempo",ylab="Residuos")

#en general varianza constante/monotona 


###Supuesto 3: Independencia, residuos no correlacionados

#H0: pk = 0 para todo k
#HA: existe k, pk=!0

Box.test(resModelo3,type=c("Box-Pierce"), lag = 4) #pvalue = 5.546e-06
Box.test(resModelo3,type=c("Ljung-Box"), lag = 4) #pvalue = 3.194e-06

# se rechaza H0, por lo que los residuos no son independientes

### Supuesto 4:Normalidad
lillie.test(resModelo3) #p-value = 9.205e-06
shapiro.test(resModelo3) #p-value = 1.117e-07
ad.test(resModelo3) #p-value = 5.829e-08

# se rechaza H0: ya que los p-values son muy pequeños
# por lo que no hay normalidad en los residuos

#Regla empírica ±2σ
supModelo3 <- mediaResModelo3 + 2 * desvResModelo3 #0.001325107
infModelo3 <- mediaResModelo3 - 2 * desvResModelo3 #-0.001436292
mean(resModelo3 > infModelo3 & resModelo3 < supModelo3) #0.9375
#Se cumple que los residuos se encuentran dentro de +-2 desv std
# Aproximadamente 94% de los residuos cae dentro de ±2σ
# Se cumple razonablemente con la regla empírica.

##### Supuesto 5: Modelo Parsimonioso
#Construir intervalos alrededor de los parámetros
#Sabemos que nuestros Parámetros(coeficientes) son theta1=0.4137609
vc <- modelo3[["var.coef"]]
varTheta3 <- as.numeric(vc["ma1","ma1"])
IC_theta3 <- c(theta1M3 - 2*sqrt(varTheta3), theta1M3 + 2*sqrt(varTheta3))
IC_theta3
#(0.9340623 , 1.0657963)
#No contiene al cero

#### Supuesto 6: Modelo admisible
#Verificar estacionariedad e invertibilidad
#Para MA(1) invertible: |ma1(param de R)| < 1
abs(as.numeric(coef(modelo3)["ma1"])) < 1
#Al ser un MA(1) es estacionario y comprobamos que es invertible.

#### Supuesto 7: Estabilidad (correlaciones entre parámetros)
# Solo hay 1 parámetro MA1 en este ajuste sin drift;
#si incluimos al drift podemos revisar correlaciones y cambia

#### Supuesto 8: Valores atípicos ±3σ
sup_m3 <- mediaResModelo3 + 3 * desvResModelo3 #0.002015457
inf_m3 <- mediaResModelo3 - 3 * desvResModelo3 #-0.002126641
sum(!(resModelo3 > inf_m3 & resModelo3 < sup_m3))
#Solo 2 residuos quedan fuera de ±3σ; se detectan pocos posibles valores atípicos.

##### Pronósticos Modelo 1
#Residuales
# Y_t = BoxCoxTs1[t]
# ΔY_t = Y_t - Y_{t-1} = muM1 + a_t + theta1M1 * a_{t-1}
theta1M1
aHatM1[1]=0
aHatM1[2]=BoxCoxTs1[2]-BoxCoxTs1[1]- muM1-theta1M1*aHatM1[1]
aHatM1[3]=BoxCoxTs1[3]-BoxCoxTs1[2]- muM1-theta1M1*aHatM1[2]
aHatM1[4]=BoxCoxTs1[4]-BoxCoxTs1[3]- muM1-theta1M1*aHatM1[3]
aHatM1[5]=BoxCoxTs1[5]-BoxCoxTs1[4]- muM1-theta1M1*aHatM1[4]
aHatM1[6]=BoxCoxTs1[6]-BoxCoxTs1[5]- muM1-theta1M1*aHatM1[5]
aHatM1[7]=BoxCoxTs1[7]-BoxCoxTs1[6]- muM1-theta1M1*aHatM1[6]
aHatM1[8]=BoxCoxTs1[8]-BoxCoxTs1[7]- muM1-theta1M1*aHatM1[7]
aHatM1[9]=BoxCoxTs1[9]-BoxCoxTs1[8]- muM1-theta1M1*aHatM1[8]
aHatM1[10]=BoxCoxTs1[10]-BoxCoxTs1[9]- muM1-theta1M1*aHatM1[9]
aHatM1[11]=BoxCoxTs1[11]-BoxCoxTs1[10]- muM1-theta1M1*aHatM1[10]
aHatM1[12]=BoxCoxTs1[12]-BoxCoxTs1[11]- muM1-theta1M1*aHatM1[11]
aHatM1[13]=BoxCoxTs1[13]-BoxCoxTs1[12]- muM1-theta1M1*aHatM1[12]
aHatM1[14]=BoxCoxTs1[14]-BoxCoxTs1[13]- muM1-theta1M1*aHatM1[13]
aHatM1[15]=BoxCoxTs1[15]-BoxCoxTs1[14]- muM1-theta1M1*aHatM1[14]
aHatM1[16]=BoxCoxTs1[16]-BoxCoxTs1[15]- muM1-theta1M1*aHatM1[15]
aHatM1[17]=BoxCoxTs1[17]-BoxCoxTs1[16]- muM1-theta1M1*aHatM1[16]
aHatM1[18]=BoxCoxTs1[18]-BoxCoxTs1[17]- muM1-theta1M1*aHatM1[17]
aHatM1[19]=BoxCoxTs1[19]-BoxCoxTs1[18]- muM1-theta1M1*aHatM1[18]
aHatM1[20]=BoxCoxTs1[20]-BoxCoxTs1[19]- muM1-theta1M1*aHatM1[19]
aHatM1[21]=BoxCoxTs1[21]-BoxCoxTs1[20]- muM1-theta1M1*aHatM1[20]
aHatM1[22]=BoxCoxTs1[22]-BoxCoxTs1[21]- muM1-theta1M1*aHatM1[21]
aHatM1[23]=BoxCoxTs1[23]-BoxCoxTs1[22]- muM1-theta1M1*aHatM1[22]
aHatM1[24]=BoxCoxTs1[24]-BoxCoxTs1[23]- muM1-theta1M1*aHatM1[23]
aHatM1[25]=BoxCoxTs1[25]-BoxCoxTs1[24]- muM1-theta1M1*aHatM1[24]
aHatM1[26]=BoxCoxTs1[26]-BoxCoxTs1[25]- muM1-theta1M1*aHatM1[25]
aHatM1[27]=BoxCoxTs1[27]-BoxCoxTs1[26]- muM1-theta1M1*aHatM1[26]
aHatM1[28]=BoxCoxTs1[28]-BoxCoxTs1[27]- muM1-theta1M1*aHatM1[27]
aHatM1[29]=BoxCoxTs1[29]-BoxCoxTs1[28]- muM1-theta1M1*aHatM1[28]
aHatM1[30]=BoxCoxTs1[30]-BoxCoxTs1[29]- muM1-theta1M1*aHatM1[29]
aHatM1[31]=BoxCoxTs1[31]-BoxCoxTs1[30]- muM1-theta1M1*aHatM1[30]
aHatM1[32]=BoxCoxTs1[32]-BoxCoxTs1[31]- muM1-theta1M1*aHatM1[31]
aHatM1[33]=BoxCoxTs1[33]-BoxCoxTs1[32]- muM1-theta1M1*aHatM1[32]
aHatM1[34]=BoxCoxTs1[34]-BoxCoxTs1[33]- muM1-theta1M1*aHatM1[33]
aHatM1[35]=BoxCoxTs1[35]-BoxCoxTs1[34]- muM1-theta1M1*aHatM1[34]
aHatM1[36]=BoxCoxTs1[36]-BoxCoxTs1[35]- muM1-theta1M1*aHatM1[35]
aHatM1[37]=BoxCoxTs1[37]-BoxCoxTs1[36]- muM1-theta1M1*aHatM1[36]
aHatM1[38]=BoxCoxTs1[38]-BoxCoxTs1[37]- muM1-theta1M1*aHatM1[37]
aHatM1[39]=BoxCoxTs1[39]-BoxCoxTs1[38]- muM1-theta1M1*aHatM1[38]
aHatM1[40]=BoxCoxTs1[40]-BoxCoxTs1[39]- muM1-theta1M1*aHatM1[39]
aHatM1[41]=BoxCoxTs1[41]-BoxCoxTs1[40]- muM1-theta1M1*aHatM1[40]
aHatM1[42]=BoxCoxTs1[42]-BoxCoxTs1[41]- muM1-theta1M1*aHatM1[41]
aHatM1[43]=BoxCoxTs1[43]-BoxCoxTs1[42]- muM1-theta1M1*aHatM1[42]
aHatM1[44]=BoxCoxTs1[44]-BoxCoxTs1[43]- muM1-theta1M1*aHatM1[43]
aHatM1[45]=BoxCoxTs1[45]-BoxCoxTs1[44]- muM1-theta1M1*aHatM1[44]
aHatM1[46]=BoxCoxTs1[46]-BoxCoxTs1[45]- muM1-theta1M1*aHatM1[45]
aHatM1[47]=BoxCoxTs1[47]-BoxCoxTs1[46]- muM1-theta1M1*aHatM1[46]
aHatM1[48]=BoxCoxTs1[48]-BoxCoxTs1[47]- muM1-theta1M1*aHatM1[47]
aHatM1[49]=BoxCoxTs1[49]-BoxCoxTs1[48]- muM1-theta1M1*aHatM1[48]
aHatM1[50]=BoxCoxTs1[50]-BoxCoxTs1[49]- muM1-theta1M1*aHatM1[49]
aHatM1[51]=BoxCoxTs1[51]-BoxCoxTs1[50]- muM1-theta1M1*aHatM1[50]
aHatM1[52]=BoxCoxTs1[52]-BoxCoxTs1[51]- muM1-theta1M1*aHatM1[51]
aHatM1[53]=BoxCoxTs1[53]-BoxCoxTs1[52]- muM1-theta1M1*aHatM1[52]
aHatM1[54]=BoxCoxTs1[54]-BoxCoxTs1[53]- muM1-theta1M1*aHatM1[53]
aHatM1[55]=BoxCoxTs1[55]-BoxCoxTs1[54]- muM1-theta1M1*aHatM1[54]
aHatM1[56]=BoxCoxTs1[56]-BoxCoxTs1[55]- muM1-theta1M1*aHatM1[55]
aHatM1[57]=BoxCoxTs1[57]-BoxCoxTs1[56]- muM1-theta1M1*aHatM1[56]
aHatM1[58]=BoxCoxTs1[58]-BoxCoxTs1[57]- muM1-theta1M1*aHatM1[57]
aHatM1[59]=BoxCoxTs1[59]-BoxCoxTs1[58]- muM1-theta1M1*aHatM1[58]
aHatM1[60]=BoxCoxTs1[60]-BoxCoxTs1[59]- muM1-theta1M1*aHatM1[59]
aHatM1[61]=BoxCoxTs1[61]-BoxCoxTs1[60]- muM1-theta1M1*aHatM1[60]
aHatM1[62]=BoxCoxTs1[62]-BoxCoxTs1[61]- muM1-theta1M1*aHatM1[61]
aHatM1[63]=BoxCoxTs1[63]-BoxCoxTs1[62]- muM1-theta1M1*aHatM1[62]
aHatM1[64]=BoxCoxTs1[64]-BoxCoxTs1[63]- muM1-theta1M1*aHatM1[63]
aHatM1[65]=BoxCoxTs1[65]-BoxCoxTs1[64]- muM1-theta1M1*aHatM1[64]
aHatM1[66]=BoxCoxTs1[66]-BoxCoxTs1[65]- muM1-theta1M1*aHatM1[65]
aHatM1[67]=BoxCoxTs1[67]-BoxCoxTs1[66]- muM1-theta1M1*aHatM1[66]
aHatM1[68]=BoxCoxTs1[68]-BoxCoxTs1[67]- muM1-theta1M1*aHatM1[67]
aHatM1[69]=BoxCoxTs1[69]-BoxCoxTs1[68]- muM1-theta1M1*aHatM1[68]
aHatM1[70]=BoxCoxTs1[70]-BoxCoxTs1[69]- muM1-theta1M1*aHatM1[69]
aHatM1[71]=BoxCoxTs1[71]-BoxCoxTs1[70]- muM1-theta1M1*aHatM1[70]
aHatM1[72]=BoxCoxTs1[72]-BoxCoxTs1[71]- muM1-theta1M1*aHatM1[71]
aHatM1[73]=BoxCoxTs1[73]-BoxCoxTs1[72]- muM1-theta1M1*aHatM1[72]
aHatM1[74]=BoxCoxTs1[74]-BoxCoxTs1[73]- muM1-theta1M1*aHatM1[73]
aHatM1[75]=BoxCoxTs1[75]-BoxCoxTs1[74]- muM1-theta1M1*aHatM1[74]
aHatM1[76]=BoxCoxTs1[76]-BoxCoxTs1[75]- muM1-theta1M1*aHatM1[75]
aHatM1[77]=BoxCoxTs1[77]-BoxCoxTs1[76]- muM1-theta1M1*aHatM1[76]
aHatM1[78]=BoxCoxTs1[78]-BoxCoxTs1[77]- muM1-theta1M1*aHatM1[77]
aHatM1[79]=BoxCoxTs1[79]-BoxCoxTs1[78]- muM1-theta1M1*aHatM1[78]
aHatM1[80]=BoxCoxTs1[80]-BoxCoxTs1[79]- muM1-theta1M1*aHatM1[79]
#el número 81 es 118.3083856
#tranformar el 81
t81 <- 118.3083856
boxCox81 <- BoxCox(t81,lambdaP)
aHatM1[81]=boxCox81-BoxCoxTs1[80]- muM1-theta1M1*aHatM1[80]


nTotal <- length(BoxCoxTs1)
aHatM1 <- numeric(nTotal)
aHatM1[1] <- 0 

for (t in 2:nTotal) {
  aHatM1[t] <- BoxCoxTs1[t] - BoxCoxTs1[t-1] -muM1 - theta1M1 * aHatM1[t-1]
}
Tn <- nTotal

pronosticos10 <- numeric(10)

for (h in 1:10) {
  pronosticos10[h] <- BoxCoxTs1[Tn] + h * muM1 + theta1M1 * aHatM1[Tn]
}

pronostico1M1 <- pronosticos10[1]
pronostico2M1 <- pronosticos10[2]
pronosticos10


#Intervalos
phi1M1 <- 1
H <- 10

ep0 <- -1
ep1 <- theta1M1+(1+phi1M1)*ep0

epsis <- numeric(H)
epsis[1] <- -1
epsis[2] <- theta1M1 + (1 + phi1M1) * epsis[1]

for (j in 3:H) {
  epsis[j] <- (1 + phi1M1) * epsis[j-1] - phi1M1 * epsis[j-2]
}


MOES<-c()
intInfM1<-c()
intSupM1<-c()

for(j in 1:10){
  MOES[j]<-qnorm(0.05/2)*(sum(epsis[1:j]^2))^(1/2)*desvResModelo1
  intInfM1[j]<-pronosticos10[j]+MOES[j]
  intSupM1[j]<-pronosticos10[j]-MOES[j]
}
intInfM1
pronosticos10
intSupM1


#####Pronósticos Serie Original
varEthM1<-c()
for(j in 1:10){
  varEthM1[j]<-(sum(epsis[1:j]^2))*(desvResModelo1)^2
}
varEthM1

#sesgo
x1<-c()
sesgoM1<-c()
tInversa<-c()
pSerieOriginal<-c()

for(j in 1:10){
  x1[j]<-2*lambdaP*(lambdaP-1)*((1+lambdaP*pronosticos10[j])^(-2))*(varEthM1[j])
  sesgoM1[j]<-(1/2+sqrt(1-x1[j])/2)^(1/lambdaP)
  tInversa[j]<-(lambdaP*pronosticos10[j]+1)^(1/lambdaP)
  pSerieOriginal[j]<-tInversa[j]*sesgoM1[j]
}

pSerieOriginal

#Intervalo de predicción

intInfOriginal<-c()
intSupOriginal<-c()
MOEOriginal<-c()

for(j in 1:10){
  MOEOriginal[j]<-qnorm(0.05/2)*sqrt(varEthM1[j])
  intInfOriginal[j]<-pSerieOriginal[j]+MOEOriginal[j]
  intSupOriginal[j]<-pSerieOriginal[j]-MOEOriginal[j]
}

intInfOriginal
pSerieOriginal
intSupOriginal

##### Actualización de pronóstico

ZtNew<-118.3083856
TZtNew<-((ZtNew^lambdaP)-1)/lambdaP
aN<-TZtNew-pronosticos10[1]

actPronosticos<-c()

for(j in 1:9){
  actPronosticos[j]<-pronosticos10[j+1]-epsis[j+1]*aN  
}

actPronosticos

forecast(Arima(BoxCoxTs1,order=c(1,1,2)))
forecast(modelo1)

#Intervalos

intInfAct<-c()
intSupAct<-c()
MOESAct<-c()

for (j in 1:9) {
  MOESAct[j] <- qnorm(0.05/2) * (sum(epsis[1:j]^2))^(1/2) * desvResModelo1
  intInfAct[j] <- actPronosticos[j] + MOESAct[j]
  intSupAct[j] <- actPronosticos[j] - MOESAct[j]
}


intInfAct
intSupAct

#### Actualización serie original

pSerieOriginalAct<-c()
sesgoAct<-c()
xAct<-c()
tInvAct<-c()

for(j in 1:9){
  xAct[j]<-2*lambdaP*(lambdaP-1)*((1+lambdaP*actPronosticos[j])^(-2))*(varEthM1[j])
  sesgoAct[j]<-(1/2+sqrt(1-xAct[j])/2)^(1/lambdaP)
  tInvAct[j]<-(lambdaP*actPronosticos[j]+1)^(1/lambdaP)
  pSerieOriginalAct[j]<-tInvAct[j]*sesgoAct[j]
}

pSerieOriginalAct

#Intervalo de predicción actualización serie original

intInfOrigAct<-c()
intSupOrigAct<-c()
MOEOrigAct<-c()

for(j in 1:9){
  MOEOrigAct[j]<-qnorm(0.05/2)*sqrt(varEthM1[j])
  intInfOrigAct[j]<-pSerieOriginalAct[j]+MOEOrigAct[j]
  intSupOrigAct[j]<-pSerieOriginalAct[j]-MOEOrigAct[j]
}

intInfOrigAct
intSupOrigAct
