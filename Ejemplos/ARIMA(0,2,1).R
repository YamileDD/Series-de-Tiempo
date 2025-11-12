library(fpp2)
library(ggplot2)
library(forecast)
library(fma)
library(expsmooth)
library(tseries)
library(nortest)
library(lmtest)

#Lectura de serie IPC

Indice<-read.csv("C:\\Users\\Areliuuu\\Dropbox\\PC\\Desktop\\Seminario de Investigacion I\\IPC.csv")
write.csv(Indice, file = "serie.csv")

IPC<-ts(Indice, start = c(1969, 01), frequency = 12) #Datos que se quieren de archivo

AUTOCORR<-acf(IPC, lag.max = 100) #covarianza y autocorrelación con graficos
PAUTOCORR<-pacf(IPC)# autocorrelación parcial
VarIPC<-var(IPC)#varianza
plot(IPC)

adf.test(IPC) #Prueba de Dickey-Fuller Aumentado


#Primera diferencia de IPC
IPCdiff1<-diff(IPC, differences = 1)
plot(IPCdiff1)
adf.test(IPCdiff1)
VarIPCdiff1<-var(IPCdiff1)
Autocorrelacionessimples1<-acf(IPCdiff1)
ggtsdisplay(IPCdiff1, lag.max = 100) #ACF=autocorrelación muestral; PACF=Autocorrelación parcial
#Todo lo que se sale de las bandas se considera significativamente distinto de cero

#Segunda diferencia de IPC
IPCdiff2<-diff(IPC, differences = 2)
plot(IPCdiff2)
adf.test(IPCdiff2)
VarIPCdiff2<-var(IPCdiff2)
Autocorrelacionessimples2<-acf(IPCdiff2)
ggtsdisplay(IPCdiff2, lag.max = 100)

ggseasonplot(IPC) #Grafico para ver estacionalidad
nsdiffs(IPC) #Diferenciar en la parte estacional
ndiffs(IPC) #Estima el numero de diferencias para hacer la serie estacionaria

#Estabilización de la varianza y transformación de Box-Cox

lambda0<-BoxCox.lambda(IPC)
lambda0
BoxCoxIPC<-BoxCox(IPC, lambda0)
autoplot(BoxCoxIPC)
adf.test(BoxCoxIPC)

auto.arima(IPC) #parametros(autoregresivo de orden,integrado de orden, media movil de orden)
#No verifica supuestos

#Grafica de funciones de autocorrelacion y test para verificar estacionariedad 
ggAcf(BoxCoxIPC, lag.max = 100)
ggtsdisplay(BoxCoxIPC, lag.max = 100)

adf.test(BoxCoxIPC) ##Prueba de Dickey-Fuller Aumentado

#Primera diferenciacion

BoxCoxIPCdiff1<-diff(BoxCoxIPC, differences = 1)
ggtsdisplay(BoxCoxIPCdiff1, lag.max = 100)
adf.test(BoxCoxIPCdiff1)
auto.arima(BoxCoxIPCdiff1)
par(mfrow=c(2,1))
FACMBoxCoxIPCdiff1<-acf(BoxCoxIPCdiff1)
FACPBoxCoxIPCdiff1<-pacf(BoxCoxIPCdiff1)
ggAcf(BoxCoxIPCdiff1, lag.max = 100)
ggtsdisplay(BoxCoxIPCdiff1, lag.max = 100)

#Segunda diferenciacion

BoxCoxIPCdiff2<-diff(BoxCoxIPC, differences = 2)
ggtsdisplay(BoxCoxIPCdiff2, lag.max = 100)
adf.test(BoxCoxIPCdiff2)
auto.arima(BoxCoxIPCdiff2)
par(mfrow=c(2,1))
FACMBoxCoxIPCdiff2<-acf(BoxCoxIPCdiff2)
FACPBoxCoxIPCdiff2<-pacf(BoxCoxIPCdiff2)
ggAcf(BoxCoxIPCdiff2, lag.max = 100)
ggtsdisplay(BoxCoxIPCdiff2, lag.max = 100)

#Modelos propuestos

#ARIMA (0,2,1), ARIMA (0,2,2)

#Estimacion de parametros e intervalos de confianza


#Estimacion de ARIMA (0,1,1)
#Modelo 1
ARIMA011<-Arima(BoxCoxIPCdiff1, order = c(0,0,1))
#Modelo 2
ARIMA012<-Arima(BoxCoxIPCdiff1, order = c(0,0,2))
#Modelo 3
AutoARIMA012<-auto.arima(BoxCoxIPCdiff1)

#Estimacion de ARIMA (0,2,1)
#Modelo 1
ARIMA021<-Arima(BoxCoxIPCdiff2, order = c(0,0,1))


#Analisis de residuos para ARIMA(0,2,1)
#Resuduales
resARIMA021<-residuals(ARIMA021)
plot(resARIMA021)
checkresiduals(resARIMA021)

#Supuestos 1: Media cero, modelo: ARIMA(0,2,1) 
Media=mean(resARIMA021)
Media 
Desv=sqrt(var(resARIMA021))
Desv
N=length(IPC)
p=0
d=2
q=1
cociente<-(sqrt(N-d-p))*(Media/Desv)
cociente #deber ser menor que 2

#Supuesto 2: Residuos con varianza constante
checkresiduals(resARIMA021) #checar p-valor<0.05 autocorrelacion y no ruido blanco

#Supuesto 3: Residuos independientes: Prueba de Ljung-Box o Box Pierce
#H0: Los datos se distribuyen de forma independiente
checkresiduals(resARIMA021)

#Supuesto 4: Normalidad y observaciones aberrantes
qqnorm(resARIMA021)
qqline(resARIMA021)
checkresiduals(resARIMA021)
shapiro.test(resARIMA021) #p-valor<0.05 no hay normalidad en datos
lillie.test(x = resARIMA021)

#Intervalo de 2 desviaciones
LI2D=-2*Desv
LS2D=2*Desv
LI2D
LS2D

#Intervalo 3 desviaciones
LI3D=-3*Desv
LS3D=3*Desv
LI3D
LS3D

plot(resARIMA021)  #95% deben estar dentro del rango
abline(h=LS2D, col="red")
abline(h=LI2D, col="red")

#Supuesto 5: EL modelo es parsimonioso

confint(resARIMA021)
