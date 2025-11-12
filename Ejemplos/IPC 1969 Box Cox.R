library(fpp2)
library(ggplot2)
library(forecast)
library(fma)
library(expsmooth)
library(tseries)
library(nortest)
library(lmtest)


#Lectura de serie IPC

Indice<-read.csv("D:\\Users\\21228\\Desktop\\Carpeta UDLAP computadora Lenovo\\UDLAP\\PRIMAVERA 2018\\Series de Tiempo\\Ejemplos_Guerrero2009\\IPC.csv")

write.csv(Indice, file = "serie.csv") 

IPC<-ts(Indice, start = c(1969, 01), frequency=12)

AUTOCORR<-acf(IPC)
PAUTOCORR<-pacf(IPC)
VarIPC<-var(IPC)
VarIPC
plot(IPC)

adf.test(IPC)
#dwtest(IPC)

#Primera diferencia de IPC
IPCdiff1<-diff(IPC, differences=1)
plot(IPCdiff1)
adf.test(IPCdiff1)
VarIPCdiff1<-var(IPCdiff1)


#Segunda diferencia de IPC
IPCdiff2<-diff(IPC, differences=2)
plot(IPCdiff2)
adf.test(IPCdiff2)
VarIPCdiff2<-var(IPCdiff2)
Autocorrelacionesimples<-acf(IPCdiff2)
ggtsdisplay(IPCdiff2, lag.max=100)


ggseasonplot(IPC)
nsdiffs(IPC)

#Estabilizacion de varianza
#Transformacion Box-Cox

lambda0<-BoxCox.lambda(IPC)
lambda0
BoxCoxIPC<-BoxCox(IPC,lambda0)
autoplot(BoxCoxIPC)
adf.test(BoxCoxIPC)
auto.arima(IPC)

#Transformacion logarotmica
BoxCoxIPC<-1/IPC
autoplot(BoxCoxIPC) 
ggAcf(BoxCoxIPC)

autoplot(log(IPC))

#Grafica de funciones de autocorrelacion y test para verificar estacionariedad

ggAcf(BoxCoxIPC, lag.max=100)
ggtsdisplay(BoxCoxIPC, lag.max=100)

adf.test(BoxCoxIPC) #Prueba de Dickey Fuller

#Primera Diferenciacion

BoxCoxIPCdiff1<-diff(BoxCoxIPC, differences=1)
ggtsdisplay(BoxCoxIPCdiff1, lag.max=100)
adf.test(BoxCoxIPCdiff1) #Prueba de Dickey Fuller
auto.arima(BoxCoxIPCdiff1)

#Segunda Diferenciacion
BoxCoxIPCdiff2<-diff(BoxCoxIPC, differences=2)
adf.test(BoxCoxIPCdiff2)
par(mfrow=c(2,1))
FACMBoxCoxIPCdiff2<-acf(BoxCoxIPCdiff2)
FACPBoxCoxIPCdiff2<-pacf(BoxCoxIPCdiff2)

ggAcf(BoxCoxIPCdiff2)
ggtsdisplay(BoxCoxIPCdiff2)


plot(BoxCoxIPCdiff2)
ggtsdisplay(BoxCoxIPCdiff2, lag.max=80)
adf.test(BoxCoxIPCdiff2) #Prueba de Dickey Fuller
Var_BCIPCdiff2<-var(BoxCoxIPCdiff2)
sqrt(Var_BCIPCdiff2)
ggAcf(BoxCoxIPCdiff2, lag.max=100)
ggPacf(BoxCoxIPCdiff2, lag.max=100)

#Tercera Diferenciacion

BoxCoxIPCdiff3<-diff(BoxCoxIPC, differences=3)
ggtsdisplay(BoxCoxIPCdiff3, lag.max=100)
adf.test(BoxCoxIPCdiff3) #Prueba de Dickey Fuller
Var_BCIPCdiff3<-var(BoxCoxIPCdiff3)
sqrt(Var_BCIPCdiff3)

#Graficas de la funcion de autocorrelacion y autocorrelacion parcial de la
#segunda diferencia

ggtsdisplay(BoxCoxIPCdiff2, lag.max=200)
auto.arima(BoxCoxIPCdiff2)

#Modelos propuestos

#ARIMA(0,2,1), ARIMA (0,2,2)

#Estimacion de parametros e intervalos de confianza

#Estimacion ARIMA(0,2,1)
#Modelo 1
(ARIMA021 <- Arima(BoxCoxIPCdiff2, order=c(0,0,1)))
#Modelo 2
(ARIMA022 <- Arima(BoxCoxIPCdiff2, order=c(0,0,2)))
#Modelo 3
autoARIMA022<-auto.arima(BoxCoxIPCdiff2)

(ARIMA022 <- Arima(BoxCoxIPCdiff2, order=c(0,0,2)))

IC_ARIMA022<-confint(ARIMA022)
IC_ARIMA022


#Analisis de residuos para ARIMA(0,2,1)
#Residuales
resARIMA021 <- residuals(ARIMA021)
resaRIMA021<-residuals(aRIMA021)
resautoARIMA022<-residuals(autoARIMA022)
checkresiduals(autoARIMA022)

#Supuesto 1: Media cero, modelo: ARIMA(0,2,1)

Media=mean(resautoARIMA022)
Media
Desv=sqrt(var(resautoARIMA022))
Desv
N=length(IPC)
p=0
d=2
q=2
cociente=(sqrt(N-d-p))*(Media/Desv)
cociente

#Supuesto 2: Residuos con varianza constante:
checkresiduals(ARIMA021)

#Supuesto 3: Residuos independientes: Prueba de Ljung-Box o Box Pierce
#H0:Los datos se distribuyen de forma independiente.
checkresiduals(ARIMA022)

#Supuesto 4: Normalidad y observaciones aberrantes
qqnorm(resARIMA021)
qqline(resARIMA021)
checkresiduals(ARIMA010)
shapiro.test(resARIMA022)
lillie.test(x = resARIMA021)

#Intervalo de dos desviaciones
LI2D=-2*Desv
LS2D=2*Desv
LI2D
LS2D

#Intervalo 3 desviaciones
LI3D=-3*Desv
LS3D=3*Desv
LI3D
LS3D

plot(resARI22, ylim=c(-2*10^(-3),6*10^(-4)))
plot(resautoARIMA022)
abline(h=LS2D, col="red")
abline(h=LI2D, col="red")
abline(h=LI3D, col="blue")
abline(h=LS3D, col="blue")

#Supuesto 5: EL modelo es parsimonioso

Theta1=-0.5443
Theta2=-0.1929

Theta=-0.6817
sigma2estim<-(1-Theta^2)/90

LITheta1=Theta1-(2*0.07712)
LITheta1
LSTheta1=Theta1+(2*0.07712)
LSTheta1

LITheta2=Theta2-(2*0.1108)
LITheta2
LSTheta2=Theta2+(2*0.1108)
LSTheta2

confint(ARIMA022)

#Supuesto 6: El modelo es admisible

#Polinomio caracteristico: 1+0.5443x+0.1929x^2=0
Theta1+Theta2
Theta2-Theta1
#Discriminante:(Theta^1)+4Theta2
discriminante<-(Theta1^2)+(4*Theta2)
discriminante
#discriminante menor que cero: raices complejas conjugadas.

par(mfrow=c(2,1))
plot(BoxCoxIPCdiff2)
acf(BoxCoxIPCdiff2, lag.max=50)
pacf(BoxCoxIPCdiff2)
#Modelo (1+0.5439B+0.1928B^2)Zt=Wt

#Supuesto 7: El modelo es estable: bajas correlaciones entre pares.
#Corr(parametro 1, parametro 2)=Cov(parametro1, parametro2)/sqrt(Var(parametro1)Var(parametro2))
VarTheta1=VarTheta2=(1-Theta2^2)/(N-d)
VarTheta1
Cov=(-Theta1*(1+Theta2))/(N-d)
Cov
#Correlacion
Correlacion=(-phi1*(1+phi2))/(1-(phi2^2))
Correlacion=Cov/sqrt(VarTheta1*VarTheta2)


Correlacion

#Pronostico

PronosticoARIMA022<-forecast(ARIMA022,h=10)
autoplot(PronosticoARIMA022)

#auto.arima

autoarima<-auto.arima(IPC)
Pronosticoautoarima<-forecast(autoarima, h=10)
autoplot(Pronosticoautoarima)




#Transformaciones para lograr normalidad (para datos con distribuci�n de cola derecha)
#log(x)
#log(x+1)
#sqrt(x)
#1/x

#Transformaciones para lograr normalidad (para datos con distribuci�n de cola izquierda)
#x^2




