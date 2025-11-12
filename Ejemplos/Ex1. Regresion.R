#Ex1. Seminario de Investigación 25 Enero 2022

# Lectura de la base de datos desde Excel
combustible<-read.csv("C:\\Users\\Areliuuu\\Dropbox\\PC\\Desktop\\Seminario de Investigacion I\\Ejemplo1_Consumo_de_Combustible.csv")

#Gráfico de dispersión de Temperatura vs Combustible
plot(combustible$Temperatura, combustible$Combustible)
plot(combustible)

#Regresión lm=ajusta modelos lineales
regresion<-lm(formula=Combustible~Temperatura, data=combustible)#y(Comb)=x(Temp); dataframe=archivo
regresion

summary(regresion) #Valor min, quartil, valorp, t, etc/ ver asteriscos
#R-squares: Coeficiente de determinación 
#R= raiz cuadrada de Rsquare

plot(combustible$Temperatura, combustible$Combustible)
abline(regresion)

cor(combustible)
pairs(combustible)
anova(regresion)

#help(plot) //ejemplo de ayuda