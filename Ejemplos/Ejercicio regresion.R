grasas <- read.table('http://verso.mat.uam.es/~joser.berrendero/datos/EdadPesoGrasas.txt', header = TRUE)
names(grasas)
pairs(grasas) #Matriz de diagrama de dispersion
cor(grasas) #Matriz de coeficientes de correlacion
regresion <- lm(grasas ~ edad, data = grasas) #lm:modelo lineal, variable dependiente ~ vindep, data=nombre fichero
summary(regresion)
#Graficar
plot(grasas$edad, grasas$grasas, xlab='Edad', ylab='Grasas')
abline(regresion)
#Calcular predicciones
nuevas.edades1 <- data.frame(edad = seq(30, 50)) #Solo edades de 30 a 50 años
predict(regresion, nuevas.edades1)
#Intervalos de confianza
confint(regresion) 
confint(regresion, level = 0.90)

#Ejemplo e intervalos de confianza
nuevas.edades <- data.frame(edad = seq(20, 60))
plot(grasas$edad, grasas$grasas, xlab='Edad', ylab='Grasas')
abline(regresion)

# ic es una matriz con tres columnas: la primera es la prediccion, las otras dos son los extremos del intervalo
ic <- predict(regresion, nuevas.edades, interval = 'confidence')
lines(nuevas.edades$edad, ic[, 2], lty = 2)
lines(nuevas.edades$edad, ic[, 3], lty = 2)

# Intervalos de prediccion
ic <- predict(regresion, nuevas.edades, interval = 'prediction')
lines(nuevas.edades$edad, ic[, 2], lty = 2, col = 'red')
lines(nuevas.edades$edad, ic[, 3], lty = 2, col = 'red')

#Analisis de varianza ANOVA
anova(regresion)
residuos <- rstandard(regresion) 
valores.ajustados <- fitted(regresion) #valores ajustados
plot(valores.ajustados, residuos)
residuals(regresion)
qqnorm(residuos) #Para hipotesis de normalidad
qqline(residuos)
