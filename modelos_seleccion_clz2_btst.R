library(tidyverse);library(viridis);library(lme4); library(fitdistrplus); library(MuMIn); library(cAIC4); library(parallel);library(pbkrtest)
library(ggcorrplot); library(PerformanceAnalytics); library(reshape2); library(hier.part); library(patchwork); library(ggeffects)
library(mgcv);library(boot); library(parallel)
df <- read.csv('df_sv.csv', header = T)

df$germ <- as.factor(df$germ)
summary(df)
sub_df <- df[,c('pob','indiv', 'fruto', 'sv', 'prob_germ', 'germ', 'masa', 
                'area_p', 'l_papus', 'vol_c', 'capitulos', 'altura_planta', 'pl_m', 'pl_v')]


### modelos sin promediar por individuo

sub_df$clz <- (1/sub_df$sv)*sub_df$prob_germ
sub_df$inverse_sv <- 1/sub_df$sv

# COLONIZACION #############################################################
## estandarizacion del fitness relativo. #####################################
sub_df$W = sub_df$clz/mean(sub_df$clz)



p3 <- ggplot(sub_df, aes(x = W)) + geom_histogram(aes(x = W, y = after_stat(density), color = germ),
                                                       fill = 'grey55', bins = 30)  + 
  geom_density(aes(color = germ)) + scale_y_continuous(expand = c(0,0)) + scale_color_viridis_d()
p3 + theme_bw()
## W ~ masa ##########################################
lm1_large = lmer(W ~ scale (masa) + I(scale(masa)^2) + (1|pob/indiv), data = sub_df)
summary(lm1_large)
lm1_small = lmer(W ~ scale (masa) + (1|pob/indiv), data = sub_df)
PBmodcomp(lm1_large, lm1_small, nsim = 10000)

fe <- function(fit) {return(fixef(fit))}
bootlm1 = bootMer(lm1_large,FUN=fe,nsim=500,use.u =F)
bootlm1$t[,2]
par(mfrow=c(1,2),cex.lab=1.5, cex.axis=1.3)
## distribucion para masa  y para el componente cuadratico
plot(density(bootlm1$t[,2]), main="", xlab="pendientes masa simuladas", ylab="Densidad", 
     lwd=3, col="red", xlim=c(0,0.7))                                                     
abline(v=quantile(bootlm1$t[,2],0.025)) ## intervalo del 95 por ciento
abline(v=quantile(bootlm1$t[,2],0.975))


plot(density(bootlm1$t[,3]), main="", xlab="pendiente cuadratica simulada", ylab="Densidad", 
     lwd=3, col="red", xlim=c(-0.25,0))                                                     
abline(v=quantile(bootlm1$t[,3],0.025)) ## intervalo del 95 por ciento
abline(v=quantile(bootlm1$t[,3],0.975))
layout(1)

## la masa es significativa ahora vamos a graficar 
pred_masa <- ggpredict(lm1_large, terms = "masa [all]")

ggplot(pred_masa, aes(x = x, y = predicted)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(title = "Efecto de la masa y su efecto cuadrático",
       x = "Masa",
       y = "Fitness Predicho") +
  theme_minimal()
## W ~ area_p ###################################
lm2_large_p = lmer(W ~ scale (area_p) + I(scale(area_p)^2) + (1|pob/indiv), data = sub_df)
summary(lm2_large_p)
lm2_small_p = lmer(W ~ scale (area_p) + (1|pob/indiv), data = sub_df)
PBmodcomp(lm2_large_p, lm2_small_p, nsim = 10000)

fe <- function(fit) {return(fixef(fit))} ## ya esta mas arriba
bootlm2 = bootMer(lm2_large_p,FUN=fe,nsim=500,use.u =F)
bootlm2$t[,2]
par(mfrow=c(1,2),cex.lab=1.5, cex.axis=1.3)
## distribucion para masa  y para el componente cuadratico
plot(density(bootlm2$t[,2]), main="", xlab="pendientes  papus", ylab="Densidad", 
     lwd=3, col="red", xlim=c(-0.5,0.7))                                                     
abline(v=quantile(bootlm2$t[,2],0.025)) ## intervalo del 95 por ciento
abline(v=quantile(bootlm2$t[,2],0.975))


plot(density(bootlm2$t[,3]), main="", xlab="pendiente cuadratica simulada", ylab="Densidad", 
     lwd=3, col="red", xlim=c(-0.25,0.25))                                                     
abline(v=quantile(bootlm2$t[,3],0.025)) ## intervalo del 95 por ciento
abline(v=quantile(bootlm2$t[,3],0.975))
layout(1)
## en ambos casos parece no ser significativo, no incluye el 0, igual voy a probar el small

lm2_small_p = lmer(W ~ scale (area_p) + (1|pob/indiv), data = sub_df, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
lm2_vsmall_p = lmer(W ~ 1 + (1|pob/indiv), data = sub_df, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

PBmodcomp(lm2_small_p, lm2_vsmall_p,  nsim = 1000)
bootlm2_small = bootMer(lm2_small_p,FUN=fe,nsim=500,use.u =F)
plot(density(bootlm2_small$t[,2]), main="", xlab="pendientes  papus", ylab="Densidad", 
     lwd=3, col="red", xlim=c(-0.5,0.7))                                                     
abline(v=quantile(bootlm2$t[,2],0.025)) ## intervalo del 95 por ciento
abline(v=quantile(bootlm2$t[,2],0.975))
## el papus no es significativo en el modelo



## la masa es significativa ahora vamos a graficar 
pred_papus <- ggpredict(lm2_large_p, terms = "area_p [all]")

ggplot(pred_papus, aes(x = scale(x), y = predicted)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) + geom_point(data = sub_df, aes(x = scale(area_p), y = W), color = '#440154FF')+
  labs(title = "Efecto del papus y su efecto cuadrático",
       x = "area papus",
       y = "Fitness Predicho") +
  theme_minimal()

## se ve muy bien como esta influyendo el componente materno, el papus viene de la flor, por lo que los frutos vanos
## siguen teniendo papus.


## W ~ volumen de la cipsela #######################
lm3_large_vol = lmer(W ~ scale(vol_c) + I(scale(vol_c)^2) + (1|pob/indiv), data = sub_df)
summary(lm3_large_vol)
lm3_small_vol = lmer(W ~ scale (vol_c) + (1|pob/indiv), data = sub_df)
lm3_vsmall_vol = lmer(W ~ 1 + (1|pob/indiv), data = sub_df)
PBmodcomp(lm3_large_vol, lm3_small_vol, nsim = 10000) ## el componente cuadratico no da significativo
PBmodcomp(lm3_small_vol, lm3_vsmall_vol, nsim = 10000) ### super significativo.


bootlm3 = bootMer(lm3_small_vol,FUN=fe,nsim=500,use.u =F)

## distribucion para masa  y para el componente cuadratico
plot(density(bootlm3$t[,2]), main="", xlab="pendientes  vol", ylab="Densidad", 
     lwd=3, col="red", xlim=c(-0.5,0.1))                                                     
abline(v=quantile(bootlm3$t[,2],0.025)) ## intervalo del 95 por ciento
abline(v=quantile(bootlm3$t[,2],0.975)) ## no incluye el 0

pred_vol <- ggpredict(lm3_small_vol, terms = "vol_c")

ggplot(pred_vol, aes(x = scale(x), y = predicted)) +
  geom_line(color = "blue") + geom_point(data = sub_df, aes(x = scale(vol_c), y = W), color = viridis(1))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(title = "Efecto del volumen en la colonización (W)",
       x = "volumen de la cipsela",
       y = "Fitness Predicho") +
  theme_minimal()


## voy a ver los residuos por qu esto esta feo

residuals <- resid(lm3_small_vol)

qqplot_data <- data.frame(sample = residuals)
ggplot(qqplot_data, aes(sample = sample)) + 
  stat_qq() + 
  stat_qq_line() + 
  ggtitle("QQ Plot de los Residuos") + 
  xlab("Cuantiles Teóricos") + 
  ylab("Cuantiles Muestrales")
 ## bastante feo en las colas
fitted_values <- fitted(lm3_small_vol)

residuals_vs_fitted_data <- data.frame(fitted = fitted_values, residuals = residuals)
ggplot(residuals_vs_fitted_data, aes(x = fitted, y = residuals)) + 
  geom_point() + 
  geom_smooth(method = "loess", col = "red") + 
  ggtitle("Residuos vs Valores Ajustados") + 
  xlab("Valores Ajustados") + 
  ylab("Residuos")

## ok es feo pero claramente le faltan variables al modelo


## usando todas las variables (gradiente de seleccion) #################
sg_large_model = lmer(W ~ scale(masa)*scale(area_p) + scale(vol_c) + (1|pob/indiv), data = sub_df)
sg_small_model1 = lmer(W ~ scale(masa) + scale(area_p) + scale(vol_c) + (1|pob/indiv), data = sub_df) ## veo la significacia de la interaccion
PBmodcomp(sg_large_model, sg_small_model1, nsim = 10000) ## la interaccion es significativa
sg_small_model2 = lmer(W ~ scale(masa)*scale(area_p) + (1|pob/indiv), data = sub_df) ## veo la significacia del volumen
PBmodcomp(sg_large_model, sg_small_model2, nsim = 10000)
sg_small_model3 = lmer(W ~  scale(masa) + scale(masa):scale(area_p) + (1|pob/indiv), data = sub_df)
PBmodcomp(sg_small_model2, sg_small_model3, nsim = 10000) ## significancia del papus da re contra significativo el papus
sg_small_model4 = lmer(W ~  scale(area_p) + scale(masa):scale(area_p) + (1|pob/indiv), data = sub_df)
PBmodcomp(sg_small_model2, sg_small_model4, nsim = 10000)


## con estas variales el modelo final con significacia estadistica


sg_final_model = lmer(W ~ scale(masa)*scale(area_p) + I(scale(masa)^3)  + (1|pob/indiv), data = sub_df) 


mean(sub_df$masa)
quantile(sub_df$masa,0.75)
quantile(sub_df$masa,0.25)
## efectos 
pred_model_papus_m <- ggpredict(sg_final_model, terms = c("area_p [all]", "masa [1.572, 2.0525,0.955]" ))## efecto del papus cuando la masa esta en la media
## en el 3er cuartil y primer cuartil

ggplot(pred_model_papus_m, aes(x = x, y = predicted)) +
  geom_line(aes(group = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) + geom_point(data = sub_df, aes(x = area_p, y = clz, color = masa))+ 
  labs(title = "Efecto del área del papus según el nivel de masa",
       x = "Área del papus",
       y = "Fitness Predicho",
       color = "Nivel de Masa",
       fill = "Nivel de Masa") + scale_fill_viridis_d()  + scale_color_viridis_c() + theme_minimal()

mean(sub_df$area_p)
quantile(sub_df$area_p,0.75)
quantile(sub_df$area_p,0.1)
pred_model_masa_p <- ggpredict(sg_final_model, terms = c("masa [all]", "area_p[164.177, 183.2296, 100.5051 ]"))## efecto del papus cuando la masa esta en la media
## en el 3er cuartil y primer cuartil

ggplot(pred_model_masa_p, aes(x = x, y = predicted)) +
  geom_line(aes(group = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) + geom_point(data = sub_df, aes(x = masa, y = clz, color = area_p))+ 
  labs(title = "Efecto de la masa según el area del papus",
       x = "masa",
       y = "Fitness Predicho",
       color = "Área de papus",
       fill = "1er cuartil, media y 3er cuartil papus") + scale_fill_viridis_d()  + scale_color_viridis_c() + theme_minimal()



residuals <- resid(sg_final_model)

qqplot_data <- data.frame(sample = residuals)
r1 = ggplot(qqplot_data, aes(sample = sample)) + 
  stat_qq() + 
  stat_qq_line() + 
  ggtitle("QQ Plot de los Residuos") + 
  xlab("Cuantiles Teóricos") + 
  ylab("Cuantiles Muestrales")
## bastante feo en las colas
fitted_values <- fitted(sg_final_model)

residuals_vs_fitted_data <- data.frame(fitted = fitted_values, residuals = residuals)
r2 = ggplot(residuals_vs_fitted_data, aes(x = fitted, y = residuals)) + 
  geom_point() + 
  geom_smooth(method = "loess", col = "red") + 
  ggtitle("Residuos vs Valores Ajustados") + 
  xlab("Valores Ajustados") + 
  ylab("Residuos")


# residuos vs predictor
residuals_vs_predictor_data <- data.frame(masa = sub_df$masa, area_p = sub_df$area_p, residuals = residuals)
r3 = ggplot(residuals_vs_predictor_data, aes(x = masa, y = residuals)) + 
  geom_point() + 
  geom_smooth(method = "loess", col = "red") + 
  ggtitle("Residuos vs Predictor (masa)") + 
  xlab("masa") + 
  ylab("Residuos")

r4 = ggplot(residuals_vs_predictor_data, aes(x = area_p, y = residuals)) + 
  geom_point() + 
  geom_smooth(method = "loess", col = "red") + 
  ggtitle("Residuos vs Predictor (area_p)") + 
  xlab("area_p") + 
  ylab("Residuos")

library(influence.ME); library(HLMdiag)
influence_data <- influence(sg_final_model, obs = TRUE)
plot(influence_data, which = 'cook')

cd.lev1 = cooks.distance(sg_final_model, level=1)
plot(cd.lev1, type="S", xlab="dato", ylab="Distancia de Cook")
## no hay puntos influyentes
library(patchwork)
(r1+r2)/(r3+r4)
### claremente hay efectos que mi modelo no esta explicando bien




# PROTOCOLO SMBV ####################
sub_df_fr <- sub_df %>%
  dplyr::select(pob, indiv,fruto, sv, prob_germ, clz, pl_m, pl_v, masa, area_p, vol_c, capitulos, altura_planta) %>%
  group_by(pob, indiv,fruto) %>%
  summarise(across(c(sv, prob_germ, clz, pl_m, pl_v, masa, area_p, vol_c, capitulos, altura_planta), mean, .names = "{col}"))


mean(sub_df_fr$masa)
sub_df_fr$W = sub_df_fr$clz/mean(sub_df_fr$clz)

corr <- round(cor(na.omit(sub_df_fr[,c('masa', 'area_p', 'vol_c', 'capitulos', 'altura_planta')])), 1)

corr
ggcorrplot(corr, type = 'lower',outline.color = "white",
           ggtheme = ggplot2::theme_bw, colors = viridis(3), lab = T)




est_df = sub_df_fr
est_v = c('masa', 'area_p', 'vol_c', 'capitulos', 
                        'altura_planta', 'pl_m', 'pl_v')


for (col in est_v) {
  est_df[,col] <- as.vector(scale(est_df[,col]))
}

## OPORTUNIDAD DE SELECCION ###############
w.os = var(sub_df_fr$W) ## var 0.4349
sd(sub_df_fr$W)/mean(sub_df_fr$W) ### coeficiente de variacion alto 0.659
w.os/sqrt(2*length(sub_df_fr$W)) # error estandar de 0.0232

## DIFERENCIAL DE SELECCION LINEAL #######################
fixef_fun <- function(model) {
  return(fixef(model))
}
n_cores <- detectCores() - 1
conf_levels <- c(0.90, 0.95, 0.99, 0.999)

### masa ####
ds_masa = lmer(W ~ masa + (1|pob/indiv), data = est_df)
summary(ds_masa)
boot_results <- bootMer(ds_masa, fixef_fun, nsim = 10000, use.u = T, re.form = NULL, ncpus = n_cores)

ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ##***
### area_p #############
ds_areap = lmer(W ~ area_p + (1|pob/indiv), data = est_df)

boot_results <- bootMer(ds_areap, fixef_fun, nsim = 10000, use.u = T, re.form = NULL, ncpus = n_cores)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ##***

### vol #########
ds_volc = lmer(W ~ vol_c + (1|pob/indiv), data = est_df)

boot_results <- bootMer(ds_volc, fixef_fun, nsim = 10000, use.u = T, re.form = NULL, ncpus = n_cores)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ##**


## diferencial de seleccion no lineal ########
## masa ########
ds_masa2 = lmer(W ~ I(0.5*(masa^2)) + (1|pob/indiv), data = est_df)
summary(ds_masa2)

boot_results <- bootMer(ds_masa2, fixef_fun, nsim = 10000, use.u = T, re.form = NULL, ncpus = n_cores)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## significativo
print(ci)
## areap #######
ds_areap2 = lmer(W ~ I(0.5*(area_p^2)) + (1|pob/indiv), data = est_df)
summary(ds_areap2)

boot_results <- bootMer(ds_areap2, fixef_fun, nsim = 10000, use.u = T, re.form = NULL, ncpus = n_cores)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## no significativo
print(ci)
### vol #########
ds_volc2 = lmer(W ~  I(0.5*(vol_c^2)) + (1|pob/indiv), data = est_df)

boot_results <- bootMer(ds_volc2, fixef_fun, nsim = 10000, use.u = T, re.form = NULL, ncpus = n_cores)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## no significativo
print(ci)

beep()
## GRADIENTE DE SELECCION DIRECCIONAL ##########################
### modelo ###################
grs = lmer(W ~ area_p + masa + vol_c + (1|pob/indiv), data = est_df)
summary(grs)

### bootstrap para significancia
boot_results <- bootMer(grs, fixef_fun, nsim = 10000, use.u = T, re.form = NULL, ncpus = n_cores)

### intervalos de confianza #######

get_conf_intervals <- function(results, conf_levels) {
  conf_intervals <- list()
  for (i in 1:ncol(results$t)) {
    ci <- boot.ci(results, index = i, type = "perc", conf = conf_levels)
    conf_intervals[[names(results$t0)[i]]] <- ci ## nombro las sub listas como la variable
  }
  return(conf_intervals)
}
conf_levels <- c(0.90, 0.95, 0.99, 0.999)
conf_intervals <- get_conf_intervals(boot_results, conf_levels)


## GRADIENTE de SELECCION CUADRATICA y CORRELACIONAL ###########
### modelo #########################
grs_cc = lmer(W ~ area_p + masa + vol_c + I(0.5*(masa^2)) +  I(0.5*(area_p^2)) + 
                I(0.5*(vol_c^2)) + area_p:masa + masa:vol_c + area_p:vol_c  + 
                (1|pob/indiv), data = est_df, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
#summary(grs_cc)
#grs_cc = lm(W ~ area_p + masa + vol_c + I(0.5*(masa^2)) +  I(0.5*(area_p^2)) + 
#                I(0.5*(vol_c^2)) + area_p:masa + masa:vol_c + area_p:vol_c, data = est_df)

## creo la funcion necesaria para para remuestrear y obtener la distribucon del modelo

## otra forma es coon butmer de lme4 que es especifico para modelos mixtos, respeta la estructura del moel ## 

### si no usara efecto aleatroio
#coef_fun <- function(data, indices) {
#
#  d <- data[indices, ]
#
#  model <- lmer(W ~ area_p + masa + vol_c + I(0.5*(masa^2)) + 
#                  I(0.5*(area_p^2)) + I(0.5*(vol_c^2)) + 
#                  area_p:masa + masa:vol_c + area_p:vol_c + 
#                  (1|pob/indiv), data = d, 
#                control = lmerControl(optimizer = "bobyqa", 
#                                      optCtrl = list(maxfun = 100000)))
#  return(fixef(model))  ## obtengo los coeficientes del modelo en cada corrida
#}
#coef_fun <- function(data, indices) {
#  
#  d <- data[indices, ]
#  
#  model <- lm(W ~ area_p + masa + vol_c + I(0.5*(masa^2)) + 
#                I(0.5*(area_p^2)) + I(0.5*(vol_c^2)) + 
#                area_p:masa + masa:vol_c + area_p:vol_c, data = d)
#  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
#}

### boostrap para significancia ###################
boot_results <- bootMer(grs_cc, fixef_fun, nsim = 10000, use.u = T, re.form = NULL, ncpus = n_cores)
#n_cores <- detectCores() - 1
#cl <- makeCluster(n_cores)
#set.seed(123) 
#results <- boot(data = est_df, statistic = coef_fun, R = 10000)
#stopCluster(cl)

### intervalos de confianza #######
conf_levels <- c(0.90, 0.95, 0.99, 0.999)
conf_intervals <- get_conf_intervals(boot_results, conf_levels)




library(plotly)
##  efectos de la seleccion cuadratica y correlacional ########

grid_data <- expand.grid(
  masa = seq(min(est_df$masa), max(est_df$masa), length = 100),
  area_p = seq(min(est_df$area_p), max(est_df$area_p), length = 100)
)
grid_data$vol_c = 0 ## para el volumen en la media
grid_data$predicted <- predict(grs_cc, newdata = grid_data, re.form = NA)


p <- ggplot(grid_data, aes(x = masa, y = area_p, z = predicted)) +
  geom_tile(aes(fill = predicted)) +
  geom_contour(color = "white") +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(x = "Masa de la diáspora", y = "Área del papus", fill = "Colonización")
p
plotly::ggplotly(p)
masa_seq = seq(min(est_df$masa), max(est_df$masa), length = 100)
area_p_seq = seq(min(est_df$area_p), max(est_df$area_p), length = 100)
predicted <- matrix(grid_data$predicted, nrow = length(masa_seq), ncol = length(area_p_seq))
plot_ly(
  x = ~masa_seq,
  y = ~area_p_seq,
  z = ~predicted,
  type = "surface",
  colorscale = 'Viridis',
  colorbar = list(title = "Colonización\npredicha")
) %>%
  layout(
    scene = list(
      xaxis = list(title = "Masa"),
      yaxis = list(title = "Área del papus"),
      zaxis = list(title = "Colonización predicha")
    )
  )



## volumen y area
grid_data <- expand.grid(
  vol_c = seq(min(est_df$vol_c), max(est_df$vol_c), length = 100),
  area_p = seq(min(est_df$area_p), max(est_df$area_p), length = 100)
)
grid_data$masa = 0 ## para el volumen en la media
grid_data$predicted <- predict(grs_cc, newdata = grid_data, re.form = NA)


p <- ggplot(grid_data, aes(x = vol_c, y = area_p, z = predicted)) +
  geom_tile(aes(fill = predicted)) +
  geom_contour(color = "white") +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "gradiente de Selección correlacional y no lineal", x = "Vol cipsela", y = "Área de papus", fill = "Fitness")
p
plotly::ggplotly(p)


plot_ly(
  x = ~grid_data$vol_c,
  y = ~grid_data$area_p,
  z = ~grid_data$predicted,
  type = "scatter3d",
  mode = "lines",
  line = list(width = 3, color = ~grid_data$predicted, colorscale = 'Viridis')
) %>%
  layout(
    title = "Gradiente de seleccion correlacional y no lineal",
    scene = list(
      xaxis = list(title = "Volumen cipsela"),
      yaxis = list(title = "area del pappus"),
      zaxis = list(title = "Fitness")
    )
  )



## GAMs gradientes no lineales y superficie de seleccion ####
plot_model = function(model, x, y, xname) {
  newx = seq(min(x), max(x), length.out = 100)
  
  newdata = data.frame(newx)
  colnames(newdata) <- xname
  
  newy = predict(model, newdata = newdata, se.fit = TRUE, type = 'response')
  
  yhat = newy$fit
  yup = newy$fit + newy$se.fit
  ydown = newy$fit - newy$se.fit
  
  df = data.frame(yhat = yhat, newx = newx, yup = yup, ydown = ydown)
  
  p <- ggplot(df, aes(x = newx, y =  yhat)) 
  p <- p + geom_ribbon(aes(ymin = ydown, ymax = yup), 
                        alpha = .15, color = 'grey', fill = 'grey', linetype = 'dashed') + geom_line(aes(y = yhat), 
                        linewidth = 0.3) + xlab(xname)  + ylab('colonización')
  return(p)
}
### masa ########
fit = gam(sub_df_fr$W ~ s(masa, bs = 'cr'), data = sub_df_fr, method = 'GCV.Cp')
p = plot_model(fit, sub_df_fr$masa, sub_df_fr$W,'masa')
p_m = p + geom_point(data = sub_df_fr, aes(x = masa , y = W), color = 'grey14', alpha = 0.3) + 
  scale_x_continuous(breaks = seq(0,3,0.5)) + xlab('masa (mg)') + theme_minimal()
p_m = p_m + geom_vline(aes(xintercept= mean(sub_df_fr$masa)), 
                 linetype = 'dashed', color = 'green') 
### papus ###### 
fit = gam(W ~ s(area_p, bs = 'cr'), data = sub_df_fr, method = 'GCV.Cp')
p = plot_model(fit, sub_df_fr$area_p, sub_df_fr$W,'area_p')
p_p = p + geom_point(data = sub_df_fr, aes(x = area_p , y = W, color = prob_germ), alpha = 0.7) + 
  labs(x = 'área del papus (mm2)', color = 'prob\ngerminación') + scale_x_continuous(breaks = seq(100,250,25)) + scale_color_viridis_c(option = 'A') + theme_minimal() 
p_p = p_p + geom_vline(aes(xintercept= mean(sub_df_fr$area_p)), 
                 linetype = 'dashed', color = 'green')

wrap_plots(p_m + ggtitle("(a)") +
             p_p + ggtitle("(b)"))
### volumen
fit = gam(W ~ s(vol_c, bs = 'cr'), data = est_df, method = 'GCV.Cp')
p = plot_model(fit, est_df$vol_c, est_df$W,'vol_c')
p + geom_point(data = est_df, aes(x = vol_c , y = W), color = 'grey14', alpha = 0.3) + xlab('volumen de la cipsela') + theme_minimal()


# PROMEDIOS INDIVIDUOS  ###############################
df_sum <- sub_df %>% dplyr::select(pob, indiv, sv, prob_germ, clz, pl_m, pl_v, masa, area_p, vol_c, capitulos, altura_planta ) %>% group_by(pob,indiv)  %>% summarise_all(lst(mean,sd))# marcamos las variables que agrupan
colnames(df_sum)
sub_df_sum = df_sum[,c(1:12)]

colnames(sub_df_sum) = c('pob', 'indiv', 'sv', 'prob_germ', 'clz', 'pl_m', 'pl_v',
                         'masa', 'area_p', 'vol_c', 'capitulos', 'altura_planta')


corr <- round(cor(na.omit(sub_df_sum[,c('masa', 'area_p', 'vol_c', 'capitulos', 'altura_planta')])), 1)

corr
ggcorrplot(corr, type = 'lower',outline.color = "white",
           ggtheme = ggplot2::theme_bw, colors = viridis(3), lab = T)


est_df = sub_df_sum ## lo nombro igual para facilitar la repeticion
est_v = c('masa', 'area_p', 'vol_c', 'capitulos', 
          'altura_planta', 'pl_m', 'pl_v')

for (col in est_v) {
  est_df[,col] <- as.vector(scale(est_df[,col]))
}

## estandarizacion del fitness relativo. #####################################
sub_df_sum$W  = est_df$W = sub_df_sum$clz/mean(sub_df_sum$clz)
## OPORTUNIDAD DE SELECCION ###############
w.os = var(sub_df_sum$W) ## var 0.269, obviamente disminuye la varianza
sd(sub_df_sum$W)/mean(sub_df_sum$W) ### coeficiente de variacion alto 0.518
w.os/sqrt(2*length(sub_df_sum$W)) # error estandar de 0.0248

## DIFERENCIAL DE SELECCION LINEAL #######################
### masa ####
ds_masa = lmer(W ~ masa + (1|pob), data = est_df)
summary(ds_masa)
ds_1 = lmer(W ~ 1 + (1|pob), data = est_df)
anova(ds_masa, ds_1) ## ahora uso esto para ascelerar tiempos ##  < 2.2e-16 *** masa re contra significativa
### area_p #############
ds_areap = lmer(W ~ area_p + (1|pob), data = est_df)
anova(ds_areap, ds_1) ## significativo 
### vol #########
ds_volc = lmer(W ~ vol_c + (1|pob), data = est_df)
anova(ds_volc, ds_1) ###### significativo
### capitulos ##################
ds_cap = lmer(W ~ capitulos + (1|pob), data = est_df)
anova(ds_cap, ds_1) ## no significativo
### altura planta ######
ds_altp = lmer(W ~ altura_planta + (1|pob), data = est_df)
summary(ds_altp)
anova(ds_altp, ds_1)  ## signficativa (marginalmente 0.04949):O


## diferencial de seleccion no lineal ########
### masa ########
ds_masa2 = lmer(W ~ I(0.5*(masa^2)) + (1|pob), data = est_df)
summary(ds_masa2)
anova(ds_masa2, ds_1) ## el componente cuadratico es no significativo
### areap #######
ds_areap2 = lmer(W ~ I(0.5*(area_p^2)) + (1|pob), data = est_df)
summary(ds_areap2)
anova(ds_areap2, ds_1) ### no significativo
### vol #########
ds_volc2 = lmer(W ~  I(0.5*(vol_c^2)) + (1|pob), data = est_df)
anova(ds_volc2, ds_1) ###### no significativo (0.09)
### capitulos ##################
ds_cap = lmer(W ~  I(0.5*(capitulos^2)) + (1|pob), data = est_df)
anova(ds_cap, ds_1) ## no significativo
### altura planta ######
ds_altp = lmer(W ~ I(0.5*(altura_planta^2)) + (1|pob), data = est_df)
anova(ds_altp, ds_1)  ## no signficativa


## GRADIENTE DE SELECCION DIRECCIONAL ##########################
grs = lmer(W ~ area_p + masa + vol_c + capitulos + altura_planta + (1|pob), data = est_df)
summary(grs)
### pappus ########
grs_s1 = lmer(W ~  masa + vol_c + capitulos + altura_planta + (1|pob), data = est_df)
anova(grs,grs_s1) ### significativo 0.00122

### masa ##########
grs_s2 = lmer(W ~  area_p + vol_c + capitulos + altura_planta + (1|pob), data = est_df)
anova(grs,grs_s2) ### significativo p < 2.2e-16 ***

### vol ##########
grs_s3 = lmer(W ~  area_p + masa + capitulos + altura_planta + (1|pob), data = est_df)
anova(grs,grs_s3) ### no significativo 0.1475

### capitulos ##########
grs_s3 = lmer(W ~  area_p + masa + vol_c + altura_planta + (1|pob), data = est_df)
anova(grs,grs_s3) ### no significativo 0.4108

### altura planta ##########
grs_s3 = lmer(W ~  area_p + masa + vol_c + capitulos + (1|pob), data = est_df)
anova(grs,grs_s3) ### no significativo 0.2279

## GRADIENTE de SELECCION CUADRATICA y CORRELACIONAL ########

grs_cc = lm(W ~ area_p + masa + vol_c + I(0.5*(masa^2)) + I(0.5*(area_p^2)) + 
                I(0.5*(vol_c^2)) + I(0.5*(capitulos^2)) + I(0.5*(altura_planta^2)) + 
                area_p:masa + masa:vol_c + area_p:vol_c + area_p:capitulos + area_p:altura_planta + 
                masa:capitulos + masa:altura_planta + vol_c:capitulos + 
                vol_c:altura_planta + capitulos:altura_planta, data = est_df) ## no tiene efecto de la poblacion


summary(grs_cc)


# GERMINADAS ####################
var <- sub_df %>%
  filter(germ == 1) %>%
  dplyr::select(pob, indiv,fruto, sv, prob_germ, clz, pl_m, pl_v, masa, area_p, vol_c, capitulos, altura_planta) %>%
  group_by(pob, indiv, fruto) %>%
  summarise(across(c(sv, prob_germ, clz, pl_m, pl_v, masa, area_p, vol_c, capitulos, altura_planta), mean, .names = "{col}"))


## varianza explicada ####
m1 = lmer(masa ~ (1|pob/indiv), data = var)
summary(m1)
m2 = lmer(vol_c ~ (1|pob/indiv), data = var)
summary(m2)
m3 = lmer(area_p ~ (1|pob/indiv), data = var)
summary(m3)

## data set final ######
sub_df_sum <- sub_df %>%
  filter(germ == 1) %>%
  dplyr::select(pob, indiv,fruto, sv, prob_germ, clz, pl_m, pl_v, masa, area_p, vol_c, capitulos, altura_planta) %>%
  group_by(pob, indiv) %>%
  summarise(across(c(sv, prob_germ, clz, pl_m, pl_v, masa, area_p, vol_c, capitulos, altura_planta), mean, .names = "{col}"),
            n = n_distinct(fruto))

mean(sub_df_sum$clz)
var(sub_df_sum$clz)
corr <- round(cor(na.omit(sub_df[,c('masa', 'area_p', 'vol_c', 'capitulos', 'altura_planta')])), 1)

corr
ggcorrplot(corr, type = 'lower',outline.color = "white",
           ggtheme = ggplot2::theme_bw, colors = viridis(3), lab = T)

est_df = sub_df_sum ## lo nombro igual para facilitar la repeticion
est_v = c('masa', 'area_p', 'vol_c', 'capitulos', 
          'altura_planta', 'pl_m', 'pl_v')

for (col in est_v) {
  est_df[,col] <- as.vector(scale(est_df[,col]))
}

## estandarizacion del fitness relativo. #####################################
sub_df_sum$W  = est_df$W = sub_df_sum$clz/mean(sub_df_sum$clz)
## OPORTUNIDAD DE SELECCION ###############
w.os = var(sub_df_sum$W) ## var 0.031, obviamente disminuye la varianza
sd(sub_df_sum$W)/mean(sub_df_sum$W) ### coeficiente de variacion alto 0.17
w.os/sqrt(2*length(sub_df_sum$W)) # error estandar de 0.003065849

## efecto de la poblacion para simplificar #######
grs_cc = lmer(W ~ area_p + masa + vol_c + altura_planta + capitulos  + I(0.5*(masa^2)) + I(0.5*(area_p^2)) + 
              I(0.5*(vol_c^2)) + I(0.5*(capitulos^2)) + I(0.5*(altura_planta^2)) + 
              area_p:masa + masa:vol_c + area_p:vol_c + area_p:capitulos + area_p:altura_planta + 
              masa:capitulos + masa:altura_planta + vol_c:capitulos + 
              vol_c:altura_planta + capitulos:altura_planta + (1|pob),data = est_df)
grs_ccs = lm(W ~ area_p + masa + vol_c + altura_planta + capitulos  + I(0.5*(masa^2)) + I(0.5*(area_p^2)) + 
                I(0.5*(vol_c^2)) + I(0.5*(capitulos^2)) + I(0.5*(altura_planta^2)) + 
                area_p:masa + masa:vol_c + area_p:vol_c + area_p:capitulos + area_p:altura_planta + 
                masa:capitulos + masa:altura_planta + vol_c:capitulos + 
                vol_c:altura_planta + capitulos:altura_planta ,data = est_df)

anova(grs_cc, grs_ccs) ## no hay efecto de la poblacion


  

## DIFERENCIAL DE SELECCION LINEAL #######################
### masa ####
ds_masa = lm(W ~ masa, data = est_df)
#summary(ds_masa)
coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ masa, data = d)
  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}
cl <- makeCluster(n_cores)
set.seed(123) 
boot_results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ##  significativo
print(ci)
### area_p #############
ds_areap = lm(W ~ area_p, data = est_df)
#summary(ds_areap)## significativo

coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ area_p, data = d)
  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}
cl <- makeCluster(n_cores)
set.seed(123) 
boot_results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## significativo
print(ci)

### vol #########
ds_volc = lm(W ~ vol_c, data = est_df)
#summary(ds_volc)
coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ vol_c, data = d)
  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}
cl <- makeCluster(n_cores)
set.seed(123) 
boot_results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## no significativo
print(ci)
### capitulos ##################
ds_cap = lm(W ~ capitulos , data = est_df)
#summary(ds_cap)
coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ capitulos, data = d)
  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}
cl <- makeCluster(n_cores)
set.seed(123) 
boot_results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## no significativo
print(ci)
### altura planta ######
ds_altp = lm(W ~ altura_planta, data = est_df)
#summary(ds_altp)  ## 
coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ altura_planta, data = d)
  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}
cl <- makeCluster(n_cores)
set.seed(123) 
boot_results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## no significativo
print(ci)

beep()
## diferencial de seleccion no lineal ########
### masa ########
ds_masa2 = lm(W ~ I(0.5*(masa^2)), data = est_df)
#summary(ds_masa2)

coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ I(0.5*(masa^2)), data = d)
  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}
cl <- makeCluster(n_cores)
set.seed(123) 
boot_results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## **significativo
print(ci)
 ## el componente cuadratico es significativo 0.002393
### areap #######
ds_areap2 = lm(W ~ I(0.5*(area_p^2)), data = est_df)
#summary(ds_areap2)

coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ I(0.5*(area_p^2)), data = d)
  return(coef(model))  
}
cl <- makeCluster(n_cores)
set.seed(123) 
boot_results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## no significativo
print(ci)

### no significativo
### vol #########
ds_volc2 = lm(W ~  I(0.5*(vol_c^2)) , data = est_df)

coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ I(0.5*(vol_c^2)), data = d)
  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}
cl <- makeCluster(n_cores)
set.seed(123) 
boot_results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## no significativo
print(ci)
### capitulos ##################
ds_cap2 = lm(W ~  I(0.5*(capitulos^2)), data = est_df)
#summary(ds_cap2) ## no significativo
coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ I(0.5*(capitulos^2)), data = d)
  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}
cl <- makeCluster(n_cores)
set.seed(123) 
boot_results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## *significativo
print(ci)
### altura planta ######
ds_altp2 = lm(W ~ I(0.5*(altura_planta^2)), data = est_df)
#summary(ds_altp2)
coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ I(0.5*(altura_planta^2)), data = d)
  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}
cl <- makeCluster(n_cores)
set.seed(123) 
boot_results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
ci <- boot.ci(boot_results, index = 2, type = "perc", conf = conf_levels) ## no significativo
print(ci)

beep()
## GRADIENTE DE SELECCION DIRECCIONAL ##########################
grs = lm(W ~ area_p + masa  + altura_planta , data = est_df)
summary(grs)
### boostrap para significancia #######
coef_fun <- function(data, indices) {
    
    d <- data[indices, ]
    
    model <- lm(W ~ area_p + masa + altura_planta, 
                data = d)
    return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}
cl <- makeCluster(n_cores)
set.seed(123) 
results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)

### intervalos de confianza #######
conf_levels <- c(0.90, 0.95, 0.99, 0.999)
conf_intervals <- get_conf_intervals(results, conf_levels)

## GRADIENTE de SELECCION CUADRATICA y CORRELACIONAL ########


grs_cc = lm(W ~ area_p + masa + altura_planta + I(0.5*(masa^2)) +  I(0.5*(area_p^2)) + 
              I(0.5*(altura_planta^2)) + area_p:masa + masa:altura_planta + area_p:altura_planta, data = est_df)
summary(grs_cc)
### boostrap para significancia ######

coef_fun <- function(data, indices) {
  
  d <- data[indices, ]
  
  model <- lm(W ~ area_p + masa + altura_planta + I(0.5*(masa^2)) +  I(0.5*(area_p^2)) + 
                I(0.5*(altura_planta^2)) + area_p:masa + masa:altura_planta + area_p:altura_planta, data = d)
  
  return(coef(model))  ## obtengo los coeficientes del modelo en cada corrida
}



cl <- makeCluster(n_cores)
set.seed(123) 
results <- boot(data = est_df, statistic = coef_fun, R = 10000)
stopCluster(cl)
### intervalo sde confianza ######
conf_levels <- c(0.90, 0.95, 0.99, 0.999)
conf_intervals <- get_conf_intervals(results, conf_levels)

##  efectos de la seleccion cuadratica y correlacional ########

area_p_seq <- seq(from = min(est_df$area_p), to = max(est_df$area_p), length.out = 100)
masa_seq <- seq(from = min(est_df$masa), to = max(est_df$masa), length.out = 100)

grid_data <- expand.grid(area_p = area_p_seq, masa = masa_seq, altura_planta = 0)

grid_data$predicted <- predict(grs_cc, newdata = grid_data)

p <- ggplot(grid_data, aes(x = masa, y = area_p, z = predicted)) +
  geom_tile(aes(fill = predicted)) +
  geom_contour(color = "white") +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs( x = "Masa", y = "Área del papus", fill = "Colonización predicha")
p
#pplotly::ggplotly(p)

predicted <- matrix(grid_data$predicted, nrow = length(masa_seq), ncol = length(area_p_seq))

plot_ly(
  x = ~masa_seq,
  y = ~area_p_seq,
  z = ~predicted,
  type = "surface",
  colorscale = 'Viridis',
  colorbar = list(title = "Colonización\npredicha")
) %>%
  layout(
    scene = list(
      xaxis = list(title = "Masa"),
      yaxis = list(title = "Área del papus"),
      zaxis = list(title = "Colonización predicha")
    )
  )
### gams #####

plot_model = function(model, x, y, xname) {
  newx = seq(min(x), max(x), length.out = 100)
  
  newdata = data.frame(newx)
  colnames(newdata) <- xname
  
  newy = predict(model, newdata = newdata, se.fit = TRUE, type = 'response')
  
  yhat = newy$fit
  yup = newy$fit + newy$se.fit
  ydown = newy$fit - newy$se.fit
  
  df = data.frame(yhat = yhat, newx = newx, yup = yup, ydown = ydown)
  
  p <- ggplot(df, aes(x = newx, y =  yhat)) 
  p <- p + geom_ribbon(aes(ymin = ydown, ymax = yup), 
                       alpha = .15, color = 'grey', 
                       fill = 'grey', linetype = 'dashed') + 
            geom_line(aes(y = yhat), linewidth = 0.3) + xlab(xname)  + 
    ylab('colonización')
  return(p)
}
### masa ########
fit = gam(W ~ s(masa, bs = 'cr'), data = sub_df_sum, method = 'GCV.Cp')
p = plot_model(fit, sub_df_sum$masa, sub_df_sum$W,'masa')
p_m = p + geom_point(data = sub_df_sum, aes(x = masa , y = W, color = 1/sv), 
                     alpha = 0.7) + scale_x_continuous(breaks = seq(0,3,0.5)) + 
          xlab('masa (mg)') + theme_minimal()
p_m = p_m + geom_vline(aes(xintercept = mean(sub_df_sum$masa)), 
                 linetype = 'dashed', color = 'green')
### papus ###### 
#sub_df_sum  = sub_df_sum[sub_df_sum$area_p < max(sub_df_sum$area_p),]
fit = gam(W ~ s(area_p, bs = 'cr'), data = sub_df_sum, method = 'GCV.Cp')
p = plot_model(fit, sub_df_sum$area_p, sub_df_sum$W,'area_p')
p_p = p + geom_point(data = sub_df_sum, aes(x = area_p , y = W, color = prob_germ), alpha = 0.7) + 
  labs(x = 'área del papus (mm2)', color = 'prob\ngerminación') + scale_x_continuous(breaks = seq(100,250,25)) + scale_color_viridis_c(option = 'A') + theme_minimal() 
p_p = p_p + geom_vline(aes(xintercept = mean(sub_df_sum$area_p)), 
                       linetype = 'dashed', color = 'green')
wrap_plots(p_m + ggtitle("(a)") +
             p_p + ggtitle("(b)"))





x_eval <- seq(min(sub_df_sum$masa), max(sub_df_sum$masa), length.out = 1000)

# Predecir valores y derivadas
predicciones <- predict(fit, newdata = data.frame(masa = x_eval), type = "response")

# Usar `deriv` para calcular las derivadas primera y segunda

# Calcular la derivada primera usando diferencias finitas
h <- diff(x_eval)[1] # Asumimos que los puntos de x_eval están equiespaciados
derivada_primera <- diff(predicciones) / h

# Calcular la derivada segunda usando diferencias finitas en la derivada primera
derivada_segunda <- diff(derivada_primera) / h

# Ajustar los vectores para tener la misma longitud que x_eval
derivada_primera <- c(NA, derivada_primera) # Añadir NA al principio para igualar longitud
derivada_segunda <- c(NA, NA, derivada_segunda) # Añadir NAs al principio para igualar longitud

# Encontrar puntos de inflexión donde la derivada segunda cambia de signo
puntos_inflexion <- x_eval[which(diff(sign(derivada_segunda)) != 0) + 2] # Ajustar índice debido a los NAs

# Visualización
plot(sub_df_sum$masa, sub_df_sum$W, main = "Spline cúbica y sus derivadas", xlab = "area_p", ylab = "W")
lines(x_eval, predicciones, col = "blue", lwd = 2, lty = 1)
lines(x_eval, derivada_primera, col = "green", lwd = 2, lty = 2)
lines(x_eval, derivada_segunda, col = "red", lwd = 2, lty = 3)
abline(h = 0, col = "gray", lwd = 0.5)
points(puntos_inflexion, predict(fit, newdata = data.frame(masa = puntos_inflexion)), col = "red", pch = 19)
legend("topright", legend = c("Spline cúbica", "Derivada primera", "Derivada segunda", "Puntos de inflexión"),
       col = c("blue", "green", "red", "red"), lty = c(1, 2, 3, NA), lwd = 2, pch = c(NA, NA, NA, 19))

# Imprimir los puntos de inflexión
puntos_inflexion








x_eval <- seq(min(sub_df_sum$area_p), max(sub_df_sum$area_p), length.out = 1000)

# Predecir valores y derivadas
predicciones <- predict(fit, newdata = data.frame(area_p = x_eval), type = "response")

# Usar `deriv` para calcular las derivadas primera y segunda

# Calcular la derivada primera usando diferencias finitas
h <- diff(x_eval)[1] # Asumimos que los puntos de x_eval están equiespaciados
derivada_primera <- diff(predicciones) / h

# Calcular la derivada segunda usando diferencias finitas en la derivada primera
derivada_segunda <- diff(derivada_primera) / h

# Ajustar los vectores para tener la misma longitud que x_eval
derivada_primera <- c(NA, derivada_primera) # Añadir NA al principio para igualar longitud
derivada_segunda <- c(NA, NA, derivada_segunda) # Añadir NAs al principio para igualar longitud

# Encontrar puntos de inflexión donde la derivada segunda cambia de signo
puntos_inflexion <- x_eval[which(diff(sign(derivada_segunda)) != 0) + 2] # Ajustar índice debido a los NAs

# Visualización
plot(sub_df_sum$area_p, sub_df_sum$W, main = "Spline cúbica y sus derivadas", xlab = "area_p", ylab = "W")
lines(x_eval, predicciones, col = "blue", lwd = 2, lty = 1)
lines(x_eval, derivada_primera, col = "green", lwd = 2, lty = 2)
lines(x_eval, derivada_segunda, col = "red", lwd = 2, lty = 3)
abline(h = 0, col = "gray", lwd = 0.5)
points(puntos_inflexion, predict(fit, newdata = data.frame(area_p = puntos_inflexion)), col = "red", pch = 19)
legend("topright", legend = c("Spline cúbica", "Derivada primera", "Derivada segunda", "Puntos de inflexión"),
       col = c("blue", "green", "red", "red"), lty = c(1, 2, 3, NA), lwd = 2, pch = c(NA, NA, NA, 19))

# Imprimir los puntos de inflexión
puntos_inflexion
