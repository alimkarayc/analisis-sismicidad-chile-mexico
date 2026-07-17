# Librerias ----
library(dplyr)
library(janitor)
library(sf)
library(units)
library(lubridate)
library(rnaturalearth)
library(ggplot2)
library(tseries)
library(geosphere)
library(MASS)
library(car)
library(pROC)
library(FactoMineR)
library(factoextra)
library(tidyr)
library(forcats)
library(ggrepel)

# Bases de datos ----

# EDIT THIS PATH BEFORE RUNNING THE SCRIPT!!!
setwd("C:/path/to/analisis-sismicidad-chile-mexico/data/raw")

CN  <- read.csv(here("data", "raw", "Cocos-Norte.csv")) # Zona interacción placa de Cocos - Norteamericana
CC  <- read.csv(here("data", "raw", "Cocos-Caribe.csv")) # Zona interacción placa de Cocos - Caribe
CH1 <- read.csv(here("data", "raw", "Chile Central 1.csv")) # Chile Central Area Norte
CH2 <- read.csv(here("data", "raw", "Chile Central 2.csv")) # Chile Central Area Sur

# Base por zonas ----

# Asignar la zona correcta ANTES de unir
CN = CN %>% mutate(zona_detallada = "México (Cocos-Norteamérica)") %>%
  mutate(zona = "México-Centroamérica")
CC = CC %>% mutate(zona_detallada = "Centroamérica (Cocos-Caribe)") %>%
  mutate(zona = "México-Centroamérica")
chile = rbind(CH1, CH2)
chile = chile %>% mutate(zona_detallada = "Chile (Nazca-Sudamérica)") %>%
  mutate(zona = "Chile Central")

# Base completa ----

camer = rbind(CN, CC)
base = rbind(chile, camer)

## Homogenizar magnitudes ----

homogenizar_mw = function(df){
  
  df$mw_hom = NA
  
  # Mw ya existente
  idx = df$magType %in% c("mwb","mwc","mwr","mww","mw")
  df$mw_hom[idx] = df$mag[idx]
  
  # mb -> Mw (ecuación IV)
  idx = df$magType == "mb"
  df$mw_hom[idx] = 1.04*df$mag[idx] - 0.02
  
  # ML -> Mw
  idx = df$magType == "ml" & df$depth <= 50
  df$mw_hom[idx] = 0.80*df$mag[idx] + 1.15
  
  idx = df$magType == "ml" & df$depth > 50
  df$mw_hom[idx] = 0.94*df$mag[idx] + 0.30
  
  # Ms -> Mw
  idx = df$magType == "ms"
  df$mw_hom[idx] = 0.74*df$mag[idx] + 1.60
  
  return(df)
}

base = homogenizar_mw(base)
chile = homogenizar_mw(chile)
camer = homogenizar_mw(camer)

# Magnitud de completitud ----

mcchile = na.omit(chile$mw_hom)
mccamer = na.omit(camer$mw_hom)

# Goodness-of-fit test (GFT)
# Wiemer & Wyss (2000)
# Implementación reproducida por Mignan & Woessner (2012)

# FUNCTIONS
fmd = function(mag,mbin){
  mi = seq(min(round(mag/mbin)*mbin), max(round(mag/mbin)*mbin), mbin)
  nbm = length(mi)
  cumnbmag = numeric(nbm)
  nbmag = numeric(nbm)
  for(i in 1:nbm) cumnbmag[i] = length(which(mag > mi[i]-mbin/2))
  cumnbmagtmp = c(cumnbmag,0)
  nbmag = abs(diff(cumnbmagtmp))
  res = list(m=mi, cum=cumnbmag, noncum=nbmag)
  return(res)
}

# Maximum Curvature (MAXC) [e.g., Wiemer & Wyss, 2000]
maxc = function(mag,mbin){
  FMD = fmd(mag,mbin)
  Mc = FMD$m[which(FMD$noncum == max(FMD$noncum))[1]]
  return(list(Mc=Mc))
}

# Goodness-of-fit test (GFT) [Wiemer & Wyss, 2000]
gft = function(mag,mbin){
  FMD = fmd(mag,mbin)
  McBound = maxc(mag,mbin)$Mc
  Mco = McBound-0.4+(seq(15)-1)/10
  R = numeric(15)
  for(i in 1:15){
    indmag = which(mag > Mco[i]-mbin/2)
    b = log10(exp(1))/(mean(mag[indmag])-(Mco[i]-mbin/2))
    a = log10(length(indmag))+b*Mco[i]
    FMDcum_model = 10^(a-b*FMD$m)
    indmi = which(FMD$m >= Mco[i])
    R[i] = sum(abs(FMD$cum[indmi]-FMDcum_model[indmi]))/sum(FMD$cum[indmi])*100
    #in Wiemer&Wyss [2000]: 100-R
  }
  indGFT = which(R <= 5) #95% confidence
  if(length(indGFT) != 0){
    Mc = Mco[indGFT[1]]
    best = "95%"
  } else{
    indGFT = which(R <= 10) #90% confidence
    if(length(indGFT) != 0){
      Mc = Mco[indGFT[1]]
      best = "90%"
    } else{
      Mc = McBound
      best = "MAXC"
    }
  }
  return(list(Mc=Mc, best=best, Mco=Mco, R=R))
}

res_chile = gft(mcchile,0.1)
res_camer = gft(mccamer,0.1)

res_chile
res_camer

#plot(res_chile$Mco,100-res_chile$R,type="b")
#abline(h=95,lty=2)
#abline(h=90,lty=3)

#plot(res_camer$Mco,100-res_camer$R,type="b")
#abline(h=95,lty=2)
#abline(h=90,lty=3)

mc_chile = res_chile$Mc # 4.8
mc_camer = res_camer$Mc # 5.1

# 5.1 Magnitud de completitud

## Filtrar por Mc ----

base = base %>%
  filter(mw_hom >= 5.1)

# Definición zonas de estudio ----

sismos_placa = function(ruta_poligono, datos_sismos, nombre_placa = "la placa") {
  
  # 1. Leer el polígono espacial
  # quiet = TRUE evita que st_read imprima texto innecesario en la consola cada vez
  poligono = st_read(ruta_poligono, quiet = TRUE) 
  
  # 2. Convertir el catalogo a formato espacial
  sismos_espaciales = st_as_sf(datos_sismos, 
                                coords = c("longitude", "latitude"), 
                                crs = 4326)
  
  # 3. Aplicar el filtro espacial (Intersección)
  sismos_filtrados = st_intersection(sismos_espaciales, poligono)
  
  # 4. catalogar
  sismos_filtrados = sismos_filtrados %>%
    mutate(placa = nombre_placa)
  
  # 5. Retornar los datos filtrados para poder guardarlos
  return(sismos_filtrados)
}

# Placa Rivera
sismos_placa_rivera = sismos_placa("Placa_Rivera.kml", camer, "Rivera")

sismos_placa_rivera = sismos_placa_rivera %>%
  mutate(tipo = "Oceanica")

# Placa Panama
sismos_placa_panama = sismos_placa("Placa_Panama.kml", camer, "Panama")

sismos_placa_panama = sismos_placa_panama %>%
  mutate(tipo = "Continental")

# Placa de Cocos
sismos_placa_cocos = sismos_placa("Placa_Cocos.kml", camer, "Cocos")

sismos_placa_cocos = sismos_placa_cocos %>%
  mutate(tipo = "Oceanica")

# Placa Norteamericana
sismos_placa_norte = sismos_placa("Placa_Norteamericana.kml", camer, "Norteamericana")

sismos_placa_norte = sismos_placa_norte %>%
  mutate(tipo = "Continental")

# Placa del Caribe
sismos_placa_caribe = sismos_placa("Placa_Caribe.kml", camer, "Caribe")

sismos_placa_caribe = sismos_placa_caribe %>%
  mutate(tipo = "Continental")

# Placa del Pacifico
sismos_placa_pacifico = sismos_placa("Placa_Pacifico.kml", camer, "Pacifico")

sismos_placa_pacifico = sismos_placa_pacifico %>%
  mutate(tipo = "Oceanica")

# Placa Norandina
sismos_placa_norandina = sismos_placa("Placa_Norand.kml", camer, "Norandina")

sismos_placa_norandina = sismos_placa_norandina %>%
  mutate(tipo = "Continental")

# Placa Nazca (Centroamerica)
sismos_placa_nazca_n = sismos_placa("Placa_Nazca_N.kml", camer, "Nazca")

sismos_placa_nazca_n = sismos_placa_nazca_n %>%
  mutate(tipo = "Oceanica")

# Placa Nazca Chile
sismos_placa_nazca = sismos_placa("Placa_Nazca.kml", chile, "Nazca")

sismos_placa_nazca = sismos_placa_nazca %>%
  mutate(tipo = "Oceanica")

# Placa Sudamericana
sismos_placa_sudamericana = sismos_placa("Placa_Sudam.kml", chile, "Sudamericana")

sismos_placa_sudamericana = sismos_placa_sudamericana %>%
  mutate(tipo = "Continental")


chile_2 = rbind(sismos_placa_nazca, sismos_placa_sudamericana)
camer_2 = rbind(sismos_placa_rivera, sismos_placa_panama, sismos_placa_cocos,
              sismos_placa_norte, sismos_placa_caribe, sismos_placa_pacifico,
              sismos_placa_norandina, sismos_placa_nazca_n)

fosa_meso = st_read("Fosa_Meso.kml") %>% st_cast("MULTILINESTRING")
dorsal_pacifico = st_read("Dorsal_Pacifico.kml") %>% st_cast("MULTILINESTRING")
golfo_california = st_read("Golfo_California.kml") %>% st_cast("MULTILINESTRING")
norte_caribe = st_read("Norte_Caribe.kml") %>% st_cast("MULTILINESTRING")
frac_panama = st_read("Fractura_Panama.kml") %>% st_cast("MULTILINESTRING")
fosa_atacama = st_read("Fosa_Atacama.kml") %>% st_cast("MULTILINESTRING")

base_esp = rbind(chile_2, camer_2)
base_esp = base_esp %>% filter(mw_hom >= 5.1)

# Convertir la base de sismos a objeto espacial
sismos_espaciales = st_as_sf(base_esp, 
                              coords = c("longitude", "latitude"), 
                              crs = 4326)

# Analisis de frecuencia -----

## Chile ----

dia_chile = base %>%
  filter(zona == "Chile Central") %>%
  # Paso 1: Convertir de texto a fecha real
  mutate(time = as.Date(time)) %>% 
  # Paso 2: Ahora sí puedes aplicar floor_date
  mutate(periodo = floor_date(time, "day")) %>%
  count(periodo) %>%
  complete(
    periodo = seq(
      as.Date("2000-01-01"),
      as.Date("2025-12-31"),
      by = "day"
    ),
    fill = list(n = 0)
  )

sem_chile = base %>%
  filter(zona == "Chile Central") %>%
  # Paso 1: Convertir de texto a fecha real
  mutate(time = as.Date(time)) %>% 
  # Paso 2: Ahora sí puedes aplicar floor_date
  mutate(periodo = floor_date(time, "week")) %>%
  count(periodo) %>%
  complete(
    periodo = seq(
      floor_date(as.Date("2000-01-01"), "week"),
      floor_date(as.Date("2025-12-31"), "week"),
      by = "week"
    ),
    fill = list(n = 0)
  )

mes_chile = base %>%
  filter(zona == "Chile Central") %>%
  # Paso 1: Convertir de texto a fecha real
  mutate(time = as.Date(time)) %>% 
  # Paso 2: Ahora sí puedes aplicar floor_date
  mutate(periodo = floor_date(time, "month")) %>%
  count(periodo) %>%
  complete(
    periodo = seq(
      floor_date(as.Date("2000-01-01"), "month"),
      floor_date(as.Date("2025-12-31"), "month"),
      by = "month"
    ),
    fill = list(n = 0)
  )

tri_chile = base %>%
  filter(zona == "Chile Central") %>%
  # Paso 1: Convertir de texto a fecha real
  mutate(time = as.Date(time)) %>% 
  # Paso 2: Ahora sí puedes aplicar floor_date
  mutate(periodo = floor_date(time, "quarter")) %>%
  count(periodo) %>%
  complete(
    periodo = seq(
      floor_date(as.Date("2000-01-01"), "quarter"),
      floor_date(as.Date("2025-12-31"), "quarter"),
      by = "quarter"
    ),
    fill = list(n = 0)
  )

anio_chile = base %>%
  filter(zona == "Chile Central") %>% 
  # Paso 1: Convertir de texto a fecha real
  mutate(time = as.Date(time)) %>% 
  # Paso 2: Ahora sí puedes aplicar floor_date
  mutate(periodo = floor_date(time, "year")) %>%
  count(periodo) %>%
  complete(
    periodo = seq(
      floor_date(as.Date("2000-01-01"), "year"),
      floor_date(as.Date("2025-12-31"), "year"),
      by = "year"
    ),
    fill = list(n = 0)
  )

compchile = data.frame(
  
  Resolucion = c("Diario",
                 "Semanal",
                 "Mensual",
                 "Trimestral",
                 "Anio"),
  
  Intervalos = c(nrow(dia_chile),
                 nrow(sem_chile),
                 nrow(mes_chile),
                 nrow(tri_chile),
                 nrow(anio_chile)),
  
  Promedio = c(mean(dia_chile$n),
               mean(sem_chile$n),
               mean(mes_chile$n),
               mean(tri_chile$n),
               mean(anio_chile$n)),
  
  Mediana = c(median(dia_chile$n),
              median(sem_chile$n),
              median(mes_chile$n),
              median(tri_chile$n),
              median(anio_chile$n)),
  
  SD = c(sd(dia_chile$n),
         sd(sem_chile$n),
         sd(mes_chile$n),
         sd(tri_chile$n),
         sd(anio_chile$n)),
  
  Ceros = c(mean(dia_chile$n == 0),
            mean(sem_chile$n == 0),
            mean(mes_chile$n == 0),
            mean(tri_chile$n == 0),
            mean(anio_chile$n==0)),
  
  Maximo = c(max(dia_chile$n),
             max(sem_chile$n),
             max(mes_chile$n),
             max(tri_chile$n),
             max(anio_chile$n))
)

compchile

## CAmer ----

dia_camer = base %>%
  filter(zona == "México-Centroamérica") %>%
  # Paso 1: Convertir de texto a fecha real
  mutate(time = as.Date(time)) %>% 
  # Paso 2: Ahora sí puedes aplicar floor_date
  mutate(periodo = floor_date(time, "day")) %>%
  count(periodo) %>%
  complete(
    periodo = seq(
      as.Date("2000-01-01"),
      as.Date("2025-12-31"),
      by = "day"
    ),
    fill = list(n = 0)
  )

sem_camer = base %>%
  filter(zona == "México-Centroamérica") %>%
  # Paso 1: Convertir de texto a fecha real
  mutate(time = as.Date(time)) %>% 
  # Paso 2: Ahora sí puedes aplicar floor_date
  mutate(periodo = floor_date(time, "week")) %>%
  count(periodo) %>%
  complete(
    periodo = seq(
      floor_date(as.Date("2000-01-01"), "week"),
      floor_date(as.Date("2025-12-31"), "week"),
      by = "week"
    ),
    fill = list(n = 0)
  )

mes_camer = base %>%
  filter(zona == "México-Centroamérica") %>%
  # Paso 1: Convertir de texto a fecha real
  mutate(time = as.Date(time)) %>% 
  # Paso 2: Ahora sí puedes aplicar floor_date
  mutate(periodo = floor_date(time, "month")) %>%
  count(periodo) %>%
  complete(
    periodo = seq(
      floor_date(as.Date("2000-01-01"), "month"),
      floor_date(as.Date("2025-12-31"), "month"),
      by = "month"
    ),
    fill = list(n = 0)
  )

tri_camer = base %>%
  filter(zona == "México-Centroamérica") %>%
  # Paso 1: Convertir de texto a fecha real
  mutate(time = as.Date(time)) %>% 
  # Paso 2: Ahora sí puedes aplicar floor_date
  mutate(periodo = floor_date(time, "quarter")) %>%
  count(periodo) %>%
  complete(
    periodo = seq(
      floor_date(as.Date("2000-01-01"), "quarter"),
      floor_date(as.Date("2025-12-31"), "quarter"),
      by = "quarter"
    ),
    fill = list(n = 0)
  )

anio_camer = base %>%
  filter(zona == "México-Centroamérica") %>% 
  # Paso 1: Convertir de texto a fecha real
  mutate(time = as.Date(time)) %>% 
  # Paso 2: Ahora sí puedes aplicar floor_date
  mutate(periodo = floor_date(time, "year")) %>%
  count(periodo) %>%
  complete(
    periodo = seq(
      floor_date(as.Date("2000-01-01"), "year"),
      floor_date(as.Date("2025-12-31"), "year"),
      by = "year"
    ),
    fill = list(n = 0)
  )

compcamer = data.frame(
  
  Resolucion = c("Diario",
                 "Semanal",
                 "Mensual",
                 "Trimestral",
                 "Anio"),
  
  Intervalos = c(nrow(dia_camer),
                 nrow(sem_camer),
                 nrow(mes_camer),
                 nrow(tri_camer),
                 nrow(anio_camer)),
  
  Promedio = c(mean(dia_camer$n),
               mean(sem_camer$n),
               mean(mes_camer$n),
               mean(tri_camer$n),
               mean(anio_camer$n)),
  
  Mediana = c(median(dia_camer$n),
              median(sem_camer$n),
              median(mes_camer$n),
              median(tri_camer$n),
              median(anio_camer$n)),
  
  SD = c(sd(dia_camer$n),
         sd(sem_camer$n),
         sd(mes_camer$n),
         sd(tri_camer$n),
         sd(anio_camer$n)),
  
  Ceros = c(mean(dia_camer$n == 0),
            mean(sem_camer$n == 0),
            mean(mes_camer$n == 0),
            mean(tri_camer$n == 0),
            mean(anio_camer$n==0)),
  
  Maximo = c(max(dia_camer$n),
             max(sem_camer$n),
             max(mes_camer$n),
             max(tri_camer$n),
             max(anio_camer$n))
)

compcamer

# Analisis por series de tiempo ----

## Chile ----

ts_tri_chile = ts(
  tri_chile$n,
  start = c(2000, 1),
  frequency = 4
)

acf(
  ts_tri_chile,
  lag.max = 12,
  main = "Función de Autocorrelación"
)

pacf(
  ts_tri_chile,
  lag.max = 12,
  main = "Función de Autocorrelación"
)

Box.test(
  ts_tri_chile,
  lag = 4,
  type = "Ljung-Box"
)

Box.test(
  ts_tri_chile,
  lag = 8,
  type = "Ljung-Box"
)

Box.test(
  ts_tri_chile,
  lag = 12,
  type = "Ljung-Box"
)

adf.test(ts_tri_chile)
forecast::auto.arima(ts_tri_chile, method="ML", max.order= 8)

## CAmer ----

ts_tri_camer = ts(
  tri_camer$n,
  start = c(2000, 1),
  frequency = 4
)

acf(
  ts_tri_camer,
  lag.max = 12,
  main = "Función de Autocorrelación"
)

pacf(
  ts_tri_camer,
  lag.max = 12,
  main = "Función de Autocorrelación"
)

Box.test(
  ts_tri_camer,
  lag = 4,
  type = "Ljung-Box"
)

Box.test(
  ts_tri_camer,
  lag = 8,
  type = "Ljung-Box"
)

Box.test(
  ts_tri_camer,
  lag = 12,
  type = "Ljung-Box"
)

adf.test(ts_tri_camer)
forecast::auto.arima(ts_tri_camer, method="ML", max.order= 8)

# Modelos lineales generalizados ----

## Chile ----

chile$time = as.POSIXct(
  chile$time,
  format = "%Y-%m-%dT%H:%M:%OSZ",
  tz = "UTC"
)

chile = chile %>%
  arrange(time)

dist_media_consecutiva = function(geom, depth){
  if(length(geom) < 2)
    return(0)
  
  coords = st_coordinates(geom)
  d_h = distHaversine(
    coords[-nrow(coords), ],
    coords[-1, ]
  ) / 1000
  
  d_v = diff(depth)
  mean(sqrt(d_h^2 + d_v^2))
}

dist_media_total = function(geom, depth){
  if(length(geom) < 2)
    return(0)
  
  coords = st_coordinates(geom)
  d_h = geosphere::distm(
    coords,
    fun = distHaversine
  ) / 1000
  
  d_v = abs(outer(depth, depth, "-"))
  d = sqrt(d_h^2 + d_v^2)
  mean(d[upper.tri(d)])
}

# Si sabes que tus datos son espaciales, reactívalos así:

base_glm = base_esp %>%
  filter(zona == "Chile Central") %>%
  # Paso 1: Convertir de texto a fecha real
  mutate(time = as.Date(time)) %>% 
  # Paso 2: Ahora sí puedes aplicar floor_date
  mutate(
    periodo = floor_date(time,"month"),
    dia = day(time),
    dias_mes = days_in_month(time),
    momento = dia/dias_mes
  ) %>%
  group_by(periodo) %>%
  arrange(time,.by_group = TRUE) %>%
  summarise(
    n = n(),
    mag_max = max(mw_hom),
    depth_media = mean(depth),
    dist_consecutiva =
      dist_media_consecutiva(
        geometry,
        depth
      ),
    dist_total =
      dist_media_total(
        geometry,
        depth
      ),
    momento_max =
      momento[which.max(mw_hom)],
    .groups="drop"
  )

base_glm = base_glm %>%
  arrange(periodo) %>%
  mutate(
    mag_max_lag = lag(mag_max),
    depth_media_lag = lag(depth_media),
    dist_consecutiva_lag = lag(dist_consecutiva),
    dist_total_lag = lag(dist_total),
    momento_max_lag = lag(momento_max)
  )

### Supuestos ----

mean(na.omit(base_glm$n)) #4.168421
var(na.omit(base_glm$n)) #221.0191

# Presenta sobredispersión

# Modelo Propuesto

modelo_nb = glm.nb(
  n ~
    mag_max +
    depth_media +
    momento_max +
    dist_total,
  data = na.omit(base_glm)
)

summary(modelo_nb)
pscl::pR2(modelo_nb)

### STEPWISE ----

mod_max = glm.nb(n ~ mag_max +
                   mag_max_lag +
                   depth_media +
                   depth_media_lag +
                   dist_consecutiva +
                   dist_total +
                   momento_max +
                   momento_max_lag,
                 data = na.omit(base_glm))

mod_step=stepAIC(
  modelo_nb,
  scope = list(
    lower = formula(modelo_nb),
    upper = mod_max
  ),
  direction = "both",
  trace = TRUE)

summary(mod_step)
pscl::pR2(mod_step)

# Diagnostico propuesto

acf(residuals(modelo_nb, type = "pearson"))
acf(residuals(modelo_nb, type = "deviance"))
Box.test(residuals(modelo_nb, type="pearson"),
         type="Ljung")

plot(fitted(modelo_nb),
     residuals(modelo_nb, type="pearson"))
abline(h=0,lty=2)

vif(modelo_nb)

### Diagnostico Step ----

acf(residuals(mod_step, type = "pearson"))
acf(residuals(mod_step, type = "deviance"))
Box.test(residuals(mod_step, type="pearson"),
         type="Ljung")

plot(fitted(mod_step),
     residuals(mod_step, type="pearson"))
abline(h=0,lty=2)

vif(mod_step)

### Graficos ----

ggplot(na.omit(base_glm), aes(x = periodo)) +
  geom_line(aes(y = n, color = "Observada"), linewidth = 0.8) +
  geom_point(aes(y = n), color = "black", size = 1.6) +
  geom_line(aes(y = predict(modelo_nb, type = "response"),
                color = "Ajustada"),
            linewidth = 0.8) +
  scale_color_manual(values = c("Observada" = "black",
                                "Ajustada" = "red"),
                     name = NULL) +
  labs(x = "Periodo",
       y = "Frecuencia de sismos",
       title = "Frecuencia observada y ajustada por el modelo binomial negativo") +
  theme_bw() +
  theme(legend.position = c(0.90, 0.82),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = "black"),
        legend.key = element_blank(),
        legend.text = element_text(size = 13),
        legend.key.width = unit(1.8, "cm"))


ggplot(na.omit(base_glm), aes(x = periodo)) +
  geom_line(aes(y = n, color = "Observada"), linewidth = 0.8) +
  geom_point(aes(y = n), color = "black", size = 1.6) +
  geom_line(aes(y = predict(mod_step, type = "response"),
                color = "Ajustada"),
            linewidth = 0.8) +
  scale_color_manual(values = c("Observada" = "black",
                                "Ajustada" = "red"),
                     name = NULL) +
  labs(x = "Periodo",
       y = "Frecuencia de sismos",
       title = "Frecuencia observada y ajustada por el modelo binomial negativo") +
  theme_bw() +
  theme(legend.position = c(0.90, 0.82),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = "black"),
        legend.key = element_blank(),
        legend.text = element_text(size = 13),
        legend.key.width = unit(1.8, "cm"))

## CAmer ----

camer$time = as.POSIXct(
  camer$time,
  format = "%Y-%m-%dT%H:%M:%OSZ",
  tz = "UTC"
)

camer = camer %>%
  arrange(time)

base_glm_camer = base_esp %>%
  filter(zona == "México-Centroamérica") %>%
  # Paso 1: Convertir de texto a fecha real
  mutate(time = as.Date(time)) %>% 
  # Paso 2: Ahora sí puedes aplicar floor_date
  mutate(
    periodo = floor_date(time,"month"),
    dia = day(time),
    dias_mes = days_in_month(time),
    momento = dia/dias_mes
  ) %>%
  group_by(periodo) %>%
  arrange(time,.by_group = TRUE) %>%
  summarise(
    n = n(),
    mag_max = max(mw_hom),
    depth_media = mean(depth),
    dist_consecutiva =
      dist_media_consecutiva(
        geometry,
        depth
      ),
    dist_total =
      dist_media_total(
        geometry,
        depth
      ),
    momento_max =
      momento[which.max(mw_hom)],
    .groups="drop"
  )

base_glm_camer = base_glm_camer %>%
  arrange(periodo) %>%
  mutate(
    mag_max_lag = lag(mag_max),
    depth_media_lag = lag(depth_media),
    dist_consecutiva_lag = lag(dist_consecutiva),
    dist_total_lag = lag(dist_total),
    momento_max_lag = lag(momento_max)
  )

modelo_camer = glm.nb(
  n ~
    mag_max +
    depth_media +
    momento_max +
    dist_total,
  data = na.omit(base_glm_camer)
)

summary(modelo_camer)
pscl::pR2(modelo_camer)

### STEPWISE ----

mod_max_camer = glm.nb(n ~ mag_max +
                         mag_max_lag +
                         depth_media +
                         depth_media_lag +
                         dist_consecutiva +
                         dist_total +
                         momento_max +
                         momento_max_lag,
                       data = na.omit(base_glm_camer))

mod_step_camer = stepAIC(
  modelo_camer,
  scope = list(
    lower = formula(modelo_camer),
    upper = mod_max_camer
  ),
  direction = "both",
  trace = TRUE)

summary(mod_step_camer)
pscl::pR2(mod_step_camer)

### Diagnostico ----

acf(residuals(modelo_camer, type = "pearson"))
acf(residuals(modelo_camer, type = "deviance"))
Box.test(residuals(modelo_camer, type="pearson"),
         type="Ljung")

plot(fitted(modelo_camer),
     residuals(modelo_camer, type="pearson"))
abline(h=0,lty=2)

vif(modelo_nb)

### Diagnostico Step ----

acf(residuals(mod_step_camer, type = "pearson"))
acf(residuals(mod_step_camer, type = "deviance"))
Box.test(residuals(mod_step_camer, type="pearson"),
         type="Ljung")

plot(fitted(mod_step_camer),
     residuals(mod_step_camer, type="pearson"))
abline(h=0,lty=2)

vif(mod_step_camer)

### Graficos ----

ggplot(na.omit(base_glm_camer), aes(x = periodo)) +
  geom_line(aes(y = n, color = "Observada"), linewidth = 0.8) +
  geom_point(aes(y = n), color = "black", size = 1.6) +
  geom_line(aes(y = predict(modelo_camer, type = "response"),
                color = "Ajustada"),
            linewidth = 0.8) +
  scale_color_manual(values = c("Observada" = "black",
                                "Ajustada" = "red"),
                     name = NULL) +
  labs(x = "Periodo",
       y = "Frecuencia de sismos",
       title = "Frecuencia observada y ajustada por el modelo binomial negativo") +
  theme_bw() +
  theme(legend.position = c(0.90, 0.82),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = "black"),
        legend.key = element_blank(),
        legend.text = element_text(size = 13),
        legend.key.width = unit(1.8, "cm"))


ggplot(na.omit(base_glm_camer), aes(x = periodo)) +
  geom_line(aes(y = n, color = "Observada"), linewidth = 0.8) +
  geom_point(aes(y = n), color = "black", size = 1.6) +
  geom_line(aes(y = predict(mod_step_camer, type = "response"),
                color = "Ajustada"),
            linewidth = 0.8) +
  scale_color_manual(values = c("Observada" = "black",
                                "Ajustada" = "red"),
                     name = NULL) +
  labs(x = "Periodo",
       y = "Frecuencia de sismos",
       title = "Frecuencia observada y ajustada por el modelo binomial negativo") +
  theme_bw() +
  theme(legend.position = c(0.90, 0.82),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = "black"),
        legend.key = element_blank(),
        legend.text = element_text(size = 13),
        legend.key.width = unit(1.8, "cm"))

# Comparación no paramétrica de magnitud y profundidad ----

mag_chile = base$mw_hom[
  base$zona == "Chile Central"
]

mag_mexico = base$mw_hom[
  base$zona == "México-Centroamérica"
]

prof_chile = base$depth[
  base$zona == "Chile Central"
]

prof_mexico = base$depth[
  base$zona == "México-Centroamérica"
]

## Magnitud homogenizada -----
mw_magnitud = wilcox.test(
  x = mag_mexico,
  y = mag_chile,
  alternative = "two.sided",
  exact = FALSE,
  correct = TRUE,
  conf.int = TRUE,
  conf.level = 0.95
)

ks_magnitud = ks.test(
  x = mag_mexico,
  y = mag_chile,
  alternative = "two.sided",
  exact = FALSE
)

mw_magnitud
ks_magnitud


## Profundidad -----

mw_profundidad = wilcox.test(
  x = prof_mexico,
  y = prof_chile,
  alternative = "two.sided",
  exact = FALSE,
  correct = TRUE,
  conf.int = TRUE,
  conf.level = 0.95
)

ks_profundidad = ks.test(
  x = prof_mexico,
  y = prof_chile,
  alternative = "two.sided",
  exact = FALSE
)

mw_profundidad
ks_profundidad

# Modelos de regresión logística ----

base = base %>%
  mutate(
    fecha_hora_utc = ymd_hms(time),
    fecha = as.Date(fecha_hora_utc),
    anio = year(fecha_hora_utc),
    mes = month(fecha_hora_utc),
    magType = tolower(magType),
    categoria_profundidad = case_when(
      depth >= 0 & depth < 70 ~ "Superficial",
      depth >= 70 & depth <= 300 ~ "Intermedio",
      depth > 300 ~ "Profundo",
      TRUE ~ NA_character_
    ),
    categoria_magnitud = case_when(
      mw_hom >= 5.0 & mw_hom < 6.0 ~ "Moderado",
      mw_hom >= 6.0 & mw_hom < 7.0 ~ "Fuerte",
      mw_hom >= 7.0 & mw_hom < 8.0 ~ "Mayor",
      mw_hom >= 8 ~ "Gran Terremoto",
      TRUE ~ NA_character_
    )
  )

## Profundidad intermedia por zona ----

# Modelo de Regresión Logística: Probabilidad de Sismos Profundos

datos_logistica = base %>%
  filter(!is.na(categoria_profundidad)) %>%
  filter(categoria_profundidad != "Profundo") %>%
  mutate(
    es_intermedia = ifelse(categoria_profundidad == "Intermedio", 1, 0),
    zona_detallada = factor(zona_detallada, levels = c("Chile (Nazca-Sudamérica)", "México (Cocos-Norteamérica)",
                                                       "Centroamérica (Cocos-Caribe)"))
  )

modelo_logistico = glm(es_intermedia ~ zona_detallada, 
                        data = datos_logistica, 
                        family = binomial)

cat("\n--- Resumen del Modelo Logístico ---\n")
summary(modelo_logistico)

cat("\n--- Odds Ratios (Exponencial de los coeficientes) ---\n")
odds_ratios = exp(coef(modelo_logistico))
print(odds_ratios)

## Magnitud fuerte ----

datos_fuertes = base %>%
  filter(!is.na(categoria_profundidad) & categoria_profundidad != "Profundo") %>%
  mutate(
    es_fuerte = ifelse(mw_hom >= 6.0, 1, 0),
    estrato_fisico = case_when(
      categoria_profundidad == "Superficial" ~ "1_Superficial",
      categoria_profundidad == "Intermedio" ~ "2_Intermedio"
    ),
    estrato_fisico = factor(estrato_fisico, levels = c("1_Superficial", "2_Intermedio")),
    zona_detallada = factor(zona_detallada, levels = c("Chile (Nazca-Sudamérica)", "México (Cocos-Norteamérica)",
                                                       "Centroamérica (Cocos-Caribe)"))
  )

modelo_riesgo_fuerte = glm(es_fuerte ~ estrato_fisico + zona_detallada, 
                            data = datos_fuertes, 
                            family = binomial)

cat("\n--- Resumen del Modelo de Sismos Fuertes ---\n")
summary(modelo_riesgo_fuerte)

cat("\n--- Odds Ratios (Impacto en la probabilidad) ---\n")
odds_ratios_fuertes = exp(coef(modelo_riesgo_fuerte))
print(odds_ratios_fuertes)

## Curvas ROC ----

### Modelo 1 (Profundidad Intermedia) ----

probabilidades_m1 = predict(modelo_logistico, type = "response")
curva_roc_m1 = roc(datos_logistica$es_intermedia, probabilidades_m1)

par(mfrow = c(1, 2))

plot(curva_roc_m1, 
     main = "Curva ROC - Prob. de Sismo Intermedio",
     col = "blue", 
     lwd = 2, 
     legacy.axes = TRUE, 
     print.auc = TRUE,
     print.thres = "best",                 # <- Agrega el umbral óptimo
     print.thres.best.method = "youden",   # <- Utiliza el método de Youden
     print.thres.pch = 19,                 # <- Diseño del punto (círculo sólido)
     print.thres.col = "darkblue",         # <- Color del texto y punto
     print.thres.cex = 0.8)                # <- Tamaño del texto del umbral

### Modelo 2 (Sismos Fuertes) ----

probabilidades_m2 = predict(modelo_riesgo_fuerte, type = "response")
curva_roc_m2 = roc(datos_fuertes$es_fuerte, probabilidades_m2)

plot(curva_roc_m2, 
     main = "Curva ROC - Prob. de Sismo Fuerte (M >= 6.0)",
     col = "red", 
     lwd = 2, 
     legacy.axes = TRUE,
     print.auc = TRUE,
     print.thres = "best",                 # <- Agrega el umbral óptimo
     print.thres.best.method = "youden",   # <- Utiliza el método de Youden
     print.thres.pch = 19,
     print.thres.col = "darkred",
     print.thres.cex = 0.8)

par(mfrow = c(1, 1))

## Matriz de confusión ----

corte_optimo_m1 <- coords(curva_roc_m1, x = "best", best.method = "youden")
print(corte_optimo_m1)
indice_youden_m1 <- corte_optimo_m1$sensitivity + corte_optimo_m1$specificity - 1
print(indice_youden_m1)

umbral_optimo_m1 <- corte_optimo_m1$threshold
clasificacion_m1 <- ifelse(probabilidades_m1 >= umbral_optimo_m1, 1, 0)
matriz_confusion_m1 <- table(
  Real = datos_logistica$es_intermedia, 
  Predicho = clasificacion_m1
)
cat("\n--- Matriz de Confusión (Modelo 1 - Umbral Youden) ---\n")
print(matriz_confusion_m1)

corte_optimo_m2 <- coords(curva_roc_m2, x = "best", best.method = "youden")
print(corte_optimo_m2)
indice_youden_m2 <- corte_optimo_m2$sensitivity + corte_optimo_m2$specificity - 1
print(indice_youden_m2)

umbral_optimo_m2 <- corte_optimo_m2$threshold
clasificacion_m2 <- ifelse(probabilidades_m2 >= umbral_optimo_m2, 1, 0)
matriz_confusion_m2 <- table(
  Real = datos_fuertes$es_fuerte, 
  Predicho = clasificacion_m2
)
cat("\n--- Matriz de Confusión (Modelo 2 - Umbral Youden) ---\n")
print(matriz_confusion_m2)

# Analisis de correspondencia ----

#Zona x perfil conjunto de magnitud y profundidad
variables_requeridas = c(
  "id",
  "zona",
  "depth",
  "mw_hom"
)

variables_faltantes = setdiff(
  variables_requeridas,
  names(base_esp)
)

if (length(variables_faltantes) > 0) {
  stop(
    paste(
      "Faltan las siguientes variables:",
      paste(variables_faltantes, collapse = ", ")
    )
  )
}

#El identificador debe existir y ser único
if (any(is.na(base$id))) {
  stop("Existen eventos sin identificador id.")
}

#2. Control del proceso de filtrado
control_base_ac = base %>%
  group_by(zona) %>%
  summarise(
    eventos_depurados = n(),
    con_mw_hom = sum(
      !is.na(mw_hom)
    ),
    sobre_mc_5_1 = sum(
      !is.na(mw_hom) &
        mw_hom >= 5.1
    ),
    con_profundidad_valida = sum(
      !is.na(mw_hom) &
        mw_hom >= 5.1 &
        !is.na(depth) &
        depth >= 0
    ),
    .groups = "drop"
  )

print(control_base_ac)

#3. Construcción de la base analítica
base_ac = base %>%
  # Protección adicional contra duplicados
  distinct(id, .keep_all = TRUE) %>%
  # Magnitud homogenizada y catalogo completo
  filter(
    !is.na(mw_hom),
    mw_hom >= 5.1,
    # Profundidad valida
    !is.na(depth),
    depth >= 0
  ) %>%
  mutate(
    # Orden explícito de las zonas
    zona = factor(
      zona,
      levels = c(
        "Chile Central",
        "México-Centroamérica"
      )
    ),
    # Clasificación operativa de profundidad
    categoria_profundidad = case_when(
      depth >= 0 & depth <= 70 ~
        "Superficial [0-70 km]",
      depth > 70 & depth <= 300 ~
        "Intermedia (70-300 km]",
      depth > 300 ~
        "Profunda (>300 km)",
      TRUE ~ NA_character_
    ),
    categoria_profundidad = factor(
      categoria_profundidad,
      levels = c(
        "Superficial [0-70 km]",
        "Intermedia (70-300 km]",
        "Profunda (>300 km)"
      ),
      ordered = TRUE
    ),
    # Categorías de magnitud para el AC
    categoria_magnitud = case_when(
      mw_hom >= 5.1 & mw_hom < 6.0 ~
        "Moderado [5,1-5,9]",
      
      mw_hom >= 6.0 & mw_hom < 7.0 ~
        "Fuerte [6,0-6,9]",
      
      mw_hom >= 7.0 ~
        "Mayor [>=7,0]",
      TRUE ~ NA_character_
    ),
    categoria_magnitud = factor(
      categoria_magnitud,
      levels = c(
        "Moderado [5,1-5,9]",
        "Fuerte [6,0-6,9]",
        "Mayor [>=7,0]"
      ),
      ordered = TRUE
    )
  ) %>%
  droplevels()

#4. Verificaciones de integridad
#No deben existir identificadores duplicados
if (anyDuplicated(base_ac$id) > 0) {
  stop("La base analítica todavía contiene eventos duplicados.")
}

#No deben existir categorías faltantes
if (
  any(is.na(base_ac$zona)) ||
  any(is.na(base_ac$categoria_profundidad)) ||
  any(is.na(base_ac$categoria_magnitud))
) {
  stop(
    "Existen observaciones que no pudieron ser categorizadas."
  )
}

#5. Resultados de control
#Tamaño analítico por zona
conteo_zona = base_ac %>%
  count(zona, name = "n_eventos")

print(conteo_zona)

#Profundidad por zona
conteo_profundidad = base_ac %>%
  count(
    zona,
    categoria_profundidad,
    name = "n_eventos",
    .drop = FALSE
  )

print(conteo_profundidad)

#Magnitud por zona
conteo_magnitud = base_ac %>%
  count(
    zona,
    categoria_magnitud,
    name = "n_eventos",
    .drop = FALSE
  )

print(conteo_magnitud)

#Combinación conjunta de magnitud y profundidad
conteo_conjunto = base_ac %>%
  count(
    zona,
    categoria_profundidad,
    categoria_magnitud,
    name = "n_eventos",
    .drop = FALSE
  )

conteo_conjunto

## Tabla contingencia ----

#Perfil conjunto de profundidad y magnitud

##perfil sísmico conjunto----
base_ac = base_ac %>%
  mutate(
    perfil_sismico = case_when(
      
      categoria_profundidad == "Superficial [0-70 km]" &
        categoria_magnitud == "Moderado [5,1-5,9]" ~
        "Superficial / Moderado",
      
      categoria_profundidad == "Superficial [0-70 km]" &
        categoria_magnitud == "Fuerte [6,0-6,9]" ~
        "Superficial / Fuerte",
      
      categoria_profundidad == "Superficial [0-70 km]" &
        categoria_magnitud == "Mayor [>=7,0]" ~
        "Superficial / Mayor",
      
      categoria_profundidad == "Intermedia (70-300 km]" &
        categoria_magnitud == "Moderado [5,1-5,9]" ~
        "Intermedia / Moderado",
      
      categoria_profundidad == "Intermedia (70-300 km]" &
        categoria_magnitud == "Fuerte [6,0-6,9]" ~
        "Intermedia / Fuerte",
      
      categoria_profundidad == "Intermedia (70-300 km]" &
        categoria_magnitud == "Mayor [>=7,0]" ~
        "Intermedia / Mayor",
      
      TRUE ~ NA_character_
    ),
    
    perfil_sismico = factor(
      perfil_sismico,
      levels = c(
        "Superficial / Moderado",
        "Superficial / Fuerte",
        "Superficial / Mayor",
        "Intermedia / Moderado",
        "Intermedia / Fuerte",
        "Intermedia / Mayor"
      )
    )
  )

## Tabla F. Absolutas ----

#filas = zonas
#columnas = perfiles sísmicos
tabla_ac = xtabs(
  ~ zona + perfil_sismico,
  data = base_ac,
  drop.unused.levels = FALSE
)

tabla_ac = tabla_ac[
  rowSums(tabla_ac) > 0,
  colSums(tabla_ac) > 0,
  drop = FALSE
]
tabla_ac

conteo_zona
#2190 evnetos en total

# Validar
stopifnot(
  sum(tabla_ac) == nrow(base_ac)
)

stopifnot(
  all(rowSums(tabla_ac) > 0),
  all(colSums(tabla_ac) > 0)
)

addmargins(tabla_ac)

## Perfiles Fila ----

perfiles_fila = prop.table(
  tabla_ac,
  margin = 1
)

perfiles_fila_pct = round(
  100 * perfiles_fila,
  2
)

perfiles_fila_pct 

## Revisión de Frec. Esperadas ----

chi_ac = chisq.test(
  tabla_ac,
  correct = FALSE
)
frecuencias_esperadas = chi_ac$expected

round(
  frecuencias_esperadas, 2
) 

## Asoc. Global ----
### Chisq.test ----

#H_0: La zona y el perfil sísmico son independientes
#H_1: La distribución de perfiles sísmicos difiere entre zonas

chi_ac = chisq.test(
  tabla_ac,
  correct = FALSE
)
chi_ac 

### Inercia total ----

# Tamaño total de la tabla

n_total = sum(tabla_ac)
n_filas = nrow(tabla_ac)
n_columnas = ncol(tabla_ac)
chi2_ac = as.numeric(chi_ac$statistic)

inercia_total = chi2_ac / n_total

# v de Cramer
v_cramer = sqrt(
  chi2_ac /
    (
      n_total *
        min(n_filas - 1, n_columnas - 1)
    )
)

# Tabla resumen

resumen_dependencia = data.frame(
  n = n_total,
  chi_cuadrado = chi2_ac,
  grados_libertad = as.numeric(chi_ac$parameter),
  p_valor = chi_ac$p.value,
  inercia_total = inercia_total,
  v_cramer = v_cramer
)
resumen_dependencia

## Residuos ajustados de Haberman----

observados = chi_ac$observed #Frecuencias observadas
esperados = chi_ac$expected #Frecuencias esperadas bajo independencia

proporcion_filas = rowSums(observados) / sum(observados)
proporcion_columnas = colSums(observados) / sum(observados)


## Matriz del factor de corrección ----

#(1 - proporción fila) * (1 - proporción columna)
factor_correccion = outer(
  1 - proporcion_filas,
  1 - proporcion_columnas,
  FUN = "*"
)

## Residuos ajustados de Haberman ----

residuos_haberman = chi_ac$stdres
round(residuos_haberman, 3) 

## Construcción de Tabla ----

tabla_residuos = as.data.frame.table(
  observados,
  responseName = "observado"
)

names(tabla_residuos)[1:2] = c(
  "zona",
  "perfil_sismico"
)

tabla_residuos = tabla_residuos %>%
  mutate(
    esperado = as.vector(esperados),
    
    residuo_haberman =
      as.vector(residuos_haberman),
    
    diferencia = observado - esperado,
    
    clasificacion = case_when(
      residuo_haberman > 1.96 ~
        "Sobrerrepresentada",
      
      residuo_haberman < -1.96 ~
        "Subrepresentada",
      
      TRUE ~
        "Sin diferencia significativa"
    )
  ) %>%
  arrange(
    zona,
    perfil_sismico
  )

tabla_residuos %>%
  mutate(
    esperado = round(esperado, 2),
    diferencia = round(diferencia, 2),
    residuo_haberman = round(
      residuo_haberman, 3
    )
  )

celdas_significativas = tabla_residuos %>%
  filter(
    abs(residuo_haberman) > 1.96
  ) %>%
  arrange(
    desc(abs(residuo_haberman))
  )

celdas_significativas %>%
  mutate(
    esperado = round(esperado, 2),
    diferencia = round(diferencia, 2),
    residuo_haberman = round(
      residuo_haberman, 3
    )
  )

## Bonferroni para comprobar ----

numero_celdas = length(
  residuos_haberman
)

# Umbral Bonferroni bilateral
umbral_bonferroni = qnorm(
  1 - 0.05 /
    (2 * numero_celdas)
)

umbral_bonferroni

# Clasificar

tabla_residuos = tabla_residuos %>%
  mutate(
    significativo_95 =
      abs(residuo_haberman) > 1.96,
    
    significativo_bonferroni =
      abs(residuo_haberman) >
      umbral_bonferroni
  )

tabla_residuos

## Mapa de calor de residuos de Haberman ----

ggplot(
  tabla_residuos,
  aes(
    x = perfil_sismico,
    y = zona,
    fill = residuo_haberman
  )
) +
  geom_tile(
    linewidth = 0.5
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        residuo_haberman
      )
    ),
    size = 3.5
  ) +
  scale_fill_gradient2(
    midpoint = 0
  ) +
  labs(
    title = "Residuos ajustados de Haberman",
    subtitle = paste0(
      "Valores absolutos superiores a 1,96 ",
      "indican desviaciones significativas"
    ),
    x = "Perfil sísmico",
    y = NULL,
    fill = "Residuo\najustado"
  ) +
  theme_minimal(
    base_size = 11
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    panel.grid = element_blank(),
    plot.title.position = "plot"
  )

## Aporte al test Chi_cuadrado ----

aporte_chi_absoluto = (
  observados - esperados
)^2 / esperados

## Contribución porcentual respecto del chi-cuadrado total ----

aporte_chi_porcentaje = 100 *
  aporte_chi_absoluto /
  chi2_ac

round(
  aporte_chi_porcentaje, 2
)

# Incorporar todo a la tabla

tabla_contribuciones = as.data.frame.table(
  aporte_chi_porcentaje,
  responseName = "contribucion_chi_pct"
)

names(tabla_contribuciones)[1:2] = c(
  "zona",
  "perfil_sismico"
)

tabla_resultados_celda = tabla_residuos %>%
  left_join(
    tabla_contribuciones,
    by = c(
      "zona",
      "perfil_sismico"
    )
  ) %>%
  mutate(
    esperado = round(esperado, 2),
    residuo_haberman = round(
      residuo_haberman,
      2
    ),
    diferencia = round(
      diferencia,
      2
    ),
    contribucion_chi_pct = round(
      contribucion_chi_pct,
      2
    )
  )
tabla_resultados_celda

## Contribución total por perfil sísmico ----

aporte_por_perfil = colSums(
  aporte_chi_porcentaje
)

tabla_aporte_perfil = data.frame(
  perfil_sismico = names(aporte_por_perfil),
  contribucion_chi_pct = as.numeric(
    aporte_por_perfil
  )
) %>%
  arrange(
    desc(contribucion_chi_pct)
  ) %>%
  mutate(
    contribucion_chi_pct = round(
      contribucion_chi_pct,
      2
    )
  )
tabla_aporte_perfil

## Grafico ----

ggplot(
  tabla_aporte_perfil,
  aes(
    x = reorder(
      perfil_sismico,
      contribucion_chi_pct
    ),
    y = contribucion_chi_pct
  )
) +
  geom_col(
    width = 0.7
  ) +
  geom_text(
    aes(
      label = paste0(
        contribucion_chi_pct,
        "%"
      )
    ),
    hjust = -0.1,
    size = 3.5
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(
      0,
      max(
        tabla_aporte_perfil$contribucion_chi_pct
      ) * 1.12
    ),
    expand = expansion(
      mult = c(0, 0.02)
    )
  ) +
  labs(
    title = "Contribución de los perfiles a la asociación",
    subtitle = "Porcentaje del estadístico chi-cuadrado explicado por cada perfil",
    x = NULL,
    y = "Contribución al chi-cuadrado (%)"
  ) +
  theme_minimal(
    base_size = 11
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title.position = "plot"
  )

## AC -----

res_ac = CA(
  tabla_ac,
  ncp = 1,
  graph = FALSE
)
res_ac
res_ac$eig 

## Comprobar la equivalencia con la inercia total ----

autovalor_dim1 = as.numeric(
  res_ac$eig[1, 1]
)

comparacion_inercia = data.frame(
  inercia_chi_cuadrado = inercia_total,
  autovalor_ac = autovalor_dim1,
  diferencia = abs(
    inercia_total - autovalor_dim1
  )
)
comparacion_inercia 

stopifnot(
  isTRUE(
    all.equal(
      autovalor_dim1,
      inercia_total,
      tolerance = 1e-10
    )
  )
)

## Extraer Resultados Factoriales ----
### Coordenadas, contribuciones y cos2 de las filas----

coord_filas = as.numeric(res_ac$row$coord)
contrib_filas = as.numeric(res_ac$row$contrib)
cos2_filas = as.numeric(res_ac$row$cos2)

tabla_filas_ac = data.frame(
  zona = rownames(tabla_ac),
  
  masa_pct = 100 *
    rowSums(tabla_ac) /
    sum(tabla_ac),
  
  coordenada_dim1 = coord_filas,
  
  contribucion_dim1 = contrib_filas,
  
  cos2_dim1 = cos2_filas,
  
  row.names = NULL
) %>%
  mutate(
    masa_pct = round(masa_pct, 2),
    coordenada_dim1 = round(coordenada_dim1, 4),
    contribucion_dim1 = round(contribucion_dim1, 2),
    cos2_dim1 = round(cos2_dim1, 4)
  )
tabla_filas_ac

### Resultados factoriales de los perfiles sísmicos -----

coord_columnas = as.numeric(res_ac$col$coord)
contrib_columnas = as.numeric(res_ac$col$contrib)
cos2_columnas = as.numeric(res_ac$col$cos2)

tabla_columnas_ac = data.frame(
  perfil_sismico = colnames(tabla_ac),
  
  masa_pct = 100 *
    colSums(tabla_ac) /
    sum(tabla_ac),
  
  coordenada_dim1 = coord_columnas,
  
  contribucion_dim1 = contrib_columnas,
  
  cos2_dim1 = cos2_columnas,
  
  row.names = NULL
) %>%
  mutate(
    masa_pct = round(masa_pct, 2),
    coordenada_dim1 = round(coordenada_dim1, 4),
    contribucion_dim1 = round(contribucion_dim1, 2),
    cos2_dim1 = round(cos2_dim1, 4)
  ) %>%
  arrange(desc(contribucion_dim1))

tabla_columnas_ac

## Graficos ----
### Grafico perfiles sísmicos ----

grafico_perfiles_ac = tabla_columnas_ac %>%
  ggplot(
    aes(
      x = coordenada_dim1,
      y = reorder(
        perfil_sismico,
        coordenada_dim1
      )
    )
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = coordenada_dim1,
      yend = reorder(
        perfil_sismico,
        coordenada_dim1
      )
    ),
    linewidth = 0.6
  ) +
  geom_point(
    aes(size = contribucion_dim1)
  ) +
  geom_text(
    aes(
      label = paste0(
        round(coordenada_dim1, 3),
        " | ",
        round(contribucion_dim1, 1),
        "%"
      )
    ),
    hjust = ifelse(
      tabla_columnas_ac$coordenada_dim1 >= 0,
      -0.1,
      1.1
    ),
    size = 3.2
  ) +
  scale_size_continuous(
    range = c(2.5, 7)
  ) +
  scale_x_continuous(
    expand = expansion(
      mult = c(0.25, 0.25)
    )
  ) +
  labs(
    title = "Perfiles sísmicos en la dimensión factorial",
    subtitle = paste0(
      "Dimensión 1: 100% de la inercia disponible; ",
      "inercia total = ",
      round(inercia_total, 4)
    ),
    x = "Coordenada en la dimensión 1",
    y = NULL,
    size = "Contribución (%)"
  ) +
  theme_minimal(
    base_size = 11
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title.position = "plot"
  )

grafico_perfiles_ac

### Grafico de las zonas ----

grafico_zonas_ac = tabla_filas_ac %>%
  ggplot(
    aes(
      x = coordenada_dim1,
      y = zona
    )
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = coordenada_dim1,
      yend = zona
    ),
    linewidth = 0.7
  ) +
  geom_point(
    aes(size = contribucion_dim1)
  ) +
  geom_text(
    aes(
      label = paste0(
        round(coordenada_dim1, 3),
        " | ",
        round(contribucion_dim1, 1),
        "%"
      )
    ),
    hjust = ifelse(
      tabla_filas_ac$coordenada_dim1 >= 0,
      -0.15,
      1.15
    ),
    size = 3.3
  ) +
  scale_size_continuous(
    range = c(4, 8)
  ) +
  scale_x_continuous(
    expand = expansion(
      mult = c(0.3, 0.3)
    )
  ) +
  labs(
    title = "Zonas en la dimensión factorial",
    x = "Coordenada en la dimensión 1",
    y = NULL,
    size = "Contribución (%)"
  ) +
  theme_minimal(
    base_size = 11
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title.position = "plot"
  )

grafico_zonas_ac

# Analisis de eventos extremos ----

resumen_extremos = base %>%
  group_by(zona) %>% 
  summarise(
    total_eventos = n(),
    n_M60 = sum(mw_hom >= 6.0, na.rm = TRUE),
    n_M65 = sum(mw_hom >= 6.5, na.rm = TRUE),
    n_M70 = sum(mw_hom >= 7.0, na.rm = TRUE)
  ) %>%
  mutate(
    prop_M60 = (n_M60 / total_eventos) * 100,
    prop_M65 = (n_M65 / total_eventos) * 100,
    prop_M70 = (n_M70 / total_eventos) * 100
  )

print("--- Resumen de Eventos Extremos ---")
print(resumen_extremos)

tiempos_fuertes = base %>%
  filter(mw_hom >= 6.0) %>% 
  mutate(time = ymd_hms(time)) %>%
  arrange(zona, time) %>%
  group_by(zona) %>%
  mutate(
    tiempo_espera_dias = as.numeric(difftime(time, lag(time), units = "days"))
  ) %>%
  filter(!is.na(tiempo_espera_dias)) %>%
  ungroup()

## Calculo de indicadores finales (Proporción, Tasa y Cv) ----

sismos_completos = base %>% mutate(time = ymd_hms(time))
anios_catalogo = as.numeric(difftime(max(sismos_completos$time, na.rm = TRUE), 
                                     min(sismos_completos$time, na.rm = TRUE), 
                                     units = "days")) / 365.25

totales_por_zona = sismos_completos %>%
  group_by(zona) %>%
  summarise(
    total_eventos = n(),
    # Calculo de Joules
    energia_total_zona = sum(10^(1.5 * mw_hom + 4.8), na.rm = TRUE)
  )

resultados_extremos = tiempos_fuertes %>%
  group_by(zona) %>%
  summarise(
    n_fuertes = n(),
    # Calculo de Joules
    energia_fuertes = sum(10^(1.5 * mw_hom + 4.8), na.rm = TRUE),
    
    # Tasa y Tiempo de retorno
    tasa_anual = n_fuertes / anios_catalogo,
    tiempo_retorno_anios = 1 / tasa_anual,
    
    # Coeficiente de Variación de tiempos de espera
    media_espera = mean(tiempo_espera_dias, na.rm = TRUE),
    desviacion_espera = sd(tiempo_espera_dias, na.rm = TRUE),
    Cv = desviacion_espera / media_espera
  ) %>%
  # 4. Cruce para calcular porcentajes
  left_join(st_drop_geometry(totales_por_zona), by = "zona") %>%
  mutate(
    proporcion_eventos_pct = (n_fuertes / total_eventos) * 100,
    proporcion_energia_pct = (energia_fuertes / energia_total_zona) * 100
  ) %>%
  # 5. Seleccionar columnas clave para el informe
  dplyr::select(zona, n_fuertes, proporcion_eventos_pct, proporcion_energia_pct, 
         tasa_anual, tiempo_retorno_anios, Cv)

# Ver los resultados
print(resultados_extremos)

## FMD ----

fmd_datos = base %>%
  filter(mw_hom >= 6.0) %>%
  mutate(mag_bin = round(mw_hom, 1)) %>%
  group_by(zona, mag_bin) %>%
  summarise(conteo = n(), .groups = "drop") %>%
  arrange(zona, desc(mag_bin)) %>%
  group_by(zona) %>%
  mutate(
    cum_n = cumsum(conteo),
    log10_cum_n = log10(cum_n)
  ) %>%
  ungroup()

grafico_gr_extremo = ggplot(fmd_datos, aes(x = mag_bin, y = log10_cum_n, color = zona)) +
  
  # A. Los puntos reales del catalogo
  geom_point(size = 3, alpha = 0.8) +
  
  # B. La línea teórica (Ajuste lineal SOLO para magnitudes de 5.1 a 7.0)
  # fullrange = TRUE permite que la línea se extienda hacia la derecha para ver la falla
  #geom_smooth(data = fmd_datos %>% filter(mag_bin >= 5.1 & mag_bin <= 7.0),
  #            method = "lm", se = FALSE, fullrange = TRUE, 
  #            linetype = "dashed", linewidth = 1.2) +
  
  # C. Estética y colores
  scale_color_manual(values = c("Chile Central" = "darkred", "México-Centroamérica" = "darkgreen")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribución Frecuencia-Magnitud (M_w >= 6.0)",
    x = "Magnitud (Mw)",
    y = expression(Log[10] * "(Número Acumulado de Sismos)"),
    color = "Región"
  ) +
  
  # D. Forzar el eje X hasta M 9.0 para visualizar bien la "cola derecha"
  scale_x_continuous(limits = c(6.0, 9.0), breaks = seq(6.0, 9.0, by = 0.5)) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Mostrar el grafico
print(grafico_gr_extremo)

## K-S Magnitudes eventos mayores a 6.0 ----

chile_2 = base %>%
  filter(zona == "Chile Central", mw_hom >= 6.0) %>%
  pull(mw_hom)

mexico_2 = base %>%
  filter(zona == "México-Centroamérica", mw_hom >= 6.0) %>%
  pull(mw_hom)

ks.test(chile_2, mexico_2)

## K-S Profundidades eventos mayores a 6.0 ----

chile_3 = base %>%
  filter(zona == "Chile Central", mw_hom >= 6.0) %>%
  pull(depth)

mexico_3 = base %>%
  filter(zona == "México-Centroamérica", mw_hom >= 6.0) %>%
  pull(depth)

ks.test(chile_3, mexico_3)

# Analisis de sensibilidad ----

base_sensibilidad = base %>%
  distinct(id, .keep_all = TRUE) %>%
  filter(
    !is.na(id),
    !is.na(zona),
    !is.na(mw_hom),
    !is.na(depth),
    depth >= 0
  ) %>%
  mutate(
    zona = factor(
      zona,
      levels = c(
        "Chile Central",
        "México-Centroamérica"
      )
    )
  ) %>%
  droplevels()


base_sensibilidad %>%
  count(zona)

# Filtrar la base

analizar_umbral = function(datos, umbral) {
  
  #Usamos .env$umbral para indicar explícitamente que umbral es un argumento de la función
  datos_umbral = datos %>%
    filter(
      .data$mw_hom >= .env$umbral
    ) %>%
    droplevels()
  
  # Vectores de magnitud
  mag_chile = datos_umbral$mw_hom[
    datos_umbral$zona == "Chile Central"
  ]
  
  mag_mexico = datos_umbral$mw_hom[
    datos_umbral$zona == "México-Centroamérica"
  ]
  
  # Vectores de profundidad
  prof_chile = datos_umbral$depth[
    datos_umbral$zona == "Chile Central"
  ]
  
  prof_mexico = datos_umbral$depth[
    datos_umbral$zona == "México-Centroamérica"
  ]
  
  # Comprobar tamaños
  if (
    length(mag_chile) < 2 ||
    length(mag_mexico) < 2
  ) {
    stop(
      paste(
        "El umbral",
        umbral,
        "no deja suficientes observaciones."
      )
    )
  }
  
  
  
  #Magnitud----
  mw_magnitud = suppressWarnings(
    wilcox.test(
      x = mag_mexico,
      y = mag_chile,
      alternative = "two.sided",
      exact = FALSE,
      correct = TRUE,
      conf.int = TRUE,
      conf.level = 0.95
    )
  )
  
  ks_magnitud = suppressWarnings(
    ks.test(
      x = mag_mexico,
      y = mag_chile,
      alternative = "two.sided",
      exact = FALSE
    )
  )
  
  
  
  
  #Profundidad----
  mw_profundidad = suppressWarnings(
    wilcox.test(
      x = prof_mexico,
      y = prof_chile,
      alternative = "two.sided",
      exact = FALSE,
      correct = TRUE,
      conf.int = TRUE,
      conf.level = 0.95
    )
  )
  
  ks_profundidad = suppressWarnings(
    ks.test(
      x = prof_mexico,
      y = prof_chile,
      alternative = "two.sided",
      exact = FALSE
    )
  )
  
  
  
  
  #Tabla de resultados----
  resultados_magnitud = data.frame(
    umbral = umbral,
    variable = "Magnitud",
    
    n_chile = length(mag_chile),
    n_mexico = length(mag_mexico),
    
    q1_chile = as.numeric(
      quantile(mag_chile, 0.25)
    ),
    mediana_chile = median(mag_chile),
    q3_chile = as.numeric(
      quantile(mag_chile, 0.75)
    ),
    
    q1_mexico = as.numeric(
      quantile(mag_mexico, 0.25)
    ),
    mediana_mexico = median(mag_mexico),
    q3_mexico = as.numeric(
      quantile(mag_mexico, 0.75)
    ),
    
    diferencia_medianas =
      median(mag_mexico) -
      median(mag_chile),
    
    W = as.numeric(
      mw_magnitud$statistic
    ),
    p_mann_whitney =
      mw_magnitud$p.value,
    
    estimador_HL =
      as.numeric(mw_magnitud$estimate),
    
    IC95_HL_inferior =
      as.numeric(mw_magnitud$conf.int[1]),
    
    IC95_HL_superior =
      as.numeric(mw_magnitud$conf.int[2]),
    
    D_KS = as.numeric(
      ks_magnitud$statistic
    ),
    p_KS = ks_magnitud$p.value
  )
  
  resultados_profundidad = data.frame(
    umbral = umbral,
    variable = "Profundidad",
    
    n_chile = length(prof_chile),
    n_mexico = length(prof_mexico),
    
    q1_chile = as.numeric(
      quantile(prof_chile, 0.25)
    ),
    mediana_chile = median(prof_chile),
    q3_chile = as.numeric(
      quantile(prof_chile, 0.75)
    ),
    
    q1_mexico = as.numeric(
      quantile(prof_mexico, 0.25)
    ),
    mediana_mexico = median(prof_mexico),
    q3_mexico = as.numeric(
      quantile(prof_mexico, 0.75)
    ),
    
    diferencia_medianas =
      median(prof_mexico) -
      median(prof_chile),
    
    W = as.numeric(
      mw_profundidad$statistic
    ),
    p_mann_whitney =
      mw_profundidad$p.value,
    
    estimador_HL =
      as.numeric(mw_profundidad$estimate),
    
    IC95_HL_inferior =
      as.numeric(mw_profundidad$conf.int[1]),
    
    IC95_HL_superior =
      as.numeric(mw_profundidad$conf.int[2]),
    
    D_KS = as.numeric(
      ks_profundidad$statistic
    ),
    p_KS = ks_profundidad$p.value
  )
  
  bind_rows(
    resultados_magnitud,
    resultados_profundidad
  )
}

## Ejecutar ----

umbrales_sensibilidad = c(5.1, 5.5, 6.0)

resultados_lista = vector(
  mode = "list",
  length = length(umbrales_sensibilidad)
)

for (
  i in seq_along(umbrales_sensibilidad)
) {
  
  umbral_actual =
    umbrales_sensibilidad[i]
  
  resultados_lista[[i]] =
    analizar_umbral(
      datos = base_sensibilidad,
      umbral = umbral_actual
    )
}

resultados_sensibilidad = bind_rows(
  resultados_lista
)

resultados_sensibilidad

## Tabla ----

tabla_sensibilidad = resultados_sensibilidad %>%
  dplyr::select(
    variable,
    umbral,
    n_chile,
    n_mexico,
    mediana_chile,
    mediana_mexico,
    diferencia_medianas,
    estimador_HL,
    IC95_HL_inferior,
    IC95_HL_superior,
    D_KS,
    p_KS,
    p_mann_whitney
  )

tabla_sensibilidad

## Graficos ----

### Medianas ----

datos_medianas = resultados_sensibilidad %>%
  dplyr::select(umbral, variable, mediana_chile, mediana_mexico) %>%
  pivot_longer(
    cols = c(mediana_chile, mediana_mexico),
    names_to = "zona",
    values_to = "mediana"
  ) %>%
  mutate(
    zona = case_match(
      zona,
      "mediana_chile"   ~ "Chile Central",
      "mediana_mexico" ~ "México-Centroamérica"
    )
  )

grafico_sensibilidad_medianas = ggplot(
  datos_medianas,
  aes(
    x = factor(umbral),
    y = mediana,
    group = zona,
    shape = zona,
    linetype = zona
  )
) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_point(
    size = 3
  ) +
  facet_wrap(
    ~ variable,
    scales = "free_y"
  ) +
  labs(
    title = "Sensibilidad al umbral mínimo de magnitud",
    subtitle = "Evolución de las medianas por zona",
    x = "Umbral mínimo de magnitud",
    y = "Mediana",
    shape = "Zona",
    linetype = "Zona"
  ) +
  theme_minimal(
    base_size = 11
  ) +
  theme(
    legend.position = "bottom",
    plot.title.position = "plot"
  )

grafico_sensibilidad_medianas

## Estadístico del KS ----

grafico_sensibilidad_ks = ggplot(
  resultados_sensibilidad,
  aes(
    x = factor(umbral),
    y = D_KS,
    group = variable
  )
) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_point(
    size = 3
  ) +
  facet_wrap(
    ~ variable,
    scales = "free_y"
  ) +
  labs(
    title = "Sensibilidad de la separación entre distribuciones",
    subtitle = "Estadístico D de Kolmogorov-Smirnov",
    x = "Umbral mínimo de magnitud",
    y = "Estadístico D"
  ) +
  theme_minimal(
    base_size = 11
  ) +
  theme(
    plot.title.position = "plot"
  )

grafico_sensibilidad_ks

# Análisis de sensibilidad espacial ----

escenarios_espaciales = c(0, 0.25, 0.50) 
resultados_espaciales = data.frame()

base_sensibilidad_esp = sismos_completos %>% filter(mw_hom >= 5.1)

for (margen in escenarios_espaciales) {
  
  # Reducción de las cajas geográficas
  sismos_reducidos = base_sensibilidad_esp %>%
    group_by(zona) %>%
    mutate(
      lat_min = min(latitude) + margen,
      lat_max = max(latitude) - margen,
      lon_min = min(longitude) + margen,
      lon_max = max(longitude) - margen
    ) %>%
    filter(
      latitude >= lat_min & latitude <= lat_max,
      longitude >= lon_min & longitude <= lon_max
    ) %>%
    ungroup()
  
  n_chile = sum(sismos_reducidos$zona == "Chile Central")
  n_mex = sum(sismos_reducidos$zona == "México-Centroamérica")
  
  # Pruebas de hipótesis no paramétricas
  test_prof = wilcox.test(depth ~ zona, data = sismos_reducidos, exact = FALSE)
  test_mag = wilcox.test(mw_hom ~ zona, data = sismos_reducidos, exact = FALSE)
  
  # Extracción de descriptivos (Min, Max, Mediana)
  stats = sismos_reducidos %>%
    group_by(zona) %>%
    summarise(
      min_prof = min(depth, na.rm = TRUE),
      max_prof = max(depth, na.rm = TRUE),
      med_prof = median(depth, na.rm = TRUE),
      min_mag = min(mw_hom, na.rm = TRUE),
      max_mag = max(mw_hom, na.rm = TRUE),
      med_mag = median(mw_hom, na.rm = TRUE),
      .groups = "drop"
    )
  
  stats_chile = stats %>% filter(zona == "Chile Central")
  stats_mex = stats %>% filter(zona == "México-Centroamérica")
  
  # Compilación de la fila de resultados
  resultados_espaciales = rbind(resultados_espaciales, data.frame(
    Reduccion_Grados = margen,
    N_Chile = n_chile,
    N_Mex = n_mex,
    
    # Bloque Magnitud
    Min_Mag_CL = stats_chile$min_mag,
    Max_Mag_CL = stats_chile$max_mag,
    Med_Mag_CL = stats_chile$med_mag,
    Min_Mag_MX = stats_mex$min_mag,
    Max_Mag_MX = stats_mex$max_mag,
    Med_Mag_MX = stats_mex$med_mag,
    P_Valor_Mag = test_mag$p.value,
    
    # Bloque Profundidad
    Min_Prof_CL = stats_chile$min_prof,
    Max_Prof_CL = stats_chile$max_prof,
    Med_Prof_CL = stats_chile$med_prof,
    Min_Prof_MX = stats_mex$min_prof,
    Max_Prof_MX = stats_mex$max_prof,
    Med_Prof_MX = stats_mex$med_prof,
    P_Valor_Prof = test_prof$p.value
  ))
}

print("--- RESULTADOS SENSIBILIDAD ESPACIAL ---")
options(scipen = 999) 

print(round(t(resultados_espaciales), 4))
