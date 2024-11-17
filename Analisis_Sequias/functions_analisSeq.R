################################################################################
################################################################################
#                               ACTIVIDAD 2                                    #
################################################################################
packages = c("data.table", "dplyr", "purrr", "terra", "ggplot2")

# Función para verificar e instalar los paquetes faltantes
check_install_packages = function(pkgs) {
  for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    } else {
      library(pkg, character.only = TRUE)
    }
  }
}

check_install_packages(packages)

################################################################################
prep_data = function(data, SPEI = NULL, Fecha_inicio, Fecha_fin) {
  if (!is.null(SPEI)) {
    names(data)[1] = "Fecha"
    data[data == Inf | data == -Inf] = NA
    data = data[Fecha >= Fecha_inicio & Fecha <= Fecha_fin, ]
    data = data[, .(Fecha, Pixel_1, Pixel_2, Pixel_3, Pixel_4, Pixel_5, Pixel_6, Pixel_7, Pixel_8, Pixel_9,
                    Pixel_10, Pixel_11, Pixel_12, Pixel_13, Pixel_14, Pixel_15, Pixel_16, Pixel_17,
                    Pixel_18,Pixel_19, Pixel_20, Pixel_21, Pixel_23, Pixel_24, Pixel_25, Pixel_26,
                    Pixel_27, Pixel_28, Pixel_29, Pixel_30, Pixel_38), ]
    
    names(data)[2:length(data)] = paste0("Pixel_", 1:length(data))
  } else {
    data = data[Fecha >= Fecha_inicio & Fecha <= Fecha_fin, ]
  }
  return(data)
}
################################################################################
categorizar_Sequias = function(data, type) {
  if (type == "SPEI") {
    cats_spei = function(x) {
      fcase(
        x <= -2, "Sequia Extrema",
        x > -2 & x <= -1.5, "Sequia Severa",
        x > -1.5 & x <= -1, "Sequia Moderada",
        x > -1 & x <= -0.5, "Sequia Leve",
        x > -0.5 & x <= 0.5, "Casi Normal",
        x > 0.5 & x <= 1, "Inundación Leve",
        x > 1 & x <= 1.5, "Inundación Moderada",
        x > 1.5 & x <= 2, "Inundación Severa",
        x > 2, "Extremadamente inundable"
      )
    }
    Clasification = list()
    for (i in setdiff(names(data), "Fecha")) {
      temp_dt = data[, .(i = get(i), Clasificacion = cats_spei(.SD[[i]])), by = Fecha]
      setnames(temp_dt, "i", i)
      Clasification[[i]] = temp_dt
    }
    
  } else {
    SSI_tibble = as_tibble(data)
    SSI_tibble = SSI_tibble %>%
      mutate(Categoria = case_when(
        Z >= 2.0 ~ "Extremadamente húmedo",
        Z >= 1.5 ~ "Severamente húmedo",
        Z >= 1.0 ~ "Moderadamente húmedo",
        Z >= 0.0 ~ "Ligeramente húmedo",
        Z >= -1.0 ~ "Sequía leve",
        Z >= -1.5 ~ "Sequía Moderada",
        Z >= -2.0 ~ "Sequía severa",
        TRUE ~ "Sequía extrema"
      ))
    Clasification = as.data.table(SSI_tibble)
  }
  return(Clasification)
}
################################################################################
##                 1. CATEGORIZACIÓN DE LOS EVENTOS DE SEQUÍA                 ##
categorizar_SPEI = function(data) {
  df = data %>%
    imap_dfr(~ mutate(.x, Pixel = .y))
  
  # Agrupamos por Píxel y Clasificación, luego contamos las frecuencias
  frecuencia_eventos = df %>%
    group_by(Pixel, Clasificacion) %>%
    summarize(Frecuencia = n(), .groups = 'drop') %>%
    arrange(Pixel, desc(Frecuencia))
  frecuencia_eventos = frecuencia_eventos[!is.na(frecuencia_eventos$Clasificacion), ]
  
  cords = data.frame(
    Pixel = paste0("Pixel_", 1:(length(SPEI_Tomebamba) - 1)),
    X = c(697318,702879,708439, 697310, 702870, 708430, 713991, 691742, 697302, 702862, 708421,
          713982, 719542, 691734,697293, 702853, 708412, 713972, 719532, 725093, 730653, 
          691725, 697285, 702844, 708403, 713963, 719523,725083,730643, 725073),
    Y = c(9698663,9698654,9698646, 9693133, 9693125, 9693116, 9693107, 9687612,9687604, 9687595, 9687587,
          9687578, 9687568, 9682083,9682075, 9682066, 9682057, 9682048, 9682038, 9682029, 9682019, 
          9676554, 9676545, 9676537, 9676527, 9676518, 9676508,9676498, 9676488, 9670968))
  
  
  ROI_SHP = terra::vect("C:/Users/Jonna/Desktop/Sequias_GIZ/Indices_Sequia/SPEI/shp_temp/Tomebamb_Shp.shp")
  ROI_SHP = terra::project(ROI_SHP, "EPSG:32717")
  dir.save = "C:/Users/Jonna/Desktop/Sequias_GIZ/Indices_Sequia/Categorias_Sequias/"
  Raster_Base = terra::rast("C:/Users/Jonna/OneDrive/Proyectos/Bases de Datos/Interpolaciones/Precipitacion/diaria/Prec_5km_2014_2024.nc")
  
  maps = function(Evento) {
    data = frecuencia_eventos[frecuencia_eventos$Clasificacion == Evento, ]
    datos_espaciales = inner_join(data, cords, by = "Pixel")
    puntos = vect(datos_espaciales, geom = c("X", "Y"), crs = "EPSG:32717")
    
    Base = Raster_Base[[1]]
    Base = Base * 0
    raster_frecuencia = rasterize(puntos, Base, field = "Frecuencia", fun = "first")
    plot(raster_frecuencia)
    plot(ROI_SHP, add = TRUE)
    terra::writeCDF(raster_frecuencia, filename=paste0(dir.save, Evento, ".nc"), overwrite=TRUE)
  }
  message(paste0("Se ha generado graficos en base a la frecuencia, revise en: ", dir.save))
  
  Sequia_leve = maps("Sequia Leve")
  Sequia_moderada = maps("Sequia Moderada")
  Sequia_severa = maps("Sequia Severa")
  Sequia_Extrema = maps("Sequia Extrema")
  return(frecuencia_eventos)
}


################################################################################
##                        CARACTERIZACIÓN DE SEQUÍAS                          ##
caracterizar_sequias = function(dt, rangos, pixel) {
  dt[, Clasificacion := as.factor(Clasificacion)]
  
  # Usar dplyr::lag() en lugar de shift para marcar el inicio de cada evento de sequía por tipo
  dt[, drought_event := cumsum(Clasificacion != dplyr::lag(Clasificacion, default = "Casi Normal") & Clasificacion %like% "Sequia")]
  
  # Filtrar solo los datos de sequía y mantener los tipos de sequía por separado
  droughts = dt[Clasificacion %like% "Sequia"]
  
  drought_summary = droughts[, .(
    Duracion = .N,  # Número de días (filas) en el evento de sequía
    Severidad = sum(value),  # Suma de los valores (severidad acumulada con signo)
    Magnitud = sum(value) / .N,  # Promedio de los valores (severidad promedio por día)
  
    Fecha_Inicio = min(Fecha),  # Fecha de inicio del evento de sequía
    Fecha_Fin = max(Fecha)  # Fecha de fin del evento de sequía
  ), by = .(Clasificacion, drought_event)]
  
}
################################################################################
##                  Emparejamiento, exclusión de sequías                      ##
# Teoría de las corridas
teori_run = function(data, umbral = NULL){
  setDT(data)
  names(data)[2] = "value"
  
  if (is.null(umbral)) {
    data$drought = fifelse(data$value < -0.5, 1, 
                           fifelse(data$value > 0.5, 0, 2))
  } else {
    data$drought = fifelse(data$value < umbral, 1, 0)
  }

  data[, drought_event := cumsum(drought != data.table::shift(drought, fill = 0))]
  
  # Sequia
  droughts = data[drought == 1]
  drought_seco <- droughts[, .(
    Tipo = "Sequia",
    Duracion = .N,  
    Severidad = round(sum(value),3),  
    Magnitud = round(sum(value) / .N, 3),  
    Fecha_Inicio = min(Fecha),  
    Fecha_Fin = max(Fecha)  
  ), by = drought_event]
  drought_seco = drought_seco[!is.na(drought_seco$drought_event), ]
  
  # Húmedo 
  droughts = data[drought == 0]
  drought_humedo <- droughts[, .(
    Tipo = "Humedo",
    Duracion = .N,  
    Severidad = round(sum(value),  3),
    Magnitud = round(sum(value) / .N,  3),
    Fecha_Inicio = min(Fecha),  
    Fecha_Fin = max(Fecha)  
  ), by = drought_event]
  
  res_eventos = rbind(drought_seco, drought_humedo)
  res_eventos = res_eventos[order(res_eventos$drought_event), ]
  return(res_eventos)
} 

drought_grouping = function(res_eventos, tc, pc) {
  sequias = res_eventos[Tipo == "Sequia", ]
  eventos_agrupados = list()
  K = 1
   
  for (i in 1:(nrow(sequias) - 1)) {  # Cambié el índice de inicio del bucle para evitar errores
    if (K == nrow(sequias)) {
      break  # Evitar acceder a índices fuera de rango
    }
    
    evento_actual = sequias[K, ]
    evento_siguiente = sequias[K + 1, ]
    
    Fecha_inicio = evento_actual$Fecha_Inicio
    Fecha_fin = evento_siguiente$Fecha_Fin
    
    ti = as.numeric(difftime(evento_siguiente$Fecha_Inicio, evento_actual$Fecha_Fin - 1, units = "days"))
    s = evento_actual$Severidad
    
    # Verificar si hay humedad entre Fecha_inicio y Fecha_fin
    humedad = res_eventos[Tipo == "Humedo" & Fecha_Inicio >= Fecha_inicio & Fecha_Fin <= Fecha_fin, ]
    
    if (nrow(humedad) > 0) {
      vi = sum(humedad$Severidad)
    } else {
      vi = 0
    }
    
    pi = vi / s
    
    if (ti <= tc & pi <= pc) {
      evento_agrupado = data.frame(
        Tipo = evento_actual$Tipo,
        Fecha_Inicio = Fecha_inicio,
        Fecha_Fin = Fecha_fin,
        Duracion = evento_actual$Duracion + evento_siguiente$Duracion + ti,
        Severidad = evento_actual$Severidad + evento_siguiente$Severidad - vi
      )
      
      
      z = 1
      contador = 0
      while (z != 0) {
        if ((K + z + 1) > nrow(sequias)) {
          break  # Evitar acceder a índices fuera de rango
        }
        
        Fecha_fin_2 = sequias[K + z + 1, ]$Fecha_Fin
        
        evento_actual = evento_agrupado
        evento_siguiente = sequias[K + z + 1, ]
        
        ti = as.numeric(difftime(evento_siguiente$Fecha_Inicio, evento_actual$Fecha_Fin - 1, units = "days"))
        s = evento_agrupado$Severidad
        
        Fecha_inicio_humedad = evento_actual$Fecha_Fin
        Fecha_fin_humedad = evento_siguiente$Fecha_Inicio
        
        humedad = res_eventos[Tipo == "Humedo" & Fecha_Inicio >= Fecha_inicio_humedad & Fecha_Fin <= Fecha_fin_humedad, ]
        
        if (nrow(humedad) > 0) {
          vi = sum(humedad$Severidad)
        } else {
          vi = 0
        }
        
        pi = vi / s
        
        if (ti <= tc & pi <= pc) {
          evento_agrupado = data.frame(
            Tipo = evento_actual$Tipo,
            Fecha_Inicio = Fecha_inicio,
            Fecha_Fin = Fecha_fin_2,
            Duracion = sum(evento_actual$Duracion, evento_siguiente$Duracion, ti),
            Severidad = evento_actual$Severidad + evento_siguiente$Severidad - vi
          )
          z = z + 1
          contador = contador + 1
        } else {
          z = 0
          contador = contador + 2
        }
      } # Cierro wl while
      
      K = K + contador
      eventos_agrupados[[i]] = evento_agrupado # A1qui estoy dentro del if, (dehberioa guardar aqui los datos)
    } else {
      eventos_agrupados[[i]] = evento_actual
      K = K + 1
      
    }
  }
  return(eventos_agrupados)  # Agregar retorno explícito
  
}

# Calculo del tiempo critico (tc) y ratio critico (pc)
tc_pc = function (SPEI, SSI){
  names(SSI) = c("Fecha", "value")
  umbral = quantile(SSI$value, 0.7)
  
  # Identificar sequias en base al umbral'
  sequias = teori_run(SSI, umbral = umbral)
  parametros = expand.grid(tc = c(0, 1, 5, 10,15,20), pc = c(0, 0.025, 0.1, 0.15, 0.2, 0.25,0.5, 0.75, 1))
  
  resultados_lista = lapply(1:nrow(parametros), function(i) {
    tc = parametros$tc[i]
    pc = parametros$pc[i]
    return(drought_grouping(sequias, tc, pc))
  })
  
  # saco el promedio de la duracuin
  resultados_finales = list()
  for (i in 1:length(resultados_lista)) {
    resultados_finales[[i]] = rbindlist(resultados_lista[[i]], fill = TRUE)
  }
  
  # zscore_duration = function(durations) {
  #   (durations - mean(durations)) / sd(durations)
  # }
  # 
  # rest_estandarizados = list()
  # for (i in 1:length(resultados_finales)) {
  #   rest_estandarizados[[i]] = resultados_finales[[i]][, .(
  #     tc = unique(tc),
  #     pc = unique(pc),
  #     Dur_Z = zscore_duration(Duracion),
  #     Sev_Z = zscore_duration(Severidad)
  #   )]
  # }
  
  
  resuls = list()
  for (i in 1:length(resultados_finales)) {
    resuls[[i]] = resultados_finales[[i]][, .(
      tc = unique(tc),
      pc = unique(pc),
      Duracion = mean(Duracion),
      Severidad = mean(Severidad)
    )]
  }
  
  resultados_df = do.call(rbind, lapply(resuls, as.data.frame))
  # resultados_df$Duracion = round(resultados_df$Duracion, 2)
  # resultados_df$Severidad = round(resultados_df$Severidad, 2)
  
  ##############################################################################
  # Gráfico de media_duracion vs p_c
  p1 = ggplot(resultados_df, aes(x = pc, y = Duracion, color = factor(tc))) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = seq(0, 0.25, by = 0.05)) +
    labs(title = "Gráfico de media_duracion vs p_c",
         x = "p_c",
         y = "media_duracion",
         color = "t_c") +
    theme_minimal()
  p1
  
  # Gráfico de media_deficit vs p_c
  ggplot(resultados_df, aes(x = pc, y = Severidad, color = factor(tc))) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1)) +
    labs(title = "Gráfico de media_deficit vs p_c",
         x = "p_c",
         y = "media_deficit",
         color = "t_c") +
    theme_minimal()
}



Drought_exclusion = function() {
  
}

