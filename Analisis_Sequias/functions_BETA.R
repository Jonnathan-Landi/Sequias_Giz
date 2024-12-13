################################################################################
#   Análisis Integral de la Recurrencia y Respuesta Espacial de las Sequías    #
#        Meteorológicas e Hidrológicas en la Región de Cuenca, Ecuador         #
################################################################################
# Autores: 
#'  Jonnathan Landi
#'  Marco Mogro 
# Fecha: 2024-11-15
# Versión: 2.0.0
################################################################################
# Librerías necesarias
packages = c("data.table", "dplyr", "purrr", "terra", "ggplot2", "plotly",
             "htmlwidgets", "gridExtra", "shiny", "shinythemes", "bslib", 
             "shinyWidgets", "lubridate", "tidyr", "reshape2","matrixStats",
             "changepoint", "zyp", "trend", "Kendall", "missForest", "parallel", "copula",
             "doParallel", "SPEI", "openair", "future.apply", "extRemes", "TLMoments",
             "fitdistrplus", "ismev", "LMoFit", "VineCopula", "extremeStat", "lmomco")

install_packages = function(pkgs) {
  for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    } else {
      library(pkg, character.only = TRUE)
    }
  }
}
install_packages(packages)
################################################################################
#                         Preparación de los datos                             #
Pixeles_cuenca = function(data, type) {
  names(data)[1] = "Fecha"
  if (type == "Yanuncay") {
    data[, .(Fecha, Pixel_8, Pixel_9, Pixel_14, Pixel_15,Pixel_16,Pixel_17,Pixel_18,Pixel_19,Pixel_22, Pixel_23, 
             Pixel_24, Pixel_25, Pixel_26, Pixel_27, Pixel_28,Pixel_31,
             Pixel_32, Pixel_33, Pixel_34, Pixel_35, Pixel_36, Pixel_39, Pixel_40, Pixel_41, Pixel_42, 
             Pixel_44, Pixel_45, Pixel_46),]
    
  } else if (type == "Tomebamba") {
    data = data[, .(Fecha, Pixel_1, Pixel_2, Pixel_3, Pixel_4, Pixel_5, Pixel_6, Pixel_7, Pixel_8, Pixel_9,
                    Pixel_10, Pixel_11, Pixel_12, Pixel_13, Pixel_14, Pixel_15, Pixel_16, Pixel_17,
                    Pixel_18,Pixel_19, Pixel_20, Pixel_21, Pixel_23, Pixel_24, Pixel_25, Pixel_26,
                    Pixel_27, Pixel_28, Pixel_29, Pixel_30, Pixel_38), ]
    
  } else {
    stop(paste0(type, " no es una subcuenca válida"))
  }
  names(data)[2:length(data)] = paste0("Pixel_", 1:length(data))
  return(data)
}

import_data = function(data, type, Fecha_Inicio, Fecha_Fin, Subcuenca, exp = NULL) {
  # compruebo que data sea un archivo data.table
  if (!is.data.table(data)) {
    stop("El archivo de datos debe ser un data.table, use 'fread' para cargar los datos")
  }
  
  if (type == "SPEI") {
    names(data)[1] = "Fecha"
    data[data == Inf | data == -Inf] = NA
    data = data[Fecha >= Fecha_Inicio & Fecha <= Fecha_Fin, ]
    
    if (!is.null(exp)) {
      if (Subcuenca == "Tomebamba") {
        print("Se ha importado los datos de SPEI para Tomebamba")
        data = Pixeles_cuenca(data, type = "Tomebamba")
      } else {
        print("Se ha importado los datos de SPEI para Yanuncay")
        data = Pixeles_cuenca(data, type = "Yanuncay")
      }
    }
  } else {
    data = data[Fecha >= Fecha_Inicio & Fecha <= Fecha_Fin, ]
  }
  return(data)
}

################################################################################
#                         Funciones complementarias                            #
anual_duration = function(SPEI, cords, dir.save, Subcuenca) {
  res_final = data.frame()
  
  cats_spei = function(x) {
    fcase(
      x <= -2, "Sequia Extrema",
      x > -2 & x <= -1.5, "Sequia Severa",
      x > -1.5 & x <= -1, "Sequia Moderada",
      x > -1 & x <= -0.5, "Sequia Leve",
      x > -0.5, "No Sequia"
    )
  }
  
  for (pixel in setdiff(names(SPEI), "Fecha")) {
    
    res_pixel = SPEI[, .(Fecha, value = get(pixel))]
    res_pixel$Categoria = sapply(res_pixel$value, cats_spei)
    
    anual = res_pixel %>% 
      group_by(year(Fecha), Categoria) %>% summarise(Duracion = n())
    
    anual_transformada = anual %>%
      tidyr::pivot_wider(
        names_from = `year(Fecha)`, # Los nombres de las columnas serán los años
        values_from = Duracion,    # Los valores vendrán de la columna Duracion
        values_fill = list(Duracion = 0) # Rellenar con 0 si no hay datos
      )
    
    pixel_un = data.frame((rep(pixel, 5)))
    names(pixel_un) = "Pixel"
    exportar = cbind(pixel_un, anual_transformada)
    
    res_final = rbind(res_final, exportar)
  }
  
  save_folder = paste0(dir.save, "/Shapes_Duracion_Anual")
  if (!dir.exists(save_folder)) {
    dir.create(save_folder)
  }
  
  data_save = merge(res_final, cords, by = "Pixel")
  puntos = vect(data_save, geom = c("X", "Y"), crs = "EPSG:32717")
  writeVector(puntos, paste0(save_folder, "/duracion_anual_", Subcuenca, ".shp", sep = ""), filetype = "ESRI Shapefile", overwrite = TRUE)
  
}

caracterizar_Min_Max = function(x) {
  x = abs(x)
  
  x_min = min(x)
  x_max = max(x)
  res = data.frame()
  
  for (i in 1:length(x)) {
    Xn = (x[i] - x_min) / (x_max - x_min)
    res = rbind(res, Xn)
  }
  
  names(res) = "X"
  
  # Calcular los límites de los intervalos basados en los cuartiles
  cuartiles = quantile(res$X, probs = seq(0, 1, by = 0.2))
  
  clases = function(y) {
    fcase(
      y <= cuartiles[2], "Muy bajo",
      y > cuartiles[2] & y <= cuartiles[3], "Bajo",
      y > cuartiles[3] & y <= cuartiles[4], "Medio",
      y > cuartiles[4] & y <= cuartiles[5], "Alto",
      y > cuartiles[5], "Muy alto"
    )
  }
  
  # Clasificar los datos en intervalos según los cuartiles
  res$Categoria = sapply(res$X, clases)
  res = data.table(res)
  
  final = cbind(x, res[,.(Categoria = Categoria)])
}

################################################################################
#                               Cálculo de SPEI                                #
#' balance_hidrico = archivo con los datos de balance hidrico formato df
#' vect = vector con los valores de n para el cálculo de SPEI (valores a cumular en dias)
# Función.
spei = function(balance_hidrico, vect) {
  names(balance_hidrico)[1] = "Fecha"
  
  SPEI_Pixel = function(Pixel, n) {
    data = balance_hidrico[,.(Fecha = Fecha, value = get(Pixel))]
    
    if (any(is.na(data$value))) {
      stop("Series con NA")
    }
    
    SPEI = SPEI::spei(data$value, scale = n, fit = "ub-pwm", verbose = F)
    data$SPEI = SPEI$fitted
    data[data == Inf | data == -Inf] = NA
    
    # if (sum(is.na(data$SPEI)) > n) {
    #   message(paste0("Se encontro," , sum(is.na(data$SPEI)), " Nas"))
    # }
    
    data = data[,c("Fecha", "SPEI")]
    return(data)
  }
  
  ejecutar = function(n) {
    Fechas = data.table(balance_hidrico)[,.(Fecha = Fecha)]
    for (i in setdiff(names(balance_hidrico), c("Fecha"))) {
      SPEI = SPEI_Pixel(Pixel = i, n = n)
      names(SPEI)[2] = i
      Fechas = merge(Fechas, SPEI, by = "Fecha")
    }
    return(Fechas)
  }
  
  resultados = list()
  k = 1
  for (i in vect) {
    message(paste0("Calculando SPEI con n = ", i))
    resultados[[k]] = ejecutar(n = i)
    k = k + 1
  }
  
  names(resultados) = paste0("SPEI_", vect)
  return(resultados)
}

spei_parallel = function(balance_hidrico, vect) {
  names(balance_hidrico)[1] = "Fecha"
  
  SPEI_Pixel = function(Pixel, n) {
    data = balance_hidrico[, .(Fecha = Fecha, value = get(Pixel))]
    
    if (any(is.na(data$value))) {
      stop("Series con NA")
    }
    
    SPEI = SPEI::spei(data$value, scale = n, fit = "ub-pwm", verbose = FALSE)
    data$SPEI = SPEI$fitted
    data[data == Inf | data == -Inf] = NA
    data = data[, .(Fecha, SPEI)]
    return(data)
  }
  
  ejecutar = function(n) {
    Fechas = data.table(balance_hidrico)[, .(Fecha = Fecha)]
    SPEI_results = future_lapply(
      setdiff(names(balance_hidrico), "Fecha"),
      function(i) {
        SPEI = SPEI_Pixel(Pixel = i, n = n)
        names(SPEI)[2] = i
        return(SPEI)
      }
    )
    
    # Combina resultados por columna
    for (SPEI in SPEI_results) {
      Fechas = merge(Fechas, SPEI, by = "Fecha", all = TRUE)
    }
    return(Fechas)
  }
  
  # Paralelizar el bucle principal
  resultados = future_lapply(vect, function(i) {
    message(paste0("Calculando SPEI con n = ", i))
    ejecutar(n = i)
  })
  
  names(resultados) = paste0("SPEI_", vect)
  return(resultados)
}

################################################################################
#                                Teori_Run                                     #
#' data = Base de datos con la columna Fecha y valor a clasificar 
#' umbral = umbral para considerar sequía y periodos húmedos. default = -0.5 y 0.5 
#' Vrecor = Lógico, si se desea considerar como sequía únicamente el área por debajo del umbral establecer en T, 
#' si se desea conservar todo no hacer nada. default = considera todo el valor como sequía. 
# Función. 
teori_run = function(data, umbral = NULL, Vrecor = NULL) {
  
  setDT(data)
  
  if (ncol(data) > 2) {
    stop("El data.table debe tener solo dos columnas")
  }
  
  names(data)[2] = "value"
  
  if (any(is.na(data))) {
    message("Los datos tienen NA, esto podria causar datos erroneos. A continuación se omitirá")
    data = na.omit(data)
  }
  
  if(any(is.na(data$value))) {
    stop("Vacíos encontrados, no es posible ejecutrar 'teori_run' con valores NA")
    message("Advertencia: Se han encontrado valores NA en los datos, esto podría alterar los resultados")
  }
  
  ################################################################################
  if (!is.null(Vrecor)) {
    if (is.null(umbral)) {
      message("Clasificación basado en el área bajo la curva, umbrales = [-0.5, 0.5]")
      data$drought = fifelse(data$value < -0.5, 1,
                             fifelse(data$value > 0.5, 0, 2))
      
      data$Magnitud = fifelse(data$value < -0.5, 0.5 + data$value,
                              fifelse(data$value > 0.5, data$value - 0.5, 0))
    } else {
      message(paste0("Clasificación basado en el area bajo la curva, umbrales = [", umbral, "]"))
      data$drought = fifelse(data$value < umbral, 1, 0)
      data$Magnitud = fifelse(data$value < umbral, 0.5 + data$value, data$value)
    }
  } else {
    if (is.null(umbral)) {
      message("Clasificación considerando como deficit todo el valor")
      data$drought = fifelse(data$value < -0.5, 1, 
                             fifelse(data$value > 0.5, 0, 2))
      
      data$Magnitud = data$value
      
    } else {
      data$drought = fifelse(data$value < umbral, 1, 0)
      data$Magnitud = data$value
    }
  } 
  
  data[, drought_event := cumsum(drought != data.table::shift(drought, fill = 0))]
  
  # Sequía
  droughts = data[drought == 1]
  drought_seco = droughts[, .(
    Tipo = "Sequia",
    Duracion = .N,  
    Magnitud = round(sum(Magnitud),3),   
    Severidad = round(sum(Magnitud) / .N, 3),  
    Fecha_Inicio = min(Fecha),  
    Fecha_Fin = max(Fecha)  
  ), by = drought_event]
  
  if (any(is.na(drought_seco))) {
    message("Advertencia: Se han encontrado valores NA en los datos, esto podría alterar los resultados")
    drought_seco = drought_seco[!is.na(drought_seco$drought_event), ]
  }
  
  # Húmedo 
  droughts = data[drought == 0]
  drought_humedo = droughts[, .(
    Tipo = "Humedo",
    Duracion = .N,  
    Magnitud = round(sum(Magnitud),3),   
    Severidad = round(sum(Magnitud) / .N, 3), 
    Fecha_Inicio = min(Fecha),  
    Fecha_Fin = max(Fecha)  
  ), by = drought_event]
  
  if (any(is.na(drought_seco$drought_event)) | any(is.na(drought_humedo$drought_event))) {
    message("Advertencia: Se han encontrado valores NA en los datos, esto podría alterar los resultados")
  }
  
  resultados = rbind(drought_seco, drought_humedo)
  resultados = resultados[order(resultados$drought_event), ]
  #resultados$Tipo = ifelse(resultados$Magnitud  <= -0.5, "Sequia", "Humedo")
  resultados = resultados[, drought_event := NULL]
  return(resultados)
} 

################################################################################
#                             Drought_Grouping                                 #
#' res_grouping = resultado del agrupamiento
#' data_original = datos del SSI
#' dir.save = dirección donde se va a guardar los graficos 
#' Subcuenca = Nombre de la subcuenca ("Tomebamba o Yanuncay")
# Función. 
graph_grouping = function(res_grouping, data_original, dir.save, Subcuenca) {
  
  eventos_agrupados = res_grouping
  
  data_graph = eventos_agrupados %>%
    mutate(Fecha_Inicio = as.Date(Fecha_Inicio),
           Fecha_Fin = as.Date(Fecha_Fin)) %>%
    mutate(Year = format(as.Date(Fecha_Fin), "%Y"))  # Extrae el año de la fecha
  
  data_graph = data_graph %>%
    group_by(Year) %>%
    mutate(
      year_start = as.Date(paste0(Year, "-01-01")),
      year_end = as.Date(paste0(Year, "-12-31"))
    ) %>%
    ungroup()
  
  # Fcuncionm para generar graficos
  graph = function(anio =NULL, f = F, Type = FALSE, ind = F, save = save) {
    
    if (Type == T) {
      seq = seq.Date(as.Date(paste0(anio, "-01-01")), as.Date(paste0(anio, "-12-31")), by = "day")
      seq = data.table(Fecha_Inicio = seq)
      datos_graph = data_graph[data_graph$Year == anio, ]
    } else {
      datos_graph = data_graph
    }
    
    fechas_i = as.Date(c(datos_graph$Fecha_Inicio))
    fechas_f = as.Date(c(datos_graph$Fecha_Fin))
    Severidad = c(datos_graph$Magnitud)
    
    ############################################################################
    #                               Datos insitu                               #
    names(data_original) = c("Fecha_Inicio", "value")
    if (Type == T) {
      data_C = merge(seq, data_original, by = "Fecha_Inicio", all = TRUE)
      data_C = data_C[data_C$Fecha >= as.Date(paste0(anio, "-01-01")) & data_C$Fecha <= as.Date(paste0(anio, "-12-30")),]
    } else {
      data_C = data_original
    }
    
    data_C$Tipo = rep("SSI", nrow(data_C))
    highlight_periods = data.frame(
      xmin = fechas_i,
      xmax = fechas_f,
      ymin = -0.5,
      ymax = c(datos_graph$Magnitud),
      Tipo = rep("Sequía", nrow(datos_graph))
    )
    ############################################################################
    #                               GRafico                                    #
    pl = ggplot() +
      geom_line(data = data_C, aes(
        x = Fecha_Inicio, y = value, color = Tipo,
        group = Tipo, 
        text = paste0("Fecha: ", Fecha_Inicio, "<br>",
                      "SSI: ", value, "<br>",
                      "Tipo: ", Tipo)
      ), linewidth = 1) +
      
      geom_rect(data = highlight_periods, aes(
        xmin = xmin, xmax = xmax, 
        ymin = ymin, ymax = ymax, 
        color = Tipo,
        text = paste0("Fecha inicio: ", xmin, "<br>",
                      "Fecha fin: ", xmax, "<br>",
                      "Magnitud: ", ymax, "<br>",
                      "Tipo: ", Tipo)
      ), 
      fill = "white", linewidth = 0.9, alpha = 0.3, inherit.aes = FALSE) +
      scale_color_manual(values = c("SSI" = "blue", "Sequía" = "red")) +
      geom_hline(yintercept = -0.5, linetype = "dashed", color = "black", linewidth = 1) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", linewidth = 1) +
      
      theme_bw() +
      
      labs(title = ifelse(Type == TRUE, anio, ""),
           x = ifelse(ind == TRUE, "Mes", "Año"), y = ifelse(f == TRUE, "Magnitud", "")) +
      # Ajustar temas de texto y ejes ind
      theme(
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 18, hjust = 0.9),
        legend.position = ifelse(Type == TRUE, "none", "bottom"),
        # legend.justification = c(0.5, 0.5),
        # legend.background = element_rect(fill = "white", colour = "black", linewidth = 1),
        plot.title = element_text(hjust = 0.5, size = rel(2), lineheight = .9,
                                  face = "bold", colour = "brown"),
      ) +
      
      scale_x_date(
        date_breaks = ifelse(Type == TRUE, "1 month", "1 year"),
        labels = function(x) {
          if (Type == TRUE) {
            meses = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
            meses[as.numeric(format(x, "%m"))]
          } else {
            as.numeric(format(x, "%Y"))
          }
        }
      ) +
      
      # Borde del gráfico
      geom_rect(aes(xmin = as.Date(-Inf), xmax = as.Date(Inf), ymin = -Inf, ymax = Inf),
                color = "black", fill = NA, size = 1.5)
    return(pl) 
  } # Cierro la funcion de graficos 
  
  message("Generando gráficos, porfavor espere....")
  grafico_general = graph(f = T, Type = F, ind = F)
  
  # Genero graficos por año
  P_1 = graph(2014, f = T, Type = T, ind = T)
  P_2 = graph(2015, f = F, Type = T, ind = T)
  P_3 = graph(2016, f = F, Type = T, ind = T)
  P_4 = graph(2017,f = T, Type = T, ind = T)
  P_5 = graph(2018,f = F, Type = T, ind = T)
  P_6 = graph(2019, f = F, Type = T, ind = T)
  P_7 = graph(2020, f = T, Type = T, ind = T)
  P_8 = graph(2021, f = F, Type = T, ind = T)
  P_9 = graph(2022, f = F, Type = T, ind = T)
  P_10 = graph(2023, f = T, Type = T, ind = T)
  P_11 = graph(2024, f = F, Type = T, ind = T)
  
  combined_plot = grid.arrange(
    grobs = list(P_1, P_2, P_3, P_4, P_5, P_6, P_7, P_8, P_9, P_10, P_11),
    ncol = 3,  # Número de columnas
    nrow = 4  # Número de filas
  )
  
  folder_save = paste0(dir.save, "/Graph_gruping_SSI/")
  if (!dir.exists(folder_save)) {
    dir.create(folder_save)
  }
  
  ggsave(
    paste(folder_save, Subcuenca, "_Grupal.png", sep = ""),
    plot = combined_plot,
    width = 12,
    height = 8,
    units = "in",
    dpi = 1000,
    bg = NULL,
  )
  
  ggsave(
    paste(folder_save, Subcuenca, "_General.png", sep = ""),
    plot = grafico_general,
    width = 12,
    height = 8,
    units = "in",
    dpi = 1000,
    bg = NULL,
  )
  
  
  # Genero gráfico interactivo
  pl_interactive = ggplotly(grafico_general, tooltip = "text")
  print(pl_interactive)
  name = paste0(folder_save, "Sequia_Iteractivo", Subcuenca, ".html")
  saveWidget(pl_interactive, file = name)
  
}

#' res_corridas = archivo formato data.table con resultados de teori_run
#' tc = Duración critica '(Solo para sequias hidrologicas)
#' pc = Magnitud critica (Solo para sequias hidrologicas)
#' graficar = Lógico, si desea generar graficos establezca en T
# Función. 
drought_grouping = function(res_corridas, tc, pc, graficar = FALSE, dir.save = NULL, Subcuenca = NULL, data_original = NULL) {
  
  calcular_intervalo = function(fecha_inicio, fecha_fin) {
    return(as.numeric(fecha_inicio - fecha_fin) - 1)
  }
  
  calcular_humedad = function(fecha_inicio, fecha_fin) {
    humedad = res_corridas[Tipo == "Humedo" & Fecha_Inicio >= fecha_inicio & Fecha_Fin <= fecha_fin, ]
    return(ifelse(nrow(humedad) > 0, sum(humedad$Magnitud), 0))
  }
  
  sequias = res_corridas[Tipo == "Sequia", ]
  eventos_agrupados = data.frame()
  indice_evento = 1
  
  for (i in 1:(nrow(sequias) - 1)) {
    
    if (indice_evento == nrow(sequias)) {
      break
    }
    
    evento_actual = sequias[indice_evento,]
    evento_siguiente = sequias[indice_evento + 1, ]
    
    ti = calcular_intervalo(evento_siguiente$Fecha_Inicio, evento_actual$Fecha_Fin)
    s = evento_actual$Magnitud
    # prueba
    s = abs(s)
    
    vi = calcular_humedad(evento_actual$Fecha_Inicio, evento_siguiente$Fecha_Fin)
    if (ti <= tc & (vi/s) <= pc) {
      
      evento_agrupado = data.frame(
        Tipo = evento_actual$Tipo,
        Fecha_Inicio = evento_actual$Fecha_Inicio,
        Fecha_Fin = evento_siguiente$Fecha_Fin,
        Duracion = evento_actual$Duracion + evento_siguiente$Duracion + ti,
        Magnitud = evento_actual$Magnitud + evento_siguiente$Magnitud + vi
      )
      
      # Compruebo que si hay mas eventos para agrupar
      iniciador = T
      n_agrupaciones = 1
      while (iniciador == T) {
        if (indice_evento + n_agrupaciones == nrow(sequias)) {
          break
        }
        
        evento_actual_2 = evento_agrupado
        evento_siguiente_2 = sequias[indice_evento + n_agrupaciones + 1,]
        
        # cálculos de ti
        ti = calcular_intervalo(evento_siguiente_2$Fecha_Inicio, evento_actual_2$Fecha_Fin)
        s = evento_actual_2$Magnitud
        
        # prueba
        s = abs(s)
        
        vi = calcular_humedad(evento_actual_2$Fecha_Fin, evento_siguiente_2$Fecha_Fin)
        
        if (ti <= tc & (vi/s) <= pc) {
          
          evento_agrupado = data.frame(
            Tipo = evento_actual_2$Tipo,
            Fecha_Inicio = evento_actual_2$Fecha_Inicio,
            Fecha_Fin = evento_siguiente_2$Fecha_Fin,
            Duracion = evento_actual_2$Duracion + evento_siguiente_2$Duracion + ti,
            Magnitud = evento_actual_2$Magnitud + evento_siguiente_2$Magnitud + vi
          )
          
          iniciador = T
        } else {# cierro el if dentro del while
          # indice_evento = indice_evento + 2
          iniciador = FALSE
        } # cierro el else del if que esta dentro del while
        n_agrupaciones = n_agrupaciones + 1
      } # cierro el while
      
      indice_evento = indice_evento + n_agrupaciones
    } else {
      evento_agrupado = data.frame(
        Tipo = evento_actual$Tipo,
        Fecha_Inicio = evento_actual$Fecha_Inicio,
        Fecha_Fin = evento_actual$Fecha_Fin,
        Duracion = evento_actual$Duracion,
        Magnitud = evento_actual$Magnitud
      )
      
      indice_evento = indice_evento + 1
    }
    
    eventos_agrupados = rbind(eventos_agrupados, evento_agrupado)
    
  } 
  
  if (graficar) {
    message("Se ha confifurado el modo gráfico")
    graph_grouping(res_grouping = eventos_agrupados, data_original = data_original,
                   dir.save = dir.save, Subcuenca = Subcuenca )
  }
  return(eventos_agrupados)
} # cierro función de agrupamiento de sequias

################################################################################
#                     Exclusión de eventos de sequía (EES)                     #
#' seq_grup = resultado del agrupamiento de sequias
exclude_droughts = function (seq_grup) {
  seq_grup$Severidad = seq_grup$Magnitud / seq_grup$Duracion
  
  eventos = data.frame()
  d_avg = abs(mean(seq_grup$Duracion))
  s_avg = abs(mean(seq_grup$Magnitud))
  
  di = d_avg * 0.3
  si = s_avg * 0.3 
  
  for (i in 1:nrow(seq_grup)) {
    d_act = seq_grup$Duracion[i]
    s_act = abs(seq_grup$Magnitud[i])
    
    if (d_act < di | s_act < si) {
      next
    } else {
      eventos = rbind(eventos, seq_grup[i,])
    }
    
  }
  return(eventos)
  #  Si la duración de la sequía dYo era menor que rd × dAv o la severidad de la sequía sYo era menor que rs × sAVE, los eventos de sequía {dYo, sYo} 
}
################################################################################
#               Identificacion, agrupacion y exclusión de sequias (IAE)        #
#' data = SSI / SPI. De ser Ser un dt de solo dos columnas, Fecha y SSI para SPEI puede tener varias columnas.
#' tc y pc son parametros que se usaran para Drought_Grouping (Parametros explicados en 'Drought_Grouping')
#' Vrecor y Umbral = Utilizado para teori_run (Parametros explicados en 'Teori_Run')
# Función. 
IAE_droughts = function(data, Type, tc = NULL, pc = NULL, Vrecor, umbral = NULL, exclude = NULL) {
  if (Type == "SSI") {
    run = teori_run(data, Vrecor = Vrecor)
    grouping = drought_grouping(res_corridas = run, tc = tc, pc = pc)
    
    if (!is.null(exclude)) {
      grouping = exclude_droughts(grouping)
    }
    
  } else {
    grouping = list()
    for (Pixel in setdiff(colnames(data), "Fecha")) {
      D = data[,.(Fecha = Fecha, SPEI = get(Pixel))]
      run = teori_run(D, umbral = umbral, Vrecor = Vrecor)
      
      # Exclusión de eventos de sequía
      agrup = run[(Tipo == "Sequia") & (Duracion >= 3),]
      grouping[[Pixel]] = agrup
    }
  }
  return(grouping)
}

################################################################################
#                 Statistical_analysis_of_droughts (SAD)                       #
# SAD es un análisis estadistico que tiene los siguientes analisis
#' 1. Análisis de tendencia de sequías hidro y meteo
#' 2. Análisis de tasas de cambio de sequías hidro y meteo
#' 3. Análisis de ruptura o puntos de cambio. 
#' balance_hídrico = archivo con los datos de balance hidrico formato df (solo SPEI)
#' SPEI = archivo con los datos de SPEI formato df (Solo SPEI)
#' SSI = archivo con los datos de SSI formato df (Solo SSI)
#' type = Tipo de análisis a realizar (SPEI, SSI)
#' cords = Coordenadas de los pixeles (Solo SPEI)
#' dir.save = Dirección donde se guardaran los resultados
#' Subcuenca = Nombre de la subcuenca ("Tomebamba o Yanuncay")
#' vect = Vector con los valores de n para el cálculo de SPEI (valores a cumular en días, Solo SPEI)

SAD = function(balance_hidrico = NULL, SPEI = NULL, SSI = NULL, type, cords = NULL, dir.save, Subcuenca, vect = NULL) {
  if (type == "SPEI") {
    
    if (Subcuenca == "Tomebamba") {
      data = Pixeles_cuenca(balance_hidrico, type = "Tomebamba")
    } else {
      data = Pixeles_cuenca(balance_hidrico, type = "Yanuncay")
    }
    
    ############################################################################
    #                       Funciones que se ejecutaran                        #
    tendencia = function(data, name, Subcuenca) {
      resultados_tendencia = data.frame()
      resultados_pendiente = data.frame()
      puntos_cambio = data.frame()
      i = 1
      
      for (Pixel in setdiff(colnames(data), "Fecha")) {
        data_analisis = data[,.(Fecha = Fecha, SPEI = get(Pixel))]
        spei_data = data_analisis$SPEI
        spei_data = na.omit(spei_data)
        
        spei_ts = ts(spei_data, start = c(2014, as.numeric(format(data_analisis$Fecha[1], "%j"))), frequency = 365)
        decomposed = stl(spei_ts, s.window = "periodic", na.action = na.omit)
        trend_component = decomposed$time.series[, "trend"]
        
        # calculo de tendencia
        result = MannKendall(trend_component)
        res_global = data.frame(
          cod = i,
          Pixel = Pixel,
          tau = as.numeric(result[1]),
          p_value = as.numeric(result[2]),
          Tipo = name
        )
        resultados_tendencia = rbind(resultados_tendencia, res_global)
        
        # Calculo de los puntos de cambio
        pettitt_result = pettitt.test(trend_component)
        change_point_date = as.Date(data$Fecha[1]) + (pettitt_result$estimate - 1)
        change_point_date = as.Date(change_point_date)
        change_point_date = change_point_date[1]
        
        res_pettitt = data.frame(
          cod = i,
          Pixel = Pixel,
          estadistico = as.numeric(pettitt_result$statistic),
          p_value = as.numeric(pettitt_result$p.value),
          fecha_cambio = change_point_date,
          Tipo = name
        )
        puntos_cambio = rbind(puntos_cambio, res_pettitt)
        
        # Calculo de la pendiente
        pdns = data.frame(
          time = 1:length(trend_component),
          trend_component = trend_component
        )
        
        m_sens = zyp.sen(trend_component ~ time, data = pdns)
        pendiente = m_sens$coefficients[[2]]
        
        res_pendiente = data.frame(
          cod = i,
          Pixel = Pixel,
          pendiente = as.numeric(pendiente),
          Tipo = "Pendiente"
        )
        
        resultados_pendiente = rbind(resultados_pendiente, res_pendiente)
        
        i = i + 1
      } # cierro el for
      
      # guardo todo
      save_tendencia = paste0(dir.save, "/Tendencia_SPEI")
      if (!dir.exists(save_tendencia)) {
        dir.create(save_tendencia)
      }
      
      # Guardo mis resultados --------------------------------------------------
      # Tendencia
      data_save_tendencia = merge(resultados_tendencia, cords, by = "Pixel")
      data_save_tendencia = data_save_tendencia %>% dplyr::arrange(cod) %>% dplyr::select(-cod)
      puntos = vect(data_save_tendencia, geom = c("X", "Y"), crs = "EPSG:32717")
      writeVector(puntos, paste0(save_tendencia, "/Tendencias_", name, "_", Subcuenca, ".shp", sep = ""), filetype = "ESRI Shapefile", overwrite = TRUE)
      
      # Pendiente
      data_save_pendiente = merge(resultados_pendiente, cords, by = "Pixel")
      data_save_pendiente = data_save_pendiente %>% dplyr::arrange(cod) %>% dplyr::select(-cod)
      puntos_m = vect(data_save_pendiente, geom = c("X", "Y"), crs = "EPSG:32717")
      writeVector(puntos_m, paste0(save_tendencia, "/Pendientes_", name, "_", Subcuenca, ".shp", sep = ""), filetype = "ESRI Shapefile", overwrite = TRUE)
      
      # puntos de cambio
      data_sava_Pcambio = merge(puntos_cambio, cords, by = "Pixel")
      data_sava_Pcambio = data_sava_Pcambio %>% dplyr::arrange(cod) %>% dplyr::select(-cod)
      data_sava_Pcambio$fecha_cambio = as.Date(data_sava_Pcambio$fecha_cambio)
      data_sava_Pcambio$fecha_caracter = as.character(data_sava_Pcambio$fecha_cambio)
      
      puntos_pc = vect(data_sava_Pcambio, geom = c("X", "Y"), crs = "EPSG:32717")
      writeVector(puntos_pc, paste0(save_tendencia, "/PuntosCambio_", name, "_", Subcuenca, ".shp", sep = ""), filetype = "ESRI Shapefile", overwrite = TRUE)
      
    } # cierro la función de tendencia
    ############################################################################
    #                                   Ejecución                              #
    analisis = spei(data, vect)
    for (i in 1:(length(analisis) + 1)) {
      message(paste0("Ejecutando análisis ", i))
      if (i == 1) {
        data = SPEI
        tendencia(data, name = paste0("dia_", i), Subcuenca = Subcuenca)
      } else {
        data = analisis[[i-1]]
        tendencia(data, name = paste0("dia_", vect[i-1]), Subcuenca = Subcuenca)
      }
      
    }
    
  } else {
    resultados_tendencia = data.frame()
    resultados_pendiente = data.frame()
    data = SSI
    names(data) = c("Fecha", "value")
    ############################################################################
    #                           Análisis de tendencia                          #  
    ############################################################################
    # Calculo la Global
    ssi_data = data$value
    ssi_ts = ts(ssi_data, start = c(2014, as.numeric(format(data$Fecha[1], "%j"))), frequency = 365)
    decomposed = stl(ssi_ts, s.window = "periodic")
    trend_component = decomposed$time.series[, "trend"]
    result = MannKendall(trend_component)
    
    # gráfico la descomposición
    decomp_df = data.frame(
      time = time(decomposed$time.series),
      seasonal = decomposed$time.series[, "seasonal"],
      trend = decomposed$time.series[, "trend"],
      remainder = decomposed$time.series[, "remainder"]
    ) %>%
      pivot_longer(cols = -time, names_to = "component", values_to = "value")
    
    
    plSTL = ggplot(decomp_df, aes(x = time, y = value)) +
      geom_line() +
      facet_wrap(~component, scales = "free_y", ncol = 1) +
      labs(title = "Descomposición STL", x = "Tiempo", y = "Valor") +
      theme_bw()
    
    # creo carpeta para guardar mis resultados
    save_tendencia = paste0(dir.save, "/Tendencia_SSI")
    if (!dir.exists(save_tendencia)) {
      dir.create(save_tendencia)
    }
    
    # Analizo los puntos de cambio en mi tendencia 
    pettitt_result = pettitt.test(trend_component)
    change_point_date = as.Date(data$Fecha[1]) + (pettitt_result$estimate - 1)
    
    plpc = ggplot() +
      geom_line(aes(x = data$Fecha, y = trend_component), color = "blue", size = 1) +
      geom_vline(xintercept = as.numeric(change_point_date), color = "red", linetype = "dashed", size = 1) +
      theme_bw() +
      labs(title = "Tendencia con punto de cambio", x = "Fecha", y = "") +
      theme(
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 18, hjust = 0.9),
        legend.position = "bottom",
        legend.justification = c(0.5, 0.5),
        legend.background = element_rect(fill = "white", colour = "black", linewidth = 1),
        plot.title = element_text(hjust = 0.5, size = rel(2), lineheight = .9,
                                  face = "bold", colour = "brown"),
      ) +
      scale_x_date(
        date_breaks = "1 year",
        labels = function(x) {
          as.numeric(format(x, "%Y"))
        }
      ) +
      
      geom_rect(aes(xmin = as.Date(-Inf), xmax = as.Date(Inf), ymin = -Inf, ymax = Inf),
                color = "black", fill = NA, size = 1.5) +
      
      # Agregar el texto
      geom_text(aes(x = change_point_date, y = max(trend_component), label = paste("Pc:", format(change_point_date, "%Y"))), 
                color = "black", size = 5, hjust = 0.5, vjust = 0.5, fontface = "bold")
    
    print(plpc)
    
    # genero mi plotly 
    pl_interactiveplpc = ggplotly(plpc)
    
    # genero mi grafico general
    plgeneral = ggplot() +
      geom_line(aes(x = data$Fecha, y = data$value), color = "blue", size = 1) +
      geom_vline(xintercept = as.numeric(change_point_date), color = "red", linetype = "dashed", size = 1) +
      theme_bw() +
      labs(title = "Punto de cambio en la serie", x = "Fecha", y = "SSI") +
      theme(
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 18, hjust = 0.9),
        legend.position = "bottom",
        legend.justification = c(0.5, 0.5),
        legend.background = element_rect(fill = "white", colour = "black", linewidth = 1),
        plot.title = element_text(hjust = 0.5, size = rel(2), lineheight = .9,
                                  face = "bold", colour = "brown"),
      ) +
      scale_x_date(
        date_breaks = "1 year",
        labels = function(x) {
          as.numeric(format(x, "%Y"))
        }
      ) +
      
      # Agregar el texto
      geom_text(aes(x = change_point_date, y = max(data$value), label = paste("Pc:", format(change_point_date, "%Y"))), 
                color = "black", size = 5, hjust = 0.5, vjust = 0.5, fontface = "bold")+
    
      
      geom_rect(aes(xmin = as.Date(-Inf), xmax = as.Date(Inf), ymin = -Inf, ymax = Inf),
                color = "black", fill = NA, size = 1.5) 
    
    plgeneral
    
    # GENERO mi iterativo 
    pl_interactiveplgeneral = ggplotly(plgeneral)
    
    # Guardo mis resultados
    # gráfico de tendencia descompuesta
    ggsave(
      paste(save_tendencia, "/", Subcuenca, "_DescomposicionSTL.png", sep = ""),
      plot = plSTL,
      width = 12,
      height = 8,
      units = "in",
      dpi = 1000,
      bg = NULL,
    )
    
    # GRAFICO DE PUNTO DE CAMBIO 
    ggsave(
      paste(save_tendencia, "/", Subcuenca, "_PuntoCambio.png", sep = ""),
      plot = plpc,
      width = 12,
      height = 8,
      units = "in",
      dpi = 1000,
      bg = NULL,
    )
    
    htmlwidgets::saveWidget(pl_interactiveplpc, paste0(save_tendencia, "/", Subcuenca, "_PuntoCambio.html", sep = ""))
    
    # GRAFICO GENERAL
    ggsave(
      paste(save_tendencia, "/", Subcuenca, "_General+PuntoCambio.png", sep = ""),
      plot = plgeneral,
      width = 12,
      height = 8,
      units = "in",
      dpi = 1000,
      bg = NULL,
    )
    
    htmlwidgets::saveWidget(pl_interactiveplgeneral, paste0(save_tendencia, "/", Subcuenca, "_General+PuntoCambio.html", sep = ""))
    # guardo mis estadísticos de tendencia
    res_global = data.frame(
      Cuenca = Subcuenca,
      tau = as.numeric(result[1]),
      p_value = as.numeric(result[2]),
      Tipo = "Tendencia Global"
    )
    
    resultados_tendencia = rbind(resultados_tendencia, res_global)
    
    # Calculo del punto de cambio
    m_sens = sens.slope(trend_component, conf.level = 0.95)
    pendiente = m_sens$estimates
    cambio_anual = as.numeric(pendiente) * 365.25
    cambio_historico = as.numeric(cambio_anual) * 14
    
    res_pendienteGlobal = data.frame(
      Cuenca = Subcuenca,
      Tipo = "Global",
      pendiente = pendiente,
      cambio_anual = cambio_anual,
      cambio_hist = cambio_historico
      # intercepto = as.numeric(sen_result_global$coefficients["Intercept"])
    )
    resultados_pendiente = rbind(resultados_pendiente, res_pendienteGlobal)
    
    rownames(resultados_tendencia) = NULL
    rownames(resultados_pendiente) = NULL
    
    # guardo en csv
    write.csv(resultados_tendencia, paste0(save_tendencia, "/", Subcuenca, "_Tendencia.csv", sep = ""), row.names = F)
    write.csv(resultados_pendiente, paste0(save_tendencia, "/", Subcuenca, "_Pendiente.csv", sep = ""), row.names = F)
    ############################################################################
    # Gráfico tendencia global
    data_graph = data
    sen_slope_value = m_sens$estimates
    data_graph$Fecha_num = as.numeric(data_graph$Fecha)
    
    # Calculamos el promedio del tiempo
    mean_date = mean(data_graph$Fecha_num)
    
    # Calculamos el intercepto usando el valor promedio del tiempo
    intercept_value = mean(data_graph$value) - sen_slope_value * mean_date
    
    pl_general = ggplot() +
      geom_line(data = data_graph, aes(
        x = Fecha, y = value), color = "blue", size = 1) + 
      
      geom_abline(intercept = intercept_value, slope = sen_slope_value, color = "red", size = 1) +
      
      theme_bw() +
      
      labs(x = "Año", y = "SSI") +
      theme(
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 18, hjust = 0.9),
        legend.position = "bottom",
        legend.justification = c(0.5, 0.5),
        legend.background = element_rect(fill = "white", colour = "black", linewidth = 1),
        plot.title = element_text(hjust = 0.5, size = rel(2), lineheight = .9,
                                  face = "bold", colour = "brown"),
      ) +
      scale_x_date(
        date_breaks = "1 year",
        labels = function(x) {
          as.numeric(format(x, "%Y"))
        }
      ) +
      geom_rect(aes(xmin = as.Date(-Inf), xmax = as.Date(Inf), ymin = -Inf, ymax = Inf),
                color = "black", fill = NA, size = 1.5) +
      
      geom_text(aes(x = as.Date("2020-01-01"), y = max(data_graph$value) * 1.1,  # Ajustar la posición de la etiqueta
                    label = paste("Pendiente:", round(sen_slope_value, 5))), 
                color = "black", size = 5, hjust = 1.1, vjust = 1.5, fontface = "bold")
    
    pl_general
    
    # guardo mi gráfico
    ggsave(
      paste(save_tendencia, "/", Subcuenca, "_cambio_pendiente.png", sep = ""),
      plot = pl_general,
      width = 12,
      height = 8,
      units = "in",
      dpi = 1000,
      bg = NULL,
    )
    
  } # cierro el else de SPEI
  
} # cierro la función de SAD

################################################################################
#                               Categorizar sequias                            #
#' Función que solo gráfica las sequías en base a las sequías agrupadas
graph_CatDroughts = function(data, type, tc = NULL, pc = NULL, Vrecor, umbral = NULL, dir.save) {
  if (type == "SSI") {
    grouping = IAE_droughts(data = data, Type = "SSI", tc = tc, pc = pc, Vrecor = Vrecor, umbral = umbral)
    caracterizar = caracterizar_Min_Max(grouping$Magnitud)
    
    data_graph = cbind(grouping, caracterizar)
    data_graph = data.frame(data_graph)
    data_graph$Fecha_Inicio = as.Date(data_graph$Fecha_Inicio)
    data_graph$Fecha_Fin = as.Date(data_graph$Fecha_Fin)
    data_graph$Año = year(data_graph$Fecha_Inicio)
    data_graph$Mes = month(data_graph$Fecha_Inicio)
    
    conteo_anual = data_graph %>%
      group_by(Año, Mes, Categoria, Duracion,  Magnitud) %>%
      summarise(conteo = n()) %>%
      ungroup()
    
    data_diaria = data_graph %>%
      rowwise() %>%
      mutate(
        Dias = list(seq(Fecha_Inicio, Fecha_Fin, by = "day"))
      ) %>%
      unnest(Dias) %>% # Expandir las fechas diarias
      mutate(
        Año = lubridate::year(Dias),
        Mes = lubridate::month(Dias)
      ) %>%
      mutate(
        # Saco el dia juliano de los dias
        Dias = as.numeric(format(Dias, "%j"))
      )
    
    data_diaria = data_diaria %>%
      mutate(
        Mes = factor(
          lubridate::month(Fecha_Inicio, label = TRUE, abbr = FALSE), 
          levels = month.name
        )
      )
    
    
    colores = c(
      "Muy bajo" = "#94d2bd",
      "Bajo" = "#ffea00",
      "Medio" = "#ff7b00",
      "Alto" = "brown",
      "Muy alto" = "red"
    )
    
    ############################################################################
    data_diaria$Categoria = factor(
      data_diaria$Categoria,
      levels = c("Muy bajo", "Bajo", "Medio", "Alto", "Muy alto")
    )
    
    grafico_diario = ggplot(data_diaria, aes(x = Dias, y = factor(Año), fill = Categoria)) +
      geom_tile(color = "white") +
      scale_fill_manual(
        values = colores,
        name = "Categoría de Sequía",
        guide = guide_legend(
          title.position = "top",
          title.hjust = 0.5,
          ncol = 1
        )
      ) +
      scale_x_continuous(
        breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
        labels = month.abb,
        expand = c(0, 0)
      ) +
      labs(
        title = "Duración de Sequías por Mes y Año",
        x = "Mes",
        y = "Año"
      ) +
      theme(
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 14),
        legend.title = element_text(colour = "brown", face = "bold", size = 12),
        legend.position = "right",
        legend.background = element_rect(fill = "white", colour = "black", linewidth = 1),
        legend.key.size = unit(1, "cm"),
        plot.title = element_text(hjust = 0.5, size = rel(2), lineheight = .9,
                                  face = "bold", colour = "brown"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    print(grafico_diario)
    
    name_folfer1 = paste0(dir.save, "/MapasCalor_SSI")
    
    if (!dir.exists(name_folfer1)) {
      dir.create(name_folfer1)
    }
    
    ggsave(
      paste(name_folfer1, "/", Subcuenca, "_Magnitud.png", sep = ""),
      plot = grafico_diario,
      width = 12,
      height = 8,
      units = "in",
      dpi = 1000,
      bg = NULL,
    )
    
    ############################################################################
  } else {
    grouping = IAE_droughts(data = data, Type = "SPEI", Vrecor = Vrecor, umbral = umbral)
    
    grup_final = list()
    for (i in 1:length(grouping)) {
      data = grouping[[i]]
      caracterizar = caracterizar_Min_Max(data$Magnitud)
      data = cbind(data, caracterizar)
      grup_final[[i]] = data
    }
    
    names(grup_final) = paste0("Pixel_", 1:length(grup_final))
    ############################################################################
    # Genero gráficos ilustrativos
    
    ############################################################################
    
    
  } # cierro el else
}

################################################################################
#   Características meteorológicas desencadenantes en sequía hidrológica       #
# Funcion que calcula los eventos de sequia meteorologica que desencadena en sequias hidrologicas
#' balance_hidrico = archivo con los datos de balance hidrico formato df (solo SPEI)
#' Subcuenca = Nombre de la subcuenca ("Tomebamba o Yanuncay")
#' data_SPEI = archivo con los datos de SPEI formato df (Solo SPEI)
#' data_SSI = archivo con los datos de SSI formato df (Solo SSI)
#' tc y pc son parametros que se usaran para Drought_Grouping (Parametros explicados en 'Drought_Grouping')
#' Vrecor y Umbral = Utilizado para teori_run (Parametros explicados en 'Teori_Run')
#' exclude = Si se desea excluir eventos de sequía (version en investigación, se excluye eventos por abajo del percentil 30)
drought_matching = function(balance_hidrico, Subcuenca, data_SPEI, data_SSI, tc, 
                            pc, Vrecor, umbral, exclude) {
  
  names(balance_hidrico)[1] = "TIMESTAMP"
  if (Subcuenca == "Tomebamba") {
    print("Tomebamba")
    balance_hidrico = balance_hidrico[, .(TIMESTAMP, Pixel_1, Pixel_2, Pixel_3, Pixel_4, Pixel_5, Pixel_6, Pixel_7, Pixel_8, Pixel_9,
                                          Pixel_10, Pixel_11, Pixel_12, Pixel_13, Pixel_14, Pixel_15, Pixel_16, Pixel_17,
                                          Pixel_18), ]
    names(balance_hidrico)[2:length(balance_hidrico)] = paste0("Pixel_", 1:length(balance_hidrico))
    
    data_SPEI = data_SPEI[, .(Fecha, Pixel_1, Pixel_2, Pixel_3, Pixel_4, Pixel_5, Pixel_6, Pixel_7, Pixel_8, Pixel_9,
                              Pixel_10, Pixel_11, Pixel_12, Pixel_13, Pixel_14, Pixel_15, Pixel_16, Pixel_17,
                              Pixel_18), ]
    names(data_SPEI)[2:length(data_SPEI)] = paste0("Pixel_", 1:length(data_SPEI))
    
  } else {
    print("Yanuncay")
    balance_hidrico = balance_hidrico[, .(TIMESTAMP, Pixel_8, Pixel_9, Pixel_14, Pixel_15,Pixel_16,Pixel_17,Pixel_18,Pixel_19,Pixel_22, Pixel_23, 
                                          Pixel_24, Pixel_25, Pixel_26, Pixel_27, Pixel_28,Pixel_31,
                                          Pixel_32, Pixel_33, Pixel_34, Pixel_35, Pixel_36, Pixel_39, Pixel_40, Pixel_41, Pixel_42, 
                                          Pixel_44, Pixel_45, Pixel_46),]
    
    names(balance_hidrico)[2:length(balance_hidrico)] = paste0("Pixel_", 1:length(balance_hidrico))
  }
  
  ##############################################################################
  # Calculo correlacion con tp de 1 a 365 días.
  plan(multisession)
  message("Espere, iniciando paralelización")
  message("---------------------------------=>")
  secuencia_lags = seq(1, 365, by = 1)
  lags = spei_parallel(balance_hidrico, secuencia_lags)
  names(lags) = paste0("lag_", secuencia_lags)
  plan(sequential)
  
  resultados = list()
  for (i in 1:length(lags)) {
    nombre_dia = names(lags)[i] 
    correlaciones_df = data.frame(Pixel = colnames(lags[[i]])[-1])  # Excluye "Fecha"
    
    # Calcular la correlación para cada pixel
    for (j in 2:ncol(lags[[i]])) {  # Desde la segunda columna porque la primera es "Fecha"
      
      datos = lags[[i]][, .(Fecha, Valor = lags[[i]][[j]])]  # Selección dinámica de pixel
      datos = merge(datos, SSI, by = "Fecha", all.x = TRUE)  # Unir con SSI por "Fecha"
      names(datos)[2:3] = c("SPEI", "SSI")  # Cambiar nombres de columnas
      
      # Calcular correlación
      correlacion = cor(datos[,-1], use = "complete.obs", method = "pearson")
      correlacion_valor = correlacion[1, 2]  # Extraer el valor específico de la matriz
      correlacion_valor = round(correlacion_valor, 3)  # Redondear a dos decimale
      correlaciones_df[j - 1, "Correlacion"] = correlacion_valor
    }
    
    correlaciones_df$Dias = nombre_dia
    
    resultados[[nombre_dia]] = correlaciones_df
    
  }
  
  resultados_df = rbindlist(resultados, idcol = "Archivo")
  
  # Extraigo días con mayor correlación (tp)
  dias_mayor_correlacion = resultados_df[, .SD[Correlacion == max(Correlacion)], by = Pixel]
  
  dias_mayor_correlacion[, Dias := as.numeric(gsub("lag_", "", Dias))]
  
  dias_mayor_correlacion = dias_mayor_correlacion[
    , .SD[which.min(Dias)], by = .(Pixel)
  ]
  
  dias_mayor_correlacion$Pixel = factor(dias_mayor_correlacion$Pixel, 
                                        levels = unique(dias_mayor_correlacion$Pixel[order(as.numeric(gsub("Pixel_", "", dias_mayor_correlacion$Pixel)))]))
  dias_mayor_correlacion$Dias = paste0(dias_mayor_correlacion$Dias, "_dias")
  
  ##############################################################################
  #                             Trigger Intervals                              #
  # Llamo a mis sequías hidrológicas emparejadas
  sequias_hidro = IAE_droughts(data_SSI, Type = "SSI", tc = tc, pc = pc, Vrecor = Vrecor, umbral = umbral, exclude = exclude) # modificar eso para emparejar con otros tc y pc
  sequias_meteo = IAE_droughts(data_SPEI, Type = "SPEI", Vrecor = Vrecor, umbral = umbral)
  
  # Tp para cada pixel
  dias_tp = dias_mayor_correlacion
  dias_tp$Pixel = as.character(dias_tp$Pixel)
  dias_tp = data.frame(dias_tp)
  
  sequias_emparejadas = list()
  for (j in 1:length(sequias_meteo)) {
    
    data_meteo = sequias_meteo[[j]]
    Pixel = names(sequias_meteo)[j]
    
    data_meteo = data.table(data_meteo)
    tp = dias_tp[dias_tp$Pixel == Pixel, "Dias"]
    tp = as.numeric(gsub("_dias", "", tp))
    res_pixel = data.frame()
    
    for (i in 2:nrow(sequias_hidro)) {
      
      evento_actual = sequias_hidro[i, ]
      evento_anterior = sequias_hidro[i - 1, ]
      
      tsi = evento_actual$Fecha_Inicio
      tei_1 = evento_anterior$Fecha_Fin
      diferencia =  as.numeric(difftime(tsi, tei_1 - 1, units = "days"))
      
      print(paste0(Pixel, " - ", i, " - ", "Sequia_n", i))
      
      # Propagación de eventos .
      if (diferencia >= tp) {
        inicio = evento_actual$Fecha_Inicio - tp
        fin = evento_actual$Fecha_Fin
        
        meteo = data_meteo[Fecha_Inicio >= inicio & Fecha_Fin <= fin, ]
        
        sequia = meteo[Tipo == "Sequia",]
        
        if (nrow(sequia) == 0) {
          sequia = data.frame(
            Magnitud = 0,
            Duracion = 0
          )
        }
        
        humedo = meteo[Tipo == "Humedo",]
        if (nrow(humedo) == 0) {
          humedo = data.frame(
            Magnitud = 0,
            Duracion = 0
          )
        }
        
        resulados = data.frame(
          Pixel = Pixel,
          Trigger_Inicio = inicio,
          Trigger_Fin = fin,
          
          # sequías hidrológica
          DateStart_Hidrologica = evento_actual$Fecha_Inicio,
          DateEnd_Hidrologica  = evento_actual$Fecha_Fin,
          Duracion_Hidrologica = evento_actual$Duracion,
          Magnitud_Hidrologica = evento_actual$Magnitud,
          
          # sequía meteorológica.
          Magnitud_Meteorologica = sum(sequia$Magnitud) + sum(humedo$Magnitud),
          Duracion_Meteorologica = sum(sequia$Duracion)
        )
        
      } else {
        
        inicio = evento_anterior$Fecha_Fin
        fin = evento_actual$Fecha_Fin
        
        meteo = data_meteo[Fecha_Inicio >= inicio & Fecha_Fin <= fin, ]
        sequia = meteo[Tipo == "Sequia",]
        
        if (nrow(sequia) == 0) {
          sequia = data.frame(
            Magnitud = 0,
            Duracion = 0
          )
        }
        
        humedo = meteo[Tipo == "Humedo",]
        if (nrow(humedo) == 0) {
          humedo = data.frame(
            Magnitud = 0,
            Duracion = 0
          )
        }
        
        resulados = data.frame(
          Pixel = Pixel,
          Trigger_Inicio = evento_anterior$Fecha_Fin,
          Trigger_Fin = evento_actual$Fecha_Inicio,
          
          # sequías hidrológica
          DateStart_Hidrologica = evento_actual$Fecha_Inicio,
          DateEnd_Hidrologica = evento_actual$Fecha_Fin,
          Duracion_Hidrologica = evento_actual$Duracion,
          Magnitud_Hidrologica = evento_actual$Magnitud,
          
          # sequía meteorológica.
          Magnitud_Meteorologica = sum(sequia$Magnitud) - sum(humedo$Magnitud),
          Duracion_Meteorologica = sum(sequia$Duracion)
        )
      }
      res_pixel = rbind(res_pixel, resulados)
    } # cierro la función del i (sequias hidro)
    sequias_emparejadas[[j]] = res_pixel
  } 
  ##############################################################################
  # Categorizacion basado en minimos y maximos\
  seq_finales = list()
  for (i in 1:length(sequias_emparejadas)) {
    data = sequias_emparejadas[[i]]
    data$Severidad_Hidrologica = data$Magnitud_Hidrologica / data$Duracion_Hidrologica
    cat = caracterizar_Min_Max(data$Magnitud_Hidrologica)
    cat = data.frame(Categoria_hidrologica = cat$Categoria)
    data = cbind(data, cat)
    
    data$Severidad_Meteorologica = data$Magnitud_Meteorologica / data$Duracion_Meteorologica
    cat_m = caracterizar_Min_Max(data$Magnitud_Meteorologica)
    cat_m = data.frame(Categoria_meteorologica = cat_m$Categoria)
    data = cbind(data, cat_m)
    
    data[data == Inf | data == -Inf] = NA
    data = na.omit(data)
    
    setDT(data)
    data = data[,.(Pixel = Pixel, Inicio_hidrologica = DateStart_Hidrologica, 
                   Fin_hidrologica = DateEnd_Hidrologica, 
                   Duracion_hidrologica = Duracion_Hidrologica, 
                   Magnitud_hidrologica = Magnitud_Hidrologica, 
                   Severidad_hidrologica = Severidad_Hidrologica,
                   Categoria_hidrologica = Categoria_hidrologica,
                   Duracion_meteorologica = Duracion_Meteorologica,
                   Magnitud_meteorologica = Magnitud_Meteorologica,
                   Severidad_meteorologica = Severidad_Meteorologica,
                   Categoria_meteorologica = Categoria_meteorologica)]
    
    seq_finales[[i]] = data
  }
  ##############################################################################
  return(seq_finales)
} # cierro mi función de emparejamiento. 

################################################################################
#                             Umbrales de propagación (antes llamado copula)                         #
#' seq_emparejadas = df con eventos de sequia meteorologica que desencadenan en sequia hidrologica 
#' dir.save = ruta para guardar archivos

# Versión 2 (Implementando método de Lmoments )
propagation_thresholds = function(seq_emparejadas, umbral) {
  # Ajuste de distribuciones marginales
  best_distribution = function(x, Pixel, nombre, tipo) {
    
    values = abs(as.numeric(x))
    ##############################################################################
    #                            L-moments                                       #
    selection = c("gev", "gum", "gam", "pe3", "ln3", "exp", "wei", "glo")
    
    ajust_1 = distLfit(
      values,
      datname = deparse(substitute(values)),
      selection = selection,
      speed = T,
      ks = T,
      quiet = T,
      order = T
    )
    
    list_vect2par = list()
    for (i in 1:length(selection)){
      list_vect2par[[i]] = vec2par(ajust_1$parameter[[selection[i]]]$para, selection[i])
    }
    
    # Bondas de ajuste de las distribuciones      
    gof = data.frame(
      Distribuciones = rownames(ajust_1$gof),
      RMSE = ajust_1$gof$RMSE,
      W1 = ajust_1$gof$weight1,
      W2 = ajust_1$gof$weight2,
      W3 = ajust_1$gof$weight3,
      ksP = ajust_1$gof$ksP,
      ksD = ajust_1$gof$ksD,
      R2 = ajust_1$gof$R2
    )
    
    gof$RMSE_norm = (max(gof$RMSE, na.rm = T) - gof$RMSE) / (max(gof$RMSE, na.rm = T) - min(gof$RMSE, na.rm = T))
    gof$ksP_norm = (gof$ksP - min(gof$ksP, na.rm = T)) / (max(gof$ksP, na.rm = T) - min(gof$ksP, na.rm = T))
    gof$ksD_norm = (max(gof$ksD, na.rm = T) - gof$ksD) / (max(gof$ksD, na.rm = T) - min(gof$ksD, na.rm = T))
    gof$R2_norm = (gof$R2 - min(gof$R2, na.rm = T)) / (max(gof$R2, na.rm = T) - min(gof$R2, na.rm = T))
    
    # Asignación de pesos
    weights = c(RMSE = 0.4, ksP = 0.3, ksD = 0.1, R2 = 0.2)
    
    # Cálculo del puntaje combinado
    gof$Score = gof$RMSE_norm * weights["RMSE"] +
      gof$ksP_norm * weights["ksP"] +
      gof$ksD_norm * weights["ksD"] +
      gof$R2_norm * weights["R2"]
    
    best_distribution = gof[which.max(gof$Score), ]
    
    cdf_functions = list(
      gev = lmomco::cdfgev,
      gum = lmomco::cdfgum,
      gam = lmomco::cdfgam,
      pe3 = lmomco::cdfpe3,
      ln3 = lmomco::cdfln3,
      exp = lmomco::cdfexp,
      wei = lmomco::cdfwei,
      glo = lmomco::cdfglo
    )
    
    # Buscar el mejor índice basado en la distribución seleccionada
    indice_mejor = which(sapply(list_vect2par, function(x) x$type) == best_distribution$Distribuciones)
    
    # Recuperar los parámetros de la mejor distribución
    params = list_vect2par[[indice_mejor]]
    
    # Seleccionar la función de CDF correspondiente de manera dinámica
    selected_cdf_function = cdf_functions[[best_distribution$Distribuciones]]
    
    # Calcular la CDF de forma eficiente
    CDF = selected_cdf_function(values, params)
    
    res_f = list(
      CDF = CDF,
      valores = values
    )
    ##############################################################################
    return(list(best_ajuste = best_distribution, distribucion = res_f))
  }
  
  distribucion_marginal = function(data, Pixel) {
    setDT(data)
    data = data[,.(H_D = Duracion_hidrologica, M_D = Duracion_meteorologica, 
                   H_S = Magnitud_hidrologica, M_S = Magnitud_meteorologica)]
    
    d_HD = best_distribution(data$H_D, Pixel = Pixel, nombre = "Duracion", tipo = "Hidrologica")
    d_MD = best_distribution(data$M_D, Pixel = Pixel, nombre = "Duracion", tipo = "Meteorologica")
    d_HS = best_distribution(data$H_S, Pixel = Pixel, nombre = "Magnitud", tipo = "Hidrologica")
    d_MS = best_distribution(data$M_S, Pixel = Pixel, nombre = "Magnitud", tipo = "Meteorologica")
    
    ############################################################################
    df_dHD = cbind(data.frame(name = rep("HD", each = nrow(d_HD$best_ajuste))), d_HD$best_ajuste)
    df_dMD = cbind(data.frame(name = rep("MD", each = nrow(d_MD$best_ajuste))), d_MD$best_ajuste)
    df_dHS = cbind(data.frame(name = rep("HS", each = nrow(d_HS$best_ajuste))), d_HS$best_ajuste)
    df_dMS = cbind(data.frame(name = rep("MS", each = nrow(d_MS$best_ajuste))), d_MS$best_ajuste)
    
    
    ta.res = rbind(df_dHD, df_dMD, df_dHS, df_dMS)
    name_pixel = rep(Pixel, nrow(ta.res))
    ta.res = cbind(name_pixel, ta.res) 

    return(list(estadisticos = ta.res, d_HD = d_HD$distribucion, d_MD = d_MD$distribucion, d_HS = d_HS$distribucion, d_MS = d_MS$distribucion))
  }
  
  distribuciones_marginales = list()
  estadisticos = data.frame()
  
  for (i in 1:length(seq_emparejadas)) {
    message(paste0("Iteración ", i, " de ", length(seq_emparejadas)))
    data = seq_emparejadas[[i]]
    Pixel = unique(data$Pixel)
    variable = distribucion_marginal(data, Pixel)
    
    distribuciones_marginales[[i]] = list(Pixel = Pixel,
      d_HD = variable$d_HD, d_MD = variable$d_MD, 
                                          d_HS = variable$d_HS, d_MS = variable$d_MS)
    
    # names[[i]] = Pixel
    
    estadisticos = rbind(estadisticos, variable$estadisticos)
  }

  ##############################################################################
  #                           Distribución conjunta                            #
  # u' = condición (meteorológicas)
  # v' = objetivo (hidrológicas)
  fmlas_copulas = function(data, u, v) {
    # Validación de datos
    if (!u %in% names(data) | !v %in% names(data)) {
      stop("Las columnas especificadas no están en los datos.")
    }
    
    u_vec = data[[u]]$CDF # meteorológicas
    v_vec = data[[v]]$CDF # hidrologicas
    
    if (any(u_vec < 0 | u_vec > 1) | any(v_vec < 0 | v_vec > 1)) {
      stop("Los valores de u y v deben estar en el rango [0, 1].")
    }
    
    # Busco la mejor copula
    fit = BiCopSelect(u_vec, v_vec, familyset = c(1, 2, 3, 4, 5, 6),
                      method = "itau")
    
    estadisticos = BiCopGofTest(u_vec, v_vec, family = fit$family)
    
    CopCDF = BiCopCDF(u_vec, v_vec, obj = fit)
    
    best_cop = data.frame(
      family = fit$familyname,
      tau = fit$tau,
      loglik = fit$logLik,
      AIC = fit$AIC,
      BIC = fit$BIC,
      Ks = estadisticos$statistic,
      p_value = estadisticos$p.value
    )
    
    return(list(
      statics_cop = best_cop,
      prob_conjunta = CopCDF
    ))
  }
  
  ##############################################################################
  dependencia_concicional = function(x_u, y_v, C_xu_yv) {
    numerador = 1 - x_u - y_v + C_xu_yv
    denominador = 1 - x_u
    Probabilidad = numerador / denominador
    
    if (Probabilidad < 0 | Probabilidad > 1) {
      stop("La probabilidad condicional debe estar en el rango [0, 1].")
    }
    
    return(Probabilidad)
  }
  
  dependencia = list()
  mejor_copula = list()
  res_copulas = list()
  
  for (i in 1:length(distribuciones_marginales)) {
    data = distribuciones_marginales[[i]]
    print(paste0("iteracion ", i, " de ", length(distribuciones_marginales)))
    
    # condición 1 # MD - HD
    cond_1 = fmlas_copulas(data = data, u = "d_MD", v = "d_HD")
    cop_1 = cond_1$prob_conjunta
    
    xu1 = distribuciones_marginales[[i]]$d_MD$CDF
    yv1 = distribuciones_marginales[[i]]$d_HD$CDF
    
    u1 = distribuciones_marginales[[i]]$d_MD$valores
    v1 = distribuciones_marginales[[i]]$d_HD$valores
    
    analisis_1 = data.frame(cbind(xu1, yv1, cop_1))
    prob_1 = apply(analisis_1, 1, function(x) dependencia_concicional(x[1], x[2], x[3]))
    
    res_prob1 = data.frame(cbind(
      u_MD = u1,
      Xu_MD = xu1,
      v_HD = v1,
      Yv_HD = yv1,
      Cop_MD_HD = cop_1,
      Prob_MD_HD = prob_1))
    
    
    # condición 2 # MD - HS
    cond_2 = fmlas_copulas(data = data, u = "d_MD", v = "d_HS")
    cop_2 = cond_2$prob_conjunta
    
    xu2 = distribuciones_marginales[[i]]$d_MD$CDF
    yv2 = distribuciones_marginales[[i]]$d_HS$CDF
    
    u2 = distribuciones_marginales[[i]]$d_MD$valores
    v2 = distribuciones_marginales[[i]]$d_HS$valores
    
    analisis_2 = data.frame(cbind(xu2, yv2, cop_2))
    prob_2 = apply(analisis_2, 1, function(x) dependencia_concicional(x[1], x[2], x[3]))
    
    res_prob2 = data.frame(cbind(
      u2_MD = u2,
      Xu2_MD = xu2,
      v2_HS = v2,
      Yv2_HS = yv2,
      Cop_MD_HS = cop_2,
      Prob_MD_HS = prob_2))
    
    # condición 3 # MS - HD
    cond_3 = fmlas_copulas(data = data, u = "d_MS", v = "d_HD")
    cop_3 = cond_3$prob_conjunta
    
    xu3 = distribuciones_marginales[[i]]$d_MS$CDF
    yv3 = distribuciones_marginales[[i]]$d_HD$CDF
    
    u3 = distribuciones_marginales[[i]]$d_MS$valores
    v3 = distribuciones_marginales[[i]]$d_HD$valores
    
    analisis_3 = data.frame(cbind(xu3, yv3, cop_3))
    prob_3 = apply(analisis_3, 1, function(x) dependencia_concicional(x[1], x[2], x[3]))
    
    res_prob3 = data.frame(cbind(
      u3_MS = u3,
      Xu3_MS = xu3,
      v3_HD = v3,
      Yv3_HD = yv3,
      Cop_MS_HD = cop_3,
      Prob_MS_HD = prob_3))
    
    
    # condición 4 # MS - HS
    cond_4 = fmlas_copulas(data = data, u = "d_MS", v = "d_HS")
    cop_4 = cond_4$prob_conjunta
    
    xu4 = distribuciones_marginales[[i]]$d_MS$CDF
    yv4 = distribuciones_marginales[[i]]$d_HS$CDF
    
    u4 = distribuciones_marginales[[i]]$d_MS$valores
    v4 = distribuciones_marginales[[i]]$d_HS$valores
    
    analisis_4 = data.frame(cbind(xu4, yv4, cop_4))
    prob_4 = apply(analisis_4, 1, function(x) dependencia_concicional(x[1], x[2], x[3]))
    
    res_prob4 = data.frame(cbind(
      u4_MS = u4,
      Xu4_MS = xu4,
      v4_HS = v4,
      Yv4_HS = yv4,
      Cop_MS_HS = cop_4,
      Prob_MS_HS = prob_4))
    
    # genero un solo df de resultados
    res = data.frame(cbind(res_prob1, res_prob2, res_prob3, res_prob4)) 
    
    
    dependencia[[i]] = res
    # genero df de las copulas que mejor se ajustan
    res_bestCopula = data.frame(
      MD_HD = cond_1$statics_cop$family,
      MD_HS = cond_2$statics_cop$family,
      MS_HD = cond_3$statics_cop$family,
      MS_HS = cond_4$statics_cop$family
    )
    
    # GENERO DF de los resultados de cada copula]
    ests_1 = cond_1$statics_cop[c("family", "AIC", "BIC", "loglik")]
    ests_2 = cond_2$statics_cop[c("family", "AIC", "BIC", "loglik")]
    ests_3 = cond_3$statics_cop[c("family", "AIC", "BIC", "loglik")]
    ests_4 = cond_4$statics_cop[c("family", "AIC", "BIC", "loglik")]
    
    res_copulas = rbind(ests_1, ests_2, ests_3, ests_4)
    res_copulas = cbind(data.frame(cond = rep(c("MD_HD", "MD_HS", "MS_HD", "MS_HS"), each = 4)), res_copulas)
    
    mejor_copula[[i]] = res_bestCopula
    res_copulas[[i]] = res_copulas
    
  }
  
  
  ##############################################################################
  # Unión de todos los datos para extraer información y clasificar umbrales.....
  clasificacion = function(x) {
    fcase(
      x >= 0 & x < 0.5, "Sequia leve",
      x >= 0.5 & x < 0.75, "Sequia moderada",
      x >= 0.75 & x < 0.9, "Sequia severa",
      x >= 0.9, "Sequia extrema"
    )
  }
  
  umbrales_propagacion = list()
  umbral = umbral
  for (i in 1:length(dependencia)) {
    data = dependencia[[i]]
    setDT(data)
    Pixel = paste0("Pixel_", i)
    # MD - HD
    cond_1 = data[,.(u1_MD = u_MD, Xu_MD = Xu_MD, v1_HD = v_HD, Yv_HD = Yv_HD,
                     Cop_MD_HD = Cop_MD_HD, Prob_MD_HD = Prob_MD_HD)]
    
    cond_1$Categoria_MD_HD = sapply(cond_1$Cop_MD_HD, clasificacion)
    cond_1 = cond_1[Prob_MD_HD >= umbral,]
    cond_1 = cond_1 %>%
      group_by(Categoria_MD_HD) %>%    
      top_n(1, Prob_MD_HD) %>%          
      ungroup()                      
    
    
    # MD - HS
    cond_2 = data[,.(u2_MD = u2_MD, Xu2_MD = Xu2_MD, v2_HS = v2_HS, Yv2_HS = Yv2_HS,
                     Cop_MD_HS = Cop_MD_HS, Prob_MD_HS = Prob_MD_HS)]

    cond_2$Categoria_MD_HS = sapply(cond_2$Cop_MD_HS, clasificacion)
    cond_2 = cond_2[Prob_MD_HS >= umbral,]
    cond_2 = cond_2 %>%
      group_by(Categoria_MD_HS) %>%    
      top_n(1, Prob_MD_HS) %>%          
      ungroup()
    
    # MS - HD
    cond_3 = data[,.(u3_MS = u3_MS, Xu3_MS = Xu3_MS, v3_HD = v3_HD, Yv3_HD = Yv3_HD,
                     Cop_MS_HD = Cop_MS_HD, Prob_MS_HD = Prob_MS_HD)]

    cond_3$Categoria_MS_HD = sapply(cond_3$Cop_MS_HD, clasificacion)
    cond_3 = cond_3[Prob_MS_HD >= umbral,]
    cond_3 = cond_3 %>%
      group_by(Categoria_MS_HD) %>%    
      top_n(1, Prob_MS_HD) %>%          
      ungroup()
    
    # MS - HS
    cond_4 = data[,.(u4_MS = u4_MS, Xu4_MS = Xu4_MS, v4_HS = v4_HS, Yv4_HS = Yv4_HS,
                     Cop_MS_HS = Cop_MS_HS, Prob_MS_HS = Prob_MS_HS)]

    cond_4$Categoria_MS_HS = sapply(cond_4$Cop_MS_HS, clasificacion)
    cond_4 = cond_4[Prob_MS_HS >= umbral,]
    cond_4 = cond_4 %>%
      group_by(Categoria_MS_HS) %>%    
      top_n(1, Prob_MS_HS) %>%          
      ungroup()
    
    umbrales_propagacion[[i]] = list(
      Pixel = Pixel,
      MD_HD = cond_1,
      MD_HS = cond_2,
      MS_HD = cond_3,
      MS_HS = cond_4)
  }
  
  return(umbrales_propagacion)
  
}

save_propThres = function(umbrales_propagacion, dir.save, Subcuenca, cords) {
  leves = data.frame()
  moderadas = data.frame()
  severas = data.frame()
  extremas = data.frame()
  
  validos = function(data) {
    setDT(data)
    data = data[,.(Categoria = names, Pixel = Pixel,
                   u1_MD = u1_MD, v1_HD = v1_HD, Prob_MD_HD = Prob_MD_HD,
                   u2_MD = u2_MD, v2_HS = v2_HS, Prob_MD_HS = Prob_MD_HS,
                   u3_MS = u3_MS, v3_HD = v3_HD, Prob_MS_HD = Prob_MS_HD,
                   u4_MS = u4_MS, v4_HS = v4_HS, Prob_MS_HS = Prob_MS_HS)]
                   

    return(data)
  }

  
  for (i in 1:length(umbrales_propagacion)) {
    print(paste0("Iteración ", i, " de ", length(umbrales_propagacion)))
    Seq_leve = data.frame(names = "sequia leve")
    Seq_moderada = data.frame(names = "sequia moderada")
    Seq_severa = data.frame(names = "sequia severa")
    Seq_extrema = data.frame(names = "sequia extrema")
    
    Pixel = umbrales_propagacion[[i]]$Pixel
    MD_HD = umbrales_propagacion[[i]]$MD_HD
    MD_HS = umbrales_propagacion[[i]]$MD_HS
    MS_HD = umbrales_propagacion[[i]]$MS_HD
    MS_HS = umbrales_propagacion[[i]]$MS_HS
    
    # Sequias leves 
    
    if (nrow(MD_HD[MD_HD$Categoria_MD_HD == "Sequia leve",]) != 0) {
      Seq_leve = cbind(Seq_leve, MD_HD[MD_HD$Categoria_MD_HD == "Sequia leve",])
    } else {
      seq_leve_na <- data.frame(matrix(NA, nrow = 1, ncol = ncol(MD_HD)))
      colnames(seq_leve_na) <- colnames(MD_HD)
      Seq_leve = cbind(Seq_leve, seq_leve_na)
      
    }

    if (nrow(MD_HS[MD_HS$Categoria_MD_HS == "Sequia leve",]) != 0) {
      Seq_leve = cbind(Seq_leve, MD_HS[MD_HS$Categoria_MD_HS == "Sequia leve",])
    } else {
      seq_leve_na <- data.frame(matrix(NA, nrow = 1, ncol = ncol(MD_HS)))
      colnames(seq_leve_na) <- colnames(MD_HS)
      Seq_leve = cbind(Seq_leve, seq_leve_na)
      
    }
    
    if (nrow(MS_HD[MS_HD$Categoria_MS_HD == "Sequia leve",]) != 0) {
      Seq_leve = cbind(Seq_leve, MS_HD[MS_HD$Categoria_MS_HD == "Sequia leve",])
    } else {
      seq_leve_na <- data.frame(matrix(NA, nrow = 1, ncol = ncol(MS_HD)))
      colnames(seq_leve_na) <- colnames(MS_HD)
      Seq_leve = cbind(Seq_leve, seq_leve_na)
      
    }
    
    if (nrow(MS_HS[MS_HS$Categoria_MS_HS == "Sequia leve",]) != 0) {
      Seq_leve = cbind(Seq_leve, MS_HS[MS_HS$Categoria_MS_HS == "Sequia leve",])
    } else {
      seq_leve_na <- data.frame(matrix(NA, nrow = 1, ncol = ncol(MS_HS)))
      colnames(seq_leve_na) <- colnames(MS_HS)
      Seq_leve = cbind(Seq_leve, seq_leve_na)
    }

    # Sequias moderadas
    if (nrow(MD_HD[MD_HD$Categoria_MD_HD == "Sequia moderada",]) != 0) {
      Seq_moderada = cbind(Seq_moderada, MD_HD[MD_HD$Categoria_MD_HD == "Sequia moderada",])
    } else {
      seq_moderada_na <- data.frame(matrix(NA, nrow = 1, ncol = ncol(MD_HD)))
      colnames(seq_moderada_na) <- colnames(MD_HD)
      Seq_moderada = cbind(Seq_moderada, seq_moderada_na)
    }
    
    if (nrow(MD_HS[MD_HS$Categoria_MD_HS == "Sequia moderada",]) != 0) {
      Seq_moderada = cbind(Seq_moderada, MD_HS[MD_HS$Categoria_MD_HS == "Sequia moderada",])
    } else {
      seq_moderada_na <- data.frame(matrix(NA, nrow = 1, ncol = ncol(MD_HS)))
      colnames(seq_moderada_na) <- colnames(MD_HS)
      Seq_moderada = cbind(Seq_moderada, seq_moderada_na)
    }
    
    if (nrow(MS_HD[MS_HD$Categoria_MS_HD == "Sequia moderada",]) != 0) {
      Seq_moderada = cbind(Seq_moderada, MS_HD[MS_HD$Categoria_MS_HD == "Sequia moderada",])
    } else {
      seq_moderada_na <- data.frame(matrix(NA, nrow = 1, ncol = ncol(MS_HD)))
      colnames(seq_moderada_na) <- colnames(MS_HD)
      Seq_moderada = cbind(Seq_moderada, seq_moderada_na)
    }
    
    if (nrow(MS_HS[MS_HS$Categoria_MS_HS == "Sequia moderada",]) != 0) {
      Seq_moderada = cbind(Seq_moderada, MS_HS[MS_HS$Categoria_MS_HS == "Sequia moderada",])
    } else {
      seq_moderada_na <- data.frame(matrix(NA, nrow = 1, ncol = ncol(MS_HS)))
      colnames(seq_moderada_na) <- colnames(MS_HS)
      Seq_moderada = cbind(Seq_moderada, seq_moderada_na)
    }
    
    
    # Sequias severas
    if (nrow(MD_HD[MD_HD$Categoria_MD_HD == "Sequia severa",]) != 0) {
      Seq_severa = cbind(Seq_severa, MD_HD[MD_HD$Categoria_MD_HD == "Sequia severa",])
    } else {
      seq_severa_na <- data.frame(matrix(NA, nrow = 1, ncol = ncol(MD_HD)))
      colnames(seq_severa_na) <- colnames(MD_HD)
      Seq_severa = cbind(Seq_severa, seq_severa_na)
    }
    
    if (nrow(MD_HS[MD_HS$Categoria_MD_HS == "Sequia severa",]) != 0) {
      Seq_severa = cbind(Seq_severa, MD_HS[MD_HS$Categoria_MD_HS == "Sequia severa",])
    } else {
      seq_severa_na <- data.frame(matrix(NA, nrow = 1, ncol = ncol(MD_HS)))
      colnames(seq_severa_na) <- colnames(MD_HS)
      Seq_severa = cbind(Seq_severa, seq_severa_na)
    }
    
    if (nrow(MS_HD[MS_HD$Categoria_MS_HD == "Sequia severa",]) != 0) {
      Seq_severa = cbind(Seq_severa, MS_HD[MS_HD$Categoria_MS_HD == "Sequia severa",])
    } else {
      seq_severa_na <- data.frame(matrix(NA, nrow = 1, ncol = ncol(MS_HD)))
      colnames(seq_severa_na) <- colnames(MS_HD)
      Seq_severa = cbind(Seq_severa, seq_severa_na)
    }
    
    if (nrow(MS_HS[MS_HS$Categoria_MS_HS == "Sequia severa",]) != 0) {
      Seq_severa = cbind(Seq_severa, MS_HS[MS_HS$Categoria_MS_HS == "Sequia severa",])
    } else {
      seq_severa_na <- data.frame(matrix(NA, nrow = 1, ncol = ncol(MS_HS)))
      colnames(seq_severa_na) <- colnames(MS_HS)
      Seq_severa = cbind(Seq_severa, seq_severa_na)
    }
    
    # Sequias extremas
    if (nrow(MD_HD[MD_HD$Categoria_MD_HD == "Sequia extrema",]) != 0) {
      Seq_extrema = cbind(Seq_extrema, MD_HD[MD_HD$Categoria_MD_HD == "Sequia extrema",])
    } else {
      seq_extrema_na <- data.frame(matrix(NA, nrow = 1, ncol = ncol(MD_HD)))
      colnames(seq_extrema_na) <- colnames(MD_HD)
      Seq_extrema = cbind(Seq_extrema, seq_extrema_na)
    }
    
    if (nrow(MD_HS[MD_HS$Categoria_MD_HS == "Sequia extrema",]) != 0) {
      Seq_extrema = cbind(Seq_extrema, MD_HS[MD_HS$Categoria_MD_HS == "Sequia extrema",])
    } else {
      seq_extrema_na <- data.frame(matrix(NA, nrow = 1, ncol = ncol(MD_HS)))
      colnames(seq_extrema_na) <- colnames(MD_HS)
      Seq_extrema = cbind(Seq_extrema, seq_extrema_na)
    }
    
    if (nrow(MS_HD[MS_HD$Categoria_MS_HD == "Sequia extrema",]) != 0) {
      Seq_extrema = cbind(Seq_extrema, MS_HD[MS_HD$Categoria_MS_HD == "Sequia extrema",])
    } else {
      seq_extrema_na <- data.frame(matrix(NA, nrow = 1, ncol = ncol(MS_HD)))
      colnames(seq_extrema_na) <- colnames(MS_HD)
      Seq_extrema = cbind(Seq_extrema, seq_extrema_na)
    }
    
    if (nrow(MS_HS[MS_HS$Categoria_MS_HS == "Sequia extrema",]) != 0) {
      Seq_extrema = cbind(Seq_extrema, MS_HS[MS_HS$Categoria_MS_HS == "Sequia extrema",])
    } else {
      seq_extrema_na <- data.frame(matrix(NA, nrow = 1, ncol = ncol(MS_HS)))
      colnames(seq_extrema_na) <- colnames(MS_HS)
      Seq_extrema = cbind(Seq_extrema, seq_extrema_na)
      
    }
    

    
    
    ## Creo el df final
    Seq_leve$Pixel = Pixel
    leves = rbind(leves, Seq_leve)
    
    Seq_moderada$Pixel = Pixel
    moderadas = rbind(moderadas, Seq_moderada)
    
    Seq_severa$Pixel = Pixel
    severas = rbind(severas, Seq_severa)
    
    Seq_extrema$Pixel = Pixel
    extremas = rbind(extremas, Seq_extrema)
  }
  
  leves = validos(leves)
  moderadas = validos(moderadas)
  severas = validos(severas)
  extremas = validos(extremas)
  
  # convierto a shape
  leves = merge(leves, cords, by = "Pixel")
  moderadas = merge(moderadas, cords, by = "Pixel")
  severas = merge(severas, cords, by = "Pixel")
  extremas = merge(extremas, cords, by = "Pixel")
  
  pnts_leves = vect(leves, geom = c("X", "Y"), crs = "EPSG:32717")
  pnts_moderadas = vect(moderadas, geom = c("X", "Y"), crs = "EPSG:32717")
  pnts_severas = vect(severas, geom = c("X", "Y"), crs = "EPSG:32717")
  pnts_extremas = vect(extremas, geom = c("X", "Y"), crs = "EPSG:32717")
  
  folder = paste0(dir.save, "/Umbrales_Propagacion")
  if (!dir.exists(folder)) {
    dir.create(folder)
  }
  
  writeVector(pnts_leves, paste0(folder, "/Leves_", Subcuenca, ".shp", sep = ""), filetype = "ESRI Shapefile", overwrite = TRUE)
  writeVector(pnts_moderadas, paste0(folder, "/Moderadas_", Subcuenca, ".shp", sep = ""), filetype = "ESRI Shapefile", overwrite = TRUE)
  writeVector(pnts_severas, paste0(folder, "/Severas_", Subcuenca, ".shp", sep = ""), filetype = "ESRI Shapefile", overwrite = TRUE)
  writeVector(pnts_extremas, paste0(folder, "/Extremas_", Subcuenca, ".shp", sep = ""), filetype = "ESRI Shapefile", overwrite = TRUE)
  
  
  
}
################################################################################
#                             Umbrales dado condiciones                        #
condition_thresolds = function() {
  
}


################################################################################
################################################################################
#             Funciones Beta (Podrian no incorporarse al analisis final)       #
################################################################################
################################################################################
copulas = function(seq_emparejadas, dir.save) {
  
  ##############################################################################
  #                           Distribución Marginal                            #
  
  distribuciones = function(x, Pixel, nombre, tipo) {
    
    # Preparo los datos para el ajuste
    values = abs(as.numeric(x))
    
    ###################### Buscar la mejor distribución ########################
    # exponencial
    fit_exp = fitdist(values, "exp", method = "mle")
    AIC_exp = fit_exp$aic
    BIC_exp = fit_exp$bic
    CV_exp = fit_exp$loglik
    
    # weibull 
    fit_weibull = fitdist(values, "weibull")
    AIC_weibull = fit_weibull$aic
    BIC_weibull = fit_weibull$bic
    CV_weibull = fit_weibull$loglik
    
    # GEV
    fit_gev = fevd(values, type = "GEV")
    AIC_gev = summary(fit_gev)$AIC
    BIC_gev = summary(fit_gev)$BIC
    CV_gev = -1 * fit_gev$results$value
    
    ############################################################################
    # Funcion para calcular el thresold de pareto
    # Función para seleccionar el threshold óptimo
    select_optimal_threshold = function(data, umin, umax, show_plot = TRUE) {
      # Ajustar el modelo GPD para un rango de thresholds
      fit_results = gpd.fitrange(data = data, umin = umin, umax = umax, show = FALSE)
      
      # Extraer parámetros relevantes
      thresholds = fit_results$thresholds
      xi_values = fit_results$mle[, 2]  # Estimaciones de xi
      ci_low = fit_results$ci.low[, 2]
      ci_up = fit_results$ci.up[, 2]
      
      # Calcular estabilidad de xi
      xi_diff = abs(diff(xi_values))  # Diferencias entre valores consecutivos
      ci_width = ci_up - ci_low       # Ancho de intervalos de confianza
      
      # Definir criterio: estabilidad de xi y confianza estrecha
      stability_scores = xi_diff / ci_width[-1]  # Ignorar el primer elemento de ci_width
      optimal_index = which.min(stability_scores)  # Índice del mejor threshold
      
      # Threshold óptimo
      optimal_threshold = thresholds[optimal_index]
      
      # Mostrar resultados
      if (show_plot) {
        library(ggplot2)
        results_df = data.frame(
          threshold = thresholds,
          xi = xi_values,
          ci.low = ci_low,
          ci.up = ci_up
        )
        
        ggplot(results_df, aes(x = threshold, y = xi)) +
          geom_line() +
          geom_point() +
          geom_errorbar(aes(ymin = ci.low, ymax = ci.up), width = 0.5) +
          geom_vline(xintercept = optimal_threshold, linetype = "dashed", color = "red") +
          labs(
            title = "Selección del Threshold Óptimo",
            x = "Threshold",
            y = "Parámetro ξ"
          ) +
          theme_minimal()
      }
      
      return(optimal_threshold)
    }
    
    umin = quantile(values, 0.1)
    umax = quantile(values, 0.95)
    sequias = select_optimal_threshold(values, umin, umax)
    # Generalizado de pareto
    
    fit_gpd = fevd(values, type = "GP", threshold = sequias)
    AIC_gpd = summary(fit_gpd)$AIC
    BIC_gpd = summary(fit_gpd)$BIC
    CV_gpd = -1 * fit_gpd$results$value
    
    ############################################################################
    #                   Ajuste de cada distribución a los datos
    ############################################################################
    CDF_exp = pexp(values, rate = fit_exp$estimate)
    CDF_weibull = pweibull(values, shape = fit_weibull$estimate[1], scale = fit_weibull$estimate[2])
    
    CDF_gev = pgev(values, loc = fit_gev$results$par[[1]], scale = fit_gev$results$par[[2]], shape = fit_gev$results$par[[3]])
    if (length(fit_gpd$results$par) > 2) {
      loc_gpd = fit_gpd$results$par[[1]]
      scale_gpd = fit_gpd$results$par[[2]]
      shape_gpd = fit_gpd$results$par[[3]]
    } else {
      loc_gpd = 0
      scale_gpd = fit_gpd$results$par[[1]]
      shape_gpd = fit_gpd$results$par[[2]]
    }
    
    CDF_gpd = pgpd(values, loc = loc_gpd, scale = scale_gpd, shape = shape_gpd)
    
    ############################################################################
    #                                    KS-TEST                               #
    ks_exp = ks.test(values, "pexp", rate = fit_exp$estimate)
    ks_weibull = ks.test(values, "pweibull", shape = fit_weibull$estimate[1], scale = fit_weibull$estimate[2])
    ks_gev = ks.test(values, "pgev", loc = fit_gev$results$par[[1]], scale = fit_gev$results$par[[2]], shape = fit_gev$results$par[[3]])
    ks_gpd = ks.test(values, "pgpd", loc = loc_gpd, scale = scale_gpd, shape = shape_gpd)
    
    ############################################################################
    #                                    QQ-PLOT                                #
    empiric_exp = ecdf(values)(values)
    empiric_weibull = ecdf(values)(values)
    empiric_gev = ecdf(values)(values)
    empiric_gpd = ecdf(values)(values)
    
    
    grafiquito = function(teorico, empirico, nombre, ks) {
      
      pl = ggplot() +
        geom_point(aes(x = empirico, y = teorico), size = 1) +
        geom_abline(intercept = 0, slope = 1, color = "red", size = 1) +
        labs(title = nombre,,
             x = "Probabilidad empírica",
             y = "Probabilidad teórica") +
        theme_bw() + 
        theme(
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 18, hjust = 0.9),
          legend.title = element_text(colour = "brown", face = "bold", size = 12),
          legend.position = "bottom",
          legend.justification = c(0.5, 0.5),
          legend.background = element_rect(fill = "white", colour = "black", linewidth = 1),
          plot.title = element_text(hjust = 0.5, size = rel(2), lineheight = .9,
                                    face = "bold", colour = "brown"),
        ) +
        
        annotate("text", x = 0.15, y = 1, label = paste("p: ", ks), size = 5, color = "black") +
        
        scale_x_continuous(breaks = seq(0, 1, by = 0.25)) +
        
        geom_rect(aes(xmin = (-Inf), xmax = (Inf), ymin = -Inf, ymax = Inf),
                  color = "black", fill = NA, size = 1.5)
      
      return(pl)
    }
    
    
    graph_exp = grafiquito(CDF_exp, empiric_exp, "Exponencial", round(ks_exp$p.value,3))
    graph_weibull = grafiquito(CDF_weibull, empiric_weibull, "Weibull", round(ks_weibull$p.value,3))
    graph_gev = grafiquito(CDF_gev, empiric_gev, "GEV", round(ks_gev$p.value,3))
    graph_gpd = grafiquito(CDF_gpd, empiric_gpd, "GPD", round(ks_gpd$p.value,3))
    
    
    combined_plot = grid.arrange(
      grobs = list(graph_exp, graph_weibull, graph_gev, graph_gpd),
      ncol = 2,  # Número de columnas
      nrow = 2  # Número de filas
    )
    
    folder_save = paste0(dir.save, "/Graficos_pp/")
    if (!dir.exists(folder_save)) {
      dir.create(folder_save)
    }
    
    name = paste(Pixel, nombre, tipo, sep = "_")
    ggsave(
      paste(folder_save, name, ".png", sep = ""),
      plot = combined_plot,
      width = 12,
      height = 8,
      units = "in",
      dpi = 1000,
      bg = NULL,
    )
    
    ############################################################################
    #                                   best distribution                      #
    ############################################################################
    
    best_disrt = data.frame(
      Distribucion = c("Exponencial", "Weibull", "GEV", "GPD"),
      AIC = c(AIC_exp, AIC_weibull, AIC_gev, AIC_gpd),
      BIC = c(BIC_exp, BIC_weibull, BIC_gev, BIC_gpd),
      Max_veros = c(CV_exp, CV_weibull, CV_gev, CV_gpd),
      Ks = c(ks_exp$statistic, ks_weibull$statistic, ks_gev$statistic, ks_gpd$statistic),
      p_Ks = c(ks_exp$p.value, ks_weibull$p.value, ks_gev$p.value, ks_gpd$p.value)
    )
    
    
    # best_disrt = data.frame(
    #   Distribucion = c("Exponencial", "Weibull", "GEV"),
    #   AIC = c(AIC_exp, AIC_weibull, AIC_gev),
    #   BIC = c(BIC_exp, BIC_weibull, BIC_gev),
    #   Max_veros = c(CV_exp, CV_weibull, CV_gev),
    #   Ks = c(ks_exp$statistic, ks_weibull$statistic, ks_gev$statistic),
    #   p_Ks = c(ks_exp$p.value, ks_weibull$p.value, ks_gev$p.value)
    # )
    
    ############################################################################
    #                             Elijo mejor distribucion                     #
    
    AIC_1 = which.min(best_disrt$AIC)
    BIC_1 = which.min(best_disrt$BIC)
    CV_1 =  which.max(best_disrt$Max_veros)
    KS_1 = which.min(best_disrt$Ks)
    p_KS_1 = which.max(best_disrt$p_Ks)
    
    res = data.frame(
      AIC = best_disrt[AIC_1, "Distribucion"],
      BIC = best_disrt[BIC_1, "Distribucion"],
      CV = best_disrt[CV_1, "Distribucion"],
      KS = best_disrt[KS_1, "Distribucion"],
      p_KS = best_disrt[p_KS_1, "Distribucion"]
    )
    
    res_vector = unlist(res)
    frecuencias = table(res_vector)
    mas_frecuente = names(frecuencias[frecuencias == max(frecuencias)])
    
    if (mas_frecuente == "Exponencial") {
      dis_f = CDF_exp
    } else if (mas_frecuente == "Weibull") {
      dis_f = CDF_weibull
    } else if (mas_frecuente == "GEV") {
      dis_f = CDF_gev
    } else if (mas_frecuente == "GPD") {
      dis_f = CDF_gev
    } else {
      stop("No se encontro la distribución")
    }
    
    return(list(best_ajuste = best_disrt, distribucion = dis_f))
  } # cierro función de la distribución 
  
  distribucion_marginal = function(data, Pixel) {
    setDT(data)
    data = data[,.(H_D = Duracion_Hidrologica, M_D = Duracion_Meteorologica, H_S = Magnitud_Hidrologica, M_S = Magnitud_Meteorologica)]
    
    d_HD = distribuciones(data$H_D, Pixel = Pixel, nombre = "Duracion", tipo = "Hidrologica")
    d_MD = distribuciones(data$M_D, Pixel = Pixel, nombre = "Duracion", tipo = "Meteorologica")
    d_HS = distribuciones(data$H_S, Pixel = Pixel, nombre = "Magnitud", tipo = "Hidrologica")
    d_MS = distribuciones(data$M_S, Pixel = Pixel, nombre = "Magnitud", tipo = "Meteorologica")
    
    
    ta.res = rbind(d_HD$best_ajuste, d_MD$best_ajuste, d_HS$best_ajuste, d_MS$best_ajuste)
    name_pixel = rep(Pixel, nrow(ta.res))
    name = rep(c("HD", "MD", "HS", "MS"), each = 4)
    
    
    ta.res = cbind(name_pixel, name, ta.res) 
    
    return(list(estadisticos = ta.res, d_HD = d_HD$distribucion, d_MD = d_MD$distribucion, d_HS = d_HS$distribucion, d_MS = d_MS$distribucion))
    
  } # Cierro la función de distribución marginal
  
  
  distribuciones_marginales = list()
  estadisticos = data.frame()
  
  for (i in 1:length(seq_emparejadas)) {
    message(paste0("Iteración", i, "de", length(seq_emparejadas)))
    data = seq_emparejadas[[i]]
    Pixel = unique(data$Pixel)
    variable = distribucion_marginal(data, Pixel)
    
    distribuciones_marginales[[i]] = list(d_HD = variable$d_HD, d_MD = variable$d_MD, 
                                          d_HS = variable$d_HS, d_MS = variable$d_MS)
    
    estadisticos = rbind(estadisticos, variable$estadisticos)
  }
  
  ##############################################################################
  #                           Distribución conjunta                            #
  # u' = condición (meteorológicas)
  # v' = objetivo (hidrológicas)
  fmlas_copulas = function(data, u, v) {
    # Validación de datos
    if (!u %in% names(data) | !v %in% names(data)) {
      stop("Las columnas especificadas no están en los datos.")
    }
    
    u_vec = data[[u]] # meteorológicas
    v_vec = data[[v]] # hidrologicas
    
    if (any(u_vec < 0 | u_vec > 1) | any(v_vec < 0 | v_vec > 1)) {
      stop("Los valores de u y v deben estar en el rango [0, 1].")
    }
    
    # Crear las cópulas
    copulas = list(
      Gumbel = gumbelCopula(dim = 2),
      Franck = frankCopula(dim = 2),
      Clayton = claytonCopula(dim = 2),
      Gaussiana = normalCopula(dim = 2)
    )
    
    # Ajustar las cópulas y calcular métricas
    resultados = lapply(copulas, function(copula) {
      tryCatch({
        ajuste = fitCopula(copula, data = cbind(u_vec, v_vec), method = "ml")
        list(
          loglik = logLik(ajuste),
          AIC = AIC(ajuste),
          BIC = BIC(ajuste),
          coef = coef(ajuste)
        )
      }, error = function(e) {
        list(loglik = NA, AIC = NA, BIC = NA, coef = NA)
      })
    })
    
    # Crear el dataframe de resultados
    resultados_df = do.call(rbind, lapply(names(resultados), function(nombre) {
      c(copula = nombre, unlist(resultados[[nombre]]))
    }))
    
    resultados_df = as.data.frame(resultados_df, stringsAsFactors = FALSE)
    resultados_df$loglik = as.numeric(resultados_df$loglik)
    resultados_df$AIC = as.numeric(resultados_df$AIC)
    resultados_df$BIC = as.numeric(resultados_df$BIC)
    
    # Selección de la mejor cópula
    mejor_AIC = resultados_df$copula[which.min(resultados_df$AIC)]
    mejor_BIC = resultados_df$copula[which.min(resultados_df$BIC)]
    mejor_loglik = resultados_df$copula[which.max(resultados_df$loglik)]
    
    # Resolver empates
    criterios = table(c(mejor_AIC, mejor_BIC, mejor_loglik))
    mas_frecuente = names(criterios[criterios == max(criterios)])
    mejor_copula = mas_frecuente[1]  # En caso de empate, toma la primera
    
    # Calcular la probabilidad conjunta
    theta = resultados[[mejor_copula]]$coef
    copula_final = switch(mejor_copula,
                          Gumbel = gumbelCopula(theta, dim = 2),
                          Franck = frankCopula(theta, dim = 2),
                          Clayton = claytonCopula(theta, dim = 2),
                          Gaussiana = normalCopula(theta, dim = 2))
    
    C_uv = pCopula(cbind(u_vec, v_vec), copula_final)
    
    # Retornar resultados
    return(list(
      resultados = resultados_df,
      mejor_copula = mejor_copula,
      prob_conjunta = C_uv
    ))
  }
  
  dependencia_concicional = function(x_u, y_v, C_xu_yv) {
    numerador = 1 - x_u - y_v + C_xu_yv
    denominador = 1 - x_u
    Probabilidad = numerador / denominador
    
    if (Probabilidad < 0 | Probabilidad > 1) {
      stop("La probabilidad condicional debe estar en el rango [0, 1].")
    }
    return(Probabilidad)
  }
  
  dependencia = list()
  mejor_copula = list()
  res_copulas = list()
  
  for (i in 1:length(distribuciones_marginales)) {
    data = distribuciones_marginales[[i]]
    print(paste0("iteracion ", i, " de ", length(distribuciones_marginales)))
    # condición 1 # MD - HD
    cond_1 = fmlas_copulas(data = data, u = "d_MD", v = "d_HD")
    cop_1 = cond_1$prob_conjunta
    xu1 = distribuciones_marginales[[i]]$d_MD
    yv1 = distribuciones_marginales[[i]]$d_HD
    analisis_1 = data.frame(cbind(xu1, yv1, cop_1))
    prob_1 = apply(analisis_1, 1, function(x) dependencia_concicional(x[1], x[2], x[3]))
    
    res_prob1 = data.frame(cbind(
      Xu_MD = xu1,
      Yv_HD = yv1,
      Cop_MD_HD = cop_1,
      Prob_MD_HD = prob_1))
    
    
    # condición 2 # MD - HS
    cond_2 = fmlas_copulas(data = data, u = "d_MD", v = "d_HS")
    cop_2 = cond_2$prob_conjunta
    xu2 = distribuciones_marginales[[i]]$d_MD
    yv2 = distribuciones_marginales[[i]]$d_HS
    analisis_2 = data.frame(cbind(xu2, yv2, cop_2))
    prob_2 = apply(analisis_2, 1, function(x) dependencia_concicional(x[1], x[2], x[3]))
    
    res_prob2 = data.frame(cbind(
      Xu2_MD = xu2,
      Yv2_HS = yv2,
      Cop_MD_HS = cop_2,
      Prob_MD_HS = prob_2))
    
    # condición 3 # MS - HD
    cond_3 = fmlas_copulas(data = data, u = "d_MS", v = "d_HD")
    cop_3 = cond_3$prob_conjunta
    xu3 = distribuciones_marginales[[i]]$d_MS
    yv3 = distribuciones_marginales[[i]]$d_HD
    analisis_3 = data.frame(cbind(xu3, yv3, cop_3))
    prob_3 = apply(analisis_3, 1, function(x) dependencia_concicional(x[1], x[2], x[3]))
    
    res_prob3 = data.frame(cbind(
      Xu3_MS = xu3,
      Yv3_HD = yv3,
      Cop_MS_HD = cop_3,
      Prob_MS_HD = prob_3))
    
    
    # condición 4 # MS - HS
    cond_4 = fmlas_copulas(data = data, u = "d_MS", v = "d_HS")
    cop_4 = cond_4$prob_conjunta
    xu4 = distribuciones_marginales[[i]]$d_MS
    yv4 = distribuciones_marginales[[i]]$d_HS
    analisis_4 = data.frame(cbind(xu4, yv4, cop_4))
    prob_4 = apply(analisis_4, 1, function(x) dependencia_concicional(x[1], x[2], x[3]))
    
    res_prob4 = data.frame(cbind(
      Xu4_MS = xu4,
      Yv4_HS = yv4,
      Cop_MS_HS = cop_4,
      Prob_MS_HS = prob_4))
    
    # genero un solo df de resultados
    res = data.frame(cbind(res_prob1, res_prob2, res_prob3, res_prob4))
    dependencia[[i]] = res
    
    # genero df de las copulas que mejor se ajustan
    res_bestCopula = data.frame(
      MD_HD = cond_1$mejor_copula,
      MD_HS = cond_2$mejor_copula,
      MS_HD = cond_3$mejor_copula,
      MS_HS = cond_4$mejor_copula
    )
    
    # GENERO DF de los resultados de cada copula]
    ests_1 = cond_1$resultados[c("copula", "AIC", "BIC", "loglik")]
    ests_2 = cond_2$resultados[c("copula", "AIC", "BIC", "loglik")]
    ests_3 = cond_3$resultados[c("copula", "AIC", "BIC", "loglik")]
    ests_4 = cond_4$resultados[c("copula", "AIC", "BIC", "loglik")]
    
    columnas_completas = unique(c(names(ests_1), 
                                  names(ests_2), 
                                  names(ests_3), 
                                  names(ests_4)))
    
    alinear_columnas = function(df, columnas) {
      faltantes = setdiff(columnas, names(df)) # Identificar columnas faltantes
      for (col in faltantes) {
        df[[col]] = NA # Agregar columnas faltantes con NA
      }
      df = df[, columnas] # Reordenar columnas en el orden correcto
      return(df)
    }
    
    # Alinear todos los data frames
    df1 = alinear_columnas(cond_1$resultados, columnas_completas)
    df2 = alinear_columnas(cond_2$resultados, columnas_completas)
    df3 = alinear_columnas(cond_3$resultados, columnas_completas)
    df4 = alinear_columnas(cond_4$resultados, columnas_completas)
    
    res_copulas = rbind(df1, df2, df3, df4)
    
    #res_copulas = rbind(cond_1$resultados, cond_2$resultados, cond_3$resultados, cond_4$resultados)
    res_copulas = cbind(data.frame(cond = rep(c("MD_HD", "MD_HS", "MS_HD", "MS_HS"), each = 4)), res_copulas)
    
    mejor_copula[[i]] = res_bestCopula
    res_copulas[[i]] = res_copulas
    
  }
  
  ##############################################################################
  # Unión de todos los datos para extraer información y clasificar umbrales.....
  
  clasificacion = function(x) {
    fcase(
      x >= 0 & x < 0.5, "Sequia leve",
      x >= 0.5 & x < 0.75, "Sequia moderada",
      x >= 0.75 & x < 0.9, "Sequia severa",
      x >= 0.9, "Sequia extrema"
    )
  }
  
  umbrales_propagacion = list()
  for (i in 1:length(dependencia)) {
    data = dependencia[[i]]
    setDT(data)
    
    # MD - HD
    cond_1 = data[,.(Xu_MD = Xu_MD, Yv_HD = Yv_HD, Cop_MD_HD = Cop_MD_HD, Prob_MD_HD = Prob_MD_HD)]
    cond_1$Categoria_MD_HD = sapply(cond_1$Cop_MD_HD, clasificacion)
    cond_1 = cond_1[Prob_MD_HD >= 0.95,]
    
    # MD - HS
    cond_2 = data[,.(Xu2_MD = Xu2_MD, Yv2_HS = Yv2_HS, Cop_MD_HS = Cop_MD_HS, Prob_MD_HS = Prob_MD_HS)]
    cond_2$Categoria_MD_HS = sapply(cond_2$Cop_MD_HS, clasificacion)
    cond_2 = cond_2[Prob_MD_HS >= 0.95,]
    
    # MS - HD
    cond_3 = data[,.(Xu3_MS = Xu3_MS, Yv3_HD = Yv3_HD, Cop_MS_HD = Cop_MS_HD, Prob_MS_HD = Prob_MS_HD)]
    cond_3$Categoria_MS_HD = sapply(cond_3$Cop_MS_HD, clasificacion)
    cond_3 = cond_3[Prob_MS_HD >= 0.95,]
    
    # MS - HS
    cond_4 = data[,.(Xu4_MS = Xu4_MS, Yv4_HS = Yv4_HS, Cop_MS_HS = Cop_MS_HS, Prob_MS_HS = Prob_MS_HS)]
    cond_4$Categoria_MS_HS = sapply(cond_4$Cop_MS_HS, clasificacion)
    cond_4 = cond_4[Prob_MS_HS >= 0.95,]
    
    umbrales_propagacion[[i]] = list(
      MD_HD = cond_1,
      MD_HS = cond_2,
      MS_HD = cond_3,
      MS_HS = cond_4)
  }
  
  ##############################################################################
  return(umbrales_propagacion)
  
} # Cierro la función copula

categorize_droughts = function(data, Type, tc = NULL, pc = NULL, Subcuenca, dir.save = NULL, 
                               save = NULL, cords = NULL, Raster_Base = NULL) {
  
  if (Type == "SSI") {
    run = teori_run(data)
    grouping = drought_grouping(res_corridas = run, tc = tc, pc = pc, graficar = T, 
                                data_original = data, dir.save = dir.save, Subcuenca = Subcuenca)
    
    grouping$Severidad = grouping$Magnitud / grouping$Duracion
    
    
    SSI_tibble = as_tibble(grouping)
    SSI_tibble = SSI_tibble %>%
      mutate(Categoria = case_when(
        Magnitud >= 2.0 ~ "No Sequia",
        Magnitud >= 1.5 ~ "No Sequia",
        Magnitud >= 1.0 ~ "No Sequia",
        Magnitud >= 0.0 ~ "No Sequia",
        Magnitud >= -1.0 ~ "Sequía leve",
        Magnitud >= -1.5 ~ "Sequía Moderada",
        Magnitud >= -2.0 ~ "Sequía severa",
        TRUE ~ "Sequía extrema"
      ))
    
    grouping = as.data.table(SSI_tibble)
    ############################################################################
    #                   Genero gráficos (Base Magnitud)                        #
    data_graph = grouping
    data_graph = data.frame(data_graph)
    data_graph$Fecha_Inicio = as.Date(data_graph$Fecha_Inicio)
    data_graph$Fecha_Fin = as.Date(data_graph$Fecha_Fin)
    data_graph$Año = year(data_graph$Fecha_Inicio)
    data_graph$Mes = month(data_graph$Fecha_Inicio)
    
    conteo_anual = data_graph %>%
      group_by(Año, Mes, Categoria, Duracion,  Magnitud) %>%
      summarise(conteo = n()) %>%
      ungroup()
    
    data_diaria = data_graph %>%
      rowwise() %>%
      mutate(
        Dias = list(seq(Fecha_Inicio, Fecha_Fin, by = "day"))
      ) %>%
      unnest(Dias) %>% # Expandir las fechas diarias
      mutate(
        Año = lubridate::year(Dias),
        Mes = lubridate::month(Dias)
      ) %>%
      mutate(
        # Saco el dia juliano de los dias
        Dias = as.numeric(format(Dias, "%j"))
      )
    
    data_diaria = data_diaria %>%
      mutate(
        Mes = factor(
          lubridate::month(Fecha_Inicio, label = TRUE, abbr = FALSE), 
          levels = month.name
        )
      )
    
    colores = c(
      "Sequía leve" = "#94d2bd",
      "Sequía Moderada" = "#ffea00",
      "Sequía severa" = "#ff7b00",
      "Sequía extrema" = "#ad2831"
    )
    
    ############################################################################
    data_diaria$Categoria = factor(
      data_diaria$Categoria,
      levels = c("Sequía leve", "Sequía Moderada", "Sequía severa", "Sequía extrema")
    )
    
    ############################################################################
    grafico_diario = ggplot(data_diaria, aes(x = Dias, y = factor(Año), fill = Categoria)) +
      geom_tile(color = "white") +
      scale_fill_manual(
        values = colores,
        name = "Categoría de Sequía",
        guide = guide_legend(
          title.position = "top",
          title.hjust = 0.5,
          ncol = 1
        )
      ) +
      scale_x_continuous(
        breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
        labels = month.abb,
        expand = c(0, 0)
      ) +
      labs(
        title = "Duración de Sequías por Mes y Año",
        x = "Mes",
        y = "Año"
      ) +
      theme(
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 14),
        legend.title = element_text(colour = "brown", face = "bold", size = 12),
        legend.position = "right",
        legend.background = element_rect(fill = "white", colour = "black", linewidth = 1),
        legend.key.size = unit(1, "cm"),
        plot.title = element_text(hjust = 0.5, size = rel(2), lineheight = .9,
                                  face = "bold", colour = "brown"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    print(grafico_diario)
    ############################################################################
    # Severidad
    SSI_tibble = as_tibble(grouping)
    SSI_tibble = SSI_tibble %>%
      mutate(Categoria = case_when(
        Severidad >= 2.0 ~ "No Sequia",
        Severidad >= 1.5 ~ "No Sequia",
        Severidad >= 1.0 ~ "No Sequia",
        Severidad >= 0.0 ~ "No Sequia",
        Severidad >= -1.0 ~ "Sequía leve",
        Severidad >= -1.5 ~ "Sequía Moderada",
        Severidad >= -2.0 ~ "Sequía severa",
        TRUE ~ "Sequía extrema"
      ))
    grouping = as.data.table(SSI_tibble)
    ############################################################################
    #                             Genero gráficos                              #
    data_graph = grouping
    data_graph = data.frame(data_graph)
    data_graph$Fecha_Inicio = as.Date(data_graph$Fecha_Inicio)
    data_graph$Fecha_Fin = as.Date(data_graph$Fecha_Fin)
    data_graph$Año = year(data_graph$Fecha_Inicio)
    data_graph$Mes = month(data_graph$Fecha_Inicio)
    
    conteo_anual = data_graph %>%
      group_by(Año, Mes, Categoria, Duracion, Severidad) %>%
      summarise(conteo = n()) %>%
      ungroup()
    
    data_diaria = data_graph %>%
      rowwise() %>%
      mutate(
        Dias = list(seq(Fecha_Inicio, Fecha_Fin, by = "day"))
      ) %>%
      unnest(Dias) %>% # Expandir las fechas diarias
      mutate(
        Año = lubridate::year(Dias),
        Mes = lubridate::month(Dias)
      ) %>%
      mutate(
        # Saco el dia juliano de los dias
        Dias = as.numeric(format(Dias, "%j"))
      )
    
    data_diaria = data_diaria %>%
      mutate(
        Mes = factor(
          lubridate::month(Fecha_Inicio, label = TRUE, abbr = FALSE), 
          levels = month.name
        )
      )
    
    colores = c(
      "Sequía leve" = "#94d2bd",
      "Sequía Moderada" = "#ffea00",
      "Sequía severa" = "#ff7b00",
      "Sequía extrema" = "#ad2831"
    )
    
    data_diaria$Categoria = factor(
      data_diaria$Categoria,
      levels = c("Sequía leve", "Sequía Moderada", "Sequía severa", "Sequía extrema")
    )
    
    
    grafico_severidad = ggplot(data_diaria, aes(x = Dias, y = factor(Año), fill = Categoria)) +
      geom_tile(color = "white") +
      scale_fill_manual(
        values = colores,
        name = "Categoría de Sequía",
        guide = guide_legend(
          title.position = "top",
          title.hjust = 0.5,
          ncol = 1
        )
      ) +
      # Aquí está el cambio principal: modificamos scale_x_continuous
      scale_x_continuous(
        breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
        labels = month.abb,  # Usa month.name para nombres completos
        expand = c(0, 0)
      ) +
      labs(
        title = "Duración de Sequías por Mes y Año",
        x = "Mes",
        y = "Año"
      ) +
      theme(
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 14),
        legend.title = element_text(colour = "brown", face = "bold", size = 12),
        legend.position = "right",
        legend.background = element_rect(fill = "white", colour = "black", linewidth = 1),
        legend.key.size = unit(1, "cm"),
        plot.title = element_text(hjust = 0.5, size = rel(2), lineheight = .9,
                                  face = "bold", colour = "brown"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    print(grafico_severidad)
    
    if (!is.null(save)) {
      name_folfer1 = paste0(dir.save, "/MapasCalor_SSI")
      
      if (!dir.exists(name_folfer1)) {
        dir.create(name_folfer1)
      }
      
      ggsave(
        paste(name_folfer1, "/", Subcuenca, "_Magnitud.png", sep = ""),
        plot = grafico_diario,
        width = 12,
        height = 8,
        units = "in",
        dpi = 1000,
        bg = NULL,
      )
      
      ggsave(
        paste(name_folfer1, "/", Subcuenca, "_Severidad.png", sep = ""),
        plot = grafico_severidad,
        width = 12,
        height = 8,
        units = "in",
        dpi = 1000,
        bg = NULL,
      )
      message(paste0("Se ha generado graficos de calor, revise en: ", name_folfer1))
    }
    
    
  } else {
    ############################################################################
    grouping = list()
    
    for (Pixel in setdiff(colnames(data), "Fecha")) {
      D = data[,.(Fecha = Fecha, SPEI = get(Pixel))]
      run = teori_run(D)
      
      # Exclusión de eventos de sequía
      agrup = run[(Tipo == "Sequia") & (Duracion >= 3),]
      cats_spei = function(x) {
        fcase(
          x <= -2, "Sequia Extrema",
          x > -2 & x <= -1.5, "Sequia Severa",
          x > -1.5 & x <= -1, "Sequia Moderada",
          x > -1 & x <= -0.5, "Sequia Leve",
          x > -0.5 & x <= 0.5, "No Sequia",
          x > 0.5 & x <= 1, "No Sequia",
          x > 1 & x <= 1.5, "No Sequia",
          x > 1.5 & x <= 2, "No Sequia",
          x > 2, "No Sequia"
        )
      }
      
      agrup$Categoria = sapply(agrup$Magnitud, cats_spei)
      grouping[[Pixel]] = agrup
    }
    
    ############################################################################
    # Genero gráficos ilustrativos
    gp_final = function(t, Raster_Base) {
      data_extrema = function(tp) {
        res = data.frame()
        for (i in 1:length(grouping)) {
          temp = grouping[[i]]
          temp = temp %>% 
            filter(Categoria == tp)
          name = names(grouping)[i]
          
          ####
          if (nrow(temp) == 0) {
            message(paste0("No se encontro eventos de tipo ", tp, " para el pixel ", name, " Se generara un evento ficticio con valores de 0"))
            temp = data.frame(
              drought_event = "Sequia",
              Tipo = "Sequia",
              Duracion = 0,
              Magnitud = 0,
              Severidad = 0,
              Fecha_Inicio = as.Date("2022-01-01"),
              Fecha_Fin = as.Date("2022-01-01"),
              Categoria = tp,
              Pixel = name
            )
          } else {
            if (t == "Magnitud") {
              max_sever = max(temp$Magnitud)
              temp = temp[temp$Magnitud == max_sever,]
              
              temp$Pixel = name
              if (nrow(temp) > 1) {
                max_dur = max(temp$Duracion)
                temp = temp[temp$Duracion == max_dur,]
                
                if (nrow(temp) > 1) {
                  temp = temp[1,]
                }
              }
            } else {
              max_sever = max(temp$Duracion)
              temp = temp[temp$Duracion == max_sever,]
              
              temp$Pixel = name
              if (nrow(temp) > 1) {
                max_dur = max(temp$Magnitud)
                temp = temp[temp$Magnitud == max_dur,]
                
                if (nrow(temp) > 1) {
                  temp = temp[1,]
                }
              }
            }
            
          }
          res = rbind(res, temp)
        }
        return(res)
      }
      Sequia_leve = data_extrema("Sequia Leve")
      Sequia_moderada = data_extrema("Sequia Moderada")
      Sequia_severa = data_extrema("Sequia Severa")
      Sequia_Extrema = data_extrema("Sequia Extrema")
      
      name.save = paste0(dir.save, "/Categorizacion_Spacial_SPEI")
      if (!dir.exists(name.save)) {
        dir.create(name.save)
      }
      
      graph_spatial = function(data, name) {
        # generar carpeta para guardar datos
        data_sat = merge(data, cords, by = "Pixel")
        data_sat = data.table(data_sat)
        
        maps = function(conds) {
          data_naliz = data_sat[,.(Pixel, X, Y, Duracion, value = get(conds), Fecha_Inicio, Fecha_Fin)]
          puntos = vect(data_naliz, geom = c("X", "Y"), crs = "EPSG:32717")
          Base = Raster_Base[[1]]
          Base = Base * 0
          raster_frecuencia = rasterize(puntos, Base, field = "value")
          plot(raster_frecuencia)
          terra::writeCDF(raster_frecuencia, filename=paste0(name.save, "/",Subcuenca, "_", name,"_", conds, ".nc", sep = ""), overwrite=TRUE)
        }
        m2 = maps("Severidad")
        m3 = maps("Magnitud")
      }
      
      Events_leve = graph_spatial(Sequia_leve, "Sequia_Leve")
      Events_moderada = graph_spatial(Sequia_moderada, "Sequia_Moderada")
      Events_severa = graph_spatial(Sequia_severa, "Sequia_Severa")
      Events_extrema = graph_spatial(Sequia_Extrema, "Sequia_Extrema")
      
      if (t == "Magnitud") {
        name.save2 = paste0(dir.save, "/Categorizacion_Magnitud_SPEI")
      } else {
        name.save2 = paste0(dir.save, "/Categorizacion_Duracion_SPEI")
      }
      
      if (!dir.exists(name.save2)) {
        dir.create(name.save2)
      }
      
      graph_insitu = function(data, Ts) {
        data$cod = seq(as.numeric(1:nrow(data),1))
        mapear = function(val) {
          
          heatmap_data = data %>%
            dcast(Pixel ~ Duracion, value.var = val)
          
          heatmap_data_long = melt(heatmap_data, id.vars = "Pixel", 
                                    variable.name = "Duracion", value.name = val)
          
          pl = ggplot(heatmap_data_long, aes(x = Duracion, y = Pixel)) +
            geom_tile(aes(fill = get(val))) +
            scale_fill_gradientn(colors = c("blue", "yellow", "red")) +
            labs(title = paste("Heatmap de", val, "por Duración y Pixel"),
                 x = "Días",
                 y = "Pixel",
                 fill = val) +
            theme_bw() +
            theme(
              axis.title = element_text(size = 20, face = "bold"),
              axis.text = element_text(size = 12, hjust = 0.9),
              legend.title = element_text(colour = "brown", face = "bold", size = 12),
              legend.position = "bottom",
              legend.justification = c(0.5, 0.5),
              legend.background = element_rect(fill = "white", colour = "black", linewidth = 1),
              plot.title = element_text(hjust = 0.5, size = rel(2), lineheight = .9,
                                        face = "bold", colour = "brown"),
            ) +
            
            geom_rect(aes(xmin = (-Inf), xmax = (Inf), ymin = -Inf, ymax = Inf),
                      color = "black", fill = NA, size = 1.5)
          
          
          print(pl)
          
          name_f = paste0(name.save2, "/", Subcuenca, "_", Ts, "_", val, ".png", sep = "")
          ggsave(
            name_f,
            plot = pl,
            width = 12,
            height = 8,
            units = "in",
            dpi = 1000,
            bg = NULL,
          )
        }
        m1 = mapear("Severidad")
        m2 = mapear("Magnitud")
      }
      
      
      heatmap_leve = graph_insitu(Sequia_leve, "Sequia_Leve")
      heatmap_moderada = graph_insitu(Sequia_moderada, "Sequia_Moderada")
      heatmap_severa = graph_insitu(Sequia_severa, "Sequia_Severa")
      heatmap_extrema = graph_insitu(Sequia_Extrema, "Sequia_Extrema")
      
      
    }
    
    Rast = terra::rast(Raster_Base)
    graph_Magnitud = gp_final(t = "Magnitud",Raster_Base = Rast)
    graph_duracion = gp_final(t = "Duracion", Raster_Base = Rast)
    
    
    ############################################################################
  } # cierro logica del SPEI
  
  return(grouping)
}

drought_MH = function(balance_hidrico, SSI, Subcuenca, dir.save, cords, data_SPEI) {
  ##############################################################################
  names(balance_hidrico)[1] = "TIMESTAMP"
  
  # Elijo pixeles por cuenca
  if (Subcuenca == "Tomebamba") {
    print("Tomebamba")
    balance_hidrico = balance_hidrico[, .(TIMESTAMP, Pixel_1, Pixel_2, Pixel_3, Pixel_4, Pixel_5, Pixel_6, Pixel_7, Pixel_8, Pixel_9,
                                          Pixel_10, Pixel_11, Pixel_12, Pixel_13, Pixel_14, Pixel_15, Pixel_16, Pixel_17,
                                          Pixel_18), ]
    names(balance_hidrico)[2:length(balance_hidrico)] = paste0("Pixel_", 1:length(balance_hidrico))
    
    data_SPEI = data_SPEI[, .(Fecha, Pixel_1, Pixel_2, Pixel_3, Pixel_4, Pixel_5, Pixel_6, Pixel_7, Pixel_8, Pixel_9,
                              Pixel_10, Pixel_11, Pixel_12, Pixel_13, Pixel_14, Pixel_15, Pixel_16, Pixel_17,
                              Pixel_18), ]
    names(data_SPEI)[2:length(data_SPEI)] = paste0("Pixel_", 1:length(data_SPEI))
    
  } else {
    print("Yanuncay")
    balance_hidrico = balance_hidrico[, .(TIMESTAMP, Pixel_8, Pixel_9, Pixel_14, Pixel_15,Pixel_16,Pixel_17,Pixel_18,Pixel_19,Pixel_22, Pixel_23, 
                                          Pixel_24, Pixel_25, Pixel_26, Pixel_27, Pixel_28,Pixel_31,
                                          Pixel_32, Pixel_33, Pixel_34, Pixel_35, Pixel_36, Pixel_39, Pixel_40, Pixel_41, Pixel_42, 
                                          Pixel_44, Pixel_45, Pixel_46),]
    
    names(balance_hidrico)[2:length(balance_hidrico)] = paste0("Pixel_", 1:length(balance_hidrico))
    
    
  }
  
  np = length(balance_hidrico) - 1
  # Fucnoin pal SPEI acumulado
  SPEI_Pixel = function(Pixel, n) {
    data = balance_hidrico[,.(Fecha = TIMESTAMP, D = get(Pixel))]
    data = data.table(data)
    
    if (any(is.na(data$D))) {
      stop("Series con NA")
    }
    
    SPEI = SPEI::spei(data$D, scale = n, fit = "ub-pwm", verbose =F)
    
    data$SPEI = SPEI$fitted
    data[data == Inf | data == -Inf] = NA
    # if (sum(is.na(data$SPEI)) > n) {
    #  # message(paste0("Se encontro," , sum(is.na(data$SPEI)), " Nas"))
    # }
    
    data = data[,c("Fecha", "SPEI")]
    
    return(data)
  }
  
  ejecutar = function(n) {
    Fechas = data.table(balance_hidrico)[,.(Fecha = TIMESTAMP)]
    for (i in setdiff(names(balance_hidrico), c("TIMESTAMP"))) {
      SPEI = SPEI_Pixel(Pixel = i, n = n)
      names(SPEI)[2] = i
      Fechas = merge(Fechas, SPEI, by = "Fecha")
    }
    return(Fechas)
  }
  
  # llamo a la funcion mediante un for
  secuencia_lags = seq(1, 365, by = 1)
  lags = list()
  for (i in secuencia_lags) {
    message(paste("Iteracion ", i , " de ", length(secuencia_lags)))
    lags[[i]] = ejecutar(i)
  }
  
  names(lags) = paste0("lag_", secuencia_lags)
  ##############################################################################
  resultados = list()
  for (i in 1:length(lags)) {
    nombre_dia = names(lags)[i] 
    correlaciones_df = data.frame(Pixel = colnames(lags[[i]])[-1])  # Excluye "Fecha"
    
    # Calcular la correlación para cada pixel
    for (j in 2:ncol(lags[[i]])) {  # Desde la segunda columna porque la primera es "Fecha"
      
      datos = lags[[i]][, .(Fecha, Valor = lags[[i]][[j]])]  # Selección dinámica de pixel
      datos = merge(datos, SSI, by = "Fecha", all.x = TRUE)  # Unir con SSI por "Fecha"
      names(datos)[2:3] = c("SPEI", "SSI")  # Cambiar nombres de columnas
      
      # Calcular correlación
      correlacion = cor(datos[,-1], use = "complete.obs", method = "pearson")
      correlacion_valor = correlacion[1, 2]  # Extraer el valor específico de la matriz
      correlacion_valor = round(correlacion_valor, 3)  # Redondear a dos decimale
      correlaciones_df[j - 1, "Correlacion"] = correlacion_valor
    }
    
    correlaciones_df$Dias = nombre_dia
    
    resultados[[nombre_dia]] = correlaciones_df
    
  }
  
  resultados_df = rbindlist(resultados, idcol = "Archivo")
  
  # Extraigo dias con mayor correlacion
  dias_mayor_correlacion = resultados_df[, .SD[Correlacion == max(Correlacion)], by = Pixel]
  
  dias_mayor_correlacion[, Dias := as.numeric(gsub("lag_", "", Dias))]
  
  dias_mayor_correlacion = dias_mayor_correlacion[
    , .SD[which.min(Dias)], by = .(Pixel)
  ]
  
  dias_mayor_correlacion$Pixel = factor(dias_mayor_correlacion$Pixel, 
                                        levels = unique(dias_mayor_correlacion$Pixel[order(as.numeric(gsub("Pixel_", "", dias_mayor_correlacion$Pixel)))]))
  
  dias_mayor_correlacion$Dias = paste0(dias_mayor_correlacion$Dias, "_dias")
  
  pl = 
    
    ggplot(dias_mayor_correlacion, aes(x = reorder(Pixel, as.numeric(gsub("Pixel_", "", Pixel))), y = Correlacion, fill = Dias)) +
    geom_bar(stat = "identity") +
    labs(title = "Máxima correlación por Pixel",
         x = "Pixel",
         y = "Correlación",
         fill = "Días") +
    theme_bw() +
    theme(
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 12, hjust = 0.9),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "bottom",
      legend.justification = c(0.5, 0.5),
      legend.background = element_rect(fill = "white", colour = "black", linewidth = 1),
      plot.title = element_text(hjust = 0.5, size = rel(2), lineheight = .9,
                                face = "bold", colour = "brown"),
    ) +
    geom_rect(aes(xmin = (-Inf), xmax = (Inf), ymin = -Inf, ymax = Inf),
              color = "black", fill = NA, size = 1.5)
  
  pl
  
  
  folder_Save = paste0(dir.save, "/Correlaciones")
  if (!dir.exists(folder_Save)) {
    dir.create(folder_Save)
  }
  
  
  ggsave(
    paste(folder_Save, "/", Subcuenca, "_Dia_MayorCorrelacion.png", sep = ""),
    plot = pl,
    width = 12,
    height = 8,
    units = "in",
    dpi = 1000,
    bg = NULL,
  )
  
  write.csv(resultados_df, paste0(folder_Save,"/", Subcuenca, "_correlaciones.csv"), row.names = F)
  
  # genero un shape 
  arc.shape = dias_mayor_correlacion[,.(Pixel = Pixel, Correlacion = Correlacion, Dias = Dias)]
  
  arc.shape = merge(arc.shape, cords, by = "Pixel")
  puntos = vect(arc.shape, geom = c("X", "Y"), crs = "EPSG:32717")
  writeVector(puntos, paste0(folder_Save, "/Mayor_correlacion", Subcuenca, ".shp", sep = ""), filetype = "ESRI Shapefile", overwrite = TRUE)
  
  # grafico final
  resultados_f = resultados_df
  resultados_f$Dias = as.numeric(sub("lag_", "", resultados_f$Dias))  # Convertir 'lag_XX' en números
  resultados_f$Pixel = as.factor(resultados_f$Pixel)
  resultados_f$Pixel = factor(resultados_f$Pixel, levels = paste0("Pixel_", 1:np))
  
  # Generar el gráfico con el orden deseado
  pl = 
    ggplot(resultados_f, aes(x = Dias, y = Correlacion)) +
    geom_line(aes(color = Pixel, group = Pixel), size = 1) +
    geom_point(size = 0.5) +
    facet_wrap(~ Pixel, scales = "free_y") +
    theme_bw() +
    labs(x = "Días", y = "Correlación") +
    theme(
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 8, hjust = 0.9),
      legend.position = "none",  # Desactiva completamente la leyenda
      plot.title = element_text(hjust = 0.5, size = rel(2), lineheight = .9,
                                face = "bold", colour = "brown")
    ) +
    #   scale_x_continuous(breaks = seq(0, 365, by = 10)) +
    geom_rect(aes(xmin = (-Inf), xmax = (Inf), ymin = -Inf, ymax = Inf),
              color = "black", fill = NA, size = 1.5)
  
  pl
  
  ggsave(
    paste(folder_Save, "/", Subcuenca, "_correlacionDiaria.png", sep = ""),
    plot = pl,
    width = 12,
    height = 8,
    units = "in",
    dpi = 1000,
    bg = NULL,
  )
  
  ##############################################################################
  #                 calculo de los intervalos de activación                    #
  ##############################################################################
  
  sequias_hidro = IAE_droughts(SSI, Type = "SSI", tc = 5, pc = 0.1)
  sequias_hidro$Severidad = sequias_hidro$Magnitud / sequias_hidro$Duracion
  
  SSI_tibble = as_tibble(sequias_hidro)
  SSI_tibble = SSI_tibble %>%
    mutate(Categoria = case_when(
      Severidad >= 2.0 ~ "No Sequia",
      Severidad >= 1.5 ~ "No Sequia",
      Severidad >= 1.0 ~ "No Sequia",
      Severidad >= 0.0 ~ "No Sequia",
      Severidad >= -1.0 ~ "Sequía leve",
      Severidad >= -1.5 ~ "Sequía Moderada",
      Severidad >= -2.0 ~ "Sequía severa",
      TRUE ~ "Sequía extrema"
    ))
  
  sequias_hidro = as.data.table(SSI_tibble)
  
  sequias_meteo = list()
  for (Pixel in setdiff(colnames(data_SPEI), "Fecha")) {
    D = data_SPEI[,.(Fecha = Fecha, SPEI = get(Pixel))]
    run = teori_run(D, umbral = - 0.5)
    cats_spei = function(x) {
      fcase(
        x <= -2, "Sequia Extrema",
        x > -2 & x <= -1.5, "Sequia Severa",
        x > -1.5 & x <= -1, "Sequia Moderada",
        x > -1 & x <= -0.5, "Sequia Leve",
        x > -0.5 & x <= 0.5, "No Sequia",
        x > 0.5 & x <= 1, "No Sequia",
        x > 1 & x <= 1.5, "No Sequia",
        x > 1.5 & x <= 2, "No Sequia",
        x > 2, "No Sequia"
      )
    }
    
    run$Categoria = sapply(run$Magnitud, cats_spei)
    sequias_meteo[[Pixel]] = run
  }
  
  dias_tp = dias_mayor_correlacion
  dias_tp$Pixel = as.character(dias_tp$Pixel)
  dias_tp = data.frame(dias_tp)
  
  sequias_emparejadas = list()
  
  for (j in 1:length(sequias_meteo)) {
    data_meteo = sequias_meteo[[j]]
    Pixel = names(sequias_meteo)[j]
    
    data_meteo = data.table(data_meteo)
    tp = dias_tp[dias_tp$Pixel == Pixel, "Dias"]
    tp = as.numeric(gsub("_dias", "", tp))
    res_pixel = data.frame()
    
    for (i in 2:(nrow(sequias_hidro) - 1)) {
      
      if (i > nrow(sequias_hidro) -1) {
        break
      }
      
      evento_actual = sequias_hidro[i, ]
      evento_anterior = sequias_hidro[i - 1, ]
      
      tsi = evento_actual$Fecha_Inicio
      tei_1 = evento_anterior$Fecha_Fin
      diferencia =  as.numeric(difftime(tsi, tei_1 - 1, units = "days"))
      
      print(paste0(Pixel, " - ", i, " - ", "Sequia_n", i))
      if (diferencia >= tp) {
        inicio = evento_actual$Fecha_Inicio - tp
        fin = evento_actual$Fecha_Fin
        
        meteo = data_meteo[Fecha_Inicio >= inicio & Fecha_Fin <= fin, ]
        sequia = meteo[Tipo == "Sequia",]
        if (nrow(sequia) == 0) {
          sequia = data.frame(
            Magnitud = 0,
            Duracion = 0
          )
        }
        
        humedo = meteo[Tipo == "Humedo",]
        if (nrow(humedo) == 0) {
          humedo = data.frame(
            Magnitud = 0,
            Duracion = 0
          )
        }
        
        resulados = data.frame(
          Pixel = Pixel,
          Trigger_Inicio = inicio,
          Trigger_Fin = fin,
          Fecha_inicioHidro = evento_actual$Fecha_Inicio,
          Fecha_finHidro = evento_actual$Fecha_Fin,
          Duracion_hidro = evento_actual$Duracion,
          Magnitud_hidro = evento_actual$Magnitud,
          Severidad_hidro = evento_actual$Severidad,
          Categoria_hidro = evento_actual$Categoria,
          Magnitud_meteo = sum(sequia$Magnitud) - sum(humedo$Magnitud),
          Duracion_meteo = sum(sequia$Duracion)
        )
        
      } else {
        
        inicio = evento_anterior$Fecha_Fin
        fin = evento_actual$Fecha_Fin
        
        meteo = data_meteo[Fecha_Inicio >= inicio & Fecha_Fin <= fin, ]
        sequia = meteo[Tipo == "Sequia",]
        
        if (nrow(sequia) == 0) {
          sequia = data.frame(
            Magnitud = 0,
            Duracion = 0
          )
        }
        
        humedo = meteo[Tipo == "Humedo",]
        if (nrow(humedo) == 0) {
          humedo = data.frame(
            Magnitud = 0,
            Duracion = 0
          )
        }
        
        resulados = data.frame(
          Pixel = Pixel,
          Trigger_Inicio = evento_anterior$Fecha_Fin,
          Trigger_Fin = evento_actual$Fecha_Inicio,
          
          Fecha_inicioHidro = evento_actual$Fecha_Inicio,
          Fecha_finHidro = evento_actual$Fecha_Fin,
          Duracion_hidro = evento_actual$Duracion,
          Magnitud_hidro = evento_actual$Magnitud,
          Severidad_hidro = evento_actual$Severidad,
          Categoria_hidro = evento_anterior$Categoria,
          Magnitud_meteo = sum(sequia$Magnitud) - sum(humedo$Magnitud),
          Duracion_meteo = sum(sequia$Duracion)
        )
      }
      
      res_pixel = rbind(res_pixel, resulados)
      cats_spei = function(x) {
        fcase(
          x <= -2, "Sequia Extrema",
          x > -2 & x <= -1.5, "Sequia Severa",
          x > -1.5 & x <= -1, "Sequia Moderada",
          x > -1 & x <= 0, "Sequia Leve",
          x > 0, "No Sequia"
        )
      }
      
      res_pixel$Categoria_hidro = sapply(res_pixel$Magnitud_hidro, cats_spei)
      
    }
    sequias_emparejadas[[j]] = res_pixel
  }
  
  # CALCULO LA SEVERIDAD Y CATEGORIZO DE NUEBO 
  seq_finales = list()
  
  cats_spei = function(x) {
    fcase(
      x >= 2, "Sequia Extrema",
      x >= 1.5 & x < 2, "Sequia Severa",
      x >= 1 & x < 1.5, "Sequia Moderada",
      x >= 0.5 & x < 1, "Sequia Leve",
      x < 0.5, "No Sequia"
    )
  }
  
  for (i in 1:length(sequias_emparejadas)) {
    data = sequias_emparejadas[[i]]
    data$Severidad_hidro = abs(data$Magnitud_hidro / data$Duracion_hidro)
    data$Categoria_hidro = sapply(data$Severidad_hidro, cats_spei)
    data$Severidad_meteo = abs(data$Magnitud_meteo / data$Duracion_meteo)
    data$Categoria_meteo = sapply(data$Severidad_meteo, cats_spei)
    
    data[data == Inf | data == -Inf] = NA
    data = na.omit(data)
    
    seq_finales[[i]] = data
  }
  
  # Por pixel, extraigo la minima cdondicion de Categoria_hidro y la maxima categoria de Categoria_hidro
  resumen = data.frame()
  for (i in 1:length(seq_finales)) {
    data = seq_finales[[i]]
    
    resultados_min = data %>%
      group_by(Categoria_hidro) %>%
      mutate(Condicion = "Minimas") %>%
      filter(Severidad_hidro == min(Severidad_hidro, na.rm = TRUE))
    
    resultados_max = data %>%
      group_by(Categoria_hidro) %>%
      mutate(Condicion = "Maximas") %>%
      filter(Severidad_hidro == max(Severidad_hidro, na.rm = TRUE))
    
    # Unir ambos resultados
    resultados = bind_rows(resultados_min, resultados_max) %>%
      arrange(Categoria_hidro, Severidad_hidro)
    
    resumen = rbind(resumen, resultados)
  } 
  
  # Guardar resultados
  folder_Save = paste0(dir.save, "/Impacto_Sequias")
  if (!dir.exists(folder_Save)) {
    dir.create(folder_Save)
  }
  
  write.csv(resumen, paste0(folder_Save, "/", Subcuenca, "_Impacto_Sequias.csv"), row.names = F)
  final = merge(resumen, cords, by = "Pixel")
  puntos = vect(final, geom = c("X", "Y"), crs = "EPSG:32717")
  writeVector(puntos, paste0(folder_Save, "/Minimas_Maximas_impactos_correlacion", Subcuenca, ".shp", sep = ""), filetype = "ESRI Shapefile", overwrite = TRUE)
}


# propagation_thresholds = function(seq_emparejadas, dir.save) {
#   # Ajuste de distribuciones marginales
#   best_distribution = function(x, Pixel, nombre, tipo) {
#     
#     values = abs(as.numeric(x))
#     
#     ###################### Buscar la mejor distribución ########################
#     # Gamma
#     fit_gamma = fitdist(values, "gamma")
#     AIC_gamma = fit_gamma$aic
#     BIC_gamma = fit_gamma$bic
#     CV_gamma = fit_gamma$loglik
#     
#     # exponencial
#     fit_exp = fitdist(values, "exp", method = "mle")
#     AIC_exp = fit_exp$aic
#     BIC_exp = fit_exp$bic
#     CV_exp = fit_exp$loglik
#     
#     # weibull 
#     # values_clean = values[values > 0]
#     fit_weibull = fitdist(values, "weibull")
#     AIC_weibull = fit_weibull$aic
#     BIC_weibull = fit_weibull$bic
#     CV_weibull = fit_weibull$loglik
#     
#     # GEV
#     fit_gev = fevd(values, type = "GEV")
#     CV_gev = -1 * fit_gev$results$value
#     k = length(fit_gev$parnames) 
#     n = length(values)
#     AIC_gev = -2 * CV_gev + 2 * k
#     BIC_gev = -2 * CV_gev + k * log(n)
#     
#     ############################################################################
#     # Función para seleccionar el threshold óptimo
#     # select_optimal_threshold = function(data, umin, umax, show_plot = TRUE) {
#     #   fit_results = gpd.fitrange(data = data, umin = umin, umax = umax, show = FALSE)
#     #   
#     #   # Extraer parámetros relevantes
#     #   thresholds = fit_results$thresholds
#     #   xi_values = fit_results$mle[, 2]  # Estimaciones de xi
#     #   ci_low = fit_results$ci.low[, 2]
#     #   ci_up = fit_results$ci.up[, 2]
#     #   
#     #   # Calcular estabilidad de xi
#     #   xi_diff = abs(diff(xi_values))  # Diferencias entre valores consecutivos
#     #   ci_width = ci_up - ci_low       # Ancho de intervalos de confianza
#     #   
#     #   # Definir criterio: estabilidad de xi y confianza estrecha
#     #   stability_scores = xi_diff / ci_width[-1]  # Ignorar el primer elemento de ci_width
#     #   optimal_index = which.min(stability_scores)  # Índice del mejor threshold
#     #   
#     #   # Threshold óptimo
#     #   optimal_threshold = thresholds[optimal_index]
#     #   
#     #   # Mostrar resultados
#     #   if (show_plot) {
#     #     library(ggplot2)
#     #     results_df = data.frame(
#     #       threshold = thresholds,
#     #       xi = xi_values,
#     #       ci.low = ci_low,
#     #       ci.up = ci_up
#     #     )
#     #     
#     #     pl = 
#     #       ggplot(results_df, aes(x = threshold, y = xi)) +
#     #       geom_line() +
#     #       geom_point() +
#     #       geom_errorbar(aes(ymin = ci.low, ymax = ci.up), width = 0.5) +
#     #       geom_vline(xintercept = optimal_threshold, linetype = "dashed", color = "red") +
#     #       labs(
#     #         title = "Selección del Threshold Óptimo",
#     #         x = "Threshold",
#     #         y = "Parámetro ξ"
#     #       ) +
#     #       theme_minimal()
#     #     print(pl)
#     #   }
#     # 
#     #   return(optimal_threshold)
#     # }
#     # 
#     # umin = quantile(values, 0.1)
#     # umax = quantile(values, 0.80)
#     # sequias = select_optimal_threshold(values, umin, umax)
#     
#     # Generalizado de pareto
#     # fit_gpd = fevd(values, type = "GP", threshold = sequias)
#     # CV_gpd = -1 * fit_gpd$results$value
#     # K = length(fit_gpd$parnames)
#     # n = length(values)
#     # AIC_gpd = -2 * CV_gpd + 2 * k
#     # BIC_gpd = -2 * CV_gpd + k * log(n)
#     
#     ############################################################################
#     ############################################################################
#     #                   Ajuste de cada distribución a los datos
#     ############################################################################
#     CDF_gamma = pgamma(values, shape = fit_gamma$estimate[1], rate = fit_gamma$estimate[2])
#     CDF_exp = pexp(values, rate = fit_exp$estimate)
#     CDF_weibull = pweibull(values, shape = fit_weibull$estimate[1], scale = fit_weibull$estimate[2])
#     CDF_gev = pgev(values, loc = fit_gev$results$par[[1]], scale = fit_gev$results$par[[2]], shape = fit_gev$results$par[[3]])
#     # if (length(fit_gpd$results$par) > 2) {
#     #   loc_gpd = fit_gpd$results$par[[1]]
#     #   scale_gpd = fit_gpd$results$par[[2]]
#     #   shape_gpd = fit_gpd$results$par[[3]]
#     # } else {
#     #   loc_gpd = 0
#     #   scale_gpd = fit_gpd$results$par[[1]]
#     #   shape_gpd = fit_gpd$results$par[[2]]
#     # }
#     # 
#     # CDF_gpd = pgpd(values, loc = loc_gpd, scale = scale_gpd, shape = shape_gpd)
#     
#     ############################################################################
#     #                                    KS-TEST                               #
#     ks_gamma = ks.test(values, "pgamma", shape = fit_gamma$estimate[1], rate = fit_gamma$estimate[2])
#     ks_exp = ks.test(values, "pexp", rate = fit_exp$estimate)
#     ks_weibull = ks.test(values, "pweibull", shape = fit_weibull$estimate[1], scale = fit_weibull$estimate[2])
#     ks_gev = ks.test(values, "pgev", loc = fit_gev$results$par[[1]], scale = fit_gev$results$par[[2]], shape = fit_gev$results$par[[3]])
#     #ks_gpd = ks.test(values, "pgpd", loc = loc_gpd, scale = scale_gpd, shape = shape_gpd)
#     
#     ############################################################################
#     #                                    QQ-PLOT                               #
#     PP_Plot = function(distr, nombre, ks) {
#       pl = ppcomp(distr, plotstyle = "ggplot") +   
#         geom_abline(intercept = 0, slope = 1, color = "black", size = 1) +   
#         geom_point(aes(color = NULL), size = 3, fill = "red", stroke = 1) +  # Desvincular color de la leyenda
#         labs(title = paste0("P-P ", nombre), 
#              x = "Probabilidades teóricas", 
#              y = "Probabilidades empíricas") +    
#         theme_bw() +    
#         guides(color = guide_legend(title = paste("p-value:", ks)))+  # Crea una leyenda personalizada
#         theme(     
#           legend.position = "none",
#           axis.title = element_text(size = 20, face = "bold"),     
#           axis.text = element_text(size = 18, hjust = 0.9),     
#           plot.title = element_text(hjust = 0.5, size = rel(2), lineheight = .9, face = "bold", colour = "brown")   
#         ) +   
#         scale_x_continuous(breaks = seq(0, 1, by = 0.25)) +   
#         geom_rect(aes(xmin = (-Inf), xmax = (Inf), ymin = -Inf, ymax = Inf), 
#                   color = "black", fill = NA, size = 1.5) +   
#         
#         annotate("text", x = 0.2 , y = 0.5, label = paste("p-value: ", ks), size = 5, color = "brown") 
#       return(pl)
#     }
#     
#     PP_Gev = function(distr, nombre, ks) {
#       
#       pp = plot(distr, type = "probprob")  # Probabilidades teóricas y empíricas
#       teo_prob = pp$model  # Probabilidades teóricas calculadas
#       emp_prob = pp$empirical  # Probabilidades empíricas calculadas
#       data_plot = data.frame(teo_prob, emp_prob)
#       
#       # Crear el gráfico ggplot
#       pl = ggplot(data_plot, aes(x = teo_prob, y = emp_prob)) +
#         geom_abline(intercept = 0, slope = 1, color = "black", size = 1) +   
#         geom_point(aes(x = teo_prob, y = emp_prob), fill = "red", color = "red", size = 3, stroke = 1) +
#         labs(title = paste0("P-P ", nombre), 
#              x = "Probabilidades teóricas", 
#              y = "Probabilidades empíricas") +    
#         theme_bw() +    
#         guides(color = guide_legend(title = paste("p-value:", ks)))+  # Crea una leyenda personalizada
#         theme(     
#           legend.position = "none",
#           axis.title = element_text(size = 20, face = "bold"),     
#           axis.text = element_text(size = 18, hjust = 0.9),     
#           plot.title = element_text(hjust = 0.5, size = rel(2), lineheight = .9, face = "bold", colour = "brown")   
#         ) +   
#         scale_x_continuous(breaks = seq(0, 1, by = 0.25)) +   
#         geom_rect(aes(xmin = (-Inf), xmax = (Inf), ymin = -Inf, ymax = Inf), 
#                   color = "black", fill = NA, size = 1.5) +   
#         
#         annotate("text", x = 0.2 , y = 0.5, label = paste("p-value: ", ks), size = 5, color = "brown") 
#       
#       return(pl)
#     }
#     
#     graph_gamma = PP_Plot(fit_gamma, "Gamma", round(ks_gamma$p.value,3))
#     graph_exp = PP_Plot(fit_exp, "Exponencial", round(ks_exp$p.value,3))
#     graph_weibull = PP_Plot(fit_weibull, "Weibull", round(ks_weibull$p.value,3))
#     graph_gev = PP_Gev(fit_gev, "GEV", round(ks_gev$p.value,3))
#     
#     combined_plot = grid.arrange(
#       grobs = list(graph_gamma, graph_exp, graph_weibull, graph_gev),
#       ncol = 2,  # Número de columnas
#       nrow = 2  # Número de filas
#     )
#     
#     folder_save = paste0(dir.save, "/Graficos_pp/")
#     if (!dir.exists(folder_save)) {
#       dir.create(folder_save)
#     }
#     
#     name = paste(Pixel, nombre, tipo, sep = "_")
#     ggsave(
#       paste(folder_save, name, ".png", sep = ""),
#       plot = combined_plot,
#       width = 12,
#       height = 8,
#       units = "in",
#       dpi = 1000,
#       bg = NULL,
#     )
#     ############################################################################
#     #                                   best distribution                      #
#     ############################################################################
#     
#     best_disrt = data.frame(
#       Distribucion = c("Exponencial", "Weibull", "GEV", "Gamma"),
#       AIC = c(AIC_exp, AIC_weibull, AIC_gev, AIC_gamma),
#       BIC = c(BIC_exp, BIC_weibull, BIC_gev, BIC_gamma),
#       Max_veros = c(CV_exp, CV_weibull, CV_gev, CV_gamma),
#       Ks = c(ks_exp$statistic, ks_weibull$statistic, ks_gev$statistic, ks_gamma$statistic),
#       p_Ks = c(ks_exp$p.value, ks_weibull$p.value, ks_gev$p.value, ks_gamma$p.value)
#     )
#     ############################################################################
#     #                             Elijo mejor distribucion                     #
#     best_disrt = best_disrt[order(best_disrt$p_Ks, decreasing = TRUE), ][1:2, ]
#     best_disrt = best_disrt[best_disrt$p_Ks > 0.05,]
#     
#     AIC_1 = which.min(best_disrt$AIC)
#     BIC_1 = which.min(best_disrt$BIC)
#     CV_1 =  which.max(best_disrt$Max_veros)
#     KS_1 = which.min(best_disrt$Ks)
#     p_KS_1 = which.max(best_disrt$p_Ks)
#     
#     res = data.frame(
#       AIC = best_disrt[AIC_1, "Distribucion"],
#       BIC = best_disrt[BIC_1, "Distribucion"],
#       CV = best_disrt[CV_1, "Distribucion"],
#       KS = best_disrt[KS_1, "Distribucion"],
#       p_KS = best_disrt[p_KS_1, "Distribucion"]
#     )
#     
#     res_vector = unlist(res)
#     frecuencias = table(res_vector)
#     mas_frecuente = names(frecuencias[frecuencias == max(frecuencias)])
#     
#     distr_elegida = best_disrt[best_disrt$Distribucion == mas_frecuente,]
#     
#     if (mas_frecuente == "Exponencial") {
#       dis_f = CDF_exp
#     } else if (mas_frecuente == "Weibull") {
#       dis_f = CDF_weibull
#     } else if (mas_frecuente == "GEV") {
#       dis_f = CDF_gev
#     } else if (mas_frecuente == "Gamma") {
#       dis_f = CDF_gamma
#     } else {
#       stop("No se encontro la distribución")
#     }
#     return(list(best_ajuste = distr_elegida, distribucion = dis_f))
#   } # cierro función de la distribución 
#   
#   distribucion_marginal = function(data, Pixel) {
#     setDT(data)
#     data = data[,.(H_D = Duracion_hidrologica, M_D = Duracion_meteorologica, 
#                    H_S = Magnitud_hidrologica, M_S = Magnitud_meteorologica)]
#     
#     d_HD = best_distribution(data$H_D, Pixel = Pixel, nombre = "Duracion", tipo = "Hidrologica")
#     d_MD = best_distribution(data$M_D, Pixel = Pixel, nombre = "Duracion", tipo = "Meteorologica")
#     d_HS = best_distribution(data$H_S, Pixel = Pixel, nombre = "Magnitud", tipo = "Hidrologica")
#     d_MS = best_distribution(data$M_S, Pixel = Pixel, nombre = "Magnitud", tipo = "Meteorologica")
#     
#     ############################################################################
#     df_dHD = cbind(data.frame(name = rep("HD", each = nrow(d_HD$best_ajuste))), d_HD$best_ajuste)
#     df_dMD = cbind(data.frame(name = rep("MD", each = nrow(d_MD$best_ajuste))), d_MD$best_ajuste)
#     df_dHS = cbind(data.frame(name = rep("HS", each = nrow(d_HS$best_ajuste))), d_HS$best_ajuste)
#     df_dMS = cbind(data.frame(name = rep("MS", each = nrow(d_MS$best_ajuste))), d_MS$best_ajuste)
#     
#     
#     ta.res = rbind(df_dHD, df_dMD, df_dHS, df_dMS)
#     name_pixel = rep(Pixel, nrow(ta.res))
#     ta.res = cbind(name_pixel, ta.res) 
#     
#     return(list(estadisticos = ta.res, d_HD = d_HD$distribucion, d_MD = d_MD$distribucion, d_HS = d_HS$distribucion, d_MS = d_MS$distribucion))
#   }
#   
#   distribuciones_marginales = list()
#   estadisticos = data.frame()
#   
#   for (i in 1:length(seq_emparejadas)) {
#     message(paste0("Iteración ", i, " de ", length(seq_emparejadas)))
#     data = seq_emparejadas[[i]]
#     Pixel = unique(data$Pixel)
#     variable = distribucion_marginal(data, Pixel)
#     
#     distribuciones_marginales[[i]] = list(d_HD = variable$d_HD, d_MD = variable$d_MD, 
#                                           d_HS = variable$d_HS, d_MS = variable$d_MS)
#     
#     estadisticos = rbind(estadisticos, variable$estadisticos)
#   }
#   
#   ##############################################################################
#   #                           Distribución conjunta                            #
#   # u' = condición (meteorológicas)
#   # v' = objetivo (hidrológicas)
#   fmlas_copulas = function(data, u, v) {
#     # Validación de datos
#     if (!u %in% names(data) | !v %in% names(data)) {
#       stop("Las columnas especificadas no están en los datos.")
#     }
#     
#     u_vec = data[[u]] # meteorológicas
#     v_vec = data[[v]] # hidrologicas
#     
#     if (any(u_vec < 0 | u_vec > 1) | any(v_vec < 0 | v_vec > 1)) {
#       stop("Los valores de u y v deben estar en el rango [0, 1].")
#     }
#     
#     # Crear las cópulas
#     copulas = list(
#       Gumbel = gumbelCopula(dim = 2),
#       Franck = frankCopula(dim = 2),
#       Clayton = claytonCopula(dim = 2),
#       Gaussiana = normalCopula(dim = 2)
#     )
#     
#     # Ajustar las cópulas y calcular métricas
#     resultados = lapply(copulas, function(copula) {
#       tryCatch({
#         ajuste = fitCopula(copula, data = cbind(u_vec, v_vec), method = "ml")
#         list(
#           loglik = logLik(ajuste),
#           AIC = AIC(ajuste),
#           BIC = BIC(ajuste),
#           coef = coef(ajuste)
#         )
#       }, error = function(e) {
#         list(loglik = NA, AIC = NA, BIC = NA, coef = NA)
#       })
#     })
#     
#     # Crear el dataframe de resultados
#     resultados_df = do.call(rbind, lapply(names(resultados), function(nombre) {
#       c(copula = nombre, unlist(resultados[[nombre]]))
#     }))
#     
#     resultados_df = as.data.frame(resultados_df, stringsAsFactors = FALSE)
#     resultados_df$loglik = as.numeric(resultados_df$loglik)
#     resultados_df$AIC = as.numeric(resultados_df$AIC)
#     resultados_df$BIC = as.numeric(resultados_df$BIC)
#     
#     # Selección de la mejor cópula
#     mejor_AIC = resultados_df$copula[which.min(resultados_df$AIC)]
#     mejor_BIC = resultados_df$copula[which.min(resultados_df$BIC)]
#     mejor_loglik = resultados_df$copula[which.max(resultados_df$loglik)]
#     
#     # Resolver empates
#     criterios = table(c(mejor_AIC, mejor_BIC, mejor_loglik))
#     mas_frecuente = names(criterios[criterios == max(criterios)])
#     mejor_copula = mas_frecuente[1]  # En caso de empate, toma la primera
#     
#     # Calcular la probabilidad conjunta
#     theta = resultados[[mejor_copula]]$coef
#     copula_final = switch(mejor_copula,
#                           Gumbel = gumbelCopula(theta, dim = 2),
#                           Franck = frankCopula(theta, dim = 2),
#                           Clayton = claytonCopula(theta, dim = 2),
#                           Gaussiana = normalCopula(theta, dim = 2))
#     
#     C_uv = pCopula(cbind(u_vec, v_vec), copula_final)
#     
#     # Retornar resultados
#     return(list(
#       resultados = resultados_df,
#       mejor_copula = mejor_copula,
#       prob_conjunta = C_uv
#     ))
#   }
#   
#   dependencia_concicional = function(x_u, y_v, C_xu_yv) {
#     numerador = 1 - x_u - y_v + C_xu_yv
#     denominador = 1 - x_u
#     Probabilidad = numerador / denominador
#     
#     if (Probabilidad < 0 | Probabilidad > 1) {
#       stop("La probabilidad condicional debe estar en el rango [0, 1].")
#     }
#     return(Probabilidad)
#   }
#   
#   dependencia = list()
#   mejor_copula = list()
#   res_copulas = list()
#   
#   for (i in 1:length(distribuciones_marginales)) {
#     data = distribuciones_marginales[[i]]
#     print(paste0("iteracion ", i, " de ", length(distribuciones_marginales)))
#     # condición 1 # MD - HD
#     cond_1 = fmlas_copulas(data = data, u = "d_MD", v = "d_HD")
#     cop_1 = cond_1$prob_conjunta
#     xu1 = distribuciones_marginales[[i]]$d_MD
#     yv1 = distribuciones_marginales[[i]]$d_HD
#     analisis_1 = data.frame(cbind(xu1, yv1, cop_1))
#     prob_1 = apply(analisis_1, 1, function(x) dependencia_concicional(x[1], x[2], x[3]))
#     
#     res_prob1 = data.frame(cbind(
#       Xu_MD = xu1,
#       Yv_HD = yv1,
#       Cop_MD_HD = cop_1,
#       Prob_MD_HD = prob_1))
#     
#     
#     # condición 2 # MD - HS
#     cond_2 = fmlas_copulas(data = data, u = "d_MD", v = "d_HS")
#     cop_2 = cond_2$prob_conjunta
#     xu2 = distribuciones_marginales[[i]]$d_MD
#     yv2 = distribuciones_marginales[[i]]$d_HS
#     analisis_2 = data.frame(cbind(xu2, yv2, cop_2))
#     prob_2 = apply(analisis_2, 1, function(x) dependencia_concicional(x[1], x[2], x[3]))
#     
#     res_prob2 = data.frame(cbind(
#       Xu2_MD = xu2,
#       Yv2_HS = yv2,
#       Cop_MD_HS = cop_2,
#       Prob_MD_HS = prob_2))
#     
#     # condición 3 # MS - HD
#     cond_3 = fmlas_copulas(data = data, u = "d_MS", v = "d_HD")
#     cop_3 = cond_3$prob_conjunta
#     xu3 = distribuciones_marginales[[i]]$d_MS
#     yv3 = distribuciones_marginales[[i]]$d_HD
#     analisis_3 = data.frame(cbind(xu3, yv3, cop_3))
#     prob_3 = apply(analisis_3, 1, function(x) dependencia_concicional(x[1], x[2], x[3]))
#     
#     res_prob3 = data.frame(cbind(
#       Xu3_MS = xu3,
#       Yv3_HD = yv3,
#       Cop_MS_HD = cop_3,
#       Prob_MS_HD = prob_3))
#     
#     
#     # condición 4 # MS - HS
#     cond_4 = fmlas_copulas(data = data, u = "d_MS", v = "d_HS")
#     cop_4 = cond_4$prob_conjunta
#     xu4 = distribuciones_marginales[[i]]$d_MS
#     yv4 = distribuciones_marginales[[i]]$d_HS
#     analisis_4 = data.frame(cbind(xu4, yv4, cop_4))
#     prob_4 = apply(analisis_4, 1, function(x) dependencia_concicional(x[1], x[2], x[3]))
#     
#     res_prob4 = data.frame(cbind(
#       Xu4_MS = xu4,
#       Yv4_HS = yv4,
#       Cop_MS_HS = cop_4,
#       Prob_MS_HS = prob_4))
#     
#     # genero un solo df de resultados
#     res = data.frame(cbind(res_prob1, res_prob2, res_prob3, res_prob4))
#     dependencia[[i]] = res
#     
#     # genero df de las copulas que mejor se ajustan
#     res_bestCopula = data.frame(
#       MD_HD = cond_1$mejor_copula,
#       MD_HS = cond_2$mejor_copula,
#       MS_HD = cond_3$mejor_copula,
#       MS_HS = cond_4$mejor_copula
#     )
#     
#     # GENERO DF de los resultados de cada copula]
#     ests_1 = cond_1$resultados[c("copula", "AIC", "BIC", "loglik")]
#     ests_2 = cond_2$resultados[c("copula", "AIC", "BIC", "loglik")]
#     ests_3 = cond_3$resultados[c("copula", "AIC", "BIC", "loglik")]
#     ests_4 = cond_4$resultados[c("copula", "AIC", "BIC", "loglik")]
#     
#     columnas_completas = unique(c(names(ests_1), 
#                                   names(ests_2), 
#                                   names(ests_3), 
#                                   names(ests_4)))
#     
#     alinear_columnas = function(df, columnas) {
#       faltantes = setdiff(columnas, names(df)) # Identificar columnas faltantes
#       for (col in faltantes) {
#         df[[col]] = NA # Agregar columnas faltantes con NA
#       }
#       df = df[, columnas] # Reordenar columnas en el orden correcto
#       return(df)
#     }
#     
#     # Alinear todos los data frames
#     df1 = alinear_columnas(cond_1$resultados, columnas_completas)
#     df2 = alinear_columnas(cond_2$resultados, columnas_completas)
#     df3 = alinear_columnas(cond_3$resultados, columnas_completas)
#     df4 = alinear_columnas(cond_4$resultados, columnas_completas)
#     
#     res_copulas = rbind(df1, df2, df3, df4)
#     
#     #res_copulas = rbind(cond_1$resultados, cond_2$resultados, cond_3$resultados, cond_4$resultados)
#     res_copulas = cbind(data.frame(cond = rep(c("MD_HD", "MD_HS", "MS_HD", "MS_HS"), each = 4)), res_copulas)
#     
#     mejor_copula[[i]] = res_bestCopula
#     res_copulas[[i]] = res_copulas
#     
#   }
#   
#   ##############################################################################
#   # Unión de todos los datos para extraer información y clasificar umbrales.....
#   clasificacion = function(x) {
#     fcase(
#       x >= 0 & x < 0.5, "Sequia leve",
#       x >= 0.5 & x < 0.75, "Sequia moderada",
#       x >= 0.75 & x < 0.9, "Sequia severa",
#       x >= 0.9, "Sequia extrema"
#     )
#   }
#   
#   umbrales_propagacion = list()
#   for (i in 1:length(dependencia)) {
#     data = dependencia[[i]]
#     setDT(data)
#     
#     # MD - HD
#     cond_1 = data[,.(Xu_MD = Xu_MD, Yv_HD = Yv_HD, Cop_MD_HD = Cop_MD_HD, Prob_MD_HD = Prob_MD_HD)]
#     cond_1$Categoria_MD_HD = sapply(cond_1$Cop_MD_HD, clasificacion)
#     cond_1 = cond_1[Prob_MD_HD >= 0.95,]
#     
#     # MD - HS
#     cond_2 = data[,.(Xu2_MD = Xu2_MD, Yv2_HS = Yv2_HS, Cop_MD_HS = Cop_MD_HS, Prob_MD_HS = Prob_MD_HS)]
#     cond_2$Categoria_MD_HS = sapply(cond_2$Cop_MD_HS, clasificacion)
#     cond_2 = cond_2[Prob_MD_HS >= 0.95,]
#     
#     # MS - HD
#     cond_3 = data[,.(Xu3_MS = Xu3_MS, Yv3_HD = Yv3_HD, Cop_MS_HD = Cop_MS_HD, Prob_MS_HD = Prob_MS_HD)]
#     cond_3$Categoria_MS_HD = sapply(cond_3$Cop_MS_HD, clasificacion)
#     cond_3 = cond_3[Prob_MS_HD >= 0.95,]
#     
#     # MS - HS
#     cond_4 = data[,.(Xu4_MS = Xu4_MS, Yv4_HS = Yv4_HS, Cop_MS_HS = Cop_MS_HS, Prob_MS_HS = Prob_MS_HS)]
#     cond_4$Categoria_MS_HS = sapply(cond_4$Cop_MS_HS, clasificacion)
#     cond_4 = cond_4[Prob_MS_HS >= 0.95,]
#     
#     umbrales_propagacion[[i]] = list(
#       MD_HD = cond_1,
#       MD_HS = cond_2,
#       MS_HD = cond_3,
#       MS_HS = cond_4)
#   }
#   
#   return(umbrales_propagacion)
#   
# }
############################################################################################################
# propagation_thresholds = function(seq_emparejadas, dir.save) {
#   # Ajuste de distribuciones marginales
#   best_distribution = function(x, Pixel, nombre, tipo) {
#     
#     values = abs(as.numeric(x))
#     ##############################################################################
#     #                            L-moments                                       #
#     selection = c("gev", "gum", "gam", "pe3", "ln3", "exp", "wei")
#     
#     ajust_1 = distLfit(
#       values,
#       datname = deparse(substitute(values)),
#       selection = selection,
#       speed = T,
#       ks = T
#     )
#     
#     # Bondas de ajuste de las distribuciones      
#     gof = data.frame(
#       Distribuciones = rownames(ajust_1$gof),
#       RMSE = ajust_1$gof$RMSE,
#       W1 = ajust_1$gof$weight1,
#       W2 = ajust_1$gof$weight2,
#       W3 = ajust_1$gof$weight3,
#       ksP = ajust_1$gof$ksP,
#       ksD = ajust_1$gof$ksD,
#       R2 = ajust_1$gof$R2
#     )
#     
#     gof$RMSE_norm = (max(gof$RMSE) - gof$RMSE) / (max(gof$RMSE) - min(gof$RMSE))
#     gof$ksP_norm = (gof$ksP - min(gof$ksP)) / (max(gof$ksP) - min(gof$ksP))
#     gof$ksD_norm = (max(gof$ksD) - gof$ksD) / (max(gof$ksD) - min(gof$ksD))
#     gof$R2_norm = (gof$R2 - min(gof$R2)) / (max(gof$R2) - min(gof$R2))
#     
#     # Asignación de pesos
#     weights = c(RMSE = 0.4, ksP = 0.3, ksD = 0.1, R2 = 0.2)
#     
#     # Cálculo del puntaje combinado
#     gof$Score = gof$RMSE_norm * weights["RMSE"] +
#       gof$ksP_norm * weights["ksP"] +
#       gof$ksD_norm * weights["ksD"] +
#       gof$R2_norm * weights["R2"]
#     
#     best_distribution = gof[which.max(gof$Score), ]
#     parametros_best = ajust_1$par[best_distribution$Distribuciones][[1]]
#     ##############################################################################
#     #                          Maximum likelihood estimation                     #
#     distr_mle = list(
#       lnorm = fitdist(values, "lnorm", method = "mle"),
#       exp = fitdist(values, "exp", method = "mle"),
#       weibull = fitdist(values, "weibull", method = "mle"),
#       gamma = fitdist(values, "gamma", method = "mle")
#     )
#     
#     fits = list(lnorm = distr_mle$lnorm, exp = distr_mle$exp , weibull = distr_mle$weibull, gamma = distr_mle$gamma)
#     gof_stats = gofstat(fits, fitnames = names(fits))
#     
#     gof_2 = data.frame(
#       Distribuciones = c("lnorm", "exp", "weibull", "gamma"),
#       KS = gof_stats$ks,
#       AD = gof_stats$ad,
#       AIC = gof_stats$aic,
#       BIC = gof_stats$bic
#     )
#     
#     weigths_2 = c(KS = 0.4, AD = 0.3, AIC = 0.2, BIC = 0.1)
#     
#     gof_2$KS_norm = (gof_2$KS - min(gof_2$KS)) / (max(gof_2$KS) - min(gof_2$KS))
#     gof_2$AD_norm = (max(gof_2$AD) - gof_2$AD) / (max(gof_2$AD) - min(gof_2$AD))
#     gof_2$AIC_norm = (max(gof_2$AIC) - gof_2$AIC) / (max(gof_2$AIC) - min(gof_2$AIC))
#     gof_2$BIC_norm = (max(gof_2$BIC) - gof_2$BIC) / (max(gof_2$BIC) - min(gof_2$BIC))
#     
#     gof_2$Score = gof_2$KS_norm * weigths_2["KS"] +
#       gof_2$AD_norm * weigths_2["AD"] +
#       gof_2$AIC_norm * weigths_2["AIC"] +
#       gof_2$BIC_norm * weigths_2["BIC"]
#     
#     best_distribution_2 = gof_2[which.max(gof_2$Score), ]
#     
#     parametros_best_mle = distr_mle[best_distribution_2$Distribuciones][[1]]
#     
#     # Parametros de la mejor distribucion 
#     
#     param_mle = data.frame(
#       loc = parametros_best_mle$estimate[1],
#       scale = ifelse (length(parametros_best_mle$estimate) > 2, parametros_best_mle$estimate[2], NA),
#       shape = ifelse (length(parametros_best_mle$estimate) > 2, parametros_best_mle$estimate[3], NA)
#     )
#     
#     param_max = data.frame(
#       loc = parametros_best$para[1],
#       scale = parametros_best$para[2],
#       shape = ifelse (length(parametros_best$para) > 2, parametros_best$para[3], NA)
#     )
#     
#     # Elijo la mejor distribución
#     merge_f = data.frame(
#       Distribuciones = c(best_distribution$Distribuciones, best_distribution_2$Distribuciones),
#       Ks = c(best_distribution$ksP, best_distribution_2$KS)
#     )
#     
#     # Mejor de las mejores
#     best_distribution_f = merge_f[which.max(merge_f$Ks), ]
#     
#     CDF_distribuciones = function(loc = NULL, scale = NULL, shape = NULL) {
#       if (best_distribution_f$Distribuciones == "gev") {
#         return(TLMoments::pgev(values, loc = loc, scale = scale, shape = shape))
#       } else if (best_distribution_f$Distribuciones == "glo") {
#         return(LMoFit::pglo(values, para = c(loc, scale, shape)))
#       } else if (best_distribution_f$Distribuciones == "gum") {
#         return(TLMoments::pgum(values, loc = loc, scale = scale))
#       } else if (best_distribution_f$Distribuciones == "lnorm") {
#         return(stats::plnorm(values, meanlog = loc, sdlog = scale))
#       } else if (best_distribution_f$Distribuciones == "exp") {
#         return(stats::pexp(values, rate = loc))
#       } else if (best_distribution_f$Distribuciones == "gamma") {
#         return(stats::pgamma(values, shape = loc, rate = scale))
#       } else if (best_distribution_f$Distribuciones == "weibull") {
#         return(stats::pweibull(values, shape = loc, scale = scale))
#       } else {
#         stop("Distribución desconocida.")
#       }
#     }
#     
#     if (best_distribution_f$Distribuciones %in% c("gev", "glo", "gum")) {
#       CDF = CDF_distribuciones(loc = param_max$loc, scale = param_max$scale, shape = param_max$shape)
#     } else {
#       CDF = CDF_distribuciones(loc = param_mle$loc, scale = param_mle$scale, shape = param_mle$shape)
#     }
#     
#     
#     return(list(best_ajuste = best_distribution_f, distribucion = CDF))
#     
#   }
#   
#   distribucion_marginal = function(data, Pixel) {
#     setDT(data)
#     data = data[,.(H_D = Duracion_hidrologica, M_D = Duracion_meteorologica, 
#                    H_S = Magnitud_hidrologica, M_S = Magnitud_meteorologica)]
#     
#     d_HD = best_distribution(data$H_D, Pixel = Pixel, nombre = "Duracion", tipo = "Hidrologica")
#     d_MD = best_distribution(data$M_D, Pixel = Pixel, nombre = "Duracion", tipo = "Meteorologica")
#     d_HS = best_distribution(data$H_S, Pixel = Pixel, nombre = "Magnitud", tipo = "Hidrologica")
#     d_MS = best_distribution(data$M_S, Pixel = Pixel, nombre = "Magnitud", tipo = "Meteorologica")
#     
#     ############################################################################
#     df_dHD = cbind(data.frame(name = rep("HD", each = nrow(d_HD$best_ajuste))), d_HD$best_ajuste)
#     df_dMD = cbind(data.frame(name = rep("MD", each = nrow(d_MD$best_ajuste))), d_MD$best_ajuste)
#     df_dHS = cbind(data.frame(name = rep("HS", each = nrow(d_HS$best_ajuste))), d_HS$best_ajuste)
#     df_dMS = cbind(data.frame(name = rep("MS", each = nrow(d_MS$best_ajuste))), d_MS$best_ajuste)
#     
#     
#     ta.res = rbind(df_dHD, df_dMD, df_dHS, df_dMS)
#     name_pixel = rep(Pixel, nrow(ta.res))
#     ta.res = cbind(name_pixel, ta.res) 
#     
#     return(list(estadisticos = ta.res, d_HD = d_HD$distribucion, d_MD = d_MD$distribucion, d_HS = d_HS$distribucion, d_MS = d_MS$distribucion))
#   }
#   
#   distribuciones_marginales = list()
#   estadisticos = data.frame()
#   
#   for (i in 1:length(seq_emparejadas)) {
#     message(paste0("Iteración ", i, " de ", length(seq_emparejadas)))
#     data = seq_emparejadas[[i]]
#     Pixel = unique(data$Pixel)
#     variable = distribucion_marginal(data, Pixel)
#     
#     distribuciones_marginales[[i]] = list(d_HD = variable$d_HD, d_MD = variable$d_MD, 
#                                           d_HS = variable$d_HS, d_MS = variable$d_MS)
#     
#     estadisticos = rbind(estadisticos, variable$estadisticos)
#   }
#   
#   ##############################################################################
#   #                           Distribución conjunta                            #
#   # u' = condición (meteorológicas)
#   # v' = objetivo (hidrológicas)
#   fmlas_copulas = function(data, u, v) {
#     # Validación de datos
#     if (!u %in% names(data) | !v %in% names(data)) {
#       stop("Las columnas especificadas no están en los datos.")
#     }
#     
#     u_vec = data[[u]] # meteorológicas
#     v_vec = data[[v]] # hidrologicas
#     
#     if (any(u_vec < 0 | u_vec > 1) | any(v_vec < 0 | v_vec > 1)) {
#       stop("Los valores de u y v deben estar en el rango [0, 1].")
#     }
#     
#     # Crear las cópulas
#     copulas = list(
#       Gumbel = gumbelCopula(dim = 2),
#       Franck = frankCopula(dim = 2),
#       Clayton = claytonCopula(dim = 2),
#       Gaussiana = normalCopula(dim = 2)
#     )
#     
#     # Ajustar las cópulas y calcular métricas
#     resultados = lapply(copulas, function(copula) {
#       tryCatch({
#         ajuste = fitCopula(copula, data = cbind(u_vec, v_vec), method = "ml")
#         list(
#           loglik = logLik(ajuste),
#           AIC = AIC(ajuste),
#           BIC = BIC(ajuste),
#           coef = coef(ajuste)
#         )
#       }, error = function(e) {
#         list(loglik = NA, AIC = NA, BIC = NA, coef = NA)
#       })
#     })
#     
#     # Crear el dataframe de resultados
#     resultados_df = do.call(rbind, lapply(names(resultados), function(nombre) {
#       c(copula = nombre, unlist(resultados[[nombre]]))
#     }))
#     
#     resultados_df = as.data.frame(resultados_df, stringsAsFactors = FALSE)
#     resultados_df$loglik = as.numeric(resultados_df$loglik)
#     resultados_df$AIC = as.numeric(resultados_df$AIC)
#     resultados_df$BIC = as.numeric(resultados_df$BIC)
#     
#     # Selección de la mejor cópula
#     mejor_AIC = resultados_df$copula[which.min(resultados_df$AIC)]
#     mejor_BIC = resultados_df$copula[which.min(resultados_df$BIC)]
#     mejor_loglik = resultados_df$copula[which.max(resultados_df$loglik)]
#     
#     # Resolver empates
#     criterios = table(c(mejor_AIC, mejor_BIC, mejor_loglik))
#     mas_frecuente = names(criterios[criterios == max(criterios)])
#     mejor_copula = mas_frecuente[1]  # En caso de empate, toma la primera
#     
#     # Calcular la probabilidad conjunta
#     theta = resultados[[mejor_copula]]$coef
#     copula_final = switch(mejor_copula,
#                           Gumbel = gumbelCopula(theta, dim = 2),
#                           Franck = frankCopula(theta, dim = 2),
#                           Clayton = claytonCopula(theta, dim = 2),
#                           Gaussiana = normalCopula(theta, dim = 2))
#     
#     C_uv = pCopula(cbind(u_vec, v_vec), copula_final)
#     
#     # Retornar resultados
#     return(list(
#       resultados = resultados_df,
#       mejor_copula = mejor_copula,
#       prob_conjunta = C_uv
#     ))
#   }
#   
#   ##############################################################################
#   dependencia_concicional = function(x_u, y_v, C_xu_yv) {
#     numerador = 1 - x_u - y_v + C_xu_yv
#     denominador = 1 - x_u
#     Probabilidad = numerador / denominador
#     
#     if (Probabilidad < 0 | Probabilidad > 1) {
#       stop("La probabilidad condicional debe estar en el rango [0, 1].")
#     }
#     return(Probabilidad)
#   }
#   
#   dependencia = list()
#   mejor_copula = list()
#   res_copulas = list()
#   
#   for (i in 1:length(distribuciones_marginales)) {
#     data = distribuciones_marginales[[i]]
#     print(paste0("iteracion ", i, " de ", length(distribuciones_marginales)))
#     # condición 1 # MD - HD
#     cond_1 = fmlas_copulas(data = data, u = "d_MD", v = "d_HD")
#     cop_1 = cond_1$prob_conjunta
#     xu1 = distribuciones_marginales[[i]]$d_MD
#     yv1 = distribuciones_marginales[[i]]$d_HD
#     analisis_1 = data.frame(cbind(xu1, yv1, cop_1))
#     prob_1 = apply(analisis_1, 1, function(x) dependencia_concicional(x[1], x[2], x[3]))
#     
#     res_prob1 = data.frame(cbind(
#       Xu_MD = xu1,
#       Yv_HD = yv1,
#       Cop_MD_HD = cop_1,
#       Prob_MD_HD = prob_1))
#     
#     
#     # condición 2 # MD - HS
#     cond_2 = fmlas_copulas(data = data, u = "d_MD", v = "d_HS")
#     cop_2 = cond_2$prob_conjunta
#     xu2 = distribuciones_marginales[[i]]$d_MD
#     yv2 = distribuciones_marginales[[i]]$d_HS
#     analisis_2 = data.frame(cbind(xu2, yv2, cop_2))
#     prob_2 = apply(analisis_2, 1, function(x) dependencia_concicional(x[1], x[2], x[3]))
#     
#     res_prob2 = data.frame(cbind(
#       Xu2_MD = xu2,
#       Yv2_HS = yv2,
#       Cop_MD_HS = cop_2,
#       Prob_MD_HS = prob_2))
#     
#     # condición 3 # MS - HD
#     cond_3 = fmlas_copulas(data = data, u = "d_MS", v = "d_HD")
#     cop_3 = cond_3$prob_conjunta
#     xu3 = distribuciones_marginales[[i]]$d_MS
#     yv3 = distribuciones_marginales[[i]]$d_HD
#     analisis_3 = data.frame(cbind(xu3, yv3, cop_3))
#     prob_3 = apply(analisis_3, 1, function(x) dependencia_concicional(x[1], x[2], x[3]))
#     
#     res_prob3 = data.frame(cbind(
#       Xu3_MS = xu3,
#       Yv3_HD = yv3,
#       Cop_MS_HD = cop_3,
#       Prob_MS_HD = prob_3))
#     
#     
#     # condición 4 # MS - HS
#     cond_4 = fmlas_copulas(data = data, u = "d_MS", v = "d_HS")
#     cop_4 = cond_4$prob_conjunta
#     xu4 = distribuciones_marginales[[i]]$d_MS
#     yv4 = distribuciones_marginales[[i]]$d_HS
#     analisis_4 = data.frame(cbind(xu4, yv4, cop_4))
#     prob_4 = apply(analisis_4, 1, function(x) dependencia_concicional(x[1], x[2], x[3]))
#     
#     res_prob4 = data.frame(cbind(
#       Xu4_MS = xu4,
#       Yv4_HS = yv4,
#       Cop_MS_HS = cop_4,
#       Prob_MS_HS = prob_4))
#     
#     # genero un solo df de resultados
#     res = data.frame(cbind(res_prob1, res_prob2, res_prob3, res_prob4))
#     dependencia[[i]] = res
#     
#     # genero df de las copulas que mejor se ajustan
#     res_bestCopula = data.frame(
#       MD_HD = cond_1$mejor_copula,
#       MD_HS = cond_2$mejor_copula,
#       MS_HD = cond_3$mejor_copula,
#       MS_HS = cond_4$mejor_copula
#     )
#     
#     # GENERO DF de los resultados de cada copula]
#     ests_1 = cond_1$resultados[c("copula", "AIC", "BIC", "loglik")]
#     ests_2 = cond_2$resultados[c("copula", "AIC", "BIC", "loglik")]
#     ests_3 = cond_3$resultados[c("copula", "AIC", "BIC", "loglik")]
#     ests_4 = cond_4$resultados[c("copula", "AIC", "BIC", "loglik")]
#     
#     columnas_completas = unique(c(names(ests_1), 
#                                   names(ests_2), 
#                                   names(ests_3), 
#                                   names(ests_4)))
#     
#     alinear_columnas = function(df, columnas) {
#       faltantes = setdiff(columnas, names(df)) # Identificar columnas faltantes
#       for (col in faltantes) {
#         df[[col]] = NA # Agregar columnas faltantes con NA
#       }
#       df = df[, columnas] # Reordenar columnas en el orden correcto
#       return(df)
#     }
#     
#     # Alinear todos los data frames
#     df1 = alinear_columnas(cond_1$resultados, columnas_completas)
#     df2 = alinear_columnas(cond_2$resultados, columnas_completas)
#     df3 = alinear_columnas(cond_3$resultados, columnas_completas)
#     df4 = alinear_columnas(cond_4$resultados, columnas_completas)
#     
#     res_copulas = rbind(df1, df2, df3, df4)
#     
#     #res_copulas = rbind(cond_1$resultados, cond_2$resultados, cond_3$resultados, cond_4$resultados)
#     res_copulas = cbind(data.frame(cond = rep(c("MD_HD", "MD_HS", "MS_HD", "MS_HS"), each = 4)), res_copulas)
#     
#     mejor_copula[[i]] = res_bestCopula
#     res_copulas[[i]] = res_copulas
#     
#   }
#   
#   ##############################################################################
#   # Unión de todos los datos para extraer información y clasificar umbrales.....
#   clasificacion = function(x) {
#     fcase(
#       x >= 0 & x < 0.5, "Sequia leve",
#       x >= 0.5 & x < 0.75, "Sequia moderada",
#       x >= 0.75 & x < 0.9, "Sequia severa",
#       x >= 0.9, "Sequia extrema"
#     )
#   }
#   
#   umbrales_propagacion = list()
#   umbral = 0.9
#   for (i in 1:length(dependencia)) {
#     data = dependencia[[i]]
#     setDT(data)
#     
#     # MD - HD
#     cond_1 = data[,.(Xu_MD = Xu_MD, Yv_HD = Yv_HD, Cop_MD_HD = Cop_MD_HD, Prob_MD_HD = Prob_MD_HD)]
#     cond_1$Categoria_MD_HD = sapply(cond_1$Cop_MD_HD, clasificacion)
#     cond_1 = cond_1[Prob_MD_HD >= umbral,]
#     
#     # MD - HS
#     cond_2 = data[,.(Xu2_MD = Xu2_MD, Yv2_HS = Yv2_HS, Cop_MD_HS = Cop_MD_HS, Prob_MD_HS = Prob_MD_HS)]
#     cond_2$Categoria_MD_HS = sapply(cond_2$Cop_MD_HS, clasificacion)
#     cond_2 = cond_2[Prob_MD_HS >= umbral,]
#     
#     # MS - HD
#     cond_3 = data[,.(Xu3_MS = Xu3_MS, Yv3_HD = Yv3_HD, Cop_MS_HD = Cop_MS_HD, Prob_MS_HD = Prob_MS_HD)]
#     cond_3$Categoria_MS_HD = sapply(cond_3$Cop_MS_HD, clasificacion)
#     cond_3 = cond_3[Prob_MS_HD >= umbral,]
#     
#     # MS - HS
#     cond_4 = data[,.(Xu4_MS = Xu4_MS, Yv4_HS = Yv4_HS, Cop_MS_HS = Cop_MS_HS, Prob_MS_HS = Prob_MS_HS)]
#     cond_4$Categoria_MS_HS = sapply(cond_4$Cop_MS_HS, clasificacion)
#     cond_4 = cond_4[Prob_MS_HS >= umbral,]
#     
#     umbrales_propagacion[[i]] = list(
#       MD_HD = cond_1,
#       MD_HS = cond_2,
#       MS_HD = cond_3,
#       MS_HS = cond_4)
#   }
#   
#   return(umbrales_propagacion)
#   
# }
############################################################################################################