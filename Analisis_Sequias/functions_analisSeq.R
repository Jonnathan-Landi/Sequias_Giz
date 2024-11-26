################################################################################
#   Análisis Integral de la Recurrencia y Respuesta Espacial de las Sequías    #
#        Meteorológicas e Hidrológicas en la Región de Cuenca, Ecuador         #
################################################################################
# Autores: 
#'  Jonnathan Landi
#'  Marco Mogro 
# Fecha: 2024-11-15
# Versión: 1.0.0
################################################################################
# Librerías necesarias
packages = c("data.table", "dplyr", "purrr", "terra", "ggplot2", "plotly",
             "htmlwidgets", "gridExtra", "shiny", "shinythemes", "bslib", 
             "shinyWidgets", "lubridate", "tidyr", "reshape2","matrixStats",
             "changepoint", "zyp", "trend", "Kendall", "missForest", "parallel",
             "doParallel")

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
import_data = function(data, type, Fecha_Inicio, Fecha_Fin, Subcuenca, dir.save = NULL) {
  # compruebo que data sea un archivo data.table
  if (!is.data.table(data)) {
    stop("El archivo de datos debe ser un data.table, use 'fread' para cargar los datos")
  }
  
  if (type == "SPEI") {
    names(data)[1] = "Fecha"
    data[data == Inf | data == -Inf] = NA
    data = data[Fecha >= Fecha_Inicio & Fecha <= Fecha_Fin, ]
    
    if (Subcuenca == "Tomebamba") {
      print("Se ha importado los datos de SPEI para Tomebamba")
      data = data[, .(Fecha, Pixel_1, Pixel_2, Pixel_3, Pixel_4, Pixel_5, Pixel_6, Pixel_7, Pixel_8, Pixel_9,
                      Pixel_10, Pixel_11, Pixel_12, Pixel_13, Pixel_14, Pixel_15, Pixel_16, Pixel_17,
                      Pixel_18,Pixel_19, Pixel_20, Pixel_21, Pixel_23, Pixel_24, Pixel_25, Pixel_26,
                      Pixel_27, Pixel_28, Pixel_29, Pixel_30, Pixel_38), ]
    } else {
      print("Se ha importado los datos de SPEI para Yanuncay")
      data = data[, .(Fecha, Pixel_8, Pixel_9, Pixel_14, Pixel_15,Pixel_16,Pixel_17,Pixel_18,Pixel_19,Pixel_22, Pixel_23, 
                      Pixel_24, Pixel_25, Pixel_26, Pixel_27, Pixel_28,Pixel_31,
                      Pixel_32, Pixel_33, Pixel_34, Pixel_35, Pixel_36, Pixel_39, Pixel_40, Pixel_41, Pixel_42, 
                      Pixel_44, Pixel_45, Pixel_46),]
    }
    
    names(data)[2:length(data)] = paste0("Pixel_", 1:length(data))
    
    if (any(is.na(data))) {
      ncom = nrow(data)
      num_cores <- detectCores() - 1 
      
      # Registrar un backend paralelo
      cl <- makeCluster(num_cores)  
      registerDoParallel(cl) 
      
      
      message("Espere un momento, estamos haciendo ajustes internamente")
      relleno = missForest::missForest(data[, -1], ntree = 2000, parallelize = "variables")
      stopCluster(cl)  # Apagar el clúster paralelo
      registerDoSEQ()  # Registrar el modo secuencial
      
      # uno con mis datos
      Fechas = data.table(data[,.(Fecha = Fecha)])
      rellenado = as.data.table(relleno$ximp)
      
      data = cbind(Fechas, rellenado)
      
      if (ncom != nrow(data)) {
        stop("Ha ocurrido un error en el proceso de imputación, revise los datos")
      }

      name = paste0(dir.save, "/relleno")
      if (!dir.exists(name)) {
        dir.create(name)
      }
      
      if (any(is.na(data))) {
        message("Se han encontrado valores NA en los datos, esto podría alterar los resultados")
      }
      
      write.csv(data, paste0(name, "/", Subcuenca, ".csv", sep = ""), row.names = F)
      
    }
    
  } else {
    data = data[Fecha >= Fecha_Inicio & Fecha <= Fecha_Fin, ]
  }
  return(data)
}
################################################################################
#                        Categorización de las sequías                         # 
categorizar_sequias = function(data, type, dir.save, Estat, save = F) {
  ##############################################################################
  # Sequias meteorologicas -----------------------------------------------------
  if (type == "SPEI") {
    
    cats_spei = function(x) {
      fcase(
        is.na(x), NA_character_,
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
    
    Clasification = list()
    for (i in setdiff(names(data), "Fecha")) {
      temp_dt = data[, .(i = get(i), Categoria = cats_spei(.SD[[i]])), by = Fecha]
      setnames(temp_dt, "i", i)
      Clasification[[i]] = temp_dt
    }
    
    # Genero gráficos ilustrativos
    valor = unique(year(data$Fecha))
    informacion_anual = data.frame(
      Año = rep(valor, 5),
      Categoria = rep(c("No Sequia", "Sequia Leve", "Sequia Moderada", "Sequia Severa", "Sequia Extrema"), length(valor))
    )
    
    informacion_anual = informacion_anual %>%
      arrange(Año)
    
    for (i in 1:length(Clasification)) {
      temp = Clasification[[i]]
      temp = na.omit(temp)
      name = colnames(temp)[2]
      
      temp = temp %>%
        mutate(Año = year(Fecha))
      
      conteo_anual <- temp %>%
        group_by(Año, Categoria) %>%
        summarise(conteo = n()) %>%
        ungroup()
      
      names(conteo_anual)[3] = name
      informacion_anual = merge(informacion_anual, conteo_anual, by = c("Año", "Categoria"), all = TRUE)
    }
    
    informacion_anual$conteo = matrixStats::rowMedians(as.matrix(informacion_anual[, 3:ncol(informacion_anual)]), na.rm = TRUE)
    informacion_anual$conteo = round(informacion_anual$conteo, 0)
    
    pl = ggplot(informacion_anual, aes(x = factor(Año), y = conteo, fill = Categoria)) +
      geom_bar(stat = "identity", position = "stack") + 
      labs(title = "Categorización de sequías anuales",
           x = "Año", y = "Frecuencia") +
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
      
      geom_text(aes(label = conteo), 
                position = position_stack(vjust = 0.5), 
                color = "black", size = 4) +  
      
      
      geom_rect(aes(xmin = (-Inf), xmax = (Inf), ymin = -Inf, ymax = Inf),
                color = "black", fill = NA, size = 1.5)
    print(pl)
    
    if (save == T) {
      name_folder = paste0(dir.save, "/Categorias_SPEI")
      message(paste0("Se ha generado gráficos, revise en: ", name_folder))
      
      if (!dir.exists(name_folder)) {
        dir.create(name_folder)
      }
      
      ggsave(
        paste(name_folder, "/", Estat, ".png", sep = ""),
        plot = pl,
        width = 12,
        height = 8,
        units = "in",
        dpi = 2000,
        bg = NULL,
      )
      
    }
  } else {
    SSI_tibble = as_tibble(data)
    SSI_tibble = SSI_tibble %>%
      mutate(Categoria = case_when(
        Z >= 2.0 ~ "No Sequia",
        Z >= 1.5 ~ "No Sequia",
        Z >= 1.0 ~ "No Sequia",
        Z >= 0.0 ~ "No Sequia",
        Z >= -1.0 ~ "Sequía leve",
        Z >= -1.5 ~ "Sequía Moderada",
        Z >= -2.0 ~ "Sequía severa",
        TRUE ~ "Sequía extrema"
      ))
    Clasification = as.data.table(SSI_tibble)
    data_graph = Clasification
    data_graph = data.frame(data_graph)
    data_graph$Fecha = as.Date(data_graph$Fecha)
    data_graph$Año = year(data_graph$Fecha)
    
    conteo_anual = data_graph %>%
      group_by(Año, Categoria) %>%
      summarise(conteo = n()) %>%
      ungroup()
    
    pl = ggplot(conteo_anual, aes(x = factor(Año), y = conteo, fill = Categoria)) +
      geom_bar(stat = "identity", position = "stack") + 
      labs(title = "Categorización de sequías anuales",
           x = "Año", y = "Frecuencia") +
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
      
      geom_text(aes(label = conteo), 
                position = position_stack(vjust = 0.5), 
                color = "black", size = 4) +  
      
      
      geom_rect(aes(xmin = (-Inf), xmax = (Inf), ymin = -Inf, ymax = Inf),
                color = "black", fill = NA, size = 1.5)
    
    print(pl)
    
    if (save == T) {
      name_folder = paste0(dir.save, "/Categorias_SSI")
      message(paste0("Se ha generado gráficos, revise en: ", name_folder))
      
      if (!dir.exists(name_folder)) {
        dir.create(name_folder)
      }
      
      ggsave(
        paste(name_folder, "/", Estat, ".png", sep = ""),
        plot = pl,
        width = 12,
        height = 8,
        units = "in",
        dpi = 2000,
        bg = NULL,
      )
      
    }
    
  } # Cierro else de seq.hidro
  return(Clasification)
}
################################################################################
#                Funciones para caracterización de las sequías                 #   
teori_run = function(data, umbral = NULL) {
  setDT(data)
  
  if (ncol(data) > 2) {
    stop("El data.table debe tener solo dos columnas")
  }
  
  names(data)[2] = "value"
  if (any(is.na(data))) {
    message("Los datos tienen NA, esto podria causar datos erroneos. A continuación se omitirá")
    data = na.omit(data)
  }
  
  # Verifico posibles an
  if(any(is.na(data$value))) {
    message("Advertencia: Se han encontrado valores NA en los datos, esto podría alterar los resultados")
  }
  
  # calculo umbrales
  if (is.null(umbral)) {
    data$drought = fifelse(data$value < -0.5, 1, 
                           fifelse(data$value > 0.5, 0, 2))
    
    data$Magnitud = fifelse(data$value < -0.5, 0.5 + data$value,
                            fifelse(data$value > 0.5, data$value - 0.5, 0))
    
    
    
  } else {
    data$drought = fifelse(data$value < umbral, 1, 0)
    data$Magnitud = fifelse(data$value < umbral, 0.5 + data$value, data$value)
  }
  
  data[, drought_event := cumsum(drought != data.table::shift(drought, fill = 0))]
  
  # Sequia
  droughts = data[drought == 1]
  drought_seco <- droughts[, .(
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
  drought_humedo <- droughts[, .(
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
  return(resultados)
}

drought_grouping = function(res_corridas, tc, pc, data_original, Estat, dir.save, save) {
  sequias = res_corridas[Tipo == "Sequia", ]
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
    s = evento_actual$Magnitud
    
    # Verificar si hay humedad entre Fecha_inicio y Fecha_fin
    humedad = res_corridas[Tipo == "Humedo" & Fecha_Inicio >= Fecha_inicio & Fecha_Fin <= Fecha_fin, ]
    
    if (nrow(humedad) > 0) {
      vi = sum(humedad$Magnitud)
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
        Magnitud = evento_actual$Magnitud + evento_siguiente$Magnitud - vi
        #Magnitud = evento_actual$Magnitud + evento_siguiente$Magnitud
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
        s = evento_agrupado$Magnitud
        
        Fecha_inicio_humedad = evento_actual$Fecha_Fin
        Fecha_fin_humedad = evento_siguiente$Fecha_Inicio
        
        humedad = res_corridas[Tipo == "Humedo" & Fecha_Inicio >= Fecha_inicio_humedad & Fecha_Fin <= Fecha_fin_humedad, ]
        
        if (nrow(humedad) > 0) {
          vi = sum(humedad$Magnitud)
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
            Magnitud = evento_actual$Magnitud + evento_siguiente$Magnitud - vi
            #Magnitud = evento_actual$Magnitud + evento_siguiente$Magnitud
          )
          z = z + 1
          contador = contador + 1
        } else {
          z = 0
          contador = contador + 2
        }
      } # Cierro wl while
      
      K = K + contador
      evento_agrupado$Severidad = evento_agrupado$Magnitud / evento_agrupado$Duracion
      eventos_agrupados[[i]] = evento_agrupado # A1qui estoy dentro del if, (dehberioa guardar aqui los datos)
    } else {
      eventos_agrupados[[i]] = evento_actual
      K = K + 1
    }
  } # Cierro el for
  resultados_df = rbindlist(eventos_agrupados, fill = TRUE)
  resultados_df = unique(resultados_df, by = "Fecha_Fin")
  resultados_df = rbindlist(eventos_agrupados, fill = TRUE)
  resultados_df = unique(resultados_df, by = "Fecha_Fin")
  resultados_df$drought_event = NULL
  ##############################################################################
  #                             Genero gráficos                                #
  data_graph = resultados_df %>%
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
      
      # Ajustar la escala de fechas
      # scale_x_date(
      #   date_breaks = "1 month",
      #   labels = function(x) {
      #     meses <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
      #     meses[as.numeric(format(x, "%m"))]
      #   }
      # ) +
      scale_x_date(
        date_breaks = ifelse(Type == TRUE, "1 month", "1 year"),
        labels = function(x) {
          if (Type == TRUE) {
            meses <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
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
  # Genero grafiuco General
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
  
  print(grafico_general)
  print(combined_plot)
  
  if (save == T) {
    message("Los graficos se estan guardando, PORFAVOR ESPERE")
    folder_save = paste0(dir.save, "/Graph_gruping_SSI/")
    if (!dir.exists(folder_save)) {
      dir.create(folder_save)
    }
    ggsave(
      paste(folder_save, Estat, "_Grupal.png", sep = ""),
      plot = combined_plot,
      width = 12,
      height = 8,
      units = "in",
      dpi = 2000,
      bg = NULL,
    )
    
    ggsave(
      paste(folder_save, Estat, "_General.png", sep = ""),
      plot = grafico_general,
      width = 12,
      height = 8,
      units = "in",
      dpi = 2000,
      bg = NULL,
    )
    
    # Genero gráfico interactivo
    pl_interactive = ggplotly(grafico_general, tooltip = "text")
    print(pl_interactive)
    name = paste0(folder_save, "Sequia_Iteractivo", Estat, ".html")
    saveWidget(pl_interactive, file = name)
  } else {
    message("Ha configurado modo visualziacion, si desea guardar: 'Save = T'")
  }
  return(resultados_df)
} # Cierro la función 

################################################################################
#                           Caracterizo mis sequias                            #
categorize_droughts = function(data, Type, tc = NULL, pc = NULL, Estat, dir.save, 
                               save = F, Raster_Base = NULL, cords = NULL) {
  # Sequías hidrológicas -------------------------------------------------------
  if (Type == "SSI") {
    
    run = teori_run(data)
    grouping = drought_grouping(run, tc, pc, data, Estat, dir.save, save)
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
    #                             Genero gráficos                              #
    data_graph = grouping
    data_graph = data.frame(data_graph)
    data_graph$Fecha_Inicio = as.Date(data_graph$Fecha_Inicio)
    data_graph$Fecha_Fin = as.Date(data_graph$Fecha_Fin)
    data_graph$Año = year(data_graph$Fecha_Inicio)
    data_graph$Mes = month(data_graph$Fecha_Inicio)
    
    conteo_anual <- data_graph %>%
      group_by(Año, Mes, Categoria, Duracion,  Magnitud) %>%
      summarise(conteo = n()) %>%
      ungroup()
    data_diaria <- data_graph %>%
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
    
    data_diaria <- data_diaria %>%
      mutate(
        Mes = factor(
          lubridate::month(Fecha_Inicio, label = TRUE, abbr = FALSE), 
          levels = month.name
        )
      )
    
    colores <- c(
      "Sequía leve" = "#94d2bd",
      "Sequía Moderada" = "#ffea00",
      "Sequía severa" = "#ff7b00",
      "Sequía extrema" = "#ad2831"
    )
    
    grafico_diario <- ggplot(data_diaria, aes(x = Dias, y = factor(Año), fill = Categoria)) +
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
    
    conteo_anual <- data_graph %>%
      group_by(Año, Mes, Categoria, Duracion, Severidad) %>%
      summarise(conteo = n()) %>%
      ungroup()
    data_diaria <- data_graph %>%
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
    
    data_diaria <- data_diaria %>%
      mutate(
        Mes = factor(
          lubridate::month(Fecha_Inicio, label = TRUE, abbr = FALSE), 
          levels = month.name
        )
      )
    
    colores <- c(
      "Sequía leve" = "#94d2bd",
      "Sequía Moderada" = "#ffea00",
      "Sequía severa" = "#ff7b00",
      "Sequía extrema" = "#ad2831"
    )
    
    grafico_severidad <- ggplot(data_diaria, aes(x = Dias, y = factor(Año), fill = Categoria)) +
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
    #######################
    if (save == T) {
      name_folfer1 = paste0(dir.save, "/MapasCalor_SSI")
      
      if (!dir.exists(name_folfer1)) {
        dir.create(name_folfer1)
      }
      
      ggsave(
        paste(name_folfer1, "/", Estat, "_Magnitud.png", sep = ""),
        plot = grafico_diario,
        width = 12,
        height = 8,
        units = "in",
        dpi = 2000,
        bg = NULL,
      )
      
      ggsave(
        paste(name_folfer1, "/", Estat, "_Severidad.png", sep = ""),
        plot = grafico_severidad,
        width = 12,
        height = 8,
        units = "in",
        dpi = 2000,
        bg = NULL,
      )
      message(paste0("Se ha generado graficos de calor, revise en: ", name_folfer1))
    }
    
  } else {
    
    grouping = list()
    for (Pixel in setdiff(colnames(data), "Fecha")) {
      D = data[,.(Fecha = Fecha, SPEI = get(Pixel))]
      run = teori_run(D, umbral = - 0.5)
      
      # Exclusión de eventos de sequía
      agrup = run[(Tipo == "Sequia") & (Duracion >= 3),]
      
      # Categorizo mis sequias #################################################
      # Espacialmente 
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
          terra::writeCDF(raster_frecuencia, filename=paste0(name.save, "/",Estat, "_", name,"_", conds, ".nc", sep = ""), overwrite=TRUE)
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
          
          heatmap_data <- data %>%
            dcast(Pixel ~ Duracion, value.var = val)
          
          heatmap_data_long <- melt(heatmap_data, id.vars = "Pixel", 
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
          
          name_f = paste0(name.save2, "/", Estat, "_", Ts, "_", val, ".png", sep = "")
          ggsave(
            name_f,
            plot = pl,
            width = 12,
            height = 8,
            units = "in",
            dpi = 2000,
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
    
  } # Cierro logica de las sequias meteorologicas
  return(grouping)
}# cierro mi funcion de categorizacion
################################################################################
#                           Análisis estadístico                               #
statistical_analysis = function(data, Typo, cords, Estat, dir.save) {
  if (Typo == "SPEI") {
    ############################################################################
    #                           Analisis de tendencia                          #
    resultados_tendencia = data.frame()
    resultados_pendiente = data.frame()
    i = 1
    for (Pixel in setdiff(colnames(data), "Fecha")) {
      data_analisis = data[,.(Fecha = Fecha, SPEI = get(Pixel))]
      spei_data = data_analisis$SPEI
      spei_ts = ts(spei_data, start = c(2014, as.numeric(format(data_analisis$Fecha[1], "%j"))), frequency = 365)
      decomposed = stl(spei_ts, s.window = "periodic")
      trend_component = decomposed$time.series[, "trend"]
      result = MannKendall(trend_component)
      res_global = data.frame(
        cod = i,
        Pixel = Pixel,
        tau = as.numeric(result[1]),
        p_value = as.numeric(result[2]),
        Tipo = "Tendencia Global"
      )
      
      # 
      resultados_tendencia = rbind(resultados_tendencia, res_global)
      
      # Escala anual
      spei_anual <- data_analisis %>%
        mutate(Año = as.numeric(format(Fecha, "%Y"))) %>%
        group_by(Año) %>%
        summarize(SPEI_anual = median(SPEI, na.rm = TRUE))
      
      # Kendall
      result_anual <- MannKendall(spei_anual$SPEI_anual)
      
      # pendiente con theil
      pendiente_anual = zyp.sen(SPEI_anual ~ Año, data = spei_anual)
      magnitud_cambio = pendiente_anual$coefficients["Año"] * 14
      
      res_anualPendiente = data.frame(
        cod = i,
        Pixel = Pixel,
        Tipo = "Anual",
        pendiente = as.numeric(pendiente_anual$coefficients["Año"]),
        cambio_acumulado = magnitud_cambio
      )
      
      resultados_pendiente = rbind(resultados_pendiente, res_anualPendiente)
      
      res_anual = data.frame(
        cod = i,
        Pixel = Pixel,
        tau = as.numeric(result_anual[1]),
        p_value = as.numeric(result_anual[2]),
        Tipo = "Tendencia Anual"
      )
      
      resultados_tendencia = rbind(resultados_tendencia, res_anual)
      
      # . Tendencias estacionales
      espei_estacional = data_analisis %>% 
        mutate(Estacion = case_when(
          format(Fecha, "%m") %in% c("03", "04", "05", "09", "10", "11") ~ "Temp_Humeda",
          format(Fecha, "%m") %in% c("01", "02", "06", "07", "08", "12") ~ "Temp_Seca"
        ))
      
      spei_estacional = espei_estacional %>%
        group_by(Año = as.numeric(format(Fecha, "%Y")), Estacion) %>%
        summarize(SPEI_estacional = median(SPEI, na.rm = TRUE))
      
      result_estacional = spei_estacional %>%
        split(.$Estacion) %>%
        map(~ MannKendall(.x$SPEI_estacional))
      
      res_humedo = data.frame(
        cod = i,
        Pixel = Pixel,
        tau = as.numeric(result_estacional[[1]][1]),
        p_value = as.numeric(result_estacional[[1]][2]),
        Tipo = "Tendencia Temporada_Humeda"
      )
      
      res_seco = data.frame(
        cod = i,
        Pixel = Pixel,
        tau = as.numeric(result_estacional[[2]][1]),
        p_value = as.numeric(result_estacional[[2]][2]),
        Tipo = "Tendencia Temporada_Seca"
      )
      resultados_tendencia = rbind(resultados_tendencia, res_humedo, res_seco)
      
      # pendientes
      result_m_humedo = zyp.sen(SPEI_estacional ~ Año, data = spei_estacional %>% filter(Estacion == "Temp_Humeda"))
      result_m_seco = zyp.sen(SPEI_estacional ~ Año, data = spei_estacional %>% filter(Estacion == "Temp_Seca"))
      
      res_humedo_pendiente = data.frame(
        cod = i,
        Pixel = Pixel,
        Tipo = "Estacional_humedo",
        pendiente = as.numeric(result_m_humedo$coefficients["Año"]),
        cambio_acumulado = result_m_humedo$coefficients["Año"] * 14
      )
      
      res_seco_pendiente = data.frame(
        cod = i,
        Pixel = Pixel,
        Tipo = "Estacional_seco",
        pendiente = as.numeric(result_m_seco$coefficients["Año"]),
        cambio_acumulado = result_m_seco$coefficients["Año"] * 14
      )
      
      resultados_pendiente = rbind(resultados_pendiente, res_humedo_pendiente, res_seco_pendiente)
      i = i + 1
    }
    
    save_tendencia = paste0(dir.save, "/Tendencia_SPEI")
    
    if (!dir.exists(save_tendencia)) {
      dir.create(save_tendencia)
    }
    
    # Guardo mis resultados 
    data_save_tendencia = merge(resultados_tendencia, cords, by = "Pixel")
    data_save_tendencia = data_save_tendencia %>% arrange(cod) %>% select(-cod)
    puntos = vect(data_save_tendencia, geom = c("X", "Y"), crs = "EPSG:32717")
    writeVector(puntos, paste0(save_tendencia, "/Tendencias_", Estat, ".shp", sep = ""), filetype = "ESRI Shapefile", overwrite = TRUE)
    
    data_save_pendiente = merge(resultados_pendiente, cords, by = "Pixel")
    data_save_pendiente = data_save_pendiente %>% arrange(cod) %>% select(-cod)
    puntos_m = vect(data_save_pendiente, geom = c("X", "Y"), crs = "EPSG:32717")
    writeVector(puntos_m, paste0(save_tendencia, "/Pendientes_", Estat, ".shp", sep = ""), filetype = "ESRI Shapefile", overwrite = TRUE)
    
    ############################################################################
    
  } else {
    ############################################################################
    resultados_tendencia = data.frame()
    resultados_pendiente = data.frame()
    names(data) = c("Fecha", "value")
    
    # Calculo la Global
    ssi_data = data$value
    ssi_ts = ts(ssi_data, start = c(2014, as.numeric(format(data$Fecha[1], "%j"))), frequency = 365)
    decomposed = stl(ssi_ts, s.window = "periodic")
    trend_component = decomposed$time.series[, "trend"]
    result = MannKendall(trend_component)
    
    
    res_global = data.frame(
      Cuenca = Estat,
      tau = as.numeric(result[1]),
      p_value = as.numeric(result[2]),
      Tipo = "Tendencia Global"
    )
    
    resultados_tendencia = rbind(resultados_tendencia, res_global)
    index <- seq_along(trend_component)
    sen_result_global = zyp.sen(trend_component ~ index)
    magnitud_cambio = sen_result_global$coefficients["index"] * 14
    
    res_pendienteGlobal = data.frame(
      Cuenca = Estat,
      Tipo = "Global",
      pendiente = as.numeric(sen_result_global$coefficients["index"]),
      cambio_acumulado = magnitud_cambio,
      intercepto = as.numeric(sen_result_global$coefficients["Intercept"])
    )
    
    resultados_pendiente = rbind(resultados_pendiente, res_pendienteGlobal)
    # Pendiente con theil
    # pendiente_global = zyp.sen(trend_component ~ time(trend_component, data = ssi_ts))
    # magnitud_cambio = pendiente_global$coefficients["time"] * 14
    
    # Escala anual
    ssi_anual = data %>% 
      mutate(Año = as.numeric(format(Fecha, "%Y"))) %>%
      group_by(Año) %>%
      summarize(SSI_anual = median(value, na.rm = TRUE))
    
    # Kendall
    result_anual <- MannKendall(ssi_anual$SSI_anual)
    
    # pendiente con theil
    pendiente_anual = zyp.sen(SSI_anual ~ Año, data = ssi_anual)
    magnitud_cambio = pendiente_anual$coefficients["Año"] * 14
    
    res_anualPendiente = data.frame(
      Cuenca = Estat,
      Tipo = "Anual",
      pendiente = as.numeric(pendiente_anual$coefficients["Año"]),
      cambio_acumulado = magnitud_cambio,
      intercepto = as.numeric(pendiente_anual$coefficients["Intercept"])
    )
    
    resultados_pendiente = rbind(resultados_pendiente, res_anualPendiente)
    
    res_anualTendencia = data.frame(
      Cuenca = Estat,
      tau = as.numeric(result_anual[1]),
      p_value = as.numeric(result_anual[2]),
      Tipo = "Tendencia Anual"
    )
    
    resultados_tendencia = rbind(resultados_tendencia, res_anualTendencia)
    
    # Tendencias estacionales
    ssi_estacional = data %>% 
      mutate(Estacion = case_when(
        format(Fecha, "%m") %in% c("03", "04", "05", "09", "10", "11") ~ "Temp_Humeda",
        format(Fecha, "%m") %in% c("01", "02", "06", "07", "08", "12") ~ "Temp_Seca"
      ))
    
    ssi_estacional = ssi_estacional %>%
      group_by(Año = as.numeric(format(Fecha, "%Y")), Estacion) %>%
      summarize(SSI_estacional = median(value, na.rm = TRUE))
    
    
    result_estacional = ssi_estacional %>%
      split(.$Estacion) %>%
      map(~ MannKendall(.x$SSI_estacional))
    
    res_humedo = data.frame(
      Cuenca = Estat,
      tau = as.numeric(result_estacional[[1]][1]),
      p_value = as.numeric(result_estacional[[1]][2]),
      Tipo = "Tendencia Temporada_Humeda"
    )
    
    res_seco = data.frame(
      Cuenca = Estat,
      tau = as.numeric(result_estacional[[2]][1]),
      p_value = as.numeric(result_estacional[[2]][2]),
      Tipo = "Tendencia Temporada_Seca"
    )
    
    resultados_tendencia = rbind(resultados_tendencia, res_humedo, res_seco)
    
    # pendientes
    result_m_humedo = zyp.sen(SSI_estacional ~ Año, data = ssi_estacional %>% filter(Estacion == "Temp_Humeda"))
    result_m_seco = zyp.sen(SSI_estacional ~ Año, data = ssi_estacional %>% filter(Estacion == "Temp_Seca"))
    
    
    res_humedo_pendiente = data.frame(
      Cuenca = Estat,
      Tipo = "Estacional_humedo",
      pendiente = as.numeric(result_m_humedo$coefficients["Año"]),
      cambio_acumulado = result_m_humedo$coefficients["Año"] * 14,
      intercepto = as.numeric(result_m_humedo$coefficients["Intercept"])
    )
    
    res_seco_pendiente = data.frame(
      Cuenca = Estat,
      Tipo = "Estacional_seco",
      pendiente = as.numeric(result_m_seco$coefficients["Año"]),
      cambio_acumulado = result_m_seco$coefficients["Año"] * 14,
      intercepto = as.numeric(result_m_seco$coefficients["Intercept"])
    )
  
    resultados_pendiente = rbind(resultados_pendiente, res_humedo_pendiente, res_seco_pendiente)
    
    rownames(resultados_tendencia) = NULL
    rownames(resultados_pendiente) = NULL
    
    # Guardo resultados
    save_tendencia = paste0(dir.save, "/Tendencia_SSI")
    if (!dir.exists(save_tendencia)) {
      dir.create(save_tendencia)
    }
    
    # guardo en csv
    write.csv(resultados_tendencia, paste0(save_tendencia, "/", Estat, "_Tendencia.csv", sep = ""), row.names = F)
    write.csv(resultados_pendiente, paste0(save_tendencia, "/", Estat, "_Pendiente.csv", sep = ""), row.names = F)
    
    ############################################################################
    # Gráfico tendencia global
    data_graph = data
    data_graph = data.table(data_graph)
    data_graph$Fecha = as.Date(data_graph$Fecha)
    data_graph[, Fecha_Num := as.numeric(Fecha - min(Fecha))]
    
    pendiente = resultados_pendiente[1, "pendiente"]
    intercepto = resultados_pendiente[1, "intercepto"]
    data_graph[, tendencia := pendiente * Fecha_Num + intercepto]
    
    
    pl = ggplot() +
      geom_line(data = data_graph, aes(
        x = Fecha, y = value), color = "blue", size = 1) + 
      
      geom_line(data = data_graph, aes(
        x = Fecha, y = tendencia), color = "black", size = 1.5) +
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
                color = "black", fill = NA, size = 1.5)
    
  pl
  
  # guardo mi gráfico
  
  ggsave(
    paste(save_tendencia, "/", Estat, "_Tendencia.png", sep = ""),
    plot = pl,
    width = 12,
    height = 8,
    units = "in",
    dpi = 2000,
    bg = NULL,
  )
  
  ################################################################################
  decomp_df <- data.frame(
    time = time(decomposed$time.series),
    seasonal = decomposed$time.series[, "seasonal"],
    trend = decomposed$time.series[, "trend"],
    remainder = decomposed$time.series[, "remainder"]
  ) %>%
    pivot_longer(cols = -time, names_to = "component", values_to = "value")
  

  pl3 <- ggplot(decomp_df, aes(x = time, y = value)) +
    geom_line() +
    facet_wrap(~component, scales = "free_y", ncol = 1) +
    labs(title = "Descomposición STL", x = "Tiempo", y = "Valor") +
    theme_minimal()
  
  print(pl3)
  ggsave(
    paste(save_tendencia, "/", Estat, "_DescomposicionSTL.png", sep = ""),
    plot = pl3,
    width = 12,
    height = 8,
    units = "in",
    dpi = 2000,
    bg = NULL,
  )
  ##############################################################################
  # Anual
  pl_anual = 
    ggplot(ssi_anual, aes(x = Año, y = SSI_anual)) +
    geom_line(color = "blue", size = 1) + # Puntos anuales
    geom_abline(intercept = coef(pendiente_anual)["Intercept"], 
                slope = coef(pendiente_anual)["Año"], color = "red") +
    labs(title = "Tendencia Anual de SSI", x = "Año", y = "SSI") +
    
    theme_minimal() +
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
    scale_x_continuous(breaks = seq(min(ssi_anual$Año), max(ssi_anual$Año), 1)) +
    geom_rect(aes(xmin = (-Inf), xmax = (Inf), ymin = -Inf, ymax = Inf),
              color = "black", fill = NA, size = 1.5)
  pl_anual

  ggsave(
    paste(save_tendencia, "/", Estat, "_TendenciaAnual.png", sep = ""),
    plot = pl_anual,
    width = 12,
    height = 8,
    units = "in",
    dpi = 2000,
    bg = NULL,
  )
  
  # Estacional
  gr = ssi_estacional
  gr$Estacion = ifelse(gr$Estacion == "Temp_Humeda", "Temporada Húmeda", "Temporada Seca")
  
  estaciones <- split(gr, gr$Estacion)
  
  graficas_estacionales <- map(estaciones, function(data_estacion) {
    pendiente_estacional <- zyp.sen(SSI_estacional ~ Año, data_estacion)
    name = unique(data_estacion$Estacion)
    pl = 
      ggplot(data_estacion, aes(x = Año, y = SSI_estacional)) +
      geom_line(color = "blue", size = 1) +
      geom_abline(intercept = coef(pendiente_estacional)["Intercept"], 
                  slope = coef(pendiente_estacional)["Año"], color = "red") +
      labs(title = paste("Tendencia Estacional -", unique(data_estacion$Estacion)),
           x = "Año", y = "SSI") +
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
      scale_x_continuous(breaks = seq(min(ssi_anual$Año), max(ssi_anual$Año), 1)) +
      geom_rect(aes(xmin = (-Inf), xmax = (Inf), ymin = -Inf, ymax = Inf),
                color = "black", fill = NA, size = 1.5)
    
    pl
    
    ggsave(
      paste(save_tendencia, "/", Estat, "_TendenciaEstacional_", name, ".png", sep = ""),
      plot = pl,
      width = 12,
      height = 8,
      units = "in",
      dpi = 2000,
      bg = NULL,
    )
  })
  
  
  graficas_estacionales[[1]] # Temporada Húmeda
  graficas_estacionales[[2]] # Temporada Seca
  
  # Guardo las gráficas
  
  
  ############################################################################
  # Puntos de cambio
  # breakpoints <- breakpoints(trend_component ~ time(trend_component))
  # 
  
  # pl2 = ggplot() +
  #   geom_line(data = data, aes(
  #     x = Fecha, y = value), color = "blue", size = 1) +
  #   
  #   geom_point(data = data[changepoints, ], aes(x = Fecha, y = value),
  #              color = "red", size = 2) +
  #   theme_bw() +
  # 
  #   labs(x = "Año", y = "SSI") +
  #   theme(
  #     axis.title = element_text(size = 20, face = "bold"),
  #     axis.text = element_text(size = 18, hjust = 0.9),
  #     legend.position = "bottom",
  #     legend.justification = c(0.5, 0.5),
  #     legend.background = element_rect(fill = "white", colour = "black", linewidth = 1),
  #     plot.title = element_text(hjust = 0.5, size = rel(2), lineheight = .9,
  #                               face = "bold", colour = "brown"),
  #   ) +
  #   scale_x_date(
  #     date_breaks = "1 year",
  #     labels = function(x) {
  #       as.numeric(format(x, "%Y"))
  #     }
  #   ) +
  # 
  #   geom_rect(aes(xmin = as.Date(-Inf), xmax = as.Date(Inf), ymin = -Inf, ymax = Inf),
  #             color = "black", fill = NA, size = 1.5)
  # pl2
  # 
  # iteractivo = ggplotly(pl2)
  # iteractivo
  # 
  # name2 = paste0(dir.save, "/PuntosCambio_SSI")
  # if (!dir.exists(name2)) {
  #   dir.create(name2)
  # }
  # 
  # ggsave(
  #   paste(name2, "/", Estat, "_PuntosCambio.png", sep = ""),
  #   plot = pl2,
  #   width = 12,
  #   height = 8,
  #   units = "in",
  #   dpi = 2000,
  #   bg = NULL,
  # )
  # 
  # # guardo el plotly
  # name3 = paste0(name2, "/", Estat, "_PuntosCambio.html")
  # saveWidget(iteractivo, file = name3)
  # lAS SEQUIAS HAN IDO AUMENTANDO ?
  } # CIERRO LOGICA DEL ELSE  
}


drought_MH = function(dir.spei, SSI) {
  
  files = list.files(dir.spei, full.names = TRUE)
  
  
  
  data = lapply(files, function(x) {
    fread(x, header = T)
  })
  
  ######### MODIFCIAR SECUENCIA DE ESTO ###### 
  if (Subcuenca == "Tomebamba") {
    print("Se ha importado los datos de SPEI para Tomebamba")
    data = data[, .(Fecha, Pixel_1, Pixel_2, Pixel_3, Pixel_4, Pixel_5, Pixel_6, Pixel_7, Pixel_8, Pixel_9,
                    Pixel_10, Pixel_11, Pixel_12, Pixel_13, Pixel_14, Pixel_15, Pixel_16, Pixel_17,
                    Pixel_18,Pixel_19, Pixel_20, Pixel_21, Pixel_23, Pixel_24, Pixel_25, Pixel_26,
                    Pixel_27, Pixel_28, Pixel_29, Pixel_30, Pixel_38), ]
  } else {
    print("Se ha importado los datos de SPEI para Yanuncay")
    data = data[, .(Fecha, Pixel_8, Pixel_9, Pixel_14, Pixel_15,Pixel_16,Pixel_17,Pixel_18,Pixel_19,Pixel_22, Pixel_23, 
                    Pixel_24, Pixel_25, Pixel_26, Pixel_27, Pixel_28,Pixel_31,
                    Pixel_32, Pixel_33, Pixel_34, Pixel_35, Pixel_36, Pixel_39, Pixel_40, Pixel_41, Pixel_42, 
                    Pixel_44, Pixel_45, Pixel_46),]
  }
  
  names(data)[2:length(data)] = paste0("Pixel_", 1:length(data))
  ###########################################
  
  names(data) = gsub(".*_(\\d+_dias)\\.csv$", "\\1", files)
  
  resultados = list()
  for (i in 1:length(data)) {
    
    nombre_dia = names(data)[i] 
    
    correlaciones_df = data.frame(Pixel = colnames(data[[i]])[-1])  # Excluye "Fecha"
    
    # Calcular la correlación para cada pixel
    for (j in 2:ncol(data[[i]])) {  # Desde la segunda columna porque la primera es "Fecha"
      
      datos = data[[i]][, .(Fecha, Valor = data[[i]][[j]])]  # Selección dinámica de pixel
      datos = merge(datos, SSI, by = "Fecha", all.x = TRUE)  # Unir con SSI por "Fecha"
      
      # Calcular correlación
      correlacion = cor(datos[,-1], use = "complete.obs", method = "pearson")
      correlacion_valor = correlacion[1, 2]  # Extraer el valor específico de la matriz
      
      correlaciones_df[j - 1, "Correlacion"] = correlacion_valor
      
    }
    
    correlaciones_df$Dias = nombre_dia
    
    resultados[[nombre_dia]] = correlaciones_df
  }
  resultados_df = rbindlist(resultados, idcol = "Archivo")
  
  dias_mayor_correlacion <- resultados_df[, .SD[Correlacion == max(Correlacion)], by = Pixel]
  dias_mayor_correlacion$Pixel <- factor(dias_mayor_correlacion$Pixel, 
                                         levels = unique(dias_mayor_correlacion$Pixel[order(as.numeric(gsub("Pixel_", "", dias_mayor_correlacion$Pixel)))]))
  
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

}



# DEBUG
PIXEL = data[['6_dias']][,.(Fecha = Fecha, Valor = get('Pixel_21'))]
res  = merge(PIXEL, SSI, by = "Fecha", all.x = TRUE)
cor(res[,-1], use = "complete.obs", method = "pearson")
