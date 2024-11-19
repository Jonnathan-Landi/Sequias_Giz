################################################################################
################################################################################
#                               ACTIVIDAD 2                                    #
################################################################################
packages = c("data.table", "dplyr", "purrr", "terra", "ggplot2", "plotly",
             "htmlwidgets", "gridExtra", "shiny", "shinythemes", "bslib", 
             "shinyWidgets", "lubridate", "tidyr")

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
categorizar_Sequias = function(data, type, dir.save, Estat) {
  
  # Logica para manejar sequias meteorologicas
  if (type == "SPEI") {
    
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
    
    Clasification = list()
    for (i in setdiff(names(data), "Fecha")) {
      temp_dt = data[, .(i = get(i), Categoria = cats_spei(.SD[[i]])), by = Fecha]
      setnames(temp_dt, "i", i)
      Clasification[[i]] = temp_dt
    }
    
    # Extrasigo informacion anual
    valor = unique(year(data$Fecha))
    informacion_anual = data.frame(
      Año = rep(valor, 5),
      Categoria = rep(c("No Sequia", "Sequia Leve", "Sequia Moderada", "Sequia Severa", "Sequia Extrema"), length(valor))
    )
    
    informacion_anual <- informacion_anual %>%
      arrange(Año)
    
    for (i in 1:length(Clasification)) {
      temp = Clasification[[i]]
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
    
    # pROMEDIO DE LOS DATOS
    informacion_anual$conteo = rowMeans(informacion_anual[, 3:length(informacion_anual)], na.rm = T)
    informacion_anual$conteo = round(informacion_anual$conteo, 0)

    # Genero mi grafico
    pl = ggplot(informacion_anual, aes(x = factor(Año), y = conteo, fill = Categoria)) +
      geom_bar(stat = "identity", position = "stack") + 
      labs(title = "Categorización de sequías anuales",
           x = "Año", y = "Frecuencia") +
      
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
    
    # Guardo mi grafico
    name_folder = paste0(dir.save, "/Categorias_SPEI")
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

    ############################################################################

    # Logica para manejar sequias hidro
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
    
    conteo_anual <- data_graph %>%
      group_by(Año, Categoria) %>%
      summarise(conteo = n()) %>%
      ungroup()
    
    pl = ggplot(conteo_anual, aes(x = factor(Año), y = conteo, fill = Categoria)) +
      geom_bar(stat = "identity", position = "stack") + 
      labs(title = "Categorización de sequías anuales",
           x = "Año", y = "Frecuencia") +
      
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
    
    # Guardo mi grafico
    name_folder = paste0(dir.save, "/Categorias_SSI")
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
teori_run = function(data, umbral = NULL){
  setDT(data)
  names(data)[2] = "value"
  data = data[!is.na(data$value), ]
  
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
    Magnitud = round(sum(value),3),   
    Severidad = round(sum(value) / .N, 3),  
    Fecha_Inicio = min(Fecha),  
    Fecha_Fin = max(Fecha)  
  ), by = drought_event]
  drought_seco = drought_seco[!is.na(drought_seco$drought_event), ]
  
  # Húmedo 
  droughts = data[drought == 0]
  drought_humedo <- droughts[, .(
    Tipo = "Humedo",
    Duracion = .N,  
    Magnitud = round(sum(value),3),   
    Severidad = round(sum(value) / .N, 3), 
    Fecha_Inicio = min(Fecha),  
    Fecha_Fin = max(Fecha)  
  ), by = drought_event]
  
  res_eventos = rbind(drought_seco, drought_humedo)
  res_eventos = res_eventos[order(res_eventos$drought_event), ]
  return(res_eventos)
} 

# Agrupamiento de sequías
drought_grouping = function(res_eventos, tc, pc, data_C, Estat, dir.save) {
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
        Severidad = evento_actual$Severidad + evento_siguiente$Severidad - vi,
        Magnitud = evento_actual$Magnitud + evento_siguiente$Magnitud
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
            Severidad = evento_actual$Severidad + evento_siguiente$Severidad - vi,
            Magnitud = evento_actual$Magnitud + evento_siguiente$Magnitud
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
  resultados_df = rbindlist(eventos_agrupados, fill = TRUE)
  resultados_df = unique(resultados_df, by = "Fecha_Fin")
  resultados_df = rbindlist(eventos_agrupados, fill = TRUE)
  resultados_df = unique(resultados_df, by = "Fecha_Fin")
  resultados_df$drought_event = NULL
  # resultados_df$Magnitud = resultados_df$Severidad / resultados_df$Duracion
  ##############################################################################
  # Genero un grafiquito
  data_graph <- resultados_df %>%
    mutate(Fecha_Inicio = as.Date(Fecha_Inicio),
           Fecha_Fin = as.Date(Fecha_Fin)) %>%
    mutate(Year = format(as.Date(Fecha_Fin), "%Y"))  # Extrae el año de la fecha
  
  data_graph <- data_graph %>%
    group_by(Year) %>%
    mutate(
      year_start = as.Date(paste0(Year, "-01-01")),
      year_end = as.Date(paste0(Year, "-12-31"))
    ) %>%
    ungroup()
  
  # Crear el gráfico con facetas por año
 graph = function(anio =NULL, f = F, Type = FALSE) {
   if (Type == T) {
     seq = seq.Date(as.Date(paste0(anio, "-01-01")), as.Date(paste0(anio, "-12-31")), by = "day")
     seq = data.table(Fecha_Inicio = seq)
     datos_graph = data_graph[data_graph$Year == anio, ]
   } else {
     datos_graph = data_graph
   }
   
   fechas_i = as.Date(c(datos_graph$Fecha_Inicio))
   fechas_f = as.Date(c(datos_graph$Fecha_Fin))
   Severidad = c(datos_graph$Severidad)
   
   # CRGO MI BASE DE DATOS DE EVENTOS DE LLUVIA
   names(data_C) = c("Fecha_Inicio", "value")
   if (Type == T) {
     data_C = merge(seq, data_C, by = "Fecha_Inicio", all = TRUE)
     data_C = data_C[data_C$Fecha >= as.Date(paste0(anio, "-01-01")) & data_C$Fecha <= as.Date(paste0(anio, "-12-30")),]
   }
   
   data_C$Tipo = rep("SSI", nrow(data_C))
   
   #data_f = merge(data_f, data_C, by = "Fecha_Inicio", all = TRUE)
   highlight_periods <- data.frame(
     xmin = fechas_i,
     xmax = fechas_f,
     ymin = -0.5,
     ymax = c(datos_graph$Severidad),
     Tipo = rep("Sequía", nrow(datos_graph))
   )
   #############################################################################
   pl = ggplot() +
     # Capa para las líneas de data_C
     geom_line(data = data_C, aes(
       x = Fecha_Inicio, y = value, color = Tipo,
       group = Tipo, # Asegura que las líneas no se combinen
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
                     "Severidad: ", ymax, "<br>",
                     "Tipo: ", Tipo)
     ), 
     fill = "white", linewidth = 0.9, alpha = 0.3, linetype = "dashed", inherit.aes = FALSE) +
     scale_color_manual(values = c("SSI" = "blue", "Sequía" = "red")) +
     geom_hline(yintercept = -0.5, linetype = "dashed", color = "black", linewidth = 1) +
     geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", linewidth = 1) +
     theme_bw() +
     labs(title = ifelse(Type == TRUE, anio, ""),
            x = "Fecha", y = ifelse(f == TRUE, "Severidad", "")) +
     # Ajustar temas de texto y ejes
     theme(
       axis.title = element_text(size = 20, face = "bold"),
       axis.text = element_text(size = 18, hjust = 0.9),
       # legend.title = element_text(colour = "brown", face = "bold", size = 12),
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
     
  # print(pl)
   return(pl)
 
   # pl_interactive
   #############################################################################
 } # CIerro funcion de grafico   pl_interactive = ggplotly(pl, tooltip = "text") anio =NULL, f = F, Type = FALSE
 P_1 = graph(2014, f = T, Type = T)
 P_2 = graph(2015, f = F, Type = T)
 P_3 = graph(2016, f = F, Type = T)
 P_4 = graph(2017,f = F, Type = T)
 P_5 = graph(2018,f = T, Type = T)
 P_6 = graph(2019, f = F, Type = T)
 P_7 = graph(2020, f = F, Type = T)
 P_8 = graph(2021, f = F, Type = T)
 P_9 = graph(2022, f = T, Type = T)
 P_10 = graph(2023, f = F, Type = T)
 P_11 = graph(2024, f = F, Type = T)
 
 combined_plot <- grid.arrange(
   grobs = list(P_1, P_2, P_3, P_4, P_5, P_6, P_7, P_8, P_9, P_10, P_11),
   ncol = 4,  # Número de columnas
   nrow = 3   # Número de filas
 )
 print(combined_plot)
 # CREO CARPETA PARA GUARDAR M,IS DATOS
 dir.save = paste0(dir.save, "/AgrupamientoSeq/")
 if (!dir.exists(dir.save)) {
   dir.create(dir.save)
 }
 
 ggsave(
   paste(dir.save, Estat, ".png", sep = ""),
   plot = combined_plot,
   width = 12,
   height = 8,
   units = "in",
   dpi = 2000,
   bg = NULL,
 )
 
 ###############################################################################
 # Crewo un grafico general de todo todo
 GENERAL = graph(f = T, Type = F)
 ggsave(
   paste(dir.save, Estat, "_General.png", sep = ""),
   plot = GENERAL,
   width = 12,
   height = 8,
   units = "in",
   dpi = 2000,
   bg = NULL,
 )
 
 # Genero el plotly
 pl_interactive = ggplotly(GENERAL, tooltip = "text")
 name = paste0(dir.save, "Sequia_Iteract", Estat, ".html")
 print(pl_interactive)
 saveWidget(pl_interactive, file = name)
 ###############################################################################
 # ui <- fluidPage(
 #   theme = bs_theme(bootswatch = "darkly"),
 #   
 #   # CSS mejorado
 #   tags$style(HTML("
 #    /* Contenedor principal */
 #    .container-fluid {
 #      padding: 0;
 #      margin: 0;
 #      width: 99%;
 #      height: 99vh;
 #    }
 #    
 #    /* Panel principal del gráfico */
 #    .main-panel {
 #      position: absolute;
 #      left: 0;
 #      width: calc(100% - 250px);
 #      height: 100vh;
 #      padding: 0;
 #      margin: 0;
 #    }
 #    
 #    /* Panel lateral (sidebar) */
 #    .sidebar-panel {
 #      position: fixed;
 #      right: 0;
 #      width: 250px;
 #      height: 100vh;
 #      padding: 15px;
 #      background-color: inherit;
 #      border-left: 1px solid rgba(255,255,255,0.1);
 #      overflow-y: auto;
 #    }
 #    
 #    /* Contenedor del gráfico */
 #    .plot-container {
 #      width: 100%;
 #      height: 100%;
 #      padding: 10px;
 #      margin: 10;
 #      
 #    }
 # 
 #    /* Estilo para el año en el título */
 #    .year-title {
 #      font-size: 24px;
 #      padding: 10px;
 #      margin-bottom: 20px;
 #      text-align: center;
 #    }
 # 
 #    /* Estilos para el modo oscuro/claro */
 #    .dark-mode {
 #      background-color: #222;
 #      color: #fff;
 #    }
 #    
 #    .light-mode {
 #      background-color: #fff;
 #      color: #222;
 #    }
 # 
 #    /* El gráfico siempre mantiene fondo blanco */
 #    .plot-background {
 #      background-color: white !important;
 #    }
 #  ")),
 #   
 #   # Layout principal usando sidebarLayout
 #   fluidRow(
 #     style = "margin: 10; height: 100vh;",
 #     
 #     # Panel principal con el gráfico (izquierda)
 #     div(
 #       class = "main-panel",
 #       div(
 #         class = "plot-container plot-background",
 #         conditionalPanel(
 #           condition = "input.interactive == false",
 #           plotOutput("general_plot", height = "100%", width = "100%")
 #         ),
 #         conditionalPanel(
 #           condition = "input.interactive == true",
 #           plotlyOutput("general_plotly", height = "100%", width = "100%")
 #         )
 #       )
 #     ),
 #     
 #     # Sidebar panel (derecha)
 #     div(
 #       class = "sidebar-panel",
 #       h4("Controles", style = "margin-bottom: 20px;"),
 #       selectInput("year_filter", 
 #                   "Selecciona el Año:",
 #                   choices = c("Todos" = "all", unique(data_graph$Year)),
 #                   selected = "all",
 #                   width = "100%"),
 #       materialSwitch(
 #         inputId = "interactive",
 #         label = "Modo interactivo",
 #         value = TRUE,
 #         status = "primary"
 #       ),
 #       materialSwitch(
 #         inputId = "toggle_theme",
 #         label = "Modo Claro/Oscuro",
 #         value = FALSE,
 #         status = "primary"
 #       )
 #     )
 #   )
 # )
 # 
 # server <- function(input, output, session) {
 #   # Tema
 #   observeEvent(input$toggle_theme, {
 #     if (input$toggle_theme) {
 #       session$setCurrentTheme(bs_theme(bootswatch = "darkly"))
 #     } else {
 #       session$setCurrentTheme(bs_theme(bootswatch = "cerulean"))
 #     }
 #   })
 #   
 #   # Datos filtrados
 #   filtered_graph <- reactive({
 #     if (input$year_filter == "all") {
 #       graph(f = TRUE, Type = FALSE)
 #     } else {
 #       graph(anio = as.numeric(input$year_filter), f = TRUE, Type = TRUE)
 #     }
 #   })
 #   
 #   # Gráfico estático
 #   output$general_plot <- renderPlot({
 #     filtered_graph() +
 #       theme(
 #         plot.margin = margin(10, 10, 10, 10),
 #         aspect.ratio = NULL,
 #         panel.background = element_rect(fill = "white"),
 #         plot.background = element_rect(fill = "white"),
 #         legend.position = "top",
 #         axis.text = element_text(size = 12),  # Aumentar tamaño del texto de los ejes
 #         axis.title = element_text(size = 14),  # Aumentar tamaño de títulos de ejes
 #         plot.title = element_text(size = 16, margin = margin(b = 20))  # Ajustar título
 #       )
 #   }, height = function() { session$clientData$output_general_plot_width * 0.6 })
 #   
 #   # Gráfico interactivo
 #   output$general_plotly <- renderPlotly({
 #     p <- ggplotly(filtered_graph()) %>%
 #       layout(
 #         autosize = TRUE,
 #         margin = list(l = 50, r = 20, b = 50, t = 20),
 #         paper_bgcolor = 'white',
 #         plot_bgcolor = 'white'
 #       )
 #     
 #     # Ajustar tamaño dinámicamente
 #     p %>% onRender("
 #      function(el) {
 #        var height = document.getElementsByClassName('plot-container')[0].offsetHeight;
 #        var width = document.getElementsByClassName('plot-container')[0].offsetWidth;
 #        Plotly.relayout(el, {
 #          'height': height,
 #          'width': width
 #        });
 #      }
 #    ")
 #   })
 # }
 # 
 # # Ejecutar la aplicación
 # shinyApp(ui = ui, server = server)
 ###############################################################################
 return(resultados_df)  # Agregar retorno explícito
}

#################################################################################
# Exclusión de eventos de sequia.
caracterizar_sequias = function(data, Type, tc = NULL, pc = NULL, Estat, dir.save) {
  ##############################################################################
  #                                SSI                                         #
  if (Type == "SSI") {
    
    run = teori_run(data)
    grouping = drought_grouping(run, tc, pc, data, Estat, dir.save)
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
    #                           SPEI                                           #
  } else {
    grouping = list()
    for (Pixel in setdiff(colnames(data), "Fecha")) {
      D = data[,.(Fecha = Fecha, SPEI = get(Pixel))]
      run = teori_run(D, umbral = - 0.5)
      
      # Exclusión de eventos de sequía
      agrup = run[(Tipo == "Sequia") & (Duracion >= 3),]
      
      # Categorizo mis sequias #################################################
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
      
       agrup$Categoria = sapply(agrup$Severidad, cats_spei)
      grouping[[Pixel]] = agrup
    }
  }
  
}
################################################################################

