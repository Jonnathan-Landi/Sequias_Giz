################################################################################
#   Análisis Integral de la Recurrencia y Respuesta Espacial de las Sequías    #
#        Meteorológicas e Hidrológicas en la Región de Cuenca, Ecuador         #
################################################################################
#' Repositorio en GitHub: https://github.com/Jonnathan-Landi/Sequias_Giz
#' Versión: 1.5.1
################################################################################
library(this.path)
dircs = this.path()
dir_Base = dirname(dircs)
dir_Base = paste0(dir_Base, "/functions_analisSeq.R")
env = new.env()
sys.source(dir_Base, envir = env)
################################################################################
# Datos generales
Fecha_Inicio = as.Date("2014-08-31", tz = "UTC")
Fecha_Fin = as.Date("2024-05-17", tz = "UTC")
dir.save = "C:/Users/Jonna/Desktop/Sequias_GIZ/Indices_Sequia/Análisis de sequias/Temp"
coordenadas_pixeles = fread("C:/Users/Jonna/Desktop/Sequias_GIZ/Indices_Sequia/Análisis de sequias/Tomebamba/data_crudo/cords.csv")
Raster_Base = "C:/Users/Jonna/OneDrive/Proyectos/Bases de Datos/Interpolaciones/Precipitacion/diaria/Prec_5km_2014_2024.nc"
Estat = "Yanuncay"
################################################################################
#                   Sequías hidrológicas                                       #
################################################################################
SSI = fread("C:/Users/Jonna/Desktop/Sequias_GIZ/Indices_Sequia/Análisis de sequias/Temp/Z-FINAL_Yanuncay.csv")
SSI = env$import_data(SSI, type = "SSI", Fecha_Inicio = Fecha_Inicio, 
                      Fecha_Fin = Fecha_Fin, Subcuenca = "Yanuncay")

# Primera categorización de mis sequías
SSI_Categorizado = env$categorizar_sequias(data = SSI, type = "SSI", dir.save = dir.save, Estat = Estat, save = T)

# Caracterizo mis sequías
SSI_Caracterizado = env$categorize_droughts(data = SSI, Type = "SSI", Estat = Estat,  tc = 5, pc = 0.1,
                                            dir.save = dir.save, save = T, Raster_Base = Raster_Base, cords = coordenadas_pixeles)

# análisis estadistico
SSI_Estadistico = env$statistical_analysis(data = SSI, Typo = "SSI",cords =  coordenadas_pixeles, Estat = Estat, dir.save = dir.save)
################################################################################
#                    Sequías meteorológicas                                    #
################################################################################
SPEI = fread("C:/Users/Jonna/Desktop/Sequias_GIZ/Indices_Sequia/Análisis de sequias/Tomebamba/data_crudo/SPEI_Tomebamba.csv")
#SPEI = env$import_data(SPEI, type = "SPEI", Fecha_Inicio = Fecha_Inicio, 
#                       Fecha_Fin = Fecha_Fin, Subcuenca = "Yanuncay", dir.save = dir.save)

# Primera categorización de mis sequías
SPEI_Categorizado = env$categorizar_sequias(data = SPEI, type = "SPEI", dir.save = dir.save, Estat = Estat, save = T)

# Caracterizo mis sequías
SPEI_Caracterizado = env$categorize_droughts(data = SPEI, Type = "SPEI", Estat = Estat, 
                                             dir.save = dir.save, save = T, Raster_Base = Raster_Base, cords = coordenadas_pixeles)

# analisis estadistico 
SPEI_Estadistico = env$statistical_analysis(data = SPEI, Typo = "SPEI",cords =  coordenadas_pixeles, Estat = Estat, dir.save = dir.save)

################################################################################
#                            Análisis mutuo                                    #
#' Solicitar dirrectorio donde estan los SPEIS diferentes intervalos
#dir.spei = "C:/Users/Jonna/Desktop/Sequias_GIZ/Indices_Sequia/Análisis de sequias/SPEI_Data"
