################################################################################
#   Análisis Integral de la Recurrencia y Respuesta Espacial de las Sequías    #
#        Meteorológicas e Hidrológicas en la Región de Cuenca, Ecuador         #
################################################################################
library(this.path)
dircs = this.path()
dir_Base = dirname(dircs)
dir_Base = paste0(dir_Base, "/functions_BETA.R")
env = new.env()
sys.source(dir_Base, envir = env)
################################################################################
# Datos generales
Fecha_Inicio = as.Date("2014-08-31", tz = "UTC")
Fecha_Fin = as.Date("2024-05-17", tz = "UTC")
dir.save = "C:/Users/Jonna/Desktop/Sequias_GIZ/Indices_Sequia/Análisis de sequias/Yanuncay"
coordenadas_pixeles = fread("C:/Users/Jonna/Desktop/Sequias_GIZ/Indices_Sequia/Yanuncay/coirds.csv")
Raster_Base = "C:/Users/Jonna/OneDrive/Proyectos/Bases de Datos/Interpolaciones/Precipitacion/diaria/Prec_5km_2014_2024.nc"
Estat = "Yanuncay"
################################################################################
#                   Sequías hidrológicas                                       #
################################################################################
SSI = fread("C:/Users/Jonna/Desktop/Sequias_GIZ/Indices_Sequia/SPEI/Tomebamba/Z-FINAL.csv")
SSI = env$import_data(SSI, type = "SSI", Fecha_Inicio = Fecha_Inicio, 
                      Fecha_Fin = Fecha_Fin, Subcuenca = "Yanuncay")

# Primera categorización de mis sequías
SSI_Categorizado = env$categorizar_sequias(data = SSI, type = "SSI", dir.save = dir.save, Estat = Estat, save = T) # rEIVSAR BUG

# Caracterizo mis sequías
SSI_Caracterizado = env$categorize_droughts(data = SSI, Type = "SSI", Estat = Estat,  tc = 5, pc = 0.1,
                                            dir.save = dir.save, save = T, Raster_Base = Raster_Base, cords = coordenadas_pixeles)

################################################################################
#                    Sequías meteorológicas                                    #
################################################################################
SPEI = fread("C:/Users/Jonna/Desktop/Sequias_GIZ/Indices_Sequia/SPEI/Results/SPEI_YanunTomeb3.csv")
SPEI = env$import_data(SPEI, type = "SPEI", Fecha_Inicio = Fecha_Inicio, 
                       Fecha_Fin = Fecha_Fin, Subcuenca = "Yanuncay")

# Primera categorización de mis sequías
SPEI_Categorizado = env$categorizar_sequias(data = SPEI, type = "SPEI", dir.save = dir.save, Estat = Estat, save = T)

# Caracterizo mis sequías
SPEI_Caracterizado = env$categorize_droughts(data = SPEI, Type = "SPEI", Estat = Estat, 
                                             dir.save = dir.save, save = T, Raster_Base = Raster_Base, cords = coordenadas_pixeles)
