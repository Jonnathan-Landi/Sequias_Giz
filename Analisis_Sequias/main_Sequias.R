########################## Análisis de Sequias #################################
# Versión:1.0.0                                                                #
# Fecha actualización: 2024-11-15                                              #
################################################################################
library(this.path)
dircs = this.path()
dir_Base = dirname(dircs)
dir_Base = paste0(dir_Base, "/functions_analisSeq.R")
env = new.env()
sys.source(dir_Base, envir = env)
################################################################################
dir.save = "C:/Users/Jonna/Desktop/Sequias_GIZ/Indices_Sequia"
Estat = "Tomebamba"
Fecha_inicio = as.Date("2014-08-31", tz = "UTC")
Fecha_fin = as.Date("2024-05-17", tz = "UTC")

################################################################################
#                           SEQUIAS HIDROLOGICAS                               #
################################################################################

# Preparo mis datos
SSI_Tomebamba = fread("C:/Users/Jonna/Desktop/Sequias_GIZ/Indices_Sequia/SPEI/Tomebamba/Z-FINAL.csv")
SSI_Tomebamba = env$prep_data(SSI_Tomebamba, Fecha_inicio = Fecha_inicio, Fecha_fin = Fecha_fin)

# Primera categorizacion de mis sequias
SSI_Cat_Tomebamba = env$categorizar_Sequias(SSI_Tomebamba, type = "SSI", dir.save, Estat)

# Caracterizo mis sequias
SSI_Carac_Tomebamba = env$caracterizar_sequias(SSI_Tomebamba, Type = "SSI", tc = 5, pc = 0.1, Estat, dir.save)

################################################################################
#                           SEQUÍAS METEOROLÓGICAS                             #
################################################################################
# Preparo mis datos
SPEI_Tomebamba = fread("C:/Users/Jonna/Desktop/Sequias_GIZ/Indices_Sequia/SPEI/Results/SPEI_YanunTomeb3.csv")
SPEI_Tomebamba = env$prep_data(SPEI_Tomebamba, SPEI = T, Fecha_inicio = Fecha_inicio, Fecha_fin = Fecha_fin)

# Primera categorizacion de mis sequias
SPEI_Cat_Tomebamba = env$categorizar_Sequias(SPEI_Tomebamba, type = "SPEI",dir.save, Estat)

# Caracterizo mis sequias
SPEI_Carac_Tomebamba = env$caracterizar_sequias(SPEI_Tomebamba, Type = "SPEI", Estat, dir.save)