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
# Cargo datos
SPEI_Tomebamba = fread("C:/Users/Jonna/Desktop/Sequias_GIZ/Indices_Sequia/SPEI/Results/SPEI_YanunTomeb.csv")
SSI_Tomebamba = fread("C:/Users/Jonna/Desktop/Sequias_GIZ/Indices_Sequia/SPEI/Tomebamba/Z-FINAL.csv")

Fecha_inicio = as.Date("2014-08-31", tz = "UTC")
Fecha_fin = as.Date("2024-05-17", tz = "UTC")

# preparo los datos
SPEI_Tomebamba = env$prep_data(SPEI_Tomebamba, SPEI = T, Fecha_inicio = Fecha_inicio, Fecha_fin = Fecha_fin)
SSI_Tomebamba = env$prep_data(SSI_Tomebamba, Fecha_inicio = Fecha_inicio, Fecha_fin = Fecha_fin)

# Categorizo mis sequias
SPEI_Cat_Tomebamba = env$categorizar_Sequias(SPEI_Tomebamba, type = "SPEI")
SSI_Cat_Tomebamba = env$categorizar_Sequias(SSI_Tomebamba, type = "SSI")

# Genero mapas de frecuencia de eventos en base al SPEI
#cat_SPEI = env$categorizar_SPEI(SPEI_Cat_Tomebamba)
#cat_SSI = env$categorizar_SSI(SSI_Cat_Tomebamba)

################################################################################
# Identificacion de sequias
run_SSI = env$teori_run(SSI_Tomebamba)
grouping = env$drought_grouping(run_SSI, tc = 5, pc = 0.1)
resultados_df <- rbindlist(grouping, fill = TRUE)

# Crear el gráfico
run_SSI$Year <- format(run_SSI$Fecha_Inicio, "%Y")
p1 = ggplot(run_SSI, aes(x = Fecha_Inicio, y = Severidad, 
                         width = Duracion, height = Severidad, 
                         fill = Tipo)) +
  
  geom_tile(color = "black", alpha = 0.5) +  # Ajustamos el color y la transparencia
  facet_wrap(~ Year, scales = "free_x") +  # Facet por año
  labs(x = "Fecha", y = "Severidad", 
       title = "Eventos de Sequías y Humedad por Año con Duración y Severidad") +
  theme_minimal() +
  scale_fill_manual(values = c("Sequia" = "red", "Humedo" = "blue")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(plotly)
p1 <- ggplotly(p1)
p1
