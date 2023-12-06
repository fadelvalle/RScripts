if (!require("deldir")) {
  install.packages("deldir")
}
library(deldir)

if (!require("MASS")) {
  install.packages("MASS")
}
library(MASS)
if (!require("dbscan")) {
  install.packages("dbscan")
}
library(dbscan)

library(tidyverse)

# Lee el archivo CSV generado por TrackMate
file_path <- file.choose()
trackmate_data <- read_csv(file_path)
trackmate_data <- trackmate_data %>%
  slice(4:n())
# Convierte todas las columnas a tipo numérico
trackmate_data <- as.data.frame(sapply(trackmate_data, as.numeric))


# grafico en X e Y de distribución de tracks
ggplot(trackmate_data, aes(x = POSITION_X, y = POSITION_Y, color = as.factor(TRACK_ID))) +
  geom_line(aes(group = TRACK_ID)) +
  labs(title = "Trayectoria de Puntos en el Tiempo",
       x = "POSITION_X",
       y = "POSITION_Y") +
  theme(legend.position = "none")  # Esto oculta la leyenda de TRACK_ID

# Calcula la MSD por trayectoria
msd_data <- trackmate_data %>%
  group_by(TRACK_ID) %>%
  arrange(FRAME) %>%
  mutate(delta_x = POSITION_X - lag(POSITION_X),
         delta_y = POSITION_Y - lag(POSITION_Y),
         msd = delta_x^2 + delta_y^2) %>%
  summarize(mean_msd = mean(msd, na.rm = TRUE))



kde_result <- kde2d(trackmate_data$POSITION_X, trackmate_data$POSITION_Y)
filled.contour(kde_result, color.palette = terrain.colors)

ggplot(trackmate_data, aes(x = POSITION_X, y = POSITION_Y)) +
  geom_bin2d() +
  labs(title = "Histograma 2D de Distribución Espacial",
       x = "POSITION_X",
       y = "POSITION_Y")



# Ejemplo con DBSCAN
dbscan_result <- dbscan(trackmate_data[, c("POSITION_X", "POSITION_Y")], eps = 0.1, minPts = 5)
ggplot(trackmate_data, aes(x = POSITION_X, y = POSITION_Y, color = as.factor(dbscan_result$cluster))) +
  geom_point() +
  labs(title = "Análisis de Clusters",
       x = "POSITION_X",
       y = "POSITION_Y",
       color = "Cluster ID") +
  theme(legend.position = "none")  # Esto oculta la leyenda de TRACK_ID

# Calcula la velocidad por frame y track
trackmate_data <- trackmate_data %>%
  group_by(TRACK_ID) %>%
  arrange(FRAME) %>%
  mutate(delta_time = FRAME - lag(FRAME),
         delta_x = POSITION_X - lag(POSITION_X),
         delta_y = POSITION_Y - lag(POSITION_Y),
         distance = sqrt(delta_x^2 + delta_y^2),
         speed = distance / delta_time) %>%
  filter(!is.na(speed))

# Calcula la velocidad media por TRACK_ID
mean_speed_data <- trackmate_data %>%
  group_by(TRACK_ID) %>%
  summarize(mean_speed = mean(speed))

# Muestra la tabla con la velocidad media por TRACK_ID
print(mean_speed_data)

trackmate_data <- trackmate_data %>%
  group_by(TRACK_ID) %>%
  arrange(FRAME) %>%
  mutate(delta_x = POSITION_X - lag(POSITION_X),
         delta_y = POSITION_Y - lag(POSITION_Y),
         msd = cumsum(delta_x^2 + delta_y^2)) %>%
  ungroup()

# Gráfico de la velocidad media por TRACK_ID
ggplot(mean_speed_data, aes(x = TRACK_ID, y = mean_speed)) +
  geom_bar(stat = "identity") +
  labs(title = "Velocidad Media por TRACK_ID",
       x = "TRACK_ID",
       y = "Velocidad Media")

trackmate_data_persistence <- trackmate_data %>%
  group_by(TRACK_ID) %>%
  summarize(persistence = sum(distance) / sum(delta_time))

# Calcula la dirección media por trayectoria en grados
directionality_data <- trackmate_data %>%
  group_by(TRACK_ID) %>%
  arrange(FRAME) %>%
  summarize(mean_direction_deg = atan2(mean(delta_y), mean(delta_x)) * (180 / pi))

# Gráfico de Direccionalidad en grados
ggplot(directionality_data, aes(x = TRACK_ID, y = mean_direction_deg)) +
  geom_bar(stat = "identity", fill = "orange", alpha = 0.7) +
  labs(title = "Direccionalidad por Trayectoria",
       x = "TRACK_ID",
       y = "Direccionalidad (grados)") +
  theme_minimal()


