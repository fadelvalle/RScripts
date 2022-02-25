library(readxl)
library(tidyverse)
library(ggplot2)
#import list of files as .csv pattern recognition and create a list of files 
file.list <- list.files(pattern='*.csv')
filelist <- as.data.frame(file.list)
df.list <- lapply(file.list, read_csv)
#put them into a DF
df <- bind_rows(df.list, .id = "id")
df_data <- as.data.frame(df)
df_data <- df_data %>% rename("Imgnum" = ...1)


grouped_by_chn<-df_data %>% group_by(Ch)
grouped_by_chn<- grouped_by_chn %>% filter(Ch==2|Ch==4)
grouped_by_chn<-as.data.frame(grouped_by_chn)

chn2 <- grouped_by_chn %>% filter(Ch==2)
chn4 <- grouped_by_chn %>% filter(Ch==4)


ratio_merge <- data.frame(id = character(length(chn2$id)), chn2 = numeric(length(chn2$id)), chn4 = numeric(length(chn2$id)), ratio = numeric(length(chn2$id)))
ratio_merge$id<-chn2$Label
ratio_merge$chn2<-chn2$Mean
ratio_merge$chn4<-chn4$Mean
ratio_merge <- ratio_merge %>% mutate(ratio = chn2/chn4)

xlsx::write.xlsx(x = grouped_by_chn, file = "resultados.xls", sheetName = "res", col.names = TRUE , row.names = FALSE, append = TRUE)

xlsx::write.xlsx(x = ratio_merge, file = "resultados_ratio.xls", sheetName = "res", col.names = TRUE , row.names = FALSE, append = TRUE)


xlsx::write.xlsx(x = chn2, file = "chn2.xls", sheetName = "res", col.names = TRUE , row.names = FALSE, append = TRUE)
xlsx::write.xlsx(x = chn4, file = "chn4.xls", sheetName = "res", col.names = TRUE , row.names = FALSE, append = TRUE)
xlsx::write.xlsx(x = filelist, file = "file_list.xls", sheetName = "res")


