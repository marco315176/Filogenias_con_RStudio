library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

genomas <- read_xlsx("Proyectos/2024/DB_2024.xlsx")

Totales <- nrow(genomas)

sum_genom <- genomas %>%
  group_by(Genoma) %>%
  summarise(Sum = n(), .groups = 'drop',
            pct=(Sum / Totales))


ggplot(sum_genom, aes(x="", y= pct, fill = Genoma)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y") +
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 25, face = "bold"),
        text = element_text(size = 15)) +
  geom_text(aes(label = percent(pct, accuracy = 0.1)), position = position_stack(vjust = 0.5), 
            color = "black", size = 10, face = "bold") 
