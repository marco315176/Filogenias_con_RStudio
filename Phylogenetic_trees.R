#!/usr/bin/Rscript

library(readxl)
library(dplyr)
library(ape)
library(phytools)
library(ggplot2)
library(ggtree)
library(ggnewscale)
library(tidytree)
library(treeio)
library(phangorn)
library(TDbook)
library(dplyr)
library(ggnewscale)
library(RColorBrewer)
library(reshape2)
library(ggExtra)
library(ggpmisc)
library(ggtree)
library(ggtreeExtra)
library(tidyr)

install.packages("tidyr")

##Instalación de paquetes##
install.packages("devtools")

install.packages("ggtreeExtra")

install.packages("Descargas/ggtreeExtra-master/", repos = NULL, type = "source")

remotes::install_local("Descargas/ggtreeExtra-master/", force = TRUE)


#Archivo para arbol filogenético
rab_cenasa <- read.tree("Filogenias/Quiropteros.fa.treefile")

#DF de datos de las muestras
sample_data <- read_xlsx("Filogenias/BASE_DE_DATOS_RABIA_CENASA_RITA_2024.xlsx")

sample_data$Año <- as.character(sample_data$Año)

#Info. para el hm
Country <- data.frame("estado" = sample_data$Estado)
rownames(Country) <- sample_data$Muestra

#1
phylo_rab <-  ggtree(rab_cenasa, 
        aes(color = Especie),
       layout = "roundrect",  
       #layout = "circular",
       branch.length = 'none',
        size = 1) %<+% sample_data +
   xlim(-10, NA) +
  geom_tiplab(color = "black", size = 4, align = TRUE, offset = 1) +
  theme(axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         plot.title = element_text(size = 15,
                                   face = "bold",
                                   hjust = 0.5,
                                   vjust = -15),
         legend.box = "vertical", legend.margin = margin()) +
         new_scale_fill() +
        # geom_tippoint(aes(fill=Año, shape = 21, size=4)) +
  #geom_text2(aes(subset=!isTip, label=label), size = 3.5, hjust = -.05, vjust = -0.5) + #Valor boostrap
  labs(title = "Análisis filogenético del virus de la rabia en mamíferos reportados en México (2020-2023)", caption = "Marco Hernández | Patricia Mora | \n| 2024 | Biología molecular y Secuenciación | CENASA ") 

  
phylo_rab
   
  hm <- gheatmap(phylo_rab, Country,
                 offset = 11,
                 width = 0.05,
                 colnames_angle = 100,
                 colnames_offset_y = .25,
                 colnames = FALSE) +
    theme_tree() +
    scale_fill_viridis_d(option = "H", name = "Estado") +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 15)) +
    geom_tiplab(aes(label = Año), color = "blue", offset = 13, size = 3.5, linetype = "blank", geom = "text") 
     
    
  hm  
