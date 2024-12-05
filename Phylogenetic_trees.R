#!/usr/bin/Rscript

library(readxl)
library(dplyr)
#install.packages("ggnewscale")
library(ape)
library(phytools)
library(ggplot2)
library(ggtree)
library(ggnewscale)
library(tidytree)
library(treeio)
library(phangorn)
library(TDbook)
library(RColorBrewer)
library(reshape2)
library(ggExtra)
library(ggpmisc)
library(ggtreeExtra)
library(tidyr)

#Rabia#####
#install.packages("tidyr")

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("ggtree")

#Cargar el archivo de arbol filogenético
rab_cenasa <- read.tree("/home/secuenciacion_cenasa/Programas_Bioinformaticos/vSNP3/vsnp3_test_dataset/AF2122_test_files/step2/Mbovis-01/Mbovis-01_2024-10-23_10-52-59.tre")

#Visualizar el arbol en primera instancia
plot(rab_cenasa, type= "p", cex=0.8, edge.width=2, font=3, label.offset=0.0005, edge.lty=1, node.pos=2)

#DF de datos de las muestras
sample_data <- read_xlsx("Ruta/al/archivo.xlsx")
sample_data <- sample_data[,c(1,2,4,5)]
colnames(sample_data) <- c("ID", "Serotipo", "RS", "Estado")
str(sample_data)

#En caso que se quiera cambiar el tipo de datos
#sample_data$Año <- as.character(sample_data$Año)

#Info. para el hm
Country <- data.frame("Estado" = sample_data$Estado)
rownames(Country) <- sample_data$ID

#Info para el hm2
Año <- data.frame("Año" = sample_data$RS)
rownames(Año) <- sample_data$ID

Año$año <- as.character(Año$año)

str(Country)

#Personaliza una paleta de colores
custom_palette <- c("#FF0000", "#00FF00", "#0000FF", "#800000", "#00FFFF", "#B0C4DE",
                    "#FFD700", "#FF1493", "#8B008B", "cadetblue4", "#8B8B00", "#4169E1",
                    "#00FA9A", "lightpink", "#FF4500")

#Personaliza una segunda paleta de colores
custom_palette2 <- c("wheat4", "black","orange2", "#8b0000", "orangered2", "#4b0082",
                     "#006400", "royalblue4", "#8b4513", "#2e8b57", "magenta4", "salmon1",
                     "magenta", "gray26")

#1
phylo_rab <-  ggtree(rab_cenasa, 
       #layout = "rect",  
       layout = "circular",
       #layout = "roundrect",
       branch.length = 'none',
        size = 1) %<+% sample_data +
  #aes(color = Estado) +
     xlim(-5, 40) +
  geom_tiplab(color = "black", size = 4, align = TRUE, offset = 0.2) +
  theme(axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         plot.title = element_text(size = 15,
                                   face = "bold",
                                   hjust = 0.5,
                                   vjust = -15),
         legend.box = "vertical", legend.margin = margin()) #+
  #geom_tippoint(mapping = aes(color = Estado), size = 2.7)       # color de la punta por continente Puedes cambiar la forma añadiendo "shape = "
    
#  guides(fill = guide_legend(title = "Genes y mutaciones"))
    #     new_scale_fill() +
        # geom_tippoint(aes(fill=Año, shape = 21, size=4)) +
  #geom_text2(aes(subset=!isTip, label=label), size = 3.5, hjust = -.05, vjust = -0.5)  #Agregar el valor de boostrap +
  #labs(title = "Comparación filogenética del gen de Nucleoproteina de rabia \n aislada de muestras del CENASA")#, caption = "Marco Hernández | Patricia Mora | \n| 2024 | Dpto. de Biología Molecular | Área de Secuenciación Masiva y Bioinformática | CENASA ") 

#Visualizar el arbol
phylo_rab

#Agregar capas de color a las ramas
 phylo_rab + geom_highlight(node = 5, fill = 'red', type  = "gradient", linetype = 3,
                            alpha = .7, extend=0.5) +
#            geom_highlight(node = 12, fill = '#7FFF00', type  = "gradient", linetype = 3,
#                  alpha = .7, extend=0.5) 

hm <- gheatmap(phylo_rab, Country,
               offset = 4,
               width = 0.05,
               colnames_angle = 100,
               colnames_offset_y = .25,
               colnames = FALSE) +
  theme_tree() +
  scale_fill_manual(values = custom_palette) +
  #scale_fill_viridis_d(option = "H", name = "Estado") +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 15)) +
  guides(fill = guide_legend(title = "Estado", hjust = 0.5)) +
  #theme(guide_legend = element_text(hjust = 0.5)) +
  theme(legend.position = c(0.95, 0.5)) +
  theme(  legend.title = element_text(size = 12, hjust = 0.5),
          legend.text = element_text(size = 12, hjust = 0.5))
#geom_tiplab(aes(label = "Especie"), color = "blue", offset = 6, size = 3.5, linetype = "blank", geom = "text") 

hm 

#Agregar un segundo hm
  h1 <- hm + new_scale_fill() 
  
  ####sEGUNDO NIVEL
  hm2 <- gheatmap(h1, Año,
                  offset = 6,
                  width = 0.05,
                  colnames = FALSE) +
    scale_fill_manual(values = custom_palette2) +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 15)) +
    guides(fill = guide_legend(title = "RS", hjust = 0.5)) +
    #theme(legend.position = c(0.95, 0.5)) +
    theme(  legend.title = element_text(size = 12, hjust = 0.5),
            legend.text = element_text(size = 12, hjust = 0.5))
  
  hm2 
  
hm2 + geom_highlight(node = c(1,2,3,4,5,60,61,62), fill= 'pink', type = "gradient", linetype = 3,
                   alpha = .7, extend=0.5, extendto = 0.5) +
  geom_highlight(node = c(29,34,26,15,16,25,27,28,39,31,32,33,24,30,36), fill= 'steelblue', type = "gradient", linetype = 3,
                 alpha = .7, extend=0.5) +
  geom_highlight(node = c(44,46,47), fill= 'darkgreen', type = "gradient", gradient.direction = 'rt',
                 alpha = .7, extend=0.5) +
  geom_highlight(node = c(52,53), fill= 'yellow', type = "gradient", linetype = 3,
                       alpha = .7, extend=0.5, extendto = 0.5) +
  geom_highlight(node = c(49,48,50,51), fill= 'black', type = "gradient", linetype = 3,
                 alpha = .7, extend=0.5, extendto = 0.5)  +
  geom_highlight(node = c(54,55,56,57,58), fill= 'red', type = "gradient", linetype = 3,
                 alpha = .7, extend=0.5, extendto = 0.5)  +
  geom_highlight(node = c(20,19,17,15,16,18,21,22,23), fill= 'orange', type = "gradient", linetype = 3,
                 alpha = .7, extend=0.5, extendto = 0.5) +
  geom_highlight(node = c(7,8,9,10), fill= 'green', type = "gradient", linetype = 3,
                 alpha = .7, extend=0.5, extendto = 0.5)  +
  geom_highlight(node = c(11,12), fill= '#838B86', type = "gradient", linetype = 3,
                 alpha = .7, extend=0.5, extendto = 0.5)  +
  geom_highlight(node = c(39,35,36), fill= 'indianred', type = "gradient", linetype = 3,
                 alpha = .7, extend=0.5, extendto = 0.5)  +
  geom_highlight(node = c(13,14), fill= 'maroon4', type = "gradient", linetype = 3,
                 alpha = .7, extend=0.5, extendto = 0.5)  +
  geom_highlight(node = c(37,38,36,39), fill= 'palegreen', type = "gradient", linetype = 3,
                 alpha = .7, extend=0.5, extendto = 0.5) +
  geom_highlight(node = c(78), fill= 'snow', type = "gradient", linetype = 3,
                 alpha = .6, extend=0.5, extendto = 0.1) 
