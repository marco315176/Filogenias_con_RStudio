#!/usr/bin/Rscript

library(readxl)
library(dplyr)
#install.packages("ggtreeExtra")
library(maps)
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

#####Rabia#####

#if (!require("BiocManager", quietly = TRUE))
 #install.packages("BiocManager")

#BiocManager::install("ggtreeExtra")

##### Archivo para arbol filogenético #####
rab_cenasa <- read.tree("/home/secuenciacion_cenasa/Analisis_corridas/kSNP4/IQ_tree/Aln_RABV_mafft.fa.treefile")

plot(rab_cenasa, type= "p", cex=0.8, edge.width=2, font=3, label.offset=0.0005, edge.lty=1, node.pos=2)

#DF de datos de las muestras
sample_data <- read_xlsx("/home/secuenciacion_cenasa/Analisis_corridas/kSNP4/Casos_de_virus.xlsx", sheet = "Rabia")

#sample_data <- sample_data[,c(1,2,4,5)]
#colnames(sample_data) <- c("ID", "Serotipo", "RS", "Estado")
str(sample_data)

#sample_data$Año <- as.character(sample_data$Año)

##### Info. para el hm #####
Country <- data.frame("Estado" = sample_data$Estado)
rownames(Country) <- sample_data$ID

Huesped <- data.frame("Huesped" = sample_data$Huesped)
rownames(Huesped) <- sample_data$ID


##### Info para el hm2 #####
Año <- data.frame("Año" = sample_data$Año)
rownames(Año) <- sample_data$ID

Año$Año <- as.character(Año$Año)

str(Año)

Majorclade <- data.frame("Clado_principal" = sample_data$Clado)
rownames(Majorclade) <- sample_data$ID

Minorclade <- data.frame("Clado_menor" = sample_data$Subclado)
rownames(Minorclade) <- sample_data$ID

#####  Personalización de paletas de colores #####

custom_palette <- c("#FF0000", "#00FF00", "#0000FF", "#800000", "#00FFFF", "#B0C4DE",
                    "#FFD700", "#FF1493", "#8B008B", "cadetblue4", "#8B8B00", "#4169E1",
                    "#00FA9A", "lightpink", "#FF4500", "#BA55D3", "#458B00", "#A0522D", "#CD6600", "#1E90FF",
                    "#00008B", "#EE3B3B", "#EEC591","#8b4513")


custom_palette2 <- c("black","orange2", "magenta4", "#1E90FF", "#2e8b57", "orangered2", "#4b0082",
                     "#006400", "royalblue4", "#8b4513",   "salmon1",
                     "magenta", "gray26")

custom_palette3 <- c("#00008B", "#EE3B3B", "#EEC591", "#009ACD", "#A2CD5A", "#2F4F4F", "#8DEEEE")

custom_palette4 <- c("#556B2F", "#9932CC")

custom_palette5 <- c("#FF8C00", "#1E90FF", "#8B7500", "#BA55D3")



######################Arbol con Raíz #########################################
rab_cenasa2 <- root(rab_cenasa, outgroup = "OU524430_ABLV", resolve.root = TRUE)

is.rooted(rab_cenasa2)

rab_cenasa2 <- ggtree(rab_cenasa, layout = "rect",  branch.length = 'none', size = .5) + 
  geom_rootedge(rootedge = 0.05, size = 2, color = "black")

phylo_rab2 <- rab_cenasa2 %<+% sample_data + 
  aes(color = Clado) +
  xlim(-8,  23) +
  #geom_tiplab(color = "black", size = 3, align = TRUE, offset = 0) + 
  geom_tiplab(aes(color = Subclado), size = 3.3, align = TRUE, offset = 0.1) +
  scale_color_manual(values = c("Bats" = "brown1",
                                "Cosmopolitan" = "cyan",
                                "Bats DR" = "black", 
                                "Bats TB1" = "green3", 
                                "Cosmopolitan Vac2" = "dodgerblue", 
                                "Desconocido" = "darkorchid2",
                                "RAC-SK" = "#B8860B",
                                "Indian-Sub" = "#8FBC8F",
                                "Cosmopolitan Vac" = "#FF7F00",
                                "ABLV" = "blue4",
                                "ABLV" = "#87CEFA",
                                "Artic" = "#2E8B57",
                                "Artic A" = "#4b0082",
                                "Africa-3" = "#009ACD",
                                "Africa-2" = "#EE4000")) + 
  #(color = Subclado, size = 3, align = TRUE, offset = 0.3) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 15,
                                  face = "bold",
                                  hjust = 0.5,
                                  vjust = -15),
        legend.box = "vertical", legend.margin = margin()) +
  geom_rootedge(rootedge = 1, size = 0.5, color = "black") + #Tamaño de raíz
  geom_text2(aes(subset=!isTip, label=label), size = 3, hjust = 0, vjust = 1, color="black")  #Valor boostrap +
#labs(title = "Análisis filogenético del gen de Glucoproteina de muestras de rabia del CENASA")#, caption = "Marco Hernández | Patricia Mora | \n| 2024 | Biología molecular y Secuenciación | CENASA ") 


phylo_rab2


######################Arbol sin raíz############################################
phylo_rab <-  ggplot(rab_cenasa, 
       #layout = "rect",  
       layout = "circular",
       #layout = "roundrect",
       branch.length = 'none',
        size = .5) %<+% sample_data +
    aes(color = Clado) +
     xlim(-15,  NA) +
  #geom_tiplab(color = "black", size = 3, align = TRUE, offset = 0) + 
  geom_tiplab(aes(color = Subclado), size = 2.8, align = TRUE, offset = 0.2) +
  scale_color_manual(values = c("Bats" = "brown1",
                                "Cosmopolitan" = "cyan",
                                "Bats DR" = "black", 
                                "Bats TB1" = "green3", 
                                "Cosmopolitan Vac2" = "dodgerblue", 
                                "Desconocido" = "darkorchid2",
                                "RAC-SK" = "#B8860B",
                                "Indian-Sub" = "#8FBC8F",
                                "Cosmopolitan Vac" = "#FF7F00",
                                "ABLV" = "blue4",
                                "ABLV" = "#87CEFA",
                                "Artic" = "#2E8B57",
                                "Artic A" = "#4b0082",
                                "Africa-3" = "#009ACD",
                                "Africa-2" = "#EE4000")) + 
#(color = Subclado, size = 3, align = TRUE, offset = 0.3) +
  theme(axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         plot.title = element_text(size = 15,
                                   face = "bold",
                                   hjust = 0.5,
                                   vjust = -15),
         legend.box = "vertical", legend.margin = margin()) +
  # geom_cladelabel(node=103, label="3 SNP's", #Anotación externa de nodos
  #                 color="red2", offset=16, align=TRUE) + 
  # geom_cladelabel(node=176, label="8 SNP's", 
  #                 color="blue", offset=2.4, align=TRUE) +
geom_text2(aes(subset=!isTip, label=label), size = 3, hjust = 0, vjust = 1, color="black")  #Valor boostrap 

phylo_rab

##### Configuración de HM ######

hm <- gheatmap(phylo_rab2, Country,
               offset = 2,
               width = 0.05,
               colnames_angle = 100,
               colnames_offset_y = .25,
               colnames = FALSE) +
  theme_tree() +
  scale_fill_manual(values = custom_palette) +
  #scale_fill_viridis_d(option = "H", name = "Estado") +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 15)) +
  guides(fill = guide_legend(title = "Estado", hjust = 0.5, order = 1)) +
  #theme(guide_legend = element_text(hjust = 0.5)) +
  theme(legend.position = c(0.9, 0.5)) +
  theme(  legend.title = element_text(size = 11, hjust = 0.5),
          legend.text = element_text(size = 10, hjust = 0.5))
#geom_tiplab(aes(label = "Especie"), color = "blue", offset = 6, size = 3.5, linetype = "blank", geom = "text") 

hm 

  h1 <- hm + new_scale_fill()  #Agregar nueva escala de colores
  
  ##### Segundo nivel #####
  hm2 <- gheatmap(h1, Huesped,
                  offset = 2.9,
                  width = 0.05,
                  colnames = FALSE) +
    scale_fill_manual(values = custom_palette2) +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 15)) +
    guides(fill = guide_legend(title = "Huesped", hjust = 0.5, order = 2)) +
    theme(legend.position = c(0.92, 0.5)) +
    theme(  legend.title = element_text(size = 10, hjust = 0.5),
            legend.text = element_text(size = 9, hjust = 0.5))
  
  hm2 

##### Rezaltar ramas de la filogenia #####

 hm2 + geom_highlight(node = c(103), fill = 'gold4', type  = "gradient", linetype = 3,
                             alpha = .7, extend=0.5) +
    geom_highlight(node = c(176), fill = 'red', type  = "gradient", linetype = 3,
                   alpha = .7, extend=0.5) +
    geom_highlight(node = c(130), fill= 'orange', type = "gradient", linetype = 3,
                   alpha = .7, extend=0.5, extendto = 0.5) +
    geom_highlight(node = c(149), fill= 'darkgreen', type = "gradient", linetype = 3,
                 alpha = .7, extend=0.5) +
  geom_highlight(node = c(150), fill= 'yellow', type = "gradient", gradient.direction = 'rt',
                 alpha = .7, extend=0.5) +
    geom_highlight(node = c(125), fill= 'maroon4', type = "gradient", linetype = 3,
                   alpha = .7, extend=0.5, extendto = 0.5)




