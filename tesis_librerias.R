# Librerías utilizadas

# Gestión de paquetes
library(BiocManager)

#Análisis de Secuencias
library(dada2)
library(DECIPHER); packageVersion("DECIPHER")

#Análisis Metagenómica
library(phyloseq)
library(metagenomeSeq)

#Estadísitica y visualización
library(DESeq2)
library(ggplot2)
library(vegan)
library(scales)

#Manipulación de data
library(dplyr)
library(data.table)

#formateo y otros
library(kableExtra)
library(knitr)
library(pander)
library(formattable)
library(gt)



# Si no funciona con la instalación por defecto usar lo siguiente:

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.16")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER")


