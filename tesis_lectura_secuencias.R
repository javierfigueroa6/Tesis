# Lectura de Secuencias (fastq)
# ADAPTADO DE https://benjjneb.github.io/dada2/tutorial.html

#Algunas pasos de este script necesitan algunos requerimientos mínimos para ejecutar.
#Recomendaría ejecutar este script con una máquina que tenga 16GB RAM mínimo.
#Funciones del modelo dada2 podrían tardar hasta 30-40 minutos en ejecutar.



#---------- Lectura de archivos y análisis de calidad ----------#



#Lecutra de archivos fastq.
path_fastq_files  <- "D:/Descargas/fastq_maquina" #Poner la ruta del archivo con las secuencias
list.files(path_fastq_files)

#Para manipular los nombres de los archivos Forward y Reverse
#En este caso, los que tenían R1 correspondian a lectura Forward y R2 a Reverse
fnFs <- sort(list.files(path_fastq_files, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path_fastq_files, pattern = "_R2_001.fastq.gz", full.names = TRUE))
#Obtenemos los ID's de cada muestra.
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

#De aquí hasta que se ejecuta la función FilterandTrim(), se debería realizar un análisis en cuanto a la calidad
# de las secuencias. En mi caso, usé el software FASTQC (descargable por Internet) para definir un recorte estándar
# en todas mis muestras y no realizar recortes específicos para cada lectura.

#De igual manera,se muestra una de las funciones de la librería dada2 para análisis de calidad de secuencias.

#Con la función plotQualityProfile() podemos inspeccionar la calidad.
# Recordar que fnFs tiene las lecturas forward y fnRs las lecturas Reverse.
plotQualityProfile(fnFs[1])
plotQualityProfile(fnRs[1])

# Dato GRAFICOS plotQualityProfile
#Línea verde -> puntuación de calidad media en cada posición (bp).
#Línea naranja -> cuartiles de la distribución de la puntuación de calidad.
#Línea roja -> proporcion escala de lecturas que se extienden al menos hasta esa posición 
#(Para Illumina no sirve mucho la línea roja).

#Las lecturas Reverse son de una calidad peor. Pero DADA2 incorpora información
#de calidad en su modelo de error, lo que hace que el algoritmo sea robusto
#para secuencias de menor calidad.

# En caso de que se hayan modificado o recortado las secuencias con una herramienta externa,
# dejarlo en el mismo formato que se estaba trabajando.

#Ya que se utilizó la función de recorte de dada2, se crean carpetas para guardar posteriormente las secuencias
# recortadas.
filtFs <- file.path("fastq_filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("fastq_filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# FUNCIÓN DE TRUNCAMIENTO
truncamiento_dada <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                                   maxN=0, maxEE= c(2,5), truncQ=20, rm.phix=TRUE, trimLeft= c(10,10),
                                   compress=TRUE, multithread=FALSE) # En Windows, multithread=FALSE
truncamiento_dada
# Esta función se iteró muchas veces ya que varía mucho dependiendo de los parámetros que se definen.
# Los que más afectan son trimLeft() que indica las bases a recortar desde la izquierda, en este caso,
# las primeras 10 bp de cada lecturas (Forward y Reverse), además se agrega un umbral de calidad en
# trunQ.





#----------      ALGORITMO DADA2 (MODELO DE ERROR Y ELIMINACIÓN DE RUIDO)       ----------#


#Se aprende el error de secuenciación en las lecturas Forward y Reverse.
errorsF <- learnErrors(filtFs, multithread=TRUE) 
errorsR <- learnErrors(filtRs, multithread=TRUE)

#Para mostrar la tasa de error en el cambio de nucleótidos.
plotErrors(errorsF, nominalQ=TRUE)
plotErrors(errorsR, nominalQ=TRUE)
#Podemos ver que la tasa de error va disminuyendo a medida que aumenta el puntaje de calidad.

#Se aplica el algoritmo dada (eliminación de ruido e inferencia de ASV's)
dadaFs <- dada(filtFs, err=errorsF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errorsR, multithread=TRUE)

#Se puede usar el siguiente comando para ver el resumen del proceso dada.
dadaFs[[1]]
dadaRs[[1]]

#Se mezclan las secuencias Forward y Reverse de cada muestra teniendo en consideración el resultado dada.
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
mergers
seqtab <- makeSequenceTable(mergers)
seqtab
table(nchar(getSequences(seqtab)))


#Por último se quitan las quimeras.
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#De esta forma, tenemos las secuencias únicas para cada muestra.


#Se crea una función para ver cuantas lecturas se eliminaron en cada paso.
#Es normal que un porcentaje grande de las lecturas se hayan eliminado en algún paso del flujo,
#sin embargo hay que tratar de que el conjunto sea representativo.

getN <- function(x) sum(getUniques(x))
track2 <- as.data.frame(cbind(truncamiento_dada, sapply(dadaFs, getN), sapply(mergers, getN),
                              rowSums(seqtab), rowSums(seqtab.nochim)))

colnames(track2) <- c('input', 'filtered', 'denoised', 'merged', 'tabled', 'nonchim')
rownames(track2) <- sample.names

if (!is.data.frame(track2)) {
  obj <- data.frame(track2)
}

track2$percent <- signif(((track2$nonchim)/(track2$input))*100, 4)
precent_desc <- track2[order(-track2$percent),]
precent_desc








