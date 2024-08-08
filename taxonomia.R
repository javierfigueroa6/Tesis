# Asignación de Taxonomía

# A continuación se mostrarán 2 formas de asignar taxonómia a las ASV's resultantes del algoritmo dada2.
# El primer método es con la librería DECIPHER, y el segundo con la librería assignTaxonomy() de la librería dada2.
# Para ejecuta ambos métodos, es necesario tener una base de datos de referencia (SILVA).



#---------- ASIGNACION CON DECIPHER ----------#

#Se almancenan las secuencias de cada ASV
dna <- DNAStringSet(getSequences(seqtab.nochim))
#Se carga la base de datos de referencia.
load("D:/Descargas/SILVA_SSU_r138_2019.RData")
#Función que realiza la asignación. 
ids_decipher <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE)
#Se definen los niveles de taxonomia. En mi caso, el algoritmo no pudo asignar especie.
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "specie")
#Se convierte el output para que sea una matriz análoga a la de assigntaxonomy().
taxid <- t(sapply(ids_decipher, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
taxid_dec <- addSpecies(taxid, path_especies)





#---------- ASIGNACION CON assignTaxonomy() ----------#

path_silva_database = "C:/Users/xavi_/silva_nr_v138_train_set.fa.gz"
path_especies = "C:/Users/xavi_/silva_species_assignment_v138.fa.gz"
#Para asignar taxonomías a las ASV.
taxa <- assignTaxonomy(seqtab.nochim, path_silva_database, multithread = TRUE)
taxa <- addSpecies(taxa, path_especies)
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)


###########################################################################
# DE AQUI EN ADELANTE SEGUIRÉ CON LOS RESULTADOS DE LA LIBRERÍA DECIPHER #
###########################################################################



#--------- CODIGO PARA ASIGNAR EL STATUS A CADA MUESTRA (CAN Y NOC) ---------#


#Se asignan 5 etiquetas que corresponden al tipo de muestra. Las de interés son 1(CAN) y 2(NOC).
sample_data1 <- data.frame(
  sample_id = sample.names,
  #dejamos cancer = 1, no cancer = 2, cn = 3, cp = 4 y Zymo = 5
  status = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,3,4,4,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,5)
)
sample_data1

sample_data1 <- sample_data1 %>% 
  tibble::column_to_rownames("sample_id")




#----------- PREPROCESAMIENTO DE DATOS -----------#

#Creamos el objeto phyloseq para usar algunas herramientas que nos pueden servir para el análisis.
#Además, es una de las herramientas estándar para análisis metagenómico.
#seqtab.nochim contiene nuestros ASV's.
ps_notprev <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                       tax_table(taxid_dec), sample_data(sample_data1))

#Ordenamos las ASV's en cuanto a la cantidad de lecturas que tienen. Lo ideal es filtrar muestras que tengan baja cantidad de lecturas.
names.ps_notprev <- as.data.frame(cbind(taxa_names(ps_notprev),
                                        paste0("ASV", seq(ntaxa(ps_notprev)))),
                                  stringsAsFactors=F)
colnames(names.ps_notprev) <- c('seq', 'asv')
taxa_names(ps_notprev) <- names.ps_notprev$asv
ps_notprev %>% sample_sums() %>% sort()
#Para ver las abundancias
taxa_sums(ps_notprev)




#Creamos ps_prev (NOS ASEGURAMOS QUE LA PREVALENCIA ESTE POR SOBRE EL 5,10 O 15%)
#Filtramos aquellas muestras que no tienen más de 10000 lecturas.
ps_notprev.pruned <- ps_notprev %>% prune_samples(sample_sums(ps_notprev) > 10000, .)
#Este filtro servía más para la taxonomía con assignTaxonomy(), pero igual se puede usar para asegurar que
# solo queremos ver bacterias (A veces se reconocen arqueas o otros microorganismos).
ps_notprev.sub <- ps_notprev.pruned %>%
  subset_taxa(
    domain %in% c('Bacteria') #&
    #!Genus %in% c('Campylobacter','Haemophilus','Halomonas','Helicobacter','Lactobacillus',
    #'Sphingomonas','Stomatobaculum','Streptococcus','TM7','Veillonella')
  )

#Definimos el umbral de prevalencia. Aquí está en 1%, pero se puede cambiar al valor que estime conveniente.
ps_notprev.sub
prevThreshold <- floor(nsamples(ps_notprev.sub) * 0.1)
prevdf <- apply(X = otu_table(ps_notprev.sub),
                MARGIN = ifelse(taxa_are_rows(ps_notprev.sub), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})

prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(ps_notprev.sub),
                     tax_table(ps_notprev.sub))

#Vemos la taxa a utilizar, y la taxa removida.
keepTaxa <- rownames(prevdf)[(prevdf$Prevalence >= prevThreshold)]
removeTaxa <- rownames(prevdf)[(prevdf$Prevalence < prevThreshold)]
keepTaxa
removeTaxa


# Creamos Phyloseq de nuevo, pero con las lecturas limpias.
ps_prev_decipher <- prune_taxa(keepTaxa, ps_notprev.sub)
ps_prev_decipher
tax_table(ps_prev_decipher)
tax_table(ps_prev)

#Filtramos para tener solo etiquetas con y sin cancer (1y2) en ps_prev.
status_decipher <- ps_prev_decipher@sam_data$status
keep_status <- c(1, 2)
ps_decipher_div <- prune_samples(status_decipher %in% keep_status, ps_prev_decipher)

#Nos aseguramos de que está todo en orden.
sam_data(ps_decipher_div)
ps_decipher_div
View(otu_table(ps_decipher_div))

#En caso de querer agrupar las bacterias por géneros (genus) utilizar:
ps_decipher_div_73 <- tax_glom(ps_decipher_div, "genus")
