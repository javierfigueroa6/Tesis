# Diversidad Beta y Log2 Fold Change

#Este script calcula la diversidad beta del estudio y Log2 Fold Change.



#---------- DIVERSIDAD BETA ----------#


# Calculamos las divergencias para cancer y no cancer
divB_dec.can <- divergence(subset_samples(ps_decipher_div, status == "1"), apply(abundances(subset_samples(ps_decipher_div, status == "1")), 1, median))
divB_dec.noc <- divergence(subset_samples(ps_decipher_div, status == "2"), apply(abundances(subset_samples(ps_decipher_div, status == "2")), 1, median))
# transformamos el resultado anterior en _dataframes_
data.frame(divB_dec.can) -> df.divB_dec.can
data.frame(divB_dec.noc) -> df.divB_dec.noc
# Gather
#library(tidyverse)
df.divB_dec.can <- gather(df.divB_dec.can, sample, divergence)
df.divB_dec.noc <- gather(df.divB_dec.noc, sample, divergence)

# Agregamos columnas a nuestros dataframes
mutate(df.divB_dec.can, status = "1") -> df.divB_dec.can
mutate(df.divB_dec.noc, status = "2") -> df.divB_dec.noc

# Cambiamos los nombres de las columans de manera que sean iguales an ambos dataframes
colnames(df.divB_dec.can) <- c("samples", "divergence", "status")
colnames(df.divB_dec.noc) <- c("samples", "divergence", "status")

df.divB_dec.can
df.divB_dec.noc

# Los combinamos en un dataframe
rbind(df.divB_dec.can, df.divB_dec.noc) -> divB_dec.boxplot
divB_dec.boxplot

#library(ggpubr)
# Y finalmente graficamos y realizamos una comparación estadística
#Boxplot 1
Plot_divB_dec1 <- ggboxplot(data = divB_dec.boxplot, x = "status", y = "divergence", fill = "status", palette = c("#a6cee3", "#fdbf6f"))
Plot_divB_dec1 + stat_compare_means(comparisons = list(c("1", "2")))

#Boxplot 2
Plot_divB_dec2 <- ggplot(data = divB_dec.boxplot, aes(x = status, y = divergence)) +
  geom_boxplot(aes(fill = status), alpha = 0.75) +
  labs(x = "Status", y = "Divergence", title = "Divergence by Status") +
  scale_fill_manual(values = c("#a6cee3", "#fdbf6f"))
Plot_divB_dec2




#---------- Log2 Fold Change ----------#

#Transformamos la data para calcular Log2 Fold Change con la librería DESeq.
pcap.ps_dec <- ps_decipher_div
pcap.data_dec <- sample.data.table(pcap.ps_dec)
pcap.asvs_dec <- asv.data.table(pcap.ps_dec)
pcap.genera.ids_dec <- taxa_names(pcap.ps_dec)
pcap.dt_dec <- pcap.data_dec[pcap.asvs_dec]
pcap.dt_dec$status <- ifelse(pcap.dt_dec$status == 1, "cancer", "no cancer")

pcap.data_dec
pcap.asvs_dec
pcap.genera.ids_dec
pcap.dt_dec

pcap.ps_dec@sam_data$status <- factor(pcap.ps_dec@sam_data$status, levels = c(1,2))
pcap.ps_dec@sam_data$status
pcap.deseq_dec <- phyloseq_to_deseq2(pcap.ps_dec, ~ status)

#Función para la media geométrica.
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#Aplicamos la función de media geométrica
geoMeans_dec = apply(counts(pcap.deseq_dec),1, gm_mean)
#Se calculan los factores de tamaño, para normalizar los tamaños de muestras.
pcap.deseq_def_dec = estimateSizeFactors(pcap.deseq_dec, geoMeans = geoMeans_dec)
#Se realiza el análisis de expresión diferencial.
pcap.deseq_def_dec = DESeq(pcap.deseq_def_dec, fitType="local")
pcap.deseq_def_dec

#Para observar los p-values en cada bacteria.
pcap.res.dt_dec <- results(pcap.deseq_def_dec) %>% as.data.table(keep.rownames = "Secuencia")
pcap.res.dt_dec
pcap.sig.dt_dec <- pcap.res.dt_dec[pvalue <= 0.1]
gt(pcap.sig.dt_dec[order(log2FoldChange)]) %>% fmt_number(columns = 2:ncol(pcap.sig.dt_dec), decimals = 5)

#Grafico log2foldchange
pcap.sig.dt_dec
ggplot(pcap.sig.dt_dec, aes(x = reorder(Secuencia, log2FoldChange), y = log2FoldChange)) + 
  geom_col() + 
  coord_flip() + 
  labs(x = "ASV", y = "Log2 Fold Change (CAN to NOC)")


#Editamos el dataset para realizar el analisis por nivel taxonómico
res = results(pcap.deseq_def_dec)
res = res[order(res$pvalue, na.last=NA), ]
alpha = 0.1
sigtab = res[(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(pcap.ps_dec)[rownames(sigtab), ], "matrix"))
sigtab

theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(genus))
# Para ordenar a nivel de phylum
x = tapply(sigtabgen$log2FoldChange, sigtabgen$phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$phylum = factor(as.character(sigtabgen$phylum), levels=names(x))
# Para ordenar a nivel de genus
x = tapply(sigtabgen$log2FoldChange, sigtabgen$genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$genus = factor(as.character(sigtabgen$genus), levels=names(x))

#Gra´fico Log2 Fold Change
ggplot(sigtabgen, aes(y=genus, x=log2FoldChange, color=phylum)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6)+
  ggtitle("Cambio en la Abundancia Relativa de Géneros\n(Log2FoldChange)") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=14),
        plot.title = element_text(size=16, hjust = 0.5))








