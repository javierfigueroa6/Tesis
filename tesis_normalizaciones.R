# Normalizaciones

# A continuación se presentan algunas normalizaciones que se pueden utilizar.
# Lo ideal es aplicar la normalización o transformación según el objetivo.
# En mi caso, la que tuvo mejor rendimiento en cuanto a visualización gráfica de las muestras
# fue la de los porcentajes relativos (NORMALIZACIÓN 6)
# Además, se muestran gráficos para visualizar la data con cada normalización.
# Recorda que ps_decipher_div contiene toda nuestra información.



# Normalización 1
total2_dec = median(sample_sums(ps_decipher_div)) 
standf2_dec = function(x, t=20000) round(total2_dec * (x / sum(x)))
ps_dec_norm1 = transform_sample_counts(ps_decipher_div, standf2_dec)
ps_dec_norm1

plot_bar(ps_dec_norm1, fill="phylum")
plot_bar(ps_dec_norm1, fill="genus")
plot_bar(ps_dec_norm1, fill="class")
plot_bar(ps_dec_norm1, fill="family")
plot_bar(ps_dec_norm1, fill="order")
plot_bar(ps_dec_norm1, fill="species")


# NORMALIZACION 2 RAREFYING Y VISUALIZACION
summary(sample_sums(ps_decipher_div))
ps_dec_norm2 <- rarefy_even_depth(ps_prev, sample.size = 10,000, rngseed = 100)
ps_dec_norm2

plot_bar(ps_dec_norm2, fill ="Phylum")
plot_bar(ps_dec_norm2, fill ="Genus")
plot_bar(ps_dec_norm2, fill ="Class")
plot_bar(ps_dec_norm2, fill ="Family")
plot_bar(ps_dec_norm2, fill ="Order")


#NORMALIZACION 3 RELATIVE ABUNDANCE
ps_dec_norm3 <- transform_sample_counts(ps_decipher_div, function(x) x / sum(x))
ps_dec_norm3

plot_bar(ps_dec_norm3, fill ="phylum")
plot_bar(ps_dec_norm3, fill ="genus")
plot_bar(ps_dec_norm3, fill ="class")
plot_bar(ps_dec_norm3, fill ="family")
plot_bar(ps_dec_norm3, fill ="order")


#NORMALIZACION 6 (SUMA 1)
#Se utiliza transform(,"compositional) para la SUMA 1. Se muestra la abundancia relativa como porcentajes 
#para cada muestra.
pseq.compositional_norm6_dec <- transform(ps_decipher_div, "compositional") 
sample_names_norm6 <- rownames(otu_table(pseq.compositional_norm6_dec))
sample_names_norm6
sample_data(pseq.compositional_norm6_dec)$SampleNames <- sample_names_norm6
sample_data(pseq.compositional_norm6_dec)

