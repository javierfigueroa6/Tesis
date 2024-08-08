# Graficos Composición Taxonómica

# Este script muestra la abundancia de los niveles taxonómicos en las muestras. Es decir, 6 gráficos
# dependiendo de cada nivel de taxonomía.
# Estos gráficos tienen como parametro de entrada pseq.compositional_norm6_dec que corresponde a la data
# que tiene la abundancia relativa de cada ASV en las muestras (porcentajes).



#Porcentaje taxa Phylum 
Phylum_prev_norm6_dec <- tax_glom(pseq.compositional_norm6_dec, taxrank = "phylum") 
tb_phyl <- psmelt(Phylum_prev_norm6_dec) 
tbphylum_prev_norm6 <- tb_phyl %>% group_by(phylum) %>% summarise(M=median(Abundance)) %>% ungroup() %>% unique() %>% top_n(M,n=6)
tbphylum_prev_norm6 <- tbphylum_prev_norm6 %>% inner_join(tb_phyl)
Plot_phylum_prev_norm6 <- ggplot(tbphylum_prev_norm6, aes(SampleNames, Abundance ,fill=phylum)) + geom_col(position="fill") # + scale_fill_manual(values=paleta5)
Plot_phylum_prev_norm6 <- Plot_phylum_prev_norm6 + theme_minimal() + xlab("SampleNames") + ylab("Relative Abundance") +
  theme(axis.text.x = element_text(color = "black", size = 7, angle = 45, hjust = .5, vjust = .5, face = "plain")) +
  theme(axis.text.y = element_text(color = "black", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain")) + 
  theme(legend.title=element_text(color= "black", size=15), legend.text=element_text(size=7, face="italic"))
Plot_phylum_prev_norm6 <- Plot_phylum_prev_norm6 +  ggtitle ("Abundancia relativa de phylums, con ajuste de prevalencia al 5% (ps_prev)") + 
  theme(plot.title = element_text(family="Arial", size=rel(2), vjust=2,face="plain", color="black", lineheight=1.5))  
Plot_phylum_prev_norm6


Genus_prev_norm6_dec <- tax_glom(pseq.compositional_norm6_dec, taxrank = "genus") 
tb_genus <- psmelt(Genus_prev_norm6_dec) 
tbgenus_prev_norm6 <- tb_genus %>% group_by(genus) %>% summarise(M=median(Abundance)) %>% ungroup() %>% unique() %>% top_n(M,n=20)
tbgenus_prev_norm6 <- tbgenus_prev_norm6 %>% inner_join(tb_genus)
Plot_genus_prev_norm6 <- ggplot(tbgenus_prev_norm6, aes(SampleNames, Abundance ,fill=genus)) + geom_col(position="fill") # + scale_fill_manual(values=paleta5)
Plot_genus_prev_norm6 <- Plot_genus_prev_norm6 + theme_minimal() + xlab("SampleNames") + ylab("Relative Abundance") +
  theme(axis.text.x = element_text(color = "black", size = 7, angle = 45, hjust = .5, vjust = .5, face = "plain")) +
  theme(axis.text.y = element_text(color = "black", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain")) + 
  theme(legend.title=element_text(color= "black", size=12), legend.text=element_text(size=7, face="italic"))
Plot_genus_prev_norm6 <- Plot_genus_prev_norm6 +  ggtitle ("Abundancia relativa de Género - Umbral de prevalencia al 5%") + 
  theme(plot.title = element_text(family="Arial", size=rel(2), vjust=2,face="plain", color="black", lineheight=1.5))  
Plot_genus_prev_norm6


Family_prev_norm6_dec <- tax_glom(pseq.compositional_norm6_dec, taxrank = "family") 
tb_family <- psmelt(Family_prev_norm6_dec) 
tbfamily_prev_norm6 <- tb_family %>% group_by(family) %>% summarise(M=median(Abundance)) %>% ungroup() %>% unique() %>% top_n(M,n=20)
tbfamily_prev_norm6 <- tbfamily_prev_norm6 %>% inner_join(tb_family)
Plot_family_prev_norm6 <- ggplot(tbfamily_prev_norm6, aes(SampleNames, Abundance ,fill=family)) + geom_col(position="fill") # + scale_fill_manual(values=paleta5)
Plot_family_prev_norm6 <- Plot_family_prev_norm6 + theme_minimal() + xlab("SampleNames") + ylab("Relative Abundance") +
  theme(axis.text.x = element_text(color = "black", size = 7, angle = 45, hjust = .5, vjust = .5, face = "plain")) +
  theme(axis.text.y = element_text(color = "black", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain")) + 
  theme(legend.title=element_text(color= "black", size=12), legend.text=element_text(size=7, face="italic"))
Plot_family_prev_norm6 <- Plot_family_prev_norm6 +  ggtitle ("Abundancia relativa de Familia - Umbral de prevalencia al 5%") + 
  theme(plot.title = element_text(family="Arial", size=rel(2), vjust=2,face="plain", color="black", lineheight=1.5))  
Plot_family_prev_norm6



Class_prev_norm6_dec <- tax_glom(pseq.compositional_norm6_dec, taxrank = "class") 
tb_class <- psmelt(Class_prev_norm6_dec) 
tbclass_prev_norm6 <- tb_class %>% group_by(class) %>% summarise(M=median(Abundance)) %>% ungroup() %>% unique() %>% top_n(M,n=20)
tbclass_prev_norm6 <- tbclass_prev_norm6 %>% inner_join(tb_class)
Plot_class_prev_norm6 <- ggplot(tbclass_prev_norm6, aes(SampleNames, Abundance ,fill=class)) + geom_col(position="fill") # + scale_fill_manual(values=paleta5)
Plot_class_prev_norm6 <- Plot_class_prev_norm6 + theme_minimal() + xlab("SampleNames") + ylab("Relative Abundance") +
  theme(axis.text.x = element_text(color = "black", size = 7, angle = 45, hjust = .5, vjust = .5, face = "plain")) +
  theme(axis.text.y = element_text(color = "black", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain")) + 
  theme(legend.title=element_text(color= "black", size=12), legend.text=element_text(size=7, face="italic"))
Plot_class_prev_norm6 <- Plot_class_prev_norm6 +  ggtitle ("Abundancia relativa de Clase - Umbral de prevalencia al 5%") + 
  theme(plot.title = element_text(family="Arial", size=rel(2), vjust=2,face="plain", color="black", lineheight=1.5))  
Plot_class_prev_norm6


Order_prev_norm6_dec <- tax_glom(pseq.compositional_norm6_dec, taxrank = "order") 
tb_order <- psmelt(Order_prev_norm6_dec) 
tborder_prev_norm6 <- tb_order %>% group_by(order) %>% summarise(M=median(Abundance)) %>% ungroup() %>% unique() %>% top_n(M,n=20)
tborder_prev_norm6 <- tborder_prev_norm6 %>% inner_join(tb_order)
Plot_order_prev_norm6 <- ggplot(tborder_prev_norm6, aes(SampleNames, Abundance ,fill=order)) + geom_col(position="fill") # + scale_fill_manual(values=paleta5)
Plot_order_prev_norm6 <- Plot_order_prev_norm6 + theme_minimal() + xlab("SampleNames") + ylab("Relative Abundance") +
  theme(axis.text.x = element_text(color = "black", size = 7, angle = 45, hjust = .5, vjust = .5, face = "plain")) +
  theme(axis.text.y = element_text(color = "black", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain")) + 
  theme(legend.title=element_text(color= "black", size=12), legend.text=element_text(size=7, face="italic"))
Plot_order_prev_norm6 <- Plot_order_prev_norm6 +  ggtitle ("Abundancia relativa de Orden - Umbral de prevalencia al 5%") + 
  theme(plot.title = element_text(family="Arial", size=rel(2), vjust=2,face="plain", color="black", lineheight=1.5))  
Plot_order_prev_norm6
