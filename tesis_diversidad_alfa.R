# Análisis de Diversidad alfa

# Este script contiene los cálculos de los índices de dievrsidad alfa Observed, Shannon, Simpson, ACE,
# chao1 y Fisher, además de sus visualizaciones correspondientes.

#Obtenemos la lista que contiene los status de cada muestra 
status_ps_decipher_div <- ps_decipher_div@sam_data$status
ps_decipher_div

# Shannon y Simpson
ER1_dec <- estimate_richness(ps_decipher_div, measures=c("Shannon",'Simpson'))
ER1_dec <- cbind(ER1_dec, status_ps_decipher_div)
ER1_dec <- data.table(ER1_dec, keep.rownames = TRUE)
ER1_dec <- melt(ER1_dec, id.vars=c('rn', 'status_ps_decipher_div'))
ER1_dec$status_ps_decipher_div <- ifelse(ER1_dec$status_ps_decipher_div == 1, "CAN", "NOC")

#Simpson solo.
ER1simpson_dec <- estimate_richness(ps_decipher_div, measures=c("Simpson"))
ER1simpson_dec <- cbind(ER1simpson_dec, status_ps_decipher_div)
ER1simpson_dec <- data.table(ER1simpson_dec, keep.rownames = TRUE)
ER1simpson_dec <- melt(ER1simpson_dec, id.vars=c('rn', 'status_ps_decipher_div'))
ER1simpson_dec$status_ps_decipher_div <- ifelse(ER1simpson_dec$status_ps_decipher_div == 1, "CAN", "NOC")

# Chao1
ER2_dec <- estimate_richness(ps_decipher_div, measures=c("Chao1"))
ER2_dec <- cbind(ER2_dec, status_ps_decipher_div)
ER2_dec <- data.table(ER2_dec, keep.rownames = TRUE)
ER2_dec <- melt(ER2_dec, id.vars=c('rn', 'status_ps_decipher_div'))
ER2_dec$status_ps_decipher_div <- ifelse(ER2_dec$status_ps_decipher_div == 1, "CAN", "NOC")

# ACE
ER3_dec <- estimate_richness(ps_decipher_div, measures=c("ACE"))
ER3_dec <- cbind(ER3_dec, status_ps_decipher_div)
ER3_dec <- data.table(ER3_dec, keep.rownames = TRUE)
ER3_dec <- melt(ER3_dec, id.vars=c('rn', 'status_ps_decipher_div'))
ER3_dec$status_ps_decipher_div <- ifelse(ER3_dec$status_ps_decipher_div == 1, "CAN", "NOC")

#Observed
ER4_dec <- estimate_richness(ps_decipher_div, measures=c("Observed"))
ER4_dec <- cbind(ER4_dec, status_ps_decipher_div)
ER4_dec <- data.table(ER4_dec, keep.rownames = TRUE)
ER4_dec <- melt(ER4_dec, id.vars=c('rn', 'status_ps_decipher_div'))
ER4_dec$status_ps_decipher_div <- ifelse(ER4_dec$status_ps_decipher_div == 1, "CAN", "NOC")

# Fisher
ER5_dec <- estimate_richness(ps_decipher_div, measures=c("Fisher"))
ER5_dec <- cbind(ER5_dec, status_ps_decipher_div)
ER5_dec <- data.table(ER5_dec, keep.rownames = TRUE)
ER5_dec <- melt(ER5_dec, id.vars=c('rn', 'status_ps_decipher_div'))
ER5_dec$status_ps_decipher_div <- ifelse(ER5_dec$status_ps_decipher_div == 1, "CAN", "NOC")


# Todos
ETotal_1_dec <- estimate_richness(ps_decipher_div, measures=c("Shannon", "Simpson","ACE","Chao1","Observed","Fisher"))
ETotal_1_dec <- cbind(ETotal_1_dec, status_ps_decipher_div)
ETotal_1_dec <- data.table(ETotal_1_dec, keep.rownames = TRUE)
ETotal_1_dec <- melt(ETotal_1_dec, id.vars=c('rn', 'status_ps_decipher_div'))
ETotal_1_dec$status_ps_decipher_div <- ifelse(ETotal_1_dec$status_ps_decipher_div == 1, "cancer", "no cancer")



# Boxplots para cada índice

#PLOT SHANNON Y SIMPSON
plot_Shannon_Simpson <- ggplot(ER1_dec, aes(x = status_ps_decipher_div, y = value, fill = status_ps_decipher_div)) +
  geom_boxplot() +
  geom_jitter(size=2, position=position_jitter(width=0.1))+
  facet_wrap(~ variable, scales = "free_y") +
  scale_fill_manual(values=c("#E69F00", "#84CF04"))+
  theme_minimal() +
  labs(title = "Distribución de valores de diversidad alfa Shannon y Simpson",
       x = "Status",
       y = "Valor",
       fill = "Status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"))

plot_Shannon_Simpson



#PLOT CHAO1
plot_Chao1 <- ggplot(ER2_dec, aes(x = status_ps_decipher_div, y = value, fill = status_ps_decipher_div)) +
  geom_boxplot() +
  geom_jitter(size=2, position=position_jitter(width=0.1))+
  facet_wrap(~ variable, scales = "free_y") +
  scale_fill_manual(values=c("#E69F00", "#84CF04"))+
  theme_minimal() +
  labs(title = "Distribución de valores de diversidad alfa Chao1",
       x = "Status",
       y = "Valor",
       fill = "Status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"))

plot_Chao1


#PLOT ACE
plot_ACE <- ggplot(ER3_dec, aes(x = status_ps_decipher_div, y = value, fill = status_ps_decipher_div)) +
  geom_boxplot() +
  geom_jitter(size=2, position=position_jitter(width=0.1))+
  facet_wrap(~ variable, scales = "free_y") +
  scale_fill_manual(values=c("#E69F00", "#84CF04"))+
  theme_minimal() +
  labs(title = "Distribución de valores de diversidad alfa ACE",
       x = "Status",
       y = "Valor",
       fill = "Status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"))

plot_ACE



#Observed
plot_Observed <- ggplot(ER4_dec, aes(x = status_ps_decipher_div, y = value, fill = status_ps_decipher_div)) +
  geom_boxplot() +
  geom_jitter(size=2, position=position_jitter(width=0.1))+
  facet_wrap(~ variable, scales = "free_y") +
  scale_fill_manual(values=c("#E69F00", "#84CF04"))+
  theme_minimal() +
  labs(title = "Distribución de valores de diversidad alfa Observed (Sobs)",
       x = "Status",
       y = "Valor",
       fill = "Status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"))

plot_Observed


#Fisher
plot_Fisher <- ggplot(ER5_dec, aes(x = status_ps_decipher_div, y = value, fill = status_ps_decipher_div)) +
  geom_boxplot() +
  geom_jitter(size=2, position=position_jitter(width=0.1))+
  facet_wrap(~ variable, scales = "free_y") +
  scale_fill_manual(values=c("#E69F00", "#84CF04"))+
  theme_minimal() +
  labs(title = "Distribución de valores de diversidad alfa Fisher",
       x = "Status",
       y = "Valor",
       fill = "Status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"))

plot_Fisher