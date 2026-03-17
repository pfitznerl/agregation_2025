library(ggplot2)
library(ggpubr) # Utilisation de ggarrange de ggpubr
library(dplyr)
library(tidyr)
library(scales)
library(grid)

dir_exe="/home/mf/dp/mcad/pfitznerl/SAVE/aggregation/"
dir_utemp="/scratch/work/pfitznerl/aggregation/"
dir_plot=paste(dir_utemp,"plot_aggregation/",sep="")
dir_rmse="rmse/"
dir_weights="weights/"
dir_mixture="mixture/"
dir_alpha="alpha/"
dir_eta="eta/"
dir_obs_pred="obs_pred/"
dir_boxplot="boxplots/"
dir_data_without="data/donnee_sans_pearp"
dir_data_with="data/donnee_avec_pearp"
dir_sorties="sorties/"

# lancer BOA et BOA sleep et enregistrer df dans la fonction plot_obs_pred_window
load("z_data_plot_papier_BOA_sleep.Rdata") # df_sleep
load("z_data_plot_papier_BOA.Rdata") # df

# il faut avoir des données du types :
# date modele temperature
#2025-06-03 BOA 5.3

levels(df_sleep$modele)=c(levels(df_sleep$modele),"BOA_sleep")
levels(df$modele)=c(levels(df$modele),"BOA_sleep")

df_sleep[which(df_sleep$modele=="BOA"),]$modele=factor("BOA_sleep",levels=levels(df_sleep$modele))
df=rbind(df,df_sleep[which(df_sleep$modele=="BOA_sleep"),])

# 1. Préparation et calcul des différences
df$date.valid <- as.Date(df$date.valid)
df$diff = df$temperature - df[which(df$modele=="obs"),]$temperature

obs <- df %>% filter(modele == "obs") %>% select(date.valid, temperature) %>% rename(temp_obs = temperature)
modeles <- df %>% filter(modele != "obs")

# couleurs d'origine
color.values = c("springgreen2","springgreen4","gold2","goldenrod4","hotpink",
                 "deeppink","blue","deepskyblue","orange","brown","red","black")

# Vecteur demandé :
# T = bleu, R = rouge, F = couleur normale BOA_sleep
select_boa_sleep <- c("R","T","R","F","T","T","T","T","T","T","T","T","T","T","T","T","R","R","R","R","T")

# On prépare le dataframe pour le plot du haut
aggregations <- modeles %>%
  filter(modele %in% c("BOA", "BOA_sleep")) %>%
  arrange(modele, date.valid)

# Vérification
if(sum(aggregations$modele == "BOA_sleep") != length(select_boa_sleep)) {
  stop("La longueur de select_boa_sleep ne correspond pas au nombre de points BOA_sleep.")
}

# Statut des points BOA_sleep
aggregations$point_status <- NA
aggregations$point_status[aggregations$modele == "BOA_sleep"] <- select_boa_sleep

# Sous-ensembles à surajouter
points_sleep_T <- aggregations %>%
  filter(modele == "BOA_sleep", point_status == "T")

points_sleep_R <- aggregations %>%
  filter(modele == "BOA_sleep", point_status == "R")

points_sleep_F <- aggregations %>%
  filter(modele == "BOA_sleep", point_status == "F")

# 2. Graphique du haut (plot.diff)
plot.diff <- ggplot(aggregations, aes(x = date.valid, y = diff, color = modele, group = modele)) +
  geom_line(size = 2) +
  geom_point(size = 6) +
  geom_point(data = points_sleep_T, aes(x = date.valid, y = diff), color = "blue", size = 6) +
  geom_point(data = points_sleep_R, aes(x = date.valid, y = diff), color = "red", size = 6) +
  geom_point(data = points_sleep_F, aes(x = date.valid, y = diff), color = "#00BFC4", size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 1) +
  scale_color_manual(
    values = c(
      "BOA" = "#8B4789",
      "BOA_sleep" = "#00BFC4"
    ),
    breaks = c("BOA", "BOA_sleep"),
    labels = parse(text = c("BOA", "BOA^s"))
  ) +
  theme_bw(base_size = 18) +
  labs(x = "Date (December 2021)", y = "Error (°C)") +
  scale_x_date(breaks = pretty_breaks(n = 15), labels = date_format("%d")) +
  guides(color = guide_legend(override.aes = list(size = 5, linewidth = 5))) +
  theme(
    plot.title = element_text(size = 30, face = "bold"),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 36),
    axis.text.x  = element_text(size = 30),
    axis.text.y  = element_text(size = 38),
    legend.title = element_blank(),
    legend.text = element_text(size = 30),
    legend.key.size = unit(2, "lines")
  )

obs_rows <- df %>% filter(modele == "obs")
df_autres <- modeles %>% filter(!modele %in% c("BOA", "BOA_sleep"))
df_autres <- bind_rows(df_autres, obs_rows)

# IMPORTANT :
# remettre ici les niveaux d'origine dans l'ordre voulu pour retrouver les bonnes couleurs
# Remplace ce vecteur par l'ordre exact historique de tes modèles si besoin
ordre_modeles_bas <- unique(as.character(df_autres$modele))

df_autres$modele <- factor(df_autres$modele, levels = ordre_modeles_bas)
obs_rows$modele   <- factor(obs_rows$modele, levels = ordre_modeles_bas)

# 3. Graphique du bas : Températures
plot.experts <- ggplot() +
  geom_line(data = df_autres, aes(x = date.valid, y = temperature, color = modele, group = modele), size = 2) +
  geom_line(data = obs_rows, aes(x = date.valid, y = temperature, color = modele, group = modele), size = 4) +
  scale_colour_manual(values = color.values) +
  scale_y_continuous(n.breaks = 6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 1) +
  theme_bw(base_size = 18) +
  labs(x = "Date (December 2021)", y = "Temperature (°C)") +
  scale_x_date(breaks = pretty_breaks(n = 15), labels = date_format("%d")) +
  guides(color = guide_legend(override.aes = list(size = 5, linewidth = 5))) +
  theme(
    plot.title = element_text(size = 30, face = "bold"),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 36),
    axis.text.x  = element_text(size = 30),
    axis.text.y  = element_text(size = 38),
    legend.title = element_blank(),
    legend.text = element_text(size = 30),
    legend.key.size = unit(2, "lines")
  )

# 4. Assemblage – on force l’alignement vertical
plot_final <- ggarrange(
  plot.diff,
  plot.experts,
  ncol  = 1,
  nrow  = 2,
  align = "v",
  widths = c(1, 1),
  common.legend = FALSE
)

# 5. Save the plot with nice dimensions
ggsave(file = paste0(dir_plot, dir_obs_pred, "obs_err_agreg_window_1253_48_cham.pdf"),
       plot = plot_final,
       width = 15,
       height = 10,
       dpi = 600,
       units = "in",
       scale = 1.2)

print(paste0(dir_plot, dir_obs_pred, "obs_err_agreg_window_1253_48_cham.pdf has been created"))