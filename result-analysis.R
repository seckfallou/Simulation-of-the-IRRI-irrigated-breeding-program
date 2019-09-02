rm(list = ls())
library(magrittr)
library(tidyverse)
library(emmeans)
library(multcomp)
library(multcompView)
# load files1  ------------------------------------------------------------
res1 <-
  dir('scheme', full.names = T) %>%
  purrr::map(readr::read_rds) %>%
  purrr::map(2)

# load files2  ------------------------------------------------------------
res2 <-
  dir('scheme-burn', full.names = T) %>%
  purrr::map(readr::read_rds) %>%
  purrr::map(2)

# TBV ---------------------------------------------------------------------
# Data handling
tbv2 <- res2 %>% purrr::map(purrr::map_dfr, 'T_BV') %>%
  purrr::map2(
    .y = rep('Option 2', 25),
    .f = function(x, y) {
      dplyr::mutate(x, Burnin = y)
    }
  )%>%
  dplyr::bind_rows() %>%
  tibble::as_tibble()%>% 
  dplyr::bind_rows() %>%
  tibble::as_tibble() 
  
tbv <- dplyr::bind_rows(tbv1, tbv2) %>% 
  dplyr::rename(val = T_BV)
tbv$nEPI %<>% as.character() %>%
  dplyr::recode('0' = 'A', '2' = 'A + E')
tbv$Interaction %<>%  stringr::str_replace("int_", "Interaction ")

#Add Ggain to graph
YGain <- c()
PGain <- c()
Burnin <- c()
nEPI <- c()
Interaction <- c()
GeneticG <- tibble()
Inter <- c("Interaction 1", "Interaction 2")
EPI <- c("A", "A + E")
Burn <- c("Option 1", "Option 2")
for (b in Burn) {
for (a in EPI) {
  for (i in Inter) {
    Burnin <- append(Burnin, b)
    nEPI <- append(nEPI, a)
    Interaction <- append(Interaction, i)
    gg <-
      tbv %>%
      dplyr::filter(nEPI == a &
                      Interaction == i, Burnin == b)
    GeneticG <- dplyr::bind_rows(GeneticG, gg)
   gg <- gg %>%
      dplyr::group_by(Year) %>%
      dplyr::summarise(val = mean(val))
    gg$Year %<>% as.numeric()
    mod <- lm(val ~ Year, data = gg[-1, ])
    GG <- coefficients(mod)[[2]]
    YGain <- append(YGain, round(GG, 2))
    PGain <- append(PGain, round(GG /coefficients(mod)[[1]], 2)*100)
  }
  
}
}
gain <-
  data.frame(
    x = rep(2, 8),
    y = c(rep(4, 4),rep(5,4)),
    YGain,
    PGain,
    nEPI,
    Interaction,
    Burnin
  ) %>%
  dplyr::as_tibble()

# Graph
tbv %>%
  ggplot(aes(
    x = Year,
    y = val,
    color = interaction(Burnin),
    group = interaction(Burnin)
  )) +
  theme_light(base_size = 11) +
  stat_summary(geom =  "line",
               fun.y = mean) +
  stat_summary(
    geom = "errorbar",
    fun.data = mean_se,
    width = 0.4,
    alpha = 0.8
  ) +
  scale_x_continuous("Cycle de selection", breaks = seq(0,14,2)) +
  scale_y_continuous('Valeur genetique rÃ©elle') +
  scale_color_manual(values = c('deepskyblue3', 'grey30'), 'Burn-in') +
  guides(shape = guide_legend(order = 1)) +
  theme(
    text = element_text(size = 17),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 17),
    strip.text = element_text(colour = "black", size = 17),
    legend.position = "bottom",
    plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "lines"),
    panel.grid.minor = element_blank()
  ) +
  facet_grid(nEPI ~  Interaction)

ggsave("Evolut_TBV.png")

#Average of genetic value 
stat <- tbv %>% 
  dplyr::group_by(Year, nEPI, Interaction, Burnin) %>% 
  dplyr::summarise(Mean = mean(val),
                   Min = min(val),
                   Max = max(val))

# -------------------------------------------------------------------------
era1 <- res1 %>% purrr::map(purrr::map_dfr, 'ERA') %>%
  purrr::map2(
    .y = rep('Option 1', 25),
    .f = function(x, y) {
      dplyr::mutate(x, Burnin = y)
    }
  )%>%
  dplyr::bind_rows() %>%
  tibble::as_tibble()%>% 
  dplyr::bind_rows() %>%
  tibble::as_tibble() %>%
  dplyr::rename(Year = cycle)
era <- dplyr::bind_rows(era1, era2) %>% 
  dplyr::rename(val = pheno)
era$nEPI %<>% as.character() %>%
  dplyr::recode('0' = 'A', '2' = 'A + E')
era$Interaction %<>%  stringr::str_replace("int_", "Interaction ")

#Add Ggain to graph
YGain <- c()
MinGG <- c()
MaxGG <- c()
PGain <- c()
Burnin <- c()
nEPI <- c()
Interaction <- c()
Inter <- c("Interaction 1", "Interaction 2")
EPI <- c("A", "A + E")
Burn <- c("Option 1", "Option 2")
for (b in Burn) {
  for (a in EPI) {
    for (i in Inter) {
      Burnin <- append(Burnin, b)
      nEPI <- append(nEPI, a)
      Interaction <- append(Interaction, i)
      gg <-
        era %>%
        dplyr::filter(nEPI == a &
                        Interaction == i  & 
                        Burnin == b) %>%
        dplyr::group_by(Year) %>%
        dplyr::summarise(pheno=mean(val),
                         GG = mean(GGain),
                         PG = mean(PercentGain),
                         MinGG = min(GGain),
                         MaxGG = max(GGain))
      
      YGain <- append(YGain, round(unique(gg$GG), 2))
      PGain <- append(PGain, round(gg$GG[1]/gg$pheno[1] ,2)*100)
      MinGG <- append(MinGG, round(unique(gg$MinGG), 2))
      MaxGG <- append(MaxGG, round(unique(gg$MaxGG), 2))
    }
    
  }
}
gain <-
  data.frame(
    x = rep(2, 8),
    y = c(rep(4, 4),rep(5,4)),
    YGain,
    MinGG,
    MaxGG,
    PGain,
    nEPI,
    Interaction,
    Burnin
  ) %>%
  dplyr::as_tibble()

# Graph
era %>% 
  ggplot(aes(x=Burnin,
             y=GGain, color=Burnin))+
  theme_light(base_size = 11) +
  geom_boxplot(position = position_dodge(width = 0.9)
  )+
  scale_x_discrete("Burn-in") +
  scale_y_continuous('Gain gÃ©nÃ©tique par cycle') +
  scale_color_manual(values=c('deepskyblue3', 'grey30'))+
  scale_fill_hue(h = c(30, 210), c = 60, l = 60) +
  theme(
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 17),
    strip.text = element_text(colour = "black", size = 17),
    legend.position = "none",
    plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "lines"),
    panel.grid.minor = element_blank()
  ) +
  
  facet_grid(nEPI ~ Interaction)
ggsave('era_boxplot.png')

era %>%
  ggplot(aes(
    x = Year,
    y = val,
    color = interaction(Burnin),
    group = interaction(Burnin)
  )) +
  theme_light(base_size = 11) +
  stat_summary(geom = "point",
               fun.y = mean) +
  stat_summary(
    geom = "errorbar",
    fun.data = mean_se,
    width = 0.4,
    alpha = 0.8
  ) +
  geom_smooth(aes(Year, val),
              method = 'lm',
              size = 0.6,
              se = F) +
  scale_x_continuous("AnnÃ©e", breaks = seq(0,14,2)) +
  scale_y_continuous('phenotype moyen') +
  scale_color_manual(values = c('deepskyblue3', 'grey30'), 'Burn-in') +
  guides(shape = guide_legend(order = 1)) +
geom_text(
  aes(
    x = x,
    y = y,
    label = paste(PGain, "% Gain", sep=''),
    group = NULL
  ),
  size = 4,
  data = gain,
  parse = F
) +
theme(
  text = element_text(size = 12),
  axis.title = element_text(size = 12),
  axis.text = element_text(size = 12),
  strip.text = element_text(colour = "black", size = 12),
  legend.position = "bottom",
  plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "lines"),
  panel.grid.minor = element_blank()
) +
  facet_grid(nEPI ~  Interaction)

ggsave("GGain_TBV.png")

# gVar --------------------------------------------------------------------

gv2 <- res2 %>% purrr::map(purrr::map_dfr, 'gvar') %>%
  purrr::map2(
    .y = rep('Option 2', 25),
    .f = function(x, y) {
      dplyr::mutate(x, Burnin = y)
    }
  )%>%
  dplyr::bind_rows() %>%
  tibble::as_tibble()
  
gv <- dplyr::bind_rows(gv1, gv2) %>% 
  dplyr::rename(val = gVar)
gv$nEPI %<>% as.character() %>%
  dplyr::recode('0' = 'A', '2' = 'A + E')
gv$Interaction %<>%  stringr::str_replace("int_", "Interaction ")
gv$Year %<>%  as.numeric()

# Graph
gv %>%
  ggplot(aes(
    x = Year,
    y = val,
    color = interaction(Burnin),
    group = interaction(Burnin)
  )) +
  theme_light(base_size = 11) +
  stat_summary(geom = "line",
               fun.y = mean) +
  stat_summary(
    geom = "errorbar",
    fun.data = mean_se,
    width = 0.4,
    alpha = 0.8
  ) +
  scale_x_continuous("Cycle de sÃ©lection", breaks = c(1:10)) +
  scale_y_continuous('Variance gÃ©nÃ©tique') +
  scale_color_manual(values = c('deepskyblue3', 'grey30'), "Burn-in") +
  theme(
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 17),
    strip.text = element_text(color = 'black', size = 17),
    legend.position = "bottom",
    plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "lines"),
    panel.grid.minor = element_blank()
  ) +
  facet_grid(nEPI ~ Interaction)
ggsave("gVar-tt.png")

#stat

stat <- gv %>% 
  dplyr::group_by(Year, nEPI, Interaction, Burnin) %>% 
  dplyr::summarise(Mean = mean(val),
                   Min = min(val),
                   Max = max(val))

gv %>% 
  dplyr::group_by(Year, Burnin) %>% 
  dplyr::summarise(Mean = mean(val),
                   Min = min(val),
                   Max = max(val)) %>% 
  dplyr::filter(Year == 10)


stat %>% 
  dplyr::filter(nEPI == 'A' & Interaction == 'Interaction 1' & Burnin == 'Option 1') %>% 
  dplyr::group_by(Year) %>% 
  dplyr::summarise(Mean = mean(Mean),
                   Min = min(Min),
                   Max = max(Max)) 

# Model
gv$Year %<>% as.factor()
gv$nEPI %<>% as.factor()
model <- lm(val ~  Burnin + nEPI + Interaction , data = gv)
hist(rstudent(model), main = paste("Histogram for the residuals"))
qqnorm(rstudent(model))
plot(rstudent(model) ~ fitted(model),
     main = paste("Residuals vs fitted values"))

anova(model)
summary(model)$r.squared
summary(model)$adj.r.squared

# h2 ----------------------------------------------------------------------

# Data handling
h22 <- res2 %>% purrr::map(purrr::map_dfr, 'h2') %>%
  purrr::map2(
    .y = rep('Option 2', 25),
    .f = function(x,y){
      dplyr::mutate(x, Burnin = y)
    }
  ) %>% 
  dplyr::bind_rows() %>%
  tibble::as_tibble()

h2 <- dplyr::bind_rows(h21,h22) %>% 
  tidyr::gather(key = "trial",
                value = "val",-c(Year, nEPI, Interaction, Burnin))
h2$nEPI %<>% as.character() %>%
  dplyr::recode('0' = 'A', '2' = "A + E")
h2$nEPI %<>% as.factor()
h2$Interaction %<>%  stringr::str_replace("int_", "Interaction ")


# Graph
h2 %>%
  ggplot(aes(x = Burnin,
             y = val,
             fill = trial)) +
  theme_light(base_size = 11) +
  stat_summary(geom = "bar",
               fun.y = mean,
               width = 0.4,
               position = position_dodge(width = 0.5)) +
  stat_summary(
    geom = "errorbar",
    fun.data = mean_se,
    width = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  scale_x_discrete("Number of QTL") +
  scale_y_continuous('Heritability') +
  scale_fill_hue(h = c(30, 210),
                 c = 60,
                 l = 60,
                 "Type d'essai :") +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    strip.text = element_text(colour = "black", size = 12),
    legend.position = "bottom",
    plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "lines"),
    panel.grid.minor = element_blank()
  ) +
  facet_grid(nEPI ~  Interaction)

ggsave("Heritability-Histogram.png")

h2 %>%
  dplyr::filter(trial == 'AYT') %>%
  ggplot(aes(
    x = Year,
    y = val,
    color = interaction(Burnin),
    group = interaction(Burnin)
  )) +
  theme_light(base_size = 11) +
  stat_summary(geom = "line",
               fun.y = mean) +
  stat_summary(
    geom = "errorbar",
    fun.data = mean_se,
    width = 0.4,
    alpha = 0.8
  ) +
  scale_x_continuous("Cycle de sÃ©lection", breaks = c(1:10)) +
  scale_y_continuous('HÃ©ritabilitÃ©') +
  scale_color_manual(values = c('deepskyblue3', 'grey30'), "Burn-in") +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    strip.text = element_text(colour = "black", size = 12),
    legend.position = "bottom",
    plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "lines"),
    panel.grid.minor = element_blank()
  ) +
  facet_grid(nEPI ~  Interaction)
ggsave("h2_AYT-Evol.png")


h2 %>%
  ggplot(aes(x = nQTL,
             y = val,
             fill = trial)) +
  theme_light(base_size = 11) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  scale_x_discrete("Number of QTL") +
  scale_y_continuous('Heritability') +
  scale_fill_hue(h = c(30, 210), c = 60, l = 60) +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "lines"),
    panel.grid.minor = element_blank()
  ) +
  facet_grid(nEPI ~ Error + Interaction)

# Stat
stat <-
  h2 %>% dplyr::group_by(Burnin, nEPI , Interaction, trial) %>%
  dplyr::summarise(
    Mean = mean(val),
    Min = min(val),
    Max = max(val),
    SD = sd(val),
    SE = sd(val) / sqrt(length(val)),
    CI = sd(val) / sqrt(length(val)) *
      qt(0.95 / 2 + 0.5, length(val) - 1)
  )


# Model
model <-
  lm(val ~ nQTL + nEPI + Error + Interaction + trial , data = h2)
hist(rstudent(model), main = paste("Histogram for the residuals"))
qqnorm(rstudent(model))
plot(rstudent(model) ~ fitted(model),
     main = paste("Residuals vs fitted values"))

anova(model)
summary(model)$r.squared
summary(model)$adj.r.squared

