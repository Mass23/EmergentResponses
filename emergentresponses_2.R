library(ggplot2)
library(tidyverse)
library(ggpubr)
library(reshape2)
library(mgcv)
library(performance)
library(maps)

setwd('~/Desktop/emergentresponses')
nomis_data = read.csv('Data/NOMIS_01_2023_FULL_preprocessed.csv')
nomis_data = rename(nomis_data, `Mountain range` = Mountain.range)

##############################################################################################################################
# 1. Chlorophyll model
##############################################################################################################################
nomis_data_chla = nomis_data %>% select(chla, latitude, longitude, turbidity, water_temp, SRP, DIN, Sample, Glacier) %>% na.omit()

mod_chla_turb1 = gamm(data = nomis_data_chla, formula = chla ~ s(latitude, longitude, bs='sos', k=-1, m=1), 
                       correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_chla_turb2 = gamm(data = nomis_data_chla, formula = chla ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(turbidity, water_temp, SRP, bs = 'ts', k=5), 
                     correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_chla_turb3 = gamm(data = nomis_data_chla, formula = chla ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(turbidity, water_temp, DIN, bs = 'ts', k=5), 
                     correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_chla_turb4 = gamm(data = nomis_data_chla, formula = chla ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(turbidity, SRP, DIN, bs = 'ts', k=5), 
                     correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_chla_turb5 = gamm(data = nomis_data_chla, formula = chla ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(water_temp, SRP, DIN, bs = 'ts', k=5), 
                     correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_chla_turb6 = gamm(data = nomis_data_chla, formula = chla ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(turbidity, water_temp, bs = 'ts', k=5), 
                     correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_chla_turb7 = gamm(data = nomis_data_chla, formula = chla ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(turbidity, DIN, bs = 'ts', k=5), 
                     correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_chla_turb8 = gamm(data = nomis_data_chla, formula = chla ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(turbidity, SRP, bs = 'ts', k=5), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_chla_turb9 = gamm(data = nomis_data_chla, formula = chla ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(water_temp, DIN, bs = 'ts', k=5), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_chla_turb10 = gamm(data = nomis_data_chla, formula = chla ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(water_temp, SRP, bs = 'ts', k=5), 
                       correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_chla_turb11 = gamm(data = nomis_data_chla, formula = chla ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(DIN, SRP, bs = 'ts', k=5), 
                       correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_chla_turb12 = gamm(data = nomis_data_chla, formula = chla ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(DIN, bs = 'ts', k=5), 
                       correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_chla_turb13 = gamm(data = nomis_data_chla, formula = chla ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(SRP, bs = 'ts', k=5), 
                       correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_chla_turb14 = gamm(data = nomis_data_chla, formula = chla ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(turbidity, bs = 'ts', k=5), 
                       correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_chla_turb15 = gamm(data = nomis_data_chla, formula = chla ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(water_temp, bs = 'ts', k=5), 
                       correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
test_bf(denominator = mod_chla_turb1$gam, mod_chla_turb2$gam, mod_chla_turb3$gam, mod_chla_turb4$gam, mod_chla_turb5$gam, 
                      mod_chla_turb6$gam, mod_chla_turb7$gam, mod_chla_turb8$gam, mod_chla_turb9$gam, mod_chla_turb10$gam, 
                      mod_chla_turb11$gam, mod_chla_turb12$gam, mod_chla_turb13$gam, mod_chla_turb14$gam, mod_chla_turb15$gam)
# Best model compared to null is 7 (BF = 1.32e+07)
capture.output(summary(mod_chla_turb7$gam), file = "Statisticsmodel_chlorophyll.txt")

################################################################
########### Plotting
################################################################

# Panel A, Turbidity
chla_turb_preds = expand_grid(sample=nomis_data$Sample, turbidity=seq(-3.5, 7.1, 0.5))
chla_turb_preds$latitude = vapply(chla_turb_preds$sample, function(x) nomis_data$latitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
chla_turb_preds$longitude = vapply(chla_turb_preds$sample, function(x) nomis_data$longitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
chla_turb_preds$DIN = vapply(chla_turb_preds$sample, function(x) nomis_data$DIN[nomis_data$Sample == x], FUN.VALUE = numeric(1))
chla_turb_preds$water_temp = vapply(chla_turb_preds$sample, function(x) nomis_data$water_temp[nomis_data$Sample == x], FUN.VALUE = numeric(1))
chla_turb_preds = na.omit(chla_turb_preds)
preds = predict.gam(mod_chla_turb7$gam, newdata=chla_turb_preds, newdata.guaranteed = T, se.fit = T)
chla_turb_preds$pred = preds$fit
chla_turb_preds$se = preds$se.fit
chla_turb_preds = chla_turb_preds %>% group_by(turbidity) %>% summarise(pred = mean(pred), se = mean(se))
p1 = ggplot() +
  geom_point(data=nomis_data, aes(x=turbidity,y=chla, colour=`Mountain range`), alpha=0.2, size=2.5) +
  labs(colour = "Mountain range") +
  geom_ribbon(data=chla_turb_preds, aes(x=turbidity, ymin=pred-se, ymax=pred+se), colour='lightgrey', alpha=.1) +
  geom_line(data=chla_turb_preds, aes(x=turbidity,y=pred), size=1.5, colour='black') +
  geom_line(data=chla_turb_preds, aes(x=turbidity,y=pred+se), size=1, colour='dimgrey', linetype='dashed') +
  geom_line(data=chla_turb_preds, aes(x=turbidity,y=pred-se), size=1, colour='dimgrey', linetype='dashed') +
  scale_y_continuous(name = bquote(""*Chlorophyll-italic(a)~(ln~mu*g~g^-1~DM)*"")) +
  scale_x_continuous(name = bquote(""*Turbidity~(ln~ NTU)*"")) +
  scale_colour_brewer(palette = "Set1") + scale_colour_brewer(palette = "Paired")+
  theme_bw() + theme(legend.position = 'none')

# Panel C, DIN
chla_din_preds = expand_grid(sample=nomis_data$Sample, DIN=seq(2, 7, 0.1))
chla_din_preds$latitude = vapply(chla_din_preds$sample, function(x) nomis_data$latitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
chla_din_preds$longitude = vapply(chla_din_preds$sample, function(x) nomis_data$longitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
chla_din_preds$turbidity = vapply(chla_din_preds$sample, function(x) nomis_data$turbidity[nomis_data$Sample == x], FUN.VALUE = numeric(1))
chla_din_preds$water_temp = vapply(chla_din_preds$sample, function(x) nomis_data$water_temp[nomis_data$Sample == x], FUN.VALUE = numeric(1))
chla_din_preds = na.omit(chla_din_preds)
preds = predict.gam(mod_chla_turb7$gam, newdata=chla_din_preds, newdata.guaranteed = T, se.fit = T, type = 'response')
chla_din_preds$pred = preds$fit
chla_din_preds$se = preds$se.fit
chla_din_preds = chla_din_preds %>% group_by(DIN) %>% summarise(pred = mean(pred), se = mean(se))
p2 = ggplot() +
  geom_point(data=nomis_data, aes(x=DIN,y=chla, colour=`Mountain range`), alpha=0.2, size=2.5) +
  labs(colour = "Mountain range") +
  geom_ribbon(data=chla_din_preds, aes(x=DIN, ymin=pred-se, ymax=pred+se), colour='lightgrey', alpha=.1) +
  geom_line(data=chla_din_preds, aes(x=DIN,y=pred), size=1.5, colour='black') +
  geom_line(data=chla_din_preds, aes(x=DIN,y=pred+se), size=1, colour='dimgrey', linetype='dashed') +
  geom_line(data=chla_din_preds, aes(x=DIN,y=pred-se), size=1, colour='dimgrey', linetype='dashed') +
  scale_y_continuous(name = '') +
  scale_x_continuous(name = bquote(""*Dissolved~inorg.~nitrogen~(ln~ N~l^-1)*"")) +
  scale_colour_brewer(palette = "Set1") + scale_colour_brewer(palette = "Paired") +
  theme_bw() + theme(legend.position = 'none')
p = ggarrange(p1, p2, ncol = 2, nrow = 1, labels = c('a', 'b'))
ggsave(p, filename = 'Plots/Fig2_chla_turb_water_temp.pdf', width=6, height = 3)


