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
# 1. Models CUE/TER
##############################################################################################################################
### CUE  ##########################
cue_ter_models_data = nomis_data %>% select(CUE, TER_cn, TER_cp, water_temp, chla, NP, latitude, longitude, Sample, Glacier, AP_cell, `Mountain range`) %>% na.omit

# Summary stats and CUE/TER distributions
cue_ter_models_data = na.omit(cue_ter_models_data)

cue_model1 = gamm(data = cue_ter_models_data, formula = CUE ~ s(latitude, longitude, bs='sos', k=-1, m=1),
                  correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
cue_model2 = gamm(data = cue_ter_models_data, formula = CUE ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(chla, water_temp, NP, bs='ts', k=5),
                  correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
cue_model3 = gamm(data = cue_ter_models_data, formula = CUE ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(chla, water_temp, bs='ts', k=5),
                 correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
cue_model4 = gamm(data = cue_ter_models_data, formula = CUE ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(chla, NP, bs='ts', k=5),
                 correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
cue_model5 = gamm(data = cue_ter_models_data, formula = CUE ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(water_temp, NP, bs='ts', k=5),
                 correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
cue_model6 = gamm(data = cue_ter_models_data, formula = CUE ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(chla, bs='ts', k=5),
                 correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
cue_model7 = gamm(data = cue_ter_models_data, formula = CUE ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(water_temp, bs='ts', k=5),
                 correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
cue_model8 = gamm(data = cue_ter_models_data, formula = CUE ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(NP, bs='ts', k=5),
                 correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
test_bf(denominator = cue_model1$gam, cue_model2$gam, cue_model3$gam, cue_model4$gam, 
                      cue_model5$gam, cue_model6$gam, cue_model7$gam, cue_model8$gam)
# Best model compared to null is 6 (BF = 5.58e+03)
capture.output(summary(cue_model6$gam), file = "Statistics/model_cue.txt")


### TER CN  ##########################
tercn_model1 = gamm(data = cue_ter_models_data, formula = TER_cn ~ s(latitude, longitude, bs='sos', k=-1, m=1),
                  correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
tercn_model2 = gamm(data = cue_ter_models_data, formula = TER_cn ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(chla, water_temp, NP, bs='ts', k=5),
                  correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
tercn_model3 = gamm(data = cue_ter_models_data, formula = TER_cn ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(chla, water_temp, bs='ts', k=5),
                  correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
tercn_model4 = gamm(data = cue_ter_models_data, formula = TER_cn ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(chla, NP, bs='ts', k=5),
                  correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
tercn_model5 = gamm(data = cue_ter_models_data, formula = TER_cn ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(water_temp, NP, bs='ts', k=5),
                  correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
tercn_model6 = gamm(data = cue_ter_models_data, formula = TER_cn ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(chla, bs='ts', k=5),
                  correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
tercn_model7 = gamm(data = cue_ter_models_data, formula = TER_cn ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(water_temp, bs='ts', k=5),
                  correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
tercn_model8 = gamm(data = cue_ter_models_data, formula = TER_cn ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(NP, bs='ts', k=5),
                  correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
test_bf(denominator = tercn_model1$gam, tercn_model2$gam, tercn_model3$gam, tercn_model4$gam, 
                      tercn_model5$gam, tercn_model6$gam, tercn_model7$gam, tercn_model8$gam)
# Best model is null
capture.output(summary(tercn_model1$gam), file = "Statistics/model_tercn.txt")

### TER CP ##########################
tercp_model1 = gamm(data = cue_ter_models_data, formula = TER_cp ~ s(latitude, longitude, bs='sos', k=-1, m=1),
                    correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
tercp_model2 = gamm(data = cue_ter_models_data, formula = TER_cp ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(chla, water_temp, NP, bs='ts', k=5),
                    correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
tercp_model3 = gamm(data = cue_ter_models_data, formula = TER_cp ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(chla, water_temp, bs='ts', k=5),
                    correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
tercp_model4 = gamm(data = cue_ter_models_data, formula = TER_cp ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(chla, NP, bs='ts', k=5),
                    correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
tercp_model5 = gamm(data = cue_ter_models_data, formula = TER_cp ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(water_temp, NP, bs='ts', k=5),
                    correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
tercp_model6 = gamm(data = cue_ter_models_data, formula = TER_cp ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(chla, bs='ts', k=5),
                    correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
tercp_model7 = gamm(data = cue_ter_models_data, formula = TER_cp ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(water_temp, bs='ts', k=5),
                    correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
tercp_model8 = gamm(data = cue_ter_models_data, formula = TER_cp ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(NP, bs='ts', k=5),
                    correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
test_bf(denominator = tercp_model1$gam, tercp_model2$gam, tercp_model3$gam, tercp_model4$gam, 
                      tercp_model5$gam, tercp_model6$gam, tercp_model7$gam, tercp_model8$gam)
# Best model is 6 (BF = 1.16e+04)
capture.output(summary(tercp_model6$gam), file = "Statistics/model_tercp.txt")

##############################################################################################################################
# 2. Plotting
##############################################################################################################################
cue_ter_pred = expand_grid(chla=seq(-14.5, 0.5, 0.5), sample=nomis_data$Sample)
cue_ter_pred$latitude = vapply(cue_ter_pred$sample, function(x) nomis_data$latitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
cue_ter_pred$longitude = vapply(cue_ter_pred$sample, function(x) nomis_data$longitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))

preds_cue = predict.gam(cue_model6$gam, newdata=cue_ter_pred, newdata.guaranteed = T, se.fit = T)
cue_ter_pred$preds_cue = preds_cue$fit
cue_ter_pred$se_cue = preds_cue$se.fit

preds_tercp = predict.gam(tercp_model6$gam, newdata=cue_ter_pred, newdata.guaranteed = T, se.fit = T)
cue_ter_pred$preds_tercp = preds_tercp$fit
cue_ter_pred$se_tercp = preds_tercp$se.fit

cue_ter_pred = cue_ter_pred %>% group_by(chla) %>% summarise(preds_cue = mean(preds_cue), 
                                                             se_cue = mean(se_cue),
                                                             preds_tercp = mean(preds_tercp), 
                                                             se_tercp = mean(se_tercp))

p1 = ggplot() +
  geom_point(data=cue_ter_models_data, aes(x=chla,y=CUE, colour=`Mountain range`), alpha=0.2, size=2.5) +
  labs(colour = "Mountain range") +
  geom_ribbon(data=cue_ter_pred, aes(x=chla, ymin=preds_cue-se_cue, ymax=preds_cue+se_cue), colour='lightgrey', alpha=.1) +
  geom_line(data=cue_ter_pred, aes(x=chla,y=preds_cue), size=1.5, colour='black') +
  geom_line(data=cue_ter_pred, aes(x=chla,y=preds_cue+se_cue), size=1, colour='dimgrey', linetype='dashed') +
  geom_line(data=cue_ter_pred, aes(x=chla,y=preds_cue-se_cue), size=1, colour='dimgrey', linetype='dashed') +
  scale_x_continuous(name = bquote(""*Chlorophyll-italic(a)~(ln~mu*g~g^-1~DM)*"")) +
  scale_y_continuous(name = bquote(""*ln~CUE*"")) +
  scale_colour_brewer(palette = "Set1") + scale_colour_brewer(palette = "Paired")+
  theme_bw() + theme(legend.position = 'none')

p2 = ggplot() +
  geom_point(data=cue_ter_models_data, aes(x=chla,y=TER_cp, colour=`Mountain range`), alpha=0.2, size=2.5) +
  labs(colour = "Mountain range") +
  geom_ribbon(data=cue_ter_pred, aes(x=chla, ymin=preds_tercp-se_tercp, ymax=preds_tercp+se_tercp), colour='lightgrey', alpha=.1) +
  geom_line(data=cue_ter_pred, aes(x=chla,y=preds_tercp), size=1.5, colour='black') +
  geom_line(data=cue_ter_pred, aes(x=chla,y=preds_tercp+se_tercp), size=1, colour='dimgrey', linetype='dashed') +
  geom_line(data=cue_ter_pred, aes(x=chla,y=preds_tercp-se_tercp), size=1, colour='dimgrey', linetype='dashed') +
  scale_x_continuous(name = bquote(""*Chlorophyll-italic(a)~(ln~mu*g~g^-1~DM)*"")) +
  scale_y_continuous(name = bquote(""*ln~TER~C:P*"")) +
  scale_colour_brewer(palette = "Set1") + scale_colour_brewer(palette = "Paired")+
  theme_bw() + theme(legend.position = 'none')

p = ggarrange(p2, p1, ncol = 2, nrow = 1, labels = c('a', 'b'), align='h')
ggsave(p, filename = 'Plots/Fig4_cue_ter.pdf', width=7, height = 3.5)


##############################################################################################################################
# 3. CUE/AP
##############################################################################################################################
cue_ap_data = nomis_data %>% select(CUE, latitude, longitude, Sample, Glacier, AP_cell, CN, `Mountain range`) %>% na.omit()
cue_ap_model = gamm(data = cue_ap_data, formula = CUE ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(AP_cell, bs='ts', k=5), 
                   correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML', select = T)
capture.output(summary(cue_ap_model$gam), file = "Statistics/model_ap_cell_cue.txt")


cue_ap_preds = expand_grid(sample=nomis_data$Sample, AP_cell=seq(-19.5,-10.4,0.5))
cue_ap_preds$latitude = vapply(cue_ap_preds$sample, function(x) nomis_data$latitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
cue_ap_preds$longitude = vapply(cue_ap_preds$sample, function(x) nomis_data$longitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
preds = predict.gam(cue_ap_model$gam, newdata=cue_ap_preds, newdata.guaranteed = T, se.fit = T)
cue_ap_preds$pred = preds$fit
cue_ap_preds$se = preds$se.fit
cue_ap_preds = cue_ap_preds %>% group_by(AP_cell) %>% summarise(pred = mean(pred), se = mean(se))

p = ggplot() + 
  geom_point(cue_ap_data, mapping=aes(x=AP_cell, y=CUE, colour=`Mountain range`)) +
  labs(colour = "Mountain range") +
  geom_ribbon(data=cue_ap_preds, aes(x=AP_cell, ymin=pred-se, ymax=pred+se), colour='lightgrey', alpha=.1) +
  geom_line(data=cue_ap_preds, aes(x=AP_cell,y=pred), size=1.5, colour='black') +
  geom_line(data=cue_ap_preds, aes(x=AP_cell,y=pred+se), size=1, colour='dimgrey', linetype='dashed') +
  geom_line(data=cue_ap_preds, aes(x=AP_cell,y=pred-se), size=1, colour='dimgrey', linetype='dashed') +
  scale_x_continuous(name = bquote(""*ln~AP~by~cell*"")) +
  scale_y_continuous(name = bquote(""*ln~CUE*"")) +
  scale_colour_brewer(palette = "Set1") + scale_colour_brewer(palette = "Paired")+
  theme_bw() + theme(legend.position = 'none')
ggsave(p, filename = 'Plots/Fig5_AP_CUE.pdf', width=5, height = 5)

