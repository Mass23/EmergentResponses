library(ggplot2)
library(tidyverse)
library(ggpubr)
library(reshape2)
library(mgcv)
library(performance)
library(maps)

setwd('~/Desktop/Tyler')
nomis_data = read.csv('NOMIS_01_2023_FULL.csv')

##############################################################################################################################
# 1. Data preparation
##############################################################################################################################
nomis_data$Sample = map_chr(nomis_data$patch, function(x) paste0(strsplit(x, '_')[[1]][1:2], collapse = '_'))

nomis_data$n6_no2..ug.l.1.[is.na(nomis_data$n6_no2..ug.l.1.)] = 0

nomis_data = nomis_data %>% group_by(Sample, Expedition) %>% summarise(latitude = mean(lat_sp..DD., na.rm=T),
                                                                       longitude = mean(lon_sp..DD., na.rm=T),
                                                                       elevation = mean(ele_sp..m., na.rm=T),
                                                                       gl_coverage = mean(gl_cov...., na.rm=T),
                                                                       gl_area = mean(gl_sa..km2., na.rm=T),
                                                                       gl_dist = mean(sn_sp_dist..m., na.rm=T),
                                                                       turbidity = mean(turb..NTU., na.rm=T),
                                                                       conductivity = mean(conductivity..uS.cm..1., na.rm=T),
                                                                       chla = mean(chla..ug.g.1., na.rm=T),
                                                                       water_temp = mean(water_temp..C., na.rm=T),
                                                                       LAP = mean(lap..nmol.g.1.h.1., na.rm=T),
                                                                       NAG = mean(nag..nmol.g.1.h.1., na.rm=T),
                                                                       BG = mean(bg..nmol.g.1.h.1., na.rm=T),
                                                                       AG = mean(ag..nmol.g.1.h.1., na.rm=T),
                                                                       AP = mean(ap..nmol.g.1.h.1., na.rm=T),
                                                                       cells = mean(sba..cells.g.1., na.rm=T),
                                                                       DOC = mean(doc..ug.l.1., na.rm=T),
                                                                       DIN = mean(n4_nh4..ug.l.1. + n5_no3..ug.l.1. + n6_no2..ug.l.1., na.rm=T),
                                                                       SRP = mean(n3_srp..ug.l.1., na.rm=T))

nomis_data = nomis_data %>% mutate(`Mountain range` = Expedition)
nomis_data$Glacier = map_chr(nomis_data$Sample, function(x) strsplit(x ,'_')[[1]][1])
nomis_data$gl_index = sqrt(nomis_data$gl_area) / (sqrt(nomis_data$gl_area) + nomis_data$gl_dist)

min_data = read.csv('Final_geology.csv')
nomis_data$Clays = map_dbl(nomis_data$Glacier, function(x) min_data$Clays[min_data$Glacier == x])
nomis_data$Quartz = map_dbl(nomis_data$Glacier, function(x) min_data$Quartz[min_data$Glacier == x])
nomis_data$Calcite = map_dbl(nomis_data$Glacier, function(x) min_data$Calcite[min_data$Glacier == x])
nomis_data$Feldspar = map_dbl(nomis_data$Glacier, function(x) min_data$Feldspar[min_data$Glacier == x])

# Compute molecular ratios
nomis_data$NP = (nomis_data$DIN/14.007) / (nomis_data$SRP/30.974)
nomis_data$CP = (nomis_data$DOC/12.011) / (nomis_data$SRP/30.974)
nomis_data$CN = (nomis_data$DOC/12.011) / (nomis_data$DIN/14.007)

nomis_data$AP_cell = nomis_data$AP / nomis_data$cells

# Compute TER and CUE (biogeochemical equilibrium = {[An * Bcn / TERcn]*[Ap*Bcp/TERcp]}^0.5)
nomis_data$LAP = nomis_data$LAP + (min(nomis_data$LAP[which(nomis_data$LAP > 0)])/2)
nomis_data$NAG = nomis_data$NAG + (min(nomis_data$NAG[which(nomis_data$NAG > 0)])/2)
nomis_data$BG = nomis_data$BG + (min(nomis_data$LAP[which(nomis_data$BG > 0)])/2)
nomis_data$AG = nomis_data$AG + (min(nomis_data$AG[which(nomis_data$AG > 0)])/2)
nomis_data$AP = nomis_data$AP + (min(nomis_data$AP[which(nomis_data$AP > 0)])/2)

##############################################################################################################################
# 2. Compute CUE and TER
##############################################################################################################################
bcn = 8.6
bcp = 60
ass_efficiency = 0.9

model_tercn <- lm(data = nomis_data %>% select(LAP, BG, NAG) %>% na.omit, formula = log(LAP + BG) ~ log(LAP + NAG))
nomis_data$TER_cn <- (((nomis_data$LAP+nomis_data$BG)/(nomis_data$LAP+nomis_data$NAG))*bcn) / as.double(coef(model_tercn)['(Intercept)'])

model_tercp <- lm(data = nomis_data %>% select(LAP, BG, AP) %>% na.omit, formula = log(LAP + BG) ~ log(AP))
nomis_data$TER_cp<-(((nomis_data$LAP+nomis_data$BG)/(nomis_data$AP))*bcp) / as.double(coef(model_tercp)['(Intercept)'])

nomis_data$CUE<-(((ass_efficiency * bcn )/ (nomis_data$TER_cn))*((ass_efficiency * bcp)/(nomis_data$TER_cp)))^0.5

# outlier analysis
nomis_data$CUE_zscore = (nomis_data$CUE-mean(nomis_data$CUE, na.rm=T))/sd(nomis_data$CUE, na.rm=T)
nomis_data$CUE_outlier = ifelse(abs(nomis_data$CUE_zscore) > 3, 'Yes', 'No')
nomis_data$TERcn_zscore = (nomis_data$TER_cn-mean(nomis_data$TER_cn, na.rm=T))/sd(nomis_data$TER_cn, na.rm=T)
nomis_data$TERcn_outlier = ifelse(abs(nomis_data$TERcn_zscore) > 3, 'Yes', 'No')
nomis_data$TERcp_zscore = (nomis_data$TER_cp-mean(nomis_data$TER_cp, na.rm=T))/sd(nomis_data$TER_cp, na.rm=T)
nomis_data$TERcp_outlier = ifelse(abs(nomis_data$TERcp_zscore) > 3, 'Yes', 'No')

nomis_data$CUE[nomis_data$CUE_outlier == 'Yes'] = NA
nomis_data$CUE[nomis_data$CUE_outlier == 'Yes'] = NA
nomis_data$CUE[nomis_data$CUE_outlier == 'Yes'] = NA

##############################################################################################################################
# 3. Summary stats
##############################################################################################################################
data_summary = nomis_data %>% ungroup() %>% select(gl_index, gl_area, gl_dist, gl_coverage, turbidity, conductivity, water_temp,
                                                   chla, DOC, DIN, SRP, NP, CN, CP, CUE, TER_cn, TER_cp, LAP, BG, AG, AP, NAG) %>% 
  summarise_all(quantile, na.rm=T, probs=c(0.25,0.5,0.75)) %>% t %>% as.data.frame

data_summary_mountain_ranges = nomis_data %>% ungroup() %>% select(gl_index, gl_area, gl_dist, gl_coverage, turbidity, conductivity, water_temp,
                                                   chla, DOC, DIN, SRP, NP, CN, CP, CUE, TER_cn, TER_cp, LAP, BG, AG, AP, NAG, `Mountain range`) %>% group_by(`Mountain range`) %>%
  summarise_all(quantile, na.rm=T, probs=c(0.25,0.5,0.75)) %>% t %>% as.data.frame

data_summary$median_difference = NA
data_summary$wilcox_p = NA
for (i in 1:nrow(data_summary)){
  var = rownames(data_summary)[i]
  if (!(var %in% c('DOC', 'CN', 'CP'))){ # DOC is only for UP sites
      vals_up = nomis_data %>% filter(grepl('UP', Sample)) %>% pull(!!sym(var))
      vals_dn = nomis_data %>% filter(grepl('DN', Sample)) %>% pull(!!sym(var))
      data_summary$median_difference[i] = median(na.omit(vals_up - vals_dn))
      data_summary$wilcox_p[i] = wilcox.test(na.omit(vals_up - vals_dn))$p.value}}

colnames(data_summary) = c('q25', 'median', 'q75', 'p_val_updn', 'median_difference_updn')
write.csv(data_summary, file = 'summary_stats.csv', quote = F, row.names = T)
write.csv(data_summary_mountain_ranges, file = 'summary_stats_mountaing_ranges.csv', quote = F, row.names = T)

# Log transforms variables if needed
log_const <- function(x){return(log(x + (min(x[which(x > 0)])/2)))}

nomis_data$AP_cell = log_const(nomis_data$AP_cell)
nomis_data$gl_index = log_const(nomis_data$gl_index)
nomis_data$gl_dist = log_const(nomis_data$gl_dist)
nomis_data$gl_area = log_const(nomis_data$gl_area)

nomis_data$turbidity = log_const(nomis_data$turbidity)
nomis_data$conductivity = log_const(nomis_data$conductivity)
nomis_data$water_temp = log_const(nomis_data$water_temp)
nomis_data$chla = log_const(nomis_data$chla)

nomis_data$DOC = log_const(nomis_data$DOC)
nomis_data$DIN = log_const(nomis_data$DIN)
nomis_data$SRP = log_const(nomis_data$SRP)

nomis_data$NP = log_const(nomis_data$NP)
nomis_data$CN = log_const(nomis_data$CN)
nomis_data$CP = log_const(nomis_data$CP)
nomis_data$CUE = log_const(nomis_data$CUE)
nomis_data$TER_cn = log_const(nomis_data$TER_cn)
nomis_data$TER_cp = log_const(nomis_data$TER_cp)

nomis_data$LAP = log(nomis_data$LAP)
nomis_data$NAG = log(nomis_data$NAG)
nomis_data$BG = log(nomis_data$BG)
nomis_data$AG = log(nomis_data$AG)
nomis_data$AP = log(nomis_data$AP)

nomis_data$Quartz = log_const(nomis_data$Quartz)
nomis_data$Calcite = log_const(nomis_data$Calcite)
nomis_data$Feldspar = log_const(nomis_data$Feldspar)
nomis_data$Clays = log_const(nomis_data$Clays)

##############################################################################################
# CUE / TER distributions
# change exp to mountain range
cue_p1 = ggplot(nomis_data %>% filter(`Mountain range` != 'Alaska Range'), aes(x=`Mountain range`, y=exp(CUE), fill=`Mountain range`)) + geom_boxplot() + scale_fill_brewer(palette = "Set1") + scale_fill_brewer(palette = "Paired") +
  theme_bw() + xlab('') + ylab('ln CUE') + theme(axis.text.x = element_blank(), legend.position = 'none') + scale_y_log10()
cue_p2 = ggplot(nomis_data, aes(y=exp(CUE))) + geom_density() + geom_hline(yintercept = median(exp(nomis_data$CUE), na.rm=T), colour='red') + theme_bw() + ylab('') + xlab('')+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_y_log10()
cue_p = ggarrange(cue_p1, cue_p2, ncol = 2, nrow = 1, widths = c(0.75, 0.25), align = 'h', labels = c('a','b'))

tercn_p1 = ggplot(nomis_data %>% filter(`Mountain range` != 'Alaska Range'), aes(x=`Mountain range`, y=exp(TER_cn), fill=`Mountain range`)) + geom_boxplot() + scale_fill_brewer(palette = "Set1") + scale_fill_brewer(palette = "Paired") +
  theme_bw() + xlab('') + ylab('ln TER C:N') + theme(axis.text.x = element_blank(), legend.position = 'none') + scale_y_log10()
tercn_p2 = ggplot(nomis_data, aes(y=exp(TER_cn))) + geom_density() + geom_hline(yintercept = median(exp(nomis_data$TER_cn), na.rm=T), colour='red') + theme_bw() + ylab('') + xlab('')+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_y_log10()
tercn_p = ggarrange(tercn_p1, tercn_p2, ncol = 2, nrow = 1, widths = c(0.75, 0.25), align = 'h', labels = c('c','d'))

tercp_p1 = ggplot(nomis_data %>% filter(`Mountain range` != 'Alaska Range'), aes(x=`Mountain range`, y=exp(TER_cp), fill=`Mountain range`)) + geom_boxplot() + scale_fill_brewer(palette = "Set1") + scale_fill_brewer(palette = "Paired") +
  theme_bw() + xlab('') + ylab('ln TER C:P') + theme(axis.text.x = element_blank(),legend.position = 'right') + scale_y_log10()
tercp_p2 = ggplot(nomis_data, aes(y=exp(TER_cp))) + geom_density() + geom_hline(yintercept = median(exp(nomis_data$TER_cp), na.rm=T), colour='red') + theme_bw() + ylab('') +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_y_log10()
leg = get_legend(tercp_p1)
tercp_p = ggarrange(tercp_p1 + theme(legend.position = 'none'), tercp_p2, ncol = 2, nrow = 1, widths = c(0.75, 0.25), align = 'h', labels = c('e','f'))

p1 = ggarrange(cue_p, tercn_p, tercp_p, nrow = 3, ncol = 1, heights = c(0.333,0.333,0.333))
p = ggarrange(p1, leg, ncol = 2, widths = c(0.73, 0.27))
ggsave(plot = p, 'Fig3_cue_ter_distributions.pdf', width = 7, height = 6)

##############################################################################################################################
# 4. Glacier coverage/dist/area vs Turbidity/Water Temperature/Water NP
##############################################################################################################################
mod_turb_glacier_data = nomis_data %>% select(`Mountain range`, latitude, longitude, turbidity, 
                      Sample, Glacier, gl_coverage, gl_index, gl_area, gl_dist) %>% na.omit()  
mod_turb_glacier1 =  gamm(data = mod_turb_glacier_data,
                        formula = turbidity ~ s(latitude, longitude, bs='sos', k=-1, m=1), 
                        correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_turb_glacier8 =  gamm(data = mod_turb_glacier_data,
                      formula = turbidity ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_coverage, gl_dist, bs = 'ts', k=5, d=c(2,1)), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML', family = gaussian())
mod_turb_glacier7 =  gamm(data = mod_turb_glacier_data,
                       formula = turbidity ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_coverage, bs = 'ts', k=5), 
                       correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML', family = gaussian())
mod_turb_glacier6 =  gamm(data = mod_turb_glacier_data,
                       formula = turbidity ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_dist, bs = 'ts', k=5), 
                       correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML', family = gaussian())
mod_turb_glacier5 =  gamm(data = mod_turb_glacier_data,
                       formula = turbidity ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_coverage, gl_dist, bs = 'ts', k=5), 
                       correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML', family = gaussian())
mod_turb_glacier4 =  gamm(data = mod_turb_glacier_data,
                       formula = turbidity ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, bs = 'ts', k=5), 
                       correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML', family = gaussian())
mod_turb_glacier3 =  gamm(data = mod_turb_glacier_data,
                       formula = turbidity ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_coverage, bs = 'ts', k=5), 
                       correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML', family = gaussian())
mod_turb_glacier2 =  gamm(data = mod_turb_glacier_data,
                       formula = turbidity ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_dist, bs = 'ts', k=5), 
                       correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML', family = gaussian())
test_bf(denominator = mod_turb_glacier1$gam, mod_turb_glacier2$gam, mod_turb_glacier3$gam, mod_turb_glacier4$gam, 
                      mod_turb_glacier5$gam, mod_turb_glacier6$gam, mod_turb_glacier7$gam, mod_turb_glacier8$gam) 
# Best model compared to null is 4: gl_area
summary(mod_turb_glacier4$gam)
capture.output(summary(mod_turb_glacier4$gam), file = "model_turbidity_glacier.txt")

# Temperature
mod_temp_glacier_data = nomis_data %>% select(`Mountain range`, latitude, longitude, water_temp, 
                                            Sample, Glacier, gl_coverage, gl_index, gl_area, gl_dist) %>% na.omit() 
mod_temp_glacier1 = gamm(data = mod_temp_glacier_data, 
                       formula = water_temp ~ s(latitude, longitude, bs='sos', k=-1, m=1), 
                       correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_temp_glacier8 = gamm(data = mod_temp_glacier_data, 
                      formula = water_temp ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_coverage, gl_dist, bs = 'ts', k=5, d=c(2,1)), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_temp_glacier7 = gamm(data = mod_temp_glacier_data, 
                      formula = water_temp ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_coverage, gl_dist, bs = 'ts', k=5), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_temp_glacier6 = gamm(data = mod_temp_glacier_data, 
                      formula = water_temp ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_dist, bs = 'ts', k=5), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_temp_glacier5 = gamm(data = mod_temp_glacier_data, 
                      formula = water_temp ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_coverage, bs = 'ts', k=5), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_temp_glacier4 = gamm(data = mod_temp_glacier_data, 
                      formula = water_temp ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, bs = 'ts', k=5), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_temp_glacier3 = gamm(data = mod_temp_glacier_data, 
                      formula = water_temp ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_coverage, bs = 'ts', k=5), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_temp_glacier2 = gamm(data = mod_temp_glacier_data, 
                      formula = water_temp ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_dist, bs = 'ts', k=5), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
test_bf(denominator = mod_temp_glacier1$gam, mod_temp_glacier2$gam, mod_temp_glacier3$gam, mod_temp_glacier4$gam, 
                      mod_temp_glacier5$gam, mod_temp_glacier6$gam, mod_temp_glacier7$gam, mod_temp_glacier8$gam)
# Best model compared to null is 6: gl_area, gl_dist
summary(mod_temp_glacier6$gam)
capture.output(summary(mod_temp_glacier6$gam), file = "model_water_temp_glacier.txt")

# DIN
mod_din_glacier_data = nomis_data %>% select(`Mountain range`, latitude, longitude, DIN, 
                                            Sample, Glacier, gl_coverage, gl_index, gl_dist, gl_area) %>% na.omit()  
mod_din_glacier1 = gamm(data = mod_din_glacier_data,
                      formula = DIN ~ s(latitude, longitude, bs='sos', k=-1, m=1), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_din_glacier8 = gamm(data = mod_din_glacier_data,
                      formula = DIN ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_coverage, gl_dist, bs = 'ts', k=5, d=c(2,1)), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_din_glacier7 = gamm(data = mod_din_glacier_data,
                     formula = DIN ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_coverage, bs = 'ts', k=5), 
                     correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_din_glacier6 = gamm(data = mod_din_glacier_data,
                     formula = DIN ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_dist, bs = 'ts', k=5), 
                     correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_din_glacier5 = gamm(data = mod_din_glacier_data,
                     formula = DIN ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_coverage, gl_dist, bs = 'ts', k=5), 
                     correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_din_glacier4 = gamm(data = mod_din_glacier_data,
                     formula = DIN ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, bs = 'ts', k=5), 
                     correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_din_glacier3 = gamm(data = mod_din_glacier_data,
                     formula = DIN ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_coverage, bs = 'ts', k=5), 
                     correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_din_glacier2 = gamm(data = mod_din_glacier_data,
                     formula = DIN ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_dist, bs = 'ts', k=5), 
                     correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
test_bf(denominator = mod_din_glacier1$gam, mod_din_glacier2$gam, mod_din_glacier3$gam, mod_din_glacier4$gam, 
                      mod_din_glacier5$gam, mod_din_glacier6$gam, mod_din_glacier7$gam, mod_din_glacier8$gam)
# Best model is null
summary(mod_din_glacier1$gam)   
capture.output(summary(mod_din_glacier1$gam), file = "model_din_glacier.txt")

# SRP
mod_srp_glacier_data = nomis_data %>% select(`Mountain range`, latitude, longitude, SRP, Feldspar,
                                           Sample, Glacier, gl_coverage, gl_index, gl_dist, gl_area) %>% na.omit() 
mod_srp_glacier1 = gamm(data = mod_srp_glacier_data,
                      formula = SRP ~ s(latitude, longitude, bs='sos', k=-1, m=1), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_srp_glacier8 = gamm(data = mod_srp_glacier_data,
                      formula = SRP ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_coverage, gl_dist, bs = 'ts', k=5, d=c(2,1)), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_srp_glacier7 = gamm(data = mod_srp_glacier_data,
                      formula = SRP ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_coverage, bs = 'ts', k=5), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_srp_glacier6 = gamm(data = mod_srp_glacier_data,
                      formula = SRP ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_dist, bs = 'ts', k=5), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_srp_glacier5 = gamm(data = mod_srp_glacier_data,
                      formula = SRP ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_coverage, gl_dist, bs = 'ts', k=5), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_srp_glacier4 = gamm(data = mod_srp_glacier_data,
                      formula = SRP ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, bs = 'ts', k=5), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_srp_glacier3 = gamm(data = mod_srp_glacier_data,
                      formula = SRP ~ s(latitude, longitude, bs='sos', k=-1, m=1) + s(gl_coverage, bs = 'ts', k=5), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_srp_glacier2 = gamm(data = mod_srp_glacier_data,
                      formula = SRP ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_dist, bs = 'ts', k=5), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
test_bf(denominator = mod_srp_glacier1$gam, mod_srp_glacier2$gam, mod_srp_glacier3$gam, mod_srp_glacier4$gam, 
                      mod_srp_glacier5$gam, mod_srp_glacier6$gam, mod_srp_glacier7$gam, mod_srp_glacier8$gam)
# Best model compared to null is 4: gl_area
summary(mod_srp_glacier4$gam)
capture.output(summary(mod_srp_glacier4$gam), file = "model_srp_glacier.txt")

nomis_data %>% filter(grepl('UP', Sample)) %$% cor.test(Feldspar, SRP, method = 'spearman')
# S = 342433, p-value = 1.919e-10
# rho = 0.4790781

# DOC
mod_doc_glacier_data = nomis_data %>% filter(grepl('UP', Sample)) %>% select(`Mountain range`, latitude, longitude, DOC, 
                                           Sample, Glacier, gl_coverage, gl_index, gl_dist, gl_area) %>% na.omit()  
mod_doc_glacier1 = gamm(data = mod_doc_glacier_data,
                      formula = DOC ~ s(latitude, longitude, bs='sos', k=-1, m=1), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_doc_glacier8 = gamm(data = mod_doc_glacier_data,
                      formula = DOC ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_coverage, gl_dist, bs = 'ts', k=5, d=c(2,1)), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_doc_glacier7 = gamm(data = mod_doc_glacier_data,
                      formula = DOC ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_coverage, bs = 'ts', k=5), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_doc_glacier6 = gamm(data = mod_doc_glacier_data,
                      formula = DOC ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_dist, bs = 'ts', k=5), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_doc_glacier5 = gamm(data = mod_doc_glacier_data,
                      formula = DOC ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_coverage, gl_dist, bs = 'ts', k=5), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_doc_glacier4 = gamm(data = mod_doc_glacier_data,
                      formula = DOC ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, bs = 'ts', k=5), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_doc_glacier3 = gamm(data = mod_doc_glacier_data,
                      formula = DOC ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_coverage, bs = 'ts', k=5), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_doc_glacier2 = gamm(data = mod_doc_glacier_data,
                      formula = DOC ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_dist, bs = 'ts', k=5), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
test_bf(denominator = mod_doc_glacier1$gam, mod_doc_glacier2$gam, mod_doc_glacier3$gam, mod_doc_glacier4$gam, 
                      mod_doc_glacier5$gam, mod_doc_glacier6$gam, mod_doc_glacier7$gam, mod_doc_glacier8$gam)
# Best model is null
summary(mod_doc_glacier1$gam)
capture.output(summary(mod_doc_glacier1$gam), file = "model_doc_glacier.txt")

################################################################
########### FIGURE 1 
################################################################

# Panel A, world map
WorldData <- ggplot2::map_data('world') %>% fortify
map = ggplot() +
  geom_map(data = WorldData, map = WorldData,
           aes(long, lat, group = group, map_id = region),
           fill = "dimgrey", colour = "dimgrey", size=0.2) +
  geom_point(nomis_data, mapping=aes(longitude, latitude, colour=`Mountain range`), size=4) + theme_bw() +
  scale_colour_brewer(palette = "Set1") + scale_colour_brewer(palette = "Paired") +
  xlab('') + ylab('') + theme(legend.title = element_blank())

# Panel B, Turbidity ~ gl_area
turb_glacier_preds = expand_grid(sample=nomis_data$Sample, gl_area=seq(-4,6,0.1))
turb_glacier_preds$latitude = vapply(turb_glacier_preds$sample, function(x) nomis_data$latitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
turb_glacier_preds$longitude = vapply(turb_glacier_preds$sample, function(x) nomis_data$longitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
preds = predict.gam(mod_turb_glacier4$gam, newdata=turb_glacier_preds, newdata.guaranteed = T, se.fit = T)
turb_glacier_preds$pred = preds$fit
turb_glacier_preds$se = preds$se.fit
turb_glacier_preds = turb_glacier_preds %>% group_by(gl_area) %>% summarise(pred = mean(pred), se = mean(se))
p1 = ggplot() +
  geom_point(data=mod_turb_glacier_data, aes(x=gl_area,y=turbidity, colour=`Mountain range`), alpha=0.2, size=2.5) +
  labs(colour = "Mountain range") +
  geom_ribbon(data=turb_glacier_preds, aes(x=gl_area, ymin=pred-se, ymax=pred+se), colour='lightgrey', alpha=.1) +
  geom_line(data=turb_glacier_preds, aes(x=gl_area,y=pred), size=1.5, colour='black') +
  geom_line(data=turb_glacier_preds, aes(x=gl_area,y=pred+se), size=1, colour='dimgrey', linetype='dashed') +
  geom_line(data=turb_glacier_preds, aes(x=gl_area,y=pred-se), size=1, colour='dimgrey', linetype='dashed') +
  scale_x_continuous(name = bquote(""*Glacier~area~(ln~ km^2)*"")) +
  scale_y_continuous(name = bquote(""*Turbidity~(ln~ NTU)*"")) +
  scale_colour_brewer(palette = "Set1") + scale_colour_brewer(palette = "Paired")+
  theme_bw() + theme(legend.position = 'none')

# Panel C, Water temp ~ gl_dist
temp_gldist_preds = expand_grid(sample=nomis_data$Sample, gl_dist=seq(-1,8,0.1))
temp_gldist_preds$latitude = vapply(temp_gldist_preds$sample, function(x) nomis_data$latitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
temp_gldist_preds$longitude = vapply(temp_gldist_preds$sample, function(x) nomis_data$longitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
temp_gldist_preds$longitude = vapply(temp_gldist_preds$sample, function(x) nomis_data$longitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
temp_gldist_preds$gl_area = vapply(temp_gldist_preds$sample, function(x) nomis_data$gl_area[nomis_data$Sample == x], FUN.VALUE = numeric(1))
preds = predict.gam(mod_temp_glacier6$gam, newdata=temp_gldist_preds, newdata.guaranteed = T, se.fit = T)
temp_gldist_preds$pred = preds$fit
temp_gldist_preds$se = preds$se.fit
temp_gldist_preds = temp_gldist_preds %>% group_by(gl_dist) %>% summarise(pred = mean(pred), se = mean(se))
p2 = ggplot() +
  geom_point(data=mod_temp_glacier_data, aes(x=gl_dist,y=water_temp, colour=`Mountain range`), alpha=0.2, size=2.5) +
  labs(colour = "Mountain range") +
  geom_ribbon(data=temp_gldist_preds, aes(x=gl_dist, ymin=pred-se, ymax=pred+se), colour='lightgrey', alpha=.1) +
  geom_line(data=temp_gldist_preds, aes(x=gl_dist,y=pred), size=1.5, colour='black') +
  geom_line(data=temp_gldist_preds, aes(x=gl_dist,y=pred+se), size=1, colour='dimgrey', linetype='dashed') +
  geom_line(data=temp_gldist_preds, aes(x=gl_dist,y=pred-se), size=1, colour='dimgrey', linetype='dashed') +
  scale_x_continuous(name = bquote(""*Distance~to~the~glacier(ln~ km)*"")) +
  scale_y_continuous(name = bquote(""*Streamwater~temperature~(ln~ degree*C)*"")) +
  scale_colour_brewer(palette = "Set1") + scale_colour_brewer(palette = "Paired")+
  theme_bw() + theme(legend.position = 'none')

# Panel D, Water temp ~ gl_area
temp_glarea_preds = expand_grid(sample=nomis_data$Sample, gl_area=seq(-4,6,0.1))
temp_glarea_preds$latitude = vapply(temp_glarea_preds$sample, function(x) nomis_data$latitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
temp_glarea_preds$longitude = vapply(temp_glarea_preds$sample, function(x) nomis_data$longitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
temp_glarea_preds$longitude = vapply(temp_glarea_preds$sample, function(x) nomis_data$longitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
temp_glarea_preds$gl_dist = vapply(temp_glarea_preds$sample, function(x) nomis_data$gl_dist[nomis_data$Sample == x], FUN.VALUE = numeric(1))
preds = predict.gam(mod_temp_glacier6$gam, newdata=temp_glarea_preds, newdata.guaranteed = T, se.fit = T)
temp_glarea_preds$pred = preds$fit
temp_glarea_preds$se = preds$se.fit
temp_glarea_preds = temp_glarea_preds %>% group_by(gl_area) %>% summarise(pred = mean(pred), se = mean(se))
p3 = ggplot() +
  geom_point(data=mod_temp_glacier_data, aes(x=gl_area,y=water_temp, colour=`Mountain range`), alpha=0.2, size=2.5) +
  labs(colour = "Mountain range") +
  geom_ribbon(data=temp_glarea_preds, aes(x=gl_area, ymin=pred-se, ymax=pred+se), colour='lightgrey', alpha=.1) +
  geom_line(data=temp_glarea_preds, aes(x=gl_area,y=pred), size=1.5, colour='black') +
  geom_line(data=temp_glarea_preds, aes(x=gl_area,y=pred+se), size=1, colour='dimgrey', linetype='dashed') +
  geom_line(data=temp_glarea_preds, aes(x=gl_area,y=pred-se), size=1, colour='dimgrey', linetype='dashed') +
  scale_x_continuous(name = bquote(""*Glacier~area(ln~ km^2)*"")) +
  scale_y_continuous(name = bquote(""*Streamwater~temperature~(ln~ degree*C)*"")) +
  scale_colour_brewer(palette = "Set1") + scale_colour_brewer(palette = "Paired")+
  theme_bw() + theme(legend.position = 'none')

# Panel E, SRP ~ gl_area
srp_glarea_preds = expand_grid(sample=nomis_data$Sample, gl_area=seq(-4,6,0.1))
srp_glarea_preds$latitude = vapply(srp_glarea_preds$sample, function(x) nomis_data$latitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
srp_glarea_preds$longitude = vapply(srp_glarea_preds$sample, function(x) nomis_data$longitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
srp_glarea_preds$longitude = vapply(srp_glarea_preds$sample, function(x) nomis_data$longitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
srp_glarea_preds$gl_coverage = vapply(srp_glarea_preds$sample, function(x) nomis_data$gl_coverage[nomis_data$Sample == x], FUN.VALUE = numeric(1))
preds = predict.gam(mod_srp_glacier4$gam, newdata=srp_glarea_preds, newdata.guaranteed = T, se.fit = T)
srp_glarea_preds$pred = preds$fit
srp_glarea_preds$se = preds$se.fit
srp_glarea_preds = srp_glarea_preds %>% group_by(gl_area) %>% summarise(pred = mean(pred), se = mean(se))
p4 = ggplot() +
  geom_point(data=mod_srp_glacier_data, aes(x=gl_area,y=SRP, colour=`Mountain range`), alpha=0.2, size=2.5) +
  labs(colour = "Mountain range") +
  geom_ribbon(data=srp_glarea_preds, aes(x=gl_area, ymin=pred-se, ymax=pred+se), colour='lightgrey', alpha=.1) +
  geom_line(data=srp_glarea_preds, aes(x=gl_area,y=pred), size=1.5, colour='black') +
  geom_line(data=srp_glarea_preds, aes(x=gl_area,y=pred+se), size=1, colour='dimgrey', linetype='dashed') +
  geom_line(data=srp_glarea_preds, aes(x=gl_area,y=pred-se), size=1, colour='dimgrey', linetype='dashed') +
  scale_x_continuous(name = bquote(""*Glacier~area(ln~ km^2)*"")) +
  scale_y_continuous(name = bquote(""*Soluble~reactive~phosphate(ln~P~l^-1)*"")) +
  scale_colour_brewer(palette = "Set1") + scale_colour_brewer(palette = "Paired") +
  theme_bw() + theme(legend.position = 'none')

sp1 = ggarrange(p1, p2, ncol = 2, nrow = 1, labels = c('b', 'c'), align = 'h')
sp2 = ggarrange(p3, p4, ncol = 2, nrow = 1, labels = c('d', 'e'), align = 'h')

p = ggarrange(map, sp1, sp2, ncol=1, nrow=3, labels=c('a','',''), align = 'v')
ggsave(p, filename = 'Fig1_params_and_glaciers.pdf', width=8.5, height = 11)


##############################################################################################################################
# 2. Chlorophyll vs Turbidity/Water temperature
##############################################################################################################################
nomis_data_chla = nomis_data %>% select(chla, latitude, longitude, turbidity, water_temp, SRP, DIN, Sample, Glacier) %>% na.omit()

mod_chla_turb1 = gamm(data = nomis_data_chla, formula = chla ~ s(latitude, longitude, bs='sos', k=-1, m=1), 
                       correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_chla_turb2 = gamm(data = nomis_data_chla, formula = chla ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(turbidity, water_temp, SRP, bs = 'ts', k=5, d=c(2,1)), 
                     correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_chla_turb3 = gamm(data = nomis_data_chla, formula = chla ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(turbidity, water_temp, DIN, bs = 'ts', k=5, d=c(2,1)), 
                     correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_chla_turb4 = gamm(data = nomis_data_chla, formula = chla ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(turbidity, SRP, DIN, bs = 'ts', k=5, d=c(2,1)), 
                     correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_chla_turb5 = gamm(data = nomis_data_chla, formula = chla ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(water_temp, SRP, DIN, bs = 'ts', k=5, d=c(2,1)), 
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
# Best model compared to null is 3
summary(mod_chla_turb3$gam)
capture.output(summary(mod_chla_turb3$gam), file = "model_chlorophyll.txt")

################################################################
########### FIGURE 2
################################################################

# Panel A, Turbidity
chla_turb_preds = expand_grid(sample=nomis_data$Sample, turbidity=seq(-3.5, 7.1, 0.5))
chla_turb_preds$latitude = vapply(chla_turb_preds$sample, function(x) nomis_data$latitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
chla_turb_preds$longitude = vapply(chla_turb_preds$sample, function(x) nomis_data$longitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
chla_turb_preds$DIN = vapply(chla_turb_preds$sample, function(x) nomis_data$DIN[nomis_data$Sample == x], FUN.VALUE = numeric(1))
chla_turb_preds$water_temp = vapply(chla_turb_preds$sample, function(x) nomis_data$water_temp[nomis_data$Sample == x], FUN.VALUE = numeric(1))
chla_turb_preds = na.omit(chla_turb_preds)
preds = predict.gam(mod_chla_turb3$gam, newdata=chla_turb_preds, newdata.guaranteed = T, se.fit = T)
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


# Panel B, Water temp
chla_wat_preds = expand_grid(sample=nomis_data$Sample, water_temp=seq(-3, 2.5, 0.1))
chla_wat_preds$latitude = vapply(chla_wat_preds$sample, function(x) nomis_data$latitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
chla_wat_preds$longitude = vapply(chla_wat_preds$sample, function(x) nomis_data$longitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
chla_wat_preds$turbidity = vapply(chla_wat_preds$sample, function(x) nomis_data$turbidity[nomis_data$Sample == x], FUN.VALUE = numeric(1))
chla_wat_preds$DIN = vapply(chla_wat_preds$sample, function(x) nomis_data$DIN[nomis_data$Sample == x], FUN.VALUE = numeric(1))
chla_wat_preds = na.omit(chla_wat_preds)
preds = predict.gam(mod_chla_turb3$gam, newdata=chla_wat_preds, newdata.guaranteed = T, se.fit = T)
chla_wat_preds$pred = preds$fit
chla_wat_preds$se = preds$se.fit
chla_wat_preds = chla_wat_preds %>% group_by(water_temp) %>% summarise(pred = mean(pred), se = mean(se))
p2 = ggplot() +
  geom_point(data=nomis_data, aes(x=water_temp,y=chla, colour=`Mountain range`), alpha=0.2, size=2.5) +
  labs(colour = "Mountain range") +
  geom_ribbon(data=chla_wat_preds, aes(x=water_temp, ymin=pred-se, ymax=pred+se), colour='lightgrey', alpha=.1) +
  geom_line(data=chla_wat_preds, aes(x=water_temp,y=pred), size=1.5, colour='black') +
  geom_line(data=chla_wat_preds, aes(x=water_temp,y=pred+se), size=1, colour='dimgrey', linetype='dashed') +
  geom_line(data=chla_wat_preds, aes(x=water_temp,y=pred-se), size=1, colour='dimgrey', linetype='dashed') +
  scale_y_continuous(name = '') +
  scale_x_continuous(name = bquote(""*Water~temperature~(ln~degree*C)*"")) +
  scale_colour_brewer(palette = "Set1") + scale_colour_brewer(palette = "Paired")+
  theme_bw() + theme(legend.position = 'none')

# Panel C, DIN
chla_din_preds = expand_grid(sample=nomis_data$Sample, DIN=seq(2, 7, 0.1))
chla_din_preds$latitude = vapply(chla_din_preds$sample, function(x) nomis_data$latitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
chla_din_preds$longitude = vapply(chla_din_preds$sample, function(x) nomis_data$longitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
chla_din_preds$turbidity = vapply(chla_din_preds$sample, function(x) nomis_data$turbidity[nomis_data$Sample == x], FUN.VALUE = numeric(1))
chla_din_preds$water_temp = vapply(chla_din_preds$sample, function(x) nomis_data$water_temp[nomis_data$Sample == x], FUN.VALUE = numeric(1))
chla_din_preds = na.omit(chla_din_preds)
preds = predict.gam(mod_chla_turb3$gam, newdata=chla_din_preds, newdata.guaranteed = T, se.fit = T, type = 'response')
chla_din_preds$pred = preds$fit
chla_din_preds$se = preds$se.fit
chla_din_preds = chla_din_preds %>% group_by(DIN) %>% summarise(pred = mean(pred), se = mean(se))
p3 = ggplot() +
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
p = ggarrange(p1, p2, p3, ncol = 3, nrow = 1, labels = c('a', 'b', 'c'))
ggsave(p, filename = 'Fig2_chla_turb_water_temp.pdf', width=8.5, height = 3)


##############################################################################################################################
# 3. CUE / TER
##############################################################################################################################
### CUE  ##########################
cue_ter_models_data = nomis_data %>% select(CUE, TER_cn, TER_cp, water_temp, chla, NP, latitude, longitude, Sample, Glacier, AP_cell, `Mountain range`) %>% na.omit

# Summary stats and CUE/TER distributions
cue_ter_models_data = na.omit(cue_ter_models_data)

cue_model1 = gamm(data = cue_ter_models_data, formula = CUE ~ s(latitude, longitude, bs='sos', k=-1, m=1),
                  correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
cue_model2 = gamm(data = cue_ter_models_data, formula = CUE ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(chla, water_temp, NP, bs='ts', k=5, d=c(2,1)),
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
# Best model compared to null is 6
summary(cue_model6$gam)
capture.output(summary(cue_model6$gam), file = "model_cue.txt")


### TER CN  ##########################
tercn_model1 = gamm(data = cue_ter_models_data, formula = TER_cn ~ s(latitude, longitude, bs='sos', k=-1, m=1),
                  correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
tercn_model2 = gamm(data = cue_ter_models_data, formula = TER_cn ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(chla, water_temp, NP, bs='ts', k=5, d=c(2,1)),
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
summary(tercn_model1$gam)
capture.output(summary(tercn_model1$gam), file = "model_tercn.txt")

### TER CP ##########################
tercp_model1 = gamm(data = cue_ter_models_data, formula = TER_cp ~ s(latitude, longitude, bs='sos', k=-1, m=1),
                    correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
tercp_model2 = gamm(data = cue_ter_models_data, formula = TER_cp ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(chla, water_temp, NP, bs='ts', k=5, d=c(2,1)),
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
# Best model is 6
summary(tercp_model6$gam)
capture.output(summary(tercp_model6$gam), file = "model_tercp.txt")

# PLOTTING
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

p = ggarrange(p1, p2, ncol = 2, nrow = 1, labels = c('a', 'b'), align='h')
ggsave(p, filename = 'Fig4_cue_ter.pdf', width=7, height = 3.5)


### CUE AP ########################################################
cue_ap_data = nomis_data %>% select(CUE, latitude, longitude, Sample, Glacier, AP_cell, CN, `Mountain range`) %>% na.omit()

cue_ap_model = gamm(data = cue_ap_data, formula = CUE ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(AP_cell, bs='ts', k=5), 
                   correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML', select = T)
capture.output(summary(cue_ap_model$gam), file = "model_ap_cell_cue.txt")


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
ggsave(p, filename = 'Fig5_AP_CUE.pdf', width=5, height = 5)

