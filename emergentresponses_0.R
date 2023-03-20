library(ggplot2)
library(tidyverse)
library(ggpubr)
library(reshape2)
library(mgcv)
library(performance)
library(maps)

setwd('~/Desktop/emergentresponses')
set.seed(23)
nomis_data = read.csv('Data/NOMIS_01_2023_FULL.csv')

dir.create('Plots')
dir.create('Statistics')

##############################################################################################################################
# 1. Data preparation
##############################################################################################################################
nomis_data$Sample = map_chr(nomis_data$patch, function(x) paste0(strsplit(x, '_')[[1]][1:2], collapse = '_'))

nomis_data$n6_no2..ug.l.1.[is.na(nomis_data$n6_no2..ug.l.1.)] = 0
nomis_data = rename(nomis_data, `Mountain range` = Mountain.range)

nomis_data = nomis_data %>% group_by(Sample, `Mountain range`) %>% summarise(latitude = mean(lat_sp..DD., na.rm=T),
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

nomis_data$Glacier = map_chr(nomis_data$Sample, function(x) strsplit(x ,'_')[[1]][1])
nomis_data$gl_index = sqrt(nomis_data$gl_area) / (sqrt(nomis_data$gl_area) + nomis_data$gl_dist)

min_data = read.csv('Data/Final_geology.csv')
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
write.csv(data_summary, file = 'Statistics/summary_stats.csv', quote = F, row.names = T)
write.csv(data_summary_mountain_ranges, file = 'Statistics/summary_stats_mountaing_ranges.csv', quote = F, row.names = T)

##############################################################################################################################
# 4. Variable transformation
##############################################################################################################################
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

write.csv(nomis_data, 'Data/NOMIS_01_2023_FULL_preprocessed.csv')

##############################################################################################################################
# 5. Plotting
##############################################################################################################################
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
ggsave(plot = p, 'Plots/Fig3_cue_ter_distributions.pdf', width = 7, height = 6)

