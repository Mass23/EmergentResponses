library(ggplot2)
library(tidyverse)
library(ggpubr)
library(reshape2)
library(mgcv)
library(performance)
library(maps)
library(ggridges)

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
                                                                             SRP = mean(n3_srp..ug.l.1., na.rm=T),
                                                                             BP = mean(bp..ngC.g.1.h.1., na.rm=T),
                                                                             respiration = mean(respiration..mg.02.g.1.h.1., na.rm=T),
                                                                             EPS = mean(eps..ugC.g.1., na.rm = T))

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

kb = 273.15
nomis_data$rel_temp = (1 / (0.0000862 * (kb + nomis_data$water_temp))) - (1 / (0.0000862 * (kb + 10)))

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
bp_dilution_factor = 2
resp_constant = 1.2

nomis_data$BP = nomis_data$BP / bp_dilution_factor # normalise with the dilution factor

model_tercn <- lm(data = nomis_data %>% select(LAP, BG, NAG) %>% na.omit, formula = log(LAP + BG) ~ log(LAP + NAG))
nomis_data$TER_cn <- (((nomis_data$LAP+nomis_data$BG)/(nomis_data$LAP+nomis_data$NAG))*bcn) / as.double(coef(model_tercn)['(Intercept)'])

model_tercp <- lm(data = nomis_data %>% select(LAP, BG, AP) %>% na.omit, formula = log(LAP + BG) ~ log(AP))
nomis_data$TER_cp<-(((nomis_data$LAP+nomis_data$BG)/(nomis_data$AP))*bcp) / as.double(coef(model_tercp)['(Intercept)'])

nomis_data$CUE<-(((ass_efficiency * bcn )/ (nomis_data$TER_cn))*((ass_efficiency * bcp)/(nomis_data$TER_cp)))^0.5
nomis_data$BGE = (nomis_data$BP/1e9) / ((nomis_data$BP/1e9) + (nomis_data$respiration*resp_constant/1e3))

quantile(nomis_data$CUE, na.rm = T)
#           0%         25%         50%         75%        100% 
#  0.005414022 0.110222971 0.150880527 0.204943998 0.471144861 
quantile(nomis_data$BGE, na.rm = T)
#           0%          25%          50%          75%         100% 
# 0.0006078487 0.0804200081 0.1898494102 0.3545086063 0.8629427993 

ggplot(nomis_data) + geom_histogram(aes(x=CUE, fill='red'), alpha=0.5) + 
  geom_histogram(aes(x=BGE, fill='blue'), alpha=0.5) + theme_bw() + xlab('Value') + ylab('Count') +
  labs(fill='') + scale_fill_discrete(labels = expression('CUE'[BGE],'CUE'[EEA])) + 
  geom_vline(xintercept = median(nomis_data$BGE, na.rm = T), colour='blue') + 
  geom_vline(xintercept = median(nomis_data$CUE, na.rm = T), colour='red')
ggsave('Plots/Fig_S1.pdf')

# outlier analysis
nomis_data$CUE_zscore = (nomis_data$CUE-mean(nomis_data$CUE, na.rm=T))/sd(nomis_data$CUE, na.rm=T)
nomis_data$CUE_outlier = ifelse(abs(nomis_data$CUE_zscore) > 3, 'Yes', 'No')
nomis_data$TERcn_zscore = (nomis_data$TER_cn-mean(nomis_data$TER_cn, na.rm=T))/sd(nomis_data$TER_cn, na.rm=T)
nomis_data$TERcn_outlier = ifelse(abs(nomis_data$TERcn_zscore) > 3, 'Yes', 'No')
nomis_data$TERcp_zscore = (nomis_data$TER_cp-mean(nomis_data$TER_cp, na.rm=T))/sd(nomis_data$TER_cp, na.rm=T)
nomis_data$TERcp_outlier = ifelse(abs(nomis_data$TERcp_zscore) > 3, 'Yes', 'No')
nomis_data$BGE_zscore = (nomis_data$BGE-mean(nomis_data$BGE, na.rm=T))/sd(nomis_data$BGE, na.rm=T)
nomis_data$BGE_outlier = ifelse(abs(nomis_data$BGE_zscore) > 3, 'Yes', 'No')

nomis_data$CUE[nomis_data$CUE_outlier == 'Yes'] = NA
nomis_data$TER_cn[nomis_data$TERcn_outlier == 'Yes'] = NA
nomis_data$TER_cp[nomis_data$TERcp_outlier == 'Yes'] = NA
nomis_data$BGE[nomis_data$BGE_outlier == 'Yes'] = NA

##############################################################################################################################
# 3. Summary stats
##############################################################################################################################
data_summary = nomis_data %>% ungroup() %>% select(gl_index, gl_area, gl_dist, gl_coverage, turbidity, conductivity, water_temp, respiration, BP,
                                                   chla, DOC, DIN, SRP, NP, CN, CP, CUE, TER_cn, TER_cp, LAP, BG, AG, AP, NAG) %>% 
  summarise_all(quantile, na.rm=T, probs=c(0.25,0.5,0.75)) %>% t %>% as.data.frame

data_summary_mountain_ranges = nomis_data %>% ungroup() %>% select(gl_index, gl_area, gl_dist, gl_coverage, turbidity, conductivity, water_temp, respiration, BP,
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
# 4. Figure 1
##############################################################################################################################
# Panel A, world map
colours = c('#2A83C7', '#52B7E0', '#61BC6A', '#1C9C31', '#E09882', '#E0422E',
            '#FFC689', '#E88B33', '#B287C4', '#9237B0', '#DECC45')
  
WorldData <- ggplot2::map_data('world') %>% fortify
map = ggplot() +
  geom_map(data = WorldData, map = WorldData,
           aes(long, lat, group = group, map_id = region),
           fill = "dimgrey", colour = "dimgrey", size=0.2) +
  geom_point(nomis_data, mapping=aes(longitude, latitude, colour=`Mountain range`), size=4) + theme_bw() +
  scale_colour_manual(values = colours) +
  xlab('') + ylab('') + theme(legend.title = element_blank())

melt_data = nomis_data %>% select(`Mountain range`, water_temp, turbidity, DIN, SRP, chla) %>% melt() %>% na.omit()
melt_data$variable = factor(melt_data$variable, labels = c("Water~temp.~(degree~C)", 
                                                           'Turbidity~(NTU)',
                                                           "DIN~(ppb)",
                                                           "SRP~(ppb)",
                                                           "Chla~(mu*g~g^-1~DM)"))
melt_data$value = map_dbl(1:nrow(melt_data), function(i) melt_data$value[i] + min(melt_data$value[(melt_data$value > 0) & (melt_data$variable == melt_data$variable[i])])/2)


expeditions = c('Alaska Range', 'Caucasus Mountains', 'Chilean Andes', 'Ecuadorian Andes', 'European Alps', 'Himalayas',
                       'Pamir & Tien Shan', 'Rwenzori Mountains', 'Scandinavian Mountains', 'Southern Alps', 'Southwest Greenland')
all_sp1 = ggplot(melt_data %>% filter(`Mountain range` != 'Rwenzori Mountains'), 
                 aes(y=factor(`Mountain range`, levels=rev(expeditions)), x=value, 
                               colour=`Mountain range`, fill=`Mountain range`, alpha=0.5)) + 
  geom_density_ridges(scale = 1.2, quantile_lines=T, quantiles = 2) + theme_bw() + 
  scale_colour_manual(values = colours[-8]) + facet_grid(~variable, scales = 'free_x', labeller = label_parsed) +
  scale_fill_manual(values = colours[-8]) + xlab('') + ylab('') + scale_x_continuous(guide = guide_axis(check.overlap = TRUE), trans = 'log10') +
  theme(legend.position = 'none', axis.text.y = element_blank(), panel.margin = unit(0.5, 'lines')) + theme(plot.margin = unit(c(0,0.5,0,0.5), 'lines'))

all_sp2 = ggplot(melt_data %>% filter(`Mountain range` == 'Rwenzori Mountains')) + 
  geom_blank(data = melt_data %>% filter(`Mountain range` != 'Rwenzori Mountains'), aes(x=value)) +
  geom_point(aes(y=factor(`Mountain range`, levels=rev(expeditions)), x=value+(min(melt_data$value[melt_data$value > 0])/2), colour=colours[8], fill=colours[8], alpha=0.5)) +
  theme_bw() + facet_grid(~variable, scales = 'free_x', labeller = label_parsed) + xlab('') + ylab('') + scale_x_log10() + 
  theme(legend.position = 'none', axis.text.y = element_blank(), 
        strip.background = element_blank(), strip.text.x = element_blank(),
        axis.ticks.y = element_blank(), panel.margin = unit(0.5, 'lines')) + theme(plot.margin = unit(c(0,0.5,0,0.5), 'lines'))

p = ggarrange(map, all_sp1, nrow = 2, heights = c(0.5, 0.5), labels = c('A', 'B'))
ggsave('Plots/Fig1_map_expe.pdf', p, width = 7.5, height = 6.5)

p = ggarrange(all_sp1 + theme(axis.text.x = element_blank(), axis.ticks = element_blank()), all_sp2, nrow = 2, heights = c(0.88, 0.12), labels = c('A', 'B'))
ggsave('Plots/SFig1_expe_dist.pdf', p, width = 8, height = 7)

##############################################################################################################################
# 5. Variable transformation
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
nomis_data$respiration = log_const(nomis_data$respiration)

nomis_data$DOC = log_const(nomis_data$DOC)
nomis_data$DIN = log_const(nomis_data$DIN)
nomis_data$SRP = log_const(nomis_data$SRP)

nomis_data$NP = log_const(nomis_data$NP)
nomis_data$CN = log_const(nomis_data$CN)
nomis_data$CP = log_const(nomis_data$CP)
nomis_data$CUE = log_const(nomis_data$CUE)
nomis_data$TER_cn = log_const(nomis_data$TER_cn)
nomis_data$TER_cp = log_const(nomis_data$TER_cp)
nomis_data$BGE = log_const(nomis_data$BGE)

nomis_data$LAP = log(nomis_data$LAP)
nomis_data$NAG = log(nomis_data$NAG)
nomis_data$BG = log(nomis_data$BG)
nomis_data$AG = log(nomis_data$AG)
nomis_data$AP = log(nomis_data$AP)
nomis_data$BP = log(nomis_data$BP)

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
colours = c('#52B7E0', '#61BC6A', '#1C9C31', '#E09882', '#E0422E',
            '#FFC689', '#E88B33', '#B287C4', '#9237B0', '#DECC45')
expeditions = c('Caucasus Mountains', 'Chilean Andes', 'Ecuadorian Andes', 'European Alps', 'Himalayas',
                'Pamir & Tien Shan', 'Rwenzori Mountains', 'Scandinavian Mountains', 'Southern Alps', 'Southwest Greenland')

cue_p1 = ggplot(nomis_data %>% filter(`Mountain range` != 'Alaska Range'), aes(x=factor(`Mountain range`, levels=rev(expeditions)), y=exp(CUE), fill=`Mountain range`)) + geom_boxplot()+ 
  theme_bw() + xlab('') + ylab('CUE') + theme(axis.text.x = element_blank(), legend.position = 'none') + scale_y_log10() + scale_fill_manual(values = colours) 
cue_p2 = ggplot(nomis_data, aes(y=exp(CUE),x=..scaled..)) + geom_density() + geom_hline(yintercept = median(exp(nomis_data$CUE), na.rm=T), colour='red') + theme_bw() + ylab('') + xlab('')+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_y_log10() + scale_x_continuous(breaks = c(0,0.5,1))
cue_p = ggarrange(cue_p1, cue_p2, ncol = 2, nrow = 1, widths = c(0.75, 0.25), align = 'h', labels = c('A',''))

tercn_p1 = ggplot(nomis_data %>% filter(`Mountain range` != 'Alaska Range'), aes(x=factor(`Mountain range`, levels=rev(expeditions)), y=exp(TER_cn), fill=`Mountain range`)) + geom_boxplot() + 
  theme_bw() + xlab('') + ylab('TER C:N') + theme(axis.text.x = element_blank(), legend.position = 'none') + scale_y_log10() + scale_fill_manual(values = colours) 
tercn_p2 = ggplot(nomis_data, aes(y=exp(TER_cn),x=..scaled..)) + geom_density() + geom_hline(yintercept = median(exp(nomis_data$TER_cn), na.rm=T), colour='red') + theme_bw() + ylab('') + xlab('')+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_y_log10() + scale_x_continuous(breaks = c(0,0.5,1))
tercn_p = ggarrange(tercn_p1, tercn_p2, ncol = 2, nrow = 1, widths = c(0.75, 0.25), align = 'h', labels = c('B',''))

tercp_p1 = ggplot(nomis_data %>% filter(`Mountain range` != 'Alaska Range'), aes(x=factor(`Mountain range`, levels=rev(expeditions)), y=exp(TER_cp), fill=`Mountain range`)) + geom_boxplot() + 
  theme_bw() + xlab('') + ylab('TER C:P') + theme(axis.text.x = element_blank(),legend.position = 'right') + scale_y_log10() + scale_fill_manual(values = colours) 
tercp_p2 = ggplot(nomis_data, aes(y=exp(TER_cp),x=..scaled..)) + geom_density() + geom_hline(yintercept = median(exp(nomis_data$TER_cp), na.rm=T), colour='red') + theme_bw() + ylab('') +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_y_log10() + scale_x_continuous(breaks = c(0,0.5,1)) + xlab('Scaled density')
leg = get_legend(tercp_p1)
tercp_p = ggarrange(tercp_p1 + theme(legend.position = 'none'), tercp_p2, ncol = 2, nrow = 1, widths = c(0.75, 0.25), align = 'h', labels = c('C',''))

p1 = ggarrange(cue_p, tercn_p, tercp_p, nrow = 3, ncol = 1, heights = c(0.333,0.333,0.333))
p = ggarrange(p1, leg, ncol = 2, widths = c(0.73, 0.27))
ggsave(plot = p, 'Plots/Fig4_cue_ter_distributions.pdf', width = 7.5, height = 6.5)


