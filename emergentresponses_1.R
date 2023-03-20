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
# 1. Model Turbidity
##############################################################################################################################
mod_turb_glacier_data = nomis_data %>% select(`Mountain range`, latitude, longitude, turbidity, 
                      Sample, Glacier, gl_coverage, gl_index, gl_area, gl_dist) %>% na.omit()  
mod_turb_glacier1 =  gamm(data = mod_turb_glacier_data,
                        formula = turbidity ~ s(latitude, longitude, bs='sos', k=-1, m=1), 
                        correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_turb_glacier8 =  gamm(data = mod_turb_glacier_data,
                      formula = turbidity ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_coverage, gl_dist, bs = 'ts', k=5), 
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
# Best model compared to null is 4: gl_area (BF = 3.06e+12)
capture.output(summary(mod_turb_glacier4$gam), file = "Statistics/model_turbidity_glacier.txt")

##############################################################################################################################
# 2. Model Temperature
##############################################################################################################################
mod_temp_glacier_data = nomis_data %>% select(`Mountain range`, latitude, longitude, water_temp, 
                                            Sample, Glacier, gl_coverage, gl_index, gl_area, gl_dist) %>% na.omit() 
mod_temp_glacier1 = gamm(data = mod_temp_glacier_data, 
                       formula = water_temp ~ s(latitude, longitude, bs='sos', k=-1, m=1), 
                       correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_temp_glacier8 = gamm(data = mod_temp_glacier_data, 
                      formula = water_temp ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_coverage, gl_dist, bs = 'ts', k=5), 
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
# Best model compared to null is 6: gl_area, gl_dist (BF = 1.84e+40)
capture.output(summary(mod_temp_glacier6$gam), file = "Statistics/model_water_temp_glacier.txt")


##############################################################################################################################
# 4. Dissolved inorganic nitrogen
##############################################################################################################################
mod_din_glacier_data = nomis_data %>% select(`Mountain range`, latitude, longitude, DIN, 
                                            Sample, Glacier, gl_coverage, gl_index, gl_dist, gl_area) %>% na.omit()  
mod_din_glacier1 = gamm(data = mod_din_glacier_data,
                      formula = DIN ~ s(latitude, longitude, bs='sos', k=-1, m=1), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_din_glacier8 = gamm(data = mod_din_glacier_data,
                      formula = DIN ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_coverage, gl_dist, bs = 'ts', k=5), 
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
capture.output(summary(mod_din_glacier1$gam), file = "Statistics/model_din_glacier.txt")

##############################################################################################################################
# 4. Soluble reactive phosphate
##############################################################################################################################
mod_srp_glacier_data = nomis_data %>% select(`Mountain range`, latitude, longitude, SRP, Feldspar,
                                           Sample, Glacier, gl_coverage, gl_index, gl_dist, gl_area) %>% na.omit() 
mod_srp_glacier1 = gamm(data = mod_srp_glacier_data,
                      formula = SRP ~ s(latitude, longitude, bs='sos', k=-1, m=1), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_srp_glacier8 = gamm(data = mod_srp_glacier_data,
                      formula = SRP ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_coverage, gl_dist, bs = 'ts', k=5), 
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
# Best model compared to null is 4: gl_area (BF = 3.46)
capture.output(summary(mod_srp_glacier4$gam), file = "Statistics/model_srp_glacier.txt")

# correlation with feldspar content of the mineralogy
feld = nomis_data %>% filter(grepl('UP', Sample)) %>% pull(Feldspar)
srp = nomis_data %>% filter(grepl('UP', Sample)) %>% pull(SRP)
cor.test(feld, srp, method = 'spearman')
# S = 342433, p-value = 1.919e-10
# rho = 0.4790781

##############################################################################################################################
# 5. Dissolved organic carbon
##############################################################################################################################
mod_doc_glacier_data = nomis_data %>% filter(grepl('UP', Sample)) %>% select(`Mountain range`, latitude, longitude, DOC, 
                                           Sample, Glacier, gl_coverage, gl_index, gl_dist, gl_area) %>% na.omit()  
mod_doc_glacier1 = gamm(data = mod_doc_glacier_data,
                      formula = DOC ~ s(latitude, longitude, bs='sos', k=-1, m=1), 
                      correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_doc_glacier8 = gamm(data = mod_doc_glacier_data,
                      formula = DOC ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_coverage, gl_dist, bs = 'ts', k=5), 
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
capture.output(summary(mod_doc_glacier1$gam), file = "Statistics/model_doc_glacier.txt")

################################################################
########### Plotting
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
ggsave(p, filename = 'Plots/Fig1_params_and_glaciers.pdf', width=8.5, height = 11)


