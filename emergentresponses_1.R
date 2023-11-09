library(ggplot2)
library(tidyverse)
library(ggpubr)
library(reshape2)
library(mgcv)
library(performance)
library(maps)
library(corrplot)

setwd('~/Desktop/emergentresponses')
nomis_data = read.csv('Data/NOMIS_01_2023_FULL_preprocessed.csv')
nomis_data = rename(nomis_data, `Mountain range` = Mountain.range)

corrplot(cor(nomis_data %>% select(turbidity,chla,water_temp,DOC,NP,CP,CN,TER_cn,TER_cp,CUE,BGE) %>% na.omit(), method = 'spearman'), 
         type = 'upper', order = 'hclust', diag = F, addCoef.col = 'black')

colours = c('#52B7E0', '#61BC6A', '#1C9C31', '#E09882', '#E0422E',
            '#FFC689', '#E88B33', '#B287C4', '#9237B0', '#DECC45')
expeditions = c('Caucasus Mountains', 'Chilean Andes', 'Ecuadorian Andes', 'European Alps', 'Himalayas',
                'Pamir & Tien Shan', 'Rwenzori Mountains', 'Scandinavian Mountains', 'Southern Alps', 'Southwest Greenland')

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
# 3. Conductivity
##############################################################################################################################
mod_cond_glacier_data = nomis_data %>% select(`Mountain range`, latitude, longitude, conductivity, 
                                              Sample, Glacier, gl_coverage, gl_index, gl_area, gl_dist) %>% na.omit() 
mod_cond_glacier1 = gamm(data = mod_cond_glacier_data, 
                         formula = conductivity ~ s(latitude, longitude, bs='sos', k=-1, m=1), 
                         correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_cond_glacier8 = gamm(data = mod_cond_glacier_data, 
                         formula = conductivity ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_coverage, gl_dist, bs = 'ts', k=5), 
                         correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_cond_glacier7 = gamm(data = mod_cond_glacier_data, 
                         formula = conductivity ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_coverage, gl_dist, bs = 'ts', k=5), 
                         correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_cond_glacier6 = gamm(data = mod_cond_glacier_data, 
                         formula = conductivity ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_dist, bs = 'ts', k=5), 
                         correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_cond_glacier5 = gamm(data = mod_cond_glacier_data, 
                         formula = conductivity ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, gl_coverage, bs = 'ts', k=5), 
                         correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_cond_glacier4 = gamm(data = mod_cond_glacier_data, 
                         formula = conductivity ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_area, bs = 'ts', k=5), 
                         correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_cond_glacier3 = gamm(data = mod_cond_glacier_data, 
                         formula = conductivity ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_coverage, bs = 'ts', k=5), 
                         correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
mod_cond_glacier2 = gamm(data = mod_cond_glacier_data, 
                         formula = conductivity ~ s(latitude, longitude, bs='sos', k=-1, m=1) + te(gl_dist, bs = 'ts', k=5), 
                         correlation = corCAR1(form = ~ 1 | Glacier), method = 'REML')
test_bf(denominator = mod_cond_glacier1$gam, mod_cond_glacier2$gam, mod_cond_glacier3$gam, mod_cond_glacier4$gam, 
        mod_cond_glacier5$gam, mod_cond_glacier6$gam, mod_cond_glacier7$gam, mod_cond_glacier8$gam)
# Best model compared to null is 3: gl_coverage 
capture.output(summary(mod_cond_glacier3$gam), file = "Statistics/model_condu_glacier.txt")

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
# 5. Soluble reactive phosphate
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

##############################################################################################################################
# 6. Dissolved organic carbon
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

# correlation with feldspar content of the mineralogy
feld = nomis_data %>% filter(grepl('UP', Sample)) %>% pull(Feldspar)
srp = nomis_data %>% filter(grepl('UP', Sample)) %>% pull(SRP)
cor.test(feld, srp, method = 'spearman')
# S = 342433, p-value = 1.919e-10
# rho = 0.4790781

################################################################
########### Plotting
################################################################
# Panel A, Turbidity
turb_gla_preds = expand_grid(sample=nomis_data$Sample, gl_area=seq(-4.5, 6.1, 0.5))
turb_gla_preds$latitude = vapply(turb_gla_preds$sample, function(x) nomis_data$latitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
turb_gla_preds$longitude = vapply(turb_gla_preds$sample, function(x) nomis_data$longitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
turb_gla_preds = na.omit(turb_gla_preds)
preds = predict.gam(mod_turb_glacier4$gam, newdata=turb_gla_preds, newdata.guaranteed = T, se.fit = T)
turb_gla_preds$pred = preds$fit
turb_gla_preds$se = preds$se.fit
turb_gla_preds = turb_gla_preds %>% group_by(gl_area) %>% summarise(pred = median(pred), se = median(se))
p1 = ggplot() +
  geom_point(data=nomis_data, aes(x=gl_area,y=turbidity, colour=factor(`Mountain range`, levels=rev(expeditions))), alpha=0.2, size=2.5) + labs(colour = "Mountain range") +
  geom_ribbon(data=turb_gla_preds, aes(x=gl_area, ymin=pred-se, ymax=pred+se), colour='lightgrey', alpha=.1) +
  geom_line(data=turb_gla_preds, aes(x=gl_area,y=pred), size=1.5, colour='black') +
  geom_line(data=turb_gla_preds, aes(x=gl_area,y=pred+se), size=1, colour='dimgrey', linetype='dashed') +
  geom_line(data=turb_gla_preds, aes(x=gl_area,y=pred-se), size=1, colour='dimgrey', linetype='dashed') +
  scale_y_continuous(name = bquote(""*Turbidity~(ln~ NTU)*"")) +
  scale_x_continuous(name = "") +
  scale_fill_manual(values = colours) + theme_bw() + theme(legend.position = 'none')

# Panel B, stream temp
temp_gld_preds = expand_grid(sample=nomis_data$Sample, gl_dist=seq(-1.5, 8.5, 0.5))
temp_gld_preds$latitude = vapply(temp_gld_preds$sample, function(x) nomis_data$latitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
temp_gld_preds$longitude = vapply(temp_gld_preds$sample, function(x) nomis_data$longitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
temp_gld_preds$gl_area = vapply(temp_gld_preds$sample, function(x) nomis_data$gl_area[nomis_data$Sample == x], FUN.VALUE = numeric(1))
temp_gld_preds = na.omit(temp_gld_preds)
preds = predict.gam(mod_temp_glacier6$gam, newdata=temp_gld_preds, newdata.guaranteed = T, se.fit = T)
temp_gld_preds$pred = preds$fit
temp_gld_preds$se = preds$se.fit
temp_gld_preds = temp_gld_preds %>% group_by(gl_dist) %>% summarise(pred = median(pred), se = median(se))
p2 = ggplot() +
  geom_point(data=nomis_data, aes(x=gl_dist,y=water_temp, colour=factor(`Mountain range`, levels=rev(expeditions))), alpha=0.2, size=2.5) + labs(colour = "Mountain range") +
  geom_ribbon(data=temp_gld_preds, aes(x=gl_dist, ymin=pred-se, ymax=pred+se), colour='lightgrey', alpha=.1) +
  geom_line(data=temp_gld_preds, aes(x=gl_dist,y=pred), size=1.5, colour='black') +
  geom_line(data=temp_gld_preds, aes(x=gl_dist,y=pred+se), size=1, colour='dimgrey', linetype='dashed') +
  geom_line(data=temp_gld_preds, aes(x=gl_dist,y=pred-se), size=1, colour='dimgrey', linetype='dashed') +
  scale_y_continuous(name = bquote(""*Water~temperature~(ln~degree~C)*"")) +
  scale_x_continuous(name = bquote(""*Distance~to~the~glacier~(ln~m)*"")) +
  scale_fill_manual(values = colours) + theme_bw() + theme(legend.position = 'none')

# Panel C, stream temp
temp_gla_preds = expand_grid(sample=nomis_data$Sample, gl_area=seq(-4.5, 6.1, 0.5))
temp_gla_preds$latitude = vapply(temp_gla_preds$sample, function(x) nomis_data$latitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
temp_gla_preds$longitude = vapply(temp_gla_preds$sample, function(x) nomis_data$longitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
temp_gla_preds$gl_dist = vapply(temp_gla_preds$sample, function(x) nomis_data$gl_dist[nomis_data$Sample == x], FUN.VALUE = numeric(1))
temp_gla_preds = na.omit(temp_gla_preds)
preds = predict.gam(mod_temp_glacier6$gam, newdata=temp_gla_preds, newdata.guaranteed = T, se.fit = T)
temp_gla_preds$pred = preds$fit
temp_gla_preds$se = preds$se.fit
temp_gla_preds = temp_gla_preds %>% group_by(gl_area) %>% summarise(pred = median(pred), se = median(se))
p3 = ggplot() +
  geom_point(data=nomis_data, aes(x=gl_area,y=water_temp, colour=factor(`Mountain range`, levels=rev(expeditions))), alpha=0.2, size=2.5) + labs(colour = "Mountain range") +
  geom_ribbon(data=temp_gla_preds, aes(x=gl_area, ymin=pred-se, ymax=pred+se), colour='lightgrey', alpha=.1) +
  geom_line(data=temp_gla_preds, aes(x=gl_area,y=pred), size=1.5, colour='black') +
  geom_line(data=temp_gla_preds, aes(x=gl_area,y=pred+se), size=1, colour='dimgrey', linetype='dashed') +
  geom_line(data=temp_gla_preds, aes(x=gl_area,y=pred-se), size=1, colour='dimgrey', linetype='dashed') +
  scale_y_continuous(name = bquote(""*Water~temperature~(ln~degree~C)*"")) +
  scale_x_continuous(name = "") +
  scale_fill_manual(values = colours) + theme_bw() + theme(legend.position = 'none')

# Panel D, srp
srp_gla_preds = expand_grid(sample=nomis_data$Sample, gl_area=seq(-4.5, 6.1, 0.5))
srp_gla_preds$latitude = vapply(srp_gla_preds$sample, function(x) nomis_data$latitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
srp_gla_preds$longitude = vapply(srp_gla_preds$sample, function(x) nomis_data$longitude[nomis_data$Sample == x], FUN.VALUE = numeric(1))
srp_gla_preds = na.omit(srp_gla_preds)
preds = predict.gam(mod_srp_glacier4$gam, newdata=srp_gla_preds, newdata.guaranteed = T, se.fit = T)
srp_gla_preds$pred = preds$fit
srp_gla_preds$se = preds$se.fit
srp_gla_preds = srp_gla_preds %>% group_by(gl_area) %>% summarise(pred = median(pred), se = median(se))
p4 = ggplot() +
  geom_point(data=nomis_data, aes(x=gl_area,y=water_temp, colour=factor(`Mountain range`, levels=rev(expeditions))), alpha=0.2, size=2.5) + labs(colour = "Mountain range") +
  geom_ribbon(data=srp_gla_preds, aes(x=gl_area, ymin=pred-se, ymax=pred+se), colour='lightgrey', alpha=.1) +
  geom_line(data=srp_gla_preds, aes(x=gl_area,y=pred), size=1.5, colour='black') +
  geom_line(data=srp_gla_preds, aes(x=gl_area,y=pred+se), size=1, colour='dimgrey', linetype='dashed') +
  geom_line(data=srp_gla_preds, aes(x=gl_area,y=pred-se), size=1, colour='dimgrey', linetype='dashed') +
  scale_y_continuous(name = bquote(""*SRP~(ln~ P~l^-1)*"")) +
  scale_x_continuous(name = bquote(""*Glacier~surface~area~(ln~km ^2)*"")) +
  scale_fill_manual(values = colours) + theme_bw() + theme(legend.position = 'none')

ggarrange(p1, p3, p4, nrow = 3, ncol = 1, labels = c('A', 'B', 'C'), align = 'v')
ggsave('Plots/Fig2_glacier_area.pdf', width = 3.5, height = 8)

p2
ggsave('Plots/Fig_S2.pdf', width = 4, height = 4)














