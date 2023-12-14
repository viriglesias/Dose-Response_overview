#This script estimates d-s curves for soy yields and makes figure 8

library(tidyverse)
library(tidymv)
library(mgcv)
# install.packages("mgcViz")
library(mgcViz)
# install.packages('directlabels')
library(directlabels)
library(ggpubr)
library(sf)
library(ggthemes)
library(foreach)
#USE ME
source('/Users/viig7608/Desktop/CCP/Agric/code/all_functions.R')

# variable <- list('summer_tmean', 'summer_prec', 'summer_vpd')

plot_gam_crop <- function(variable, win, min_yr, what){
  new <- function(variable){
    a <- paste(strsplit(variable, '_')[[1]][1], ifelse(strsplit(variable, '_')[[1]][2] %in% 'vpd',
                                                       'VPD',
                                                       ifelse(strsplit(variable, '_')[[1]][2] %in% 'prec',
                                                              'precipitation',
                                                              ifelse(strsplit(variable, '_')[[1]][2] %in% 'tmean',
                                                                     'temperature', NA))), sep = ' ')
    return(a)
  }
  x <- read.csv(paste0('/Users/viig7608/Desktop/CCP/Overview/Soybeans/Data/Processed/soy_county_diff/soy_county_diff_1980_', win, '_yr.csv'))[,-1]
  x <- filter(x, !is.na(nonirrigated_mean_diff) & !is.na(irrigated_mean_diff))
  colnames(x)[1] <- 'id'
  climate <- read.csv(paste0('/Users/viig7608/Desktop/CCP/Agric/All_seasons/', variable, '/', variable, '_diff_', min_yr, '-2020_', win, '_yr.csv'))[,-1]
  climate <- filter(climate, Year>= min_yr)
  fips <- read.csv('/Users/viig7608/Desktop/CCP/state-geocodes-v2016.csv')#fips and state names
  climate <- left_join(climate, fips) 
  climate <- climate %>% 
    mutate(STATE = toupper(State), 
           COUNTY = toupper(NAME),
           id = paste(STATE, COUNTY, sep = '_'))#create unique identifier that matches corn data
  complete <- inner_join(x, climate, by = c('id', 'Year'))
  complete <- complete %>% 
    mutate(period = ifelse(Year<1979, '1950-1980',
                           '1980-2019')) %>% 
    dplyr::select(c(irrigated_mean_diff,
             nonirrigated_mean_diff,
             sym(paste0('Mean_', sub('.*_', '', variable))), 
             sym(paste0('mean_diff_', sub('.*_', '', variable))),
             sym(paste0('relative_diff_mean_perc_', sub('.*_', '', variable))),
             irrigated_relative_diff_mean_perc,
             nonirrigated_relative_diff_mean_perc,
             period))
  complete <- complete %>% 
    pivot_longer(-c(3, 4, 5, 8))
  relative <- complete %>% 
    filter(name %in% c('irrigated_relative_diff_mean_perc', 'nonirrigated_relative_diff_mean_perc')) %>% 
    mutate(name = ifelse(name %in% 'irrigated_relative_diff_mean_perc', paste0('Irrigated'), paste0('Non-irrigated')))
  compl2 <- relative[,c(6, 3, 5)] 
  names(compl2) <- c('yield', 'clim', 'type')
  compl2.2 <- compl2 %>% 
    mutate(type = factor(type))
  y <-  new(variable)
  comp <- ggplot(compl2.2) +
    geom_smooth(aes(y = yield, x = clim, fill = type, group = type),
                col = 'black',
                method = 'gam', 
                formula = y~s(x, k = 3),
                size = .35) +
    scale_fill_manual(values = c('5ab4ac', '#d8b365'), name = '') +
    scale_color_manual(values = c('black', 'black'), name = '') +
    theme_bw(base_size = 11) +
    theme(plot.tag = element_text(size = 11, face = 'bold'),
          legend.title.align = 0.5,
          legend.key.height = unit(0.25, 'cm'),
          legend.position = 'bottom',
          strip.background = element_rect(fill = NA)) +
    geom_rug(sides = 'tr',
             length = unit(0.02, "npc"),
             data = compl2.2, aes(clim, yield),
             col = 'gray') +
    xlab(paste0('Relative change in \n ', y, ' (%)')) +
    ylab('Relative change \n in yield (%)') 
  g1 <- gam(yield~s(clim, k=3), 
            data = filter(compl2.2, grepl('Irrigated', type)), method = "REML")
  irrig_ww <- summary(g1)
  g2 <- gam(yield~s(clim, k=3), 
            data = filter(compl2.2, grepl('Non-irrigated', type)), method = "REML")
  nonirrig_ww <- summary(g2)
  df2 <- data.frame(climate = c(filter(compl2.2, grepl('Irrigated', type))$clim,
                                filter(compl2.2, grepl('Non-irrigated', type))$clim),
                    yield = c(g1$fitted.values, g2$fitted.values),
                    type = c(rep(paste('Irrigated'), length(g1$fitted.values)),
                             rep(paste('Non-irrigated'), length(g2$fitted.values))),
                    variable = variable,
                    se = c(predict(g1, newdata = filter(compl2.2, grepl('Irrigated', type)), type = 'response', se.fit = T)$se.fit, 
                           predict(g2, newdata = filter(compl2.2, grepl('Non-irrigated', type)), type = 'response', se.fit = T)$se.fit),
                    res = c(g1$residuals, g2$residuals))
  
  df <- data.frame(crop = c(paste('Irrigated'),
                            paste('Non-irrigated')),
                   variable = rep(variable, 2),
                   dev = c(irrig_ww$dev.expl*100, nonirrig_ww$dev.expl*100),
                   pval = c(irrig_ww$s.table[1, 4], nonirrig_ww$s.table[1, 4]))
  if(what %in% 'Plot'){
    return(comp)
  }else if(what %in% 'df'){
    return(df)
  }else if (what %in% 'fitted'){
    return(df2)
  }else{
    print('what not supported')
  }
}





b <- plot_gam_crop(variable = 'summer_vpd', win = 4, min_yr = 1980,  what = 'Plot')
cc <- plot_gam_crop(variable = 'summer_tmean', win = 4, min_yr = 1980,  what = 'Plot')
d <- plot_gam_crop(variable = 'summer_prec', win = 4, min_yr = 1980,  what = 'Plot')



county_st <- st_read( '/Users/viig7608/Desktop/CCP/cb_2019_us_county_500k/cb_2019_us_county_500k.shp') 

county_st <- subset(county_st, !(STATEFP %in% c('78', '72', '69', '66', '60', '15', '02'))) %>% 
  mutate(STATEFP = as.numeric(STATEFP))#keep lower 48 + DC

soy <- read.csv('/Users/viig7608/Desktop/CCP/Overview/Soybeans/Data/Processed/soy_1980_clean.csv')#clean NDVI data 
st <- function(state){
  str_to_title(str_split_fixed(state, '_', 2)[1])
}
cou <- function(county){
  str_split_fixed(county, '_', 2)[2]
}
soy <- soy %>% 
  mutate(prop = irrigated_yield/nonirrigated_yield) %>% 
  group_by(Location) %>% 
  summarize(prop = mean(prop, na.rm = T)) %>% 
  ungroup()
soy <- soy %>% 
  mutate(State = sapply(Location, FUN = st),
         County = sapply(Location, FUN = cou))
fips <- read.csv('/Users/viig7608/Desktop/CCP/state-geocodes-v2016.csv')#fips and state names
soy <- left_join(soy, fips)
soy <- soy %>% 
  mutate(County = str_to_sentence(County))

combo <- left_join(county_st, soy, by = c('STATEFP', 'NAME' = 'County'))

#Get counties with soy
ss <- read.csv('/Users/viig7608/Desktop/CCP/Overview/Soybeans/Data/Raw/soy1.csv') %>% 
  group_by(State, County) %>% 
  summarize(Total = sum(Value, na.rm = T)) %>% 
  filter(Total>0) %>% 
  ungroup()
ss2 <- read.csv('/Users/viig7608/Desktop/CCP/Overview/Soybeans/Data/Raw/soy2.csv') %>% 
  group_by(State, County) %>% 
  summarize(Total = sum(Value, na.rm = T)) %>% 
  filter(Total>0) %>% 
  ungroup()
ss <- filter(ss, !State %in% c("MARYLAND", "MICHIGAN", "MINNESOTA", "MISSISSIPPI","MISSOURI"))
ss <- rbind(ss, ss2)

ss <- ss %>% 
  mutate(State = str_to_sentence(State),
         County = str_to_sentence(County))
ss <- left_join(ss, fips)

combo <- left_join(combo, ss, by = c('STATEFP', 'NAME' = 'County'))
us_map <- st_read('/Users/viig7608/Desktop/Fast fires_Virginia/Data/Raw/cb_2018_us_state_500k')#US state boundaries
us_map <- filter(us_map, !NAME %in% c('Puerto Rico', 'Commonwealth of the Northern Mariana Islands', 'Alaska', 'Hawaii', 'American Samoa', 'Guam', 'United States Virgin Islands'))
ref <- st_read('/Users/viig7608/Desktop/Fast fires_Virginia/Data/Raw/event_poly_conus_albers_w_stats_lc_to2020.gpkg')
us_map <- st_transform(us_map, st_crs(ref))
combo <- st_transform(combo, st_crs(ref))
map <- ggplot() + 
  geom_sf(aes(fill = obs), color = NA, fill = 'gray75', 
          data = filter(combo, !is.na(Total))) +
  geom_sf(aes(fill = prop), color = NA, data = combo) +
  viridis::scale_fill_viridis(option = 'A', na.value = NA, 'Irrigated yield/ \n Non-irrigated yield') +
  geom_sf(data = us_map, fill = NA) +
  ggthemes::theme_map(base_size = 11) +
  theme(plot.tag = element_text(size = 11, face = 'bold'),
        legend.title.align = 0.5,
        legend.title = element_text(size = 9),
        legend.key.height = unit(0.25, 'cm'),
        legend.position = 'bottom',
        legend.justification = 'center',
        strip.background = element_rect(fill = NA)) +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

fig <- ggarrange(ggarrange(map, b, 
                 ncol = 2,
                 nrow = 1,
                 common.legend = F, labels = c('a', 'b')),
                 ggarrange(cc, d,
                           ncol = 2,
                           nrow = 1,
                           common.legend = T, labels = c('c', 'd'), legend = 'bottom'),
                 ncol = 1, nrow = 2)

fig

ggsave('/Users/viig7608/Desktop/CCP/Overview/Soybeans/Figures/fig_soy.jpeg',
       plot = fig, width = 6.85, height = 6.85*5.2/5.3, dpi = 750)

ggsave('/Users/viig7608/Desktop/CCP/Overview/Soybeans/Figures/fig_soy.pdf',
       plot = fig, width = 6.85, height = 6.85*5.2/5.3, dpi = 750)


b <- foreach(i = 1:length(variable))%do%
  plot_gam_crop(crop = crop, variable = variable[[i]], win = 4, Absolute = 'yes', min_yr = 1980, type = type, what = 'df')

b <- bind_rows(b)
write_csv(b, paste0('/Users/viig7608/Desktop/CCP/Agric/Processed/',crop,'_relative_gams_dev_expl_1980.csv'))
write_csv(b, paste0('/Users/viig7608/Desktop/CCP/Agric/Processed/',crop, '_absolute_gams_dev_expl_1980.csv'))

d <- foreach(i = 1:length(variable))%do%
  plot_gam_crop(crop = crop, variable = variable[[i]], win = 4, Absolute = 'no', min_yr = 1980, type = type, what = 'fitted')
d <- bind_rows(d)
write_csv(d, paste0('/Users/viig7608/Desktop/CCP/Agric/Processed/',crop,'_relative_gams_fitted_vals_mean_1980.csv'))



