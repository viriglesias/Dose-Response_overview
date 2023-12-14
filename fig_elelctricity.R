#This script estimates d-s curves for electricity and makes figure 4

library(tidyverse)
library(tidymv)
library(mgcv)
library(mgcViz)
library(foreach)
library(ggpubr)
library(viridis)
library(sf)
library(ggthemes)
library(patchwork)
files <- list.files('/Users/viig7608/Desktop/CCP/Overview/Electricity/Data/Raw/',
                    pattern = '_full_',
                    full.names = T)
electr <- foreach::foreach(i = 1:length(files), .combine = rbind)%do%{
  read.csv(files[[i]])[,-1]
}

electr <- electr %>% 
  mutate(sales_adj = sales/customers)
# tmin <-  read.csv('/Users/viig7608/Desktop/CCP/Overview/Electricity/Data/Processed/tmin.csv')[,-1]
# tmin$state.NAME[tmin$state.NAME == 'Louisiana'] <- 'LA'
# tmin$state.NAME[tmin$state.NAME == 'Missouri'] <- 'MO'
# tmin$state.NAME[tmin$state.NAME == 'Minnesota'] <- 'MN'
# names(tmin)[1] <- 'stateid'

# electr <- electr %>% 
#   mutate(sales_adj = sales/customers)
# 
# combo <- inner_join(tmin, electr) %>% 
#   mutate(stateid = as.factor(stateid))
# 
# 
plot_gam <- function(data, climate){
    q1 <- ggplot() +
      geom_smooth(data = combo,
                  method = 'gam',
                  formula = y~s(x, k = 20),
                  aes(y = sales_adj, x = var_mean, color = stateid)) +
      geom_smooth(data = combo,
                  method = 'gam',
                  formula = y~s(x, k = 20),
                  aes(y = sales_adj, x = var_mean),
                  linetype = 'dotted') +
      geom_rug(sides = 'br',
               length = unit(0.02, "npc"),
               data = combo,
               aes(var_mean)) +
      theme_bw(base_size = 11) +
      xlab(climate) +
      ylab('Sales per customer')
}

# a <- plot_gam(climate = 'Tmin', data = 'combo')

tmean <-  read.csv('/Users/viig7608/Desktop/CCP/Overview/Electricity/Data/Processed/tmean.csv')[,-1]
# tmean$state.NAME[tmean$state.NAME == 'Louisiana'] <- 'LA'
# tmean$state.NAME[tmean$state.NAME == 'Missouri'] <- 'MO'
# tmean$state.NAME[tmean$state.NAME == 'Minnesota'] <- 'MN'
# names(tmean)[1] <- 'stateid'

us <- read.csv('/Users/viig7608/Desktop/CCP/Overview/Electricity/Data/Processed/us_lat.csv')[,-1]
names(us)[1:2] <- c('stateid', 'state.NAME')
combo <- inner_join(tmean, us) 
combo <- inner_join(combo, electr)
combo <- combo %>% 
  mutate(stateid = as.factor(stateid))
combo <- filter(combo, !is.na(sales_adj))
plot_gam(climate = 'tmean', data = 'combo')

# tmax <-  read.csv('/Users/viig7608/Desktop/CCP/Overview/Electricity/Data/Processed/tmax.csv')[,-1]
# tmax$state.NAME[tmax$state.NAME == 'Louisiana'] <- 'LA'
# tmax$state.NAME[tmax$state.NAME == 'Missouri'] <- 'MO'
# tmax$state.NAME[tmax$state.NAME == 'Minnesota'] <- 'MN'
# names(tmax)[1] <- 'stateid'

# combo <- inner_join(tmax, electr) %>% 
#   mutate(stateid = as.factor(stateid))
# 
# b <- plot_gam(climate = 'Tmax', data = 'combo')
# fig <- ggarrange(a, d, b, 
#                  ncol = 3,
#                  nrow = 1,
#                  common.legend = T,
#                  legend = 'bottom', labels = 'AUTO')
# fig

# ggsave(paste0('/Users/viig7608/Desktop/CCP/Overview/Electricity/Figures/electricity_no_epoch_gam.jpeg'), 
#        plot = fig, width = 11.1, height = 7.35)
# rm(list=ls())

combo <- filter(combo, !stateid %in% c('AK', 'HI'))
cbo <- filter(combo, stateid %in% c('AK', 'HI'))

a <- gam(sales_adj~s(var_mean, k = 10) + s(lat, k = 10) + s(var_mean, lat, k = 10), data = combo)
aa <- gam(sales_adj~s(var_mean, k = 10) + s(lat, k = 10), data = combo)
aaa <- gam(sales_adj~s(var_mean, lat, k = 10), data = combo)
aaaa <- gam(sales_adj~s(var_mean, k = 10),data = combo)
summary(a)
summary(aa)
summary(aaa)
AIC(a, aa, aaa, aaaa)

comp <- ggplot() +
  geom_smooth(data = combo,
              aes(y = sales_adj*1000000, x = var_mean),
              col = 'black',
              method = 'gam', 
              formula = y~s(x, k = 3),
              linewidth = .35) +
  geom_smooth(data = filter(combo, stateid %in% 'TX'),
              aes(y = sales_adj*1000000, x = var_mean),
              fill = 'red',
              col = 'black',
              linetype = 'dashed',
              method = 'gam', 
              formula = y~s(x, k = 3),
              linewidth = .35) +
  geom_smooth(data = filter(combo, stateid %in% 'ND'),
              aes(y = sales_adj*1000000, x = var_mean),
              fill = 'blue',
              col = 'black',
              linetype = 'dashed',
              method = 'gam', 
              formula = y~s(x, k = 3),
              linewidth = .35) +
  theme_bw(base_size = 11) +
  theme(plot.tag = element_text(size = 11, face = 'bold'),
        legend.title.align = 0.5,
        legend.key.height = unit(0.25, 'cm'),
        legend.position = 'bottom',
        strip.background = element_rect(fill = NA)) +
  geom_rug(sides = 'tr',
           length = unit(0.02, "npc"),
           data = combo, aes(var_mean, sales_adj*1000000),
           col = 'gray') +
  xlab(expression('Temperature'~(degree*C))) +
  ylab(expression('Energy sales (kilowatthours/customer)')) 

combo2 <- filter(combo, var_mean>20)
cs <- combo2 %>% 
  group_by(stateid) %>% 
  summarise(ss = n())
combo_split <- split(combo2, as.factor(combo2$stateid))

cff <- purrr::possibly(function(data){
  ml <- lm(sales_adj~var_mean, data = data)
  cf <- unname(ml$coefficients[2])*1000000
  return(cf)
  }, otherwise = NA)

j <- foreach(i = 1:length(combo_split), .combine = rbind) %do%
  cff(combo_split[[i]])

stateid <- names(combo_split)

df <- data.frame(stateid = stateid, coeffs = j)
str(df)

us <- inner_join(us, df)
t_sum <- combo %>% 
  group_by(state.NAME) %>% 
  summarise(Temperature = mean(var_mean))
us <- inner_join(us, t_sum)
f1 <- ggplot(us) +
  geom_point(aes(x = lat, y = coeffs, color = Temperature), size = 2) +
  geom_text(aes(x = lat, y = coeffs, label = stateid), size = 2, check_overlap = T,
            nudge_y = 2.75) +
   theme_bw(base_size = 11) +
  scale_color_viridis(option = 'H', expression("Temperature "  (degree*C)),
                      breaks = c(18, 20, 22, 24, 26, 28)) +
  xlab(expression('Latitude'~(degree*N))) +
  ylab(expression("Regression coefficient (kwh/customer/"*degree*C*')')) +
  theme_bw(base_size = 11) +
  theme(plot.tag = element_text(size = 11, face = 'bold'),
        legend.title.align = 0.5,
        legend.key.height = unit(0.25, 'cm'),
        legend.position = 'bottom',
        strip.background = element_rect(fill = NA)) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))
# 
# ggsave(paste0('/Users/viig7608/Desktop/CCP/Overview/Electricity/Figures/electricity_linear.jpeg'), 
#        width = 11.1, height = 7.35) 

us_map <- st_read('/Users/viig7608/Desktop/Fast fires_Virginia/Data/Raw/cb_2018_us_state_500k')#US state boundaries
us_map <- filter(us_map, !STUSPS%in% c('AK', 'HI'))

us_map <- inner_join(us_map, us, by = c('NAME' = 'state.NAME'))
ref <- st_read('/Users/viig7608/Desktop/Fast fires_Virginia/Data/Raw/event_poly_conus_albers_w_stats_lc_to2020.gpkg')
us_map <- st_transform(us_map, st_crs(ref))
states <- filter(us_map, stateid %in% c('TX', 'ND'))
f2 <- ggplot() +
  geom_sf(data = us_map, aes(fill = Temperature)) +
  geom_sf(data = states, linewidth = 1.25, fill = NA, color = 'black') +
  theme_map(base_size = 11) +
  scale_fill_viridis(option = 'H', 
                     breaks = c(18, 20, 22, 24, 26, 28)) +
  theme(legend.position = 'none') 
  
  
# ggsave(paste0('/Users/viig7608/Desktop/CCP/Overview/Electricity/Figures/electricity_map.jpeg'), 
#        width = 11.1, height = 7.35) 


B <- f1 +
  inset_element(f2, 0.55, 0.7, 1, 1, align_to = 'full')

fig <- ggpubr::ggarrange(comp, B,
                 ncol = 2,
                 nrow = 1,
                 labels = 'auto')
ggsave(paste0('/Users/viig7608/Desktop/CCP/Overview/Electricity/Figures/fig_electricity.jpeg'), 
              width = 6.85, height = 7.35*6.85/11.1, plot = fig, dpi = 750)

ggsave(paste0('/Users/viig7608/Desktop/CCP/Overview/Electricity/Figures/fig_electricity.pdf'), 
       width = 6.85, height = 7.35*6.85/11.1, plot = fig, dpi = 750)
