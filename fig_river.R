#This script estimates d-s curves for streaflw and makes figure 6

library(tidyverse)
library(raster)
library(sf)
library(foreach)
library(doParallel)
library(zoo)
library(patchwork)
library(tidymv)
library(mgcv)
library(mgcViz)
library(ggthemes)
library(ggpubr)
ref <- st_read('/Users/viig7608/Desktop/Fast fires_Virginia/Data/Raw/event_poly_conus_albers_w_stats_lc_to2020.gpkg')

us_map <- st_read('/Users/viig7608/Desktop/Fast fires_Virginia/Data/Raw/cb_2018_us_state_500k') %>% 
  st_transform(crs(ref))#US state boundaries
us_map <- filter(us_map, !STUSPS%in% c('AK', 'HI', 'PR', 'GU', 'AS', 'VI', 'MP'))
basin <- st_read('/Users/viig7608/Desktop/CCP/Overview/Colorado river/Data/Raw/Colorado_River_Basin_Hydrological_Boundaries_with_Areas_served_by_Colorado_River 2') %>% 
  st_transform(crs(ref))
river <- st_read('/Users/viig7608/Desktop/CCP/Overview/Colorado river/Data/Raw/Colorado_River_Basin_Rivers/2ecafda8-df4e-48ea-beac-fd22238c26592020313-1-19nxhp4.a488.shp') %>% 
  st_transform(crs(ref))
mp <- ggplot() +
  geom_sf(data = us_map, fill = NA, color = 'black') +
  geom_sf(data = filter(basin, HU_6_Name %in% 'Lower Colorado'), fill = 'red', alpha = .5, color = 'black') +
  geom_sf(data = filter(basin, HU_6_Name %in% 'Upper Colorado'), fill = 'blue', alpha = .5, color = 'black') +
  # geom_sf(data = river, fill = NA, color = 'blue', linewidth = 0.3) +
   theme_map() 

#Discharge
files <- list.files('/Users/viig7608/Desktop/CCP/Overview/Colorado river/Data/Raw/river',
                    full.names = T)

wy <- function(month, year){
  wy <- ifelse(month<10, year, (year + 1))
  return(wy)#Calculates water year (Oct 1 - Sep 30)
}
read_flow <- function(file){
  flow <- read.table(file, comment.char = '#', header = T, fill = T)
  flow <- flow[-1,]
  flow <- flow %>% 
    mutate(year = as.numeric(substr(datetime, 1, 4)),
           month = as.numeric(substr(datetime, 6, 7)),
           water_year = wy(month, year))
  names(flow)[4] <- 'flow'
  flow_mean <- flow %>% 
    mutate(flow = as.numeric(flow)) %>% 
    group_by(water_year) %>% 
    summarise(discharge_mean = mean(flow, na.rm = T),
              discharge_var = var(flow, na.rm = T))
  return(flow_mean)
}

#Only work with upper 
flow_upper <- read_flow(files[1])
flow_upper <- filter(flow_upper, water_year>1921 & water_year<2020)
flow_lower <- read_flow(files[2])
flow_lower <- filter(flow_lower, water_year>1921 & water_year<1970)
flow_lower2 <- read_flow(files[3])
flow_lower2 <- filter(flow_lower2, water_year>1921 & water_year<2020)
upper_dams <- c(1951, 1970, 1950, 1946)
lower_dams <- c(1950, 1938, 1958, 1941, 1998)

flow_g <- ggplot() +
  geom_line(data = flow_upper,
            aes(water_year, discharge_mean), color = 'blue') +
  geom_line(data = flow_lower,
            aes(water_year, discharge_mean), color = 'red') +
  geom_line(data = flow_lower2,
            aes(water_year, discharge_mean), color = 'red', linetype = 'dashed') +
  geom_point(aes(upper_dams, rep(31600, length(upper_dams))), color = 'blue', shape = '+', size = 4) +
  geom_point(aes(lower_dams, rep(31600, length(lower_dams))), color = 'red', shape = '+', size = 4) +
  geom_rect(aes(ymin = 30000, ymax = 30600, xmin = 1931, xmax = 1936), fill = 'blue') +
  geom_rect(aes(ymin = 30000, ymax = 30600, xmin = 1956, xmax = 1966), fill = 'red') +
  geom_text(aes(x = 1933.5, y = 28700, label = 'Hoover Dam'), size = 2.5) +
  geom_text(aes(x = 1961, y = 28700, label = 'Glen Canyon Dam'), size = 2.5) +
  theme_bw(base_size = 11) +
  theme(plot.tag = element_text(size = 11, face = 'bold'))  +
  ylab('Streamflow (cfs)') +
  xlab('Water year') +
  labs(tag = 'a')

A <- flow_g +
  inset_element(mp, 0.75, 0.6, 1, .95, align_to = 'plot')
library(ggpubr)

#GAMs
flow <- read_flow(files[1])
prec <-  read.csv('/Users/viig7608/Desktop/CCP/Overview/Colorado river/Data/Processed/precip_annual.csv')[,-1]

combo <- inner_join(prec, flow)
combo <- combo %>% 
  mutate(period = as.factor(ifelse(water_year<1964, 'pre', 'post')))

plot_smooths2 <- function (model, series, comparison = NULL, facet_terms = NULL, 
                           conditions = NULL, exclude_random = TRUE, exclude_terms = NULL, 
                           series_length = 25, split = NULL, sep = "\\.", transform = NULL, 
                           ci_z = 1.96, time_series) 
{
  if (!missing(time_series)) {
    warning("The time_series argument has been deprecated and will be removed in the future. Please use `series` instead.")
    series_q = dplyr::enquo(time_series)
  }
  else {
    time_series = NULL
    series_q <- dplyr::enquo(series)
  }
  comparison_q <- dplyr::enquo(comparison)
  facet_terms_q <- dplyr::enquo(facet_terms)
  if (rlang::quo_is_null(comparison_q)) {
    comparison_q <- NULL
  }
  if (rlang::quo_is_null(facet_terms_q)) {
    facet_terms_q <- NULL
  }
  outcome_q <- model$formula[[2]]
  predicted_tbl <- get_gam_predictions(model, {
    {
      series
    }
  }, conditions, exclude_random = exclude_random, exclude_terms = exclude_terms, 
  series_length = series_length, split = split, sep = sep, 
  transform = transform, ci_z = ci_z, .comparison = {
    {
      comparison
    }
  })
  smooths_plot <- predicted_tbl %>% ggplot2::ggplot(ggplot2::aes_string(rlang::quo_name(series_q), 
                                                                        rlang::sym(rlang::quo_name(outcome_q)))) + {
                                                                          if (!is.null(comparison_q)) {
                                                                            ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "CI_lower", 
                                                                                                                     ymax = "CI_upper", fill = rlang::quo_name(comparison_q)), 
                                                                                                 alpha = 0.2)
                                                                          }
                                                                        } + {
                                                                          if (is.null(comparison_q)) {
                                                                            ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "CI_lower", 
                                                                                                                     ymax = "CI_upper"), alpha = 0.2)
                                                                          }
                                                                        } + {
                                                                          if (!is.null(comparison_q)) {
                                                                            ggplot2::geom_path(ggplot2::aes_string(colour = rlang::quo_name(comparison_q)))
                                                                          }
                                                                        } + {
                                                                          if (is.null(comparison_q)) {
                                                                            ggplot2::geom_path(ggplot2::aes_string())
                                                                          }
                                                                        } + {
                                                                          if (!is.null(facet_terms_q)) {
                                                                            ggplot2::facet_wrap(facet_terms_q)
                                                                          }
                                                                        }
  return(smooths_plot)
}


combo$period <- recode_factor(combo$period, pre = 'Before Glen Canyon Dam', post = 'After Glen Canyon Dam')
combo$period <- factor(combo$period, levels = c('Before Glen Canyon Dam', 'After Glen Canyon Dam'))  
plot_gam <- function(climate, meanvar, tag){
  if(meanvar %in% 'Mean' & climate %in% 'Mean'){
    g1 <- gam(discharge_mean~s(prec_annual,by = period), data = combo, method = "REML")
  }else if(meanvar %in% 'Var' & climate %in% 'Mean'){
    g1 <- gam(discharge_var~s(prec_annual, k=3, by = period), data = combo, method = "REML")
  }else if(meanvar %in% 'Mean' & climate %in% 'Var'){
    g1 <- gam(discharge_mean~s(prec_sd, k=3, by = period), data = combo, method = "REML")
  }else{
    g1 <- gam(discharge_var~s(prec_sd, k=3, by = period), data = combo, method = "REML")
  }
  x <- ifelse(climate %in% 'Var', 'Precipitation SD (mm)', paste0(climate, ' precipitation (mm)'))
  q1 <- if(climate %in% 'Mean'){
    plot_smooths2(model = g1,
                 series = prec_annual,
                 comparison = period) +
      scale_fill_manual(values = c('5ab4ac', '#d8b365'), name = '') +
      scale_color_manual(values = c('black', 'black'), name = '') +
      theme_bw(base_size = 11) +
      theme(legend.key.height = unit(0.25, 'cm'),
            legend.position = 'bottom',
            strip.background = element_rect(fill = NA),
            plot.tag = element_text(size = 11, face = 'bold')) +     
      geom_rug(sides = 'br',
               length = unit(0.02, "npc"),
               data = combo,
               aes(prec_annual),
               col = 'gray') +
      xlab(x) +
      ylab(paste0(meanvar, ' streamflow (csf)')) +
      labs(tag = tag)
  }else{
    plot_smooths2(model = g1,
                 series = prec_sd,
                 comparison = period) +
      scale_fill_manual(values = c('5ab4ac', '#d8b365'), name = '') +
      scale_color_manual(values = c('black', 'black'), name = '') +
      theme_bw(base_size = 11) +
      theme(legend.key.height = unit(0.25, 'cm'),
            legend.position = 'bottom',
            strip.background = element_rect(fill = NA),
            plot.tag = element_text(size = 11, face = 'bold')) +
      geom_rug(sides = 'br',
               length = unit(0.02, "npc"),
               data = combo,
               aes(prec_sd),
               col = 'gray') +
      xlab(x) +
      ylab(paste0(meanvar, ' streamflow (csf)')) +
      labs(tag = tag)
  }
}
B <- plot_gam(climate = 'Mean', meanvar = 'Mean', tag = 'b')
C <- plot_gam(climate = 'Var', meanvar = 'Mean', tag = 'c')

fig <- ggarrange(A, ggarrange(B, C, ncol = 2, common.legend = T, legend = 'bottom'), nrow = 2)
fig
ggsave(paste0('/Users/viig7608/Desktop/CCP/Overview/Colorado river/Figures/fig_river.jpeg'), 
       width = 6.85, height = 6.5, plot = fig, dpi = 750)
ggsave(paste0('/Users/viig7608/Desktop/CCP/Overview/Colorado river/Figures/fig_river.pdf'), 
       width = 6.85, height = 6.5, plot = fig, dpi = 750)
