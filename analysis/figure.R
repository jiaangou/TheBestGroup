library('readr')   # for importing csv and rds files
library('dplyr')   # for data wrangling, e.g., %>%, mutate()
library('ggplot2') # for fancy plots

# change ggplot theme
theme_set(theme_bw() +
            theme(legend.position = 'top', # place the legend above plots
                  strip.background = element_blank(), # no box for the facets
                  strip.placement = 'outside', # place the facet text outside
                  panel.grid = element_blank()))

preds <- read_rds('analysis/predictions.rds') # from our data
preds_g <- read_rds('analysis/predictions-gushulak.rds') # from Gushulak et al.

isotopes_long <-
  read_csv('data/simulated-data.csv',
           col_types = cols(.default = 'd', replicate = 'f', basin = 'f')) %>%
  mutate(perc_c = frac_c * 100, # change data to percentages
         perc_n = frac_n * 100, # change data to percentages
         c_n_ratio = frac_c / frac_n) %>% # create C:N ratio
  select(- frac_c, - frac_n) %>% # remove fractions
  tidyr::pivot_longer(-c('depth_m', 'replicate', 'basin'),
                      names_to = 'parameter') %>%
  mutate(label = case_when(parameter == 'd15n' ~ 'delta^{15}~N~\'\211\'',
                           parameter == 'perc_n' ~ 'N~\'\045\'~dry~mass',
                           parameter == 'd13c' ~ 'delta^{13}~C~\'\211\'',
                           parameter == 'perc_c' ~ 'C~\'\045\'~dry~mass',
                           parameter == 'c_n_ratio' ~ 'C:N~ratio') %>%
           # convert to factor to sort variables non-alphabetically
           factor(levels = c('delta^{13}~C~\'\211\'',
                             'delta^{15}~N~\'\211\'',
                             'C~\'\045\'~dry~mass',
                             'N~\'\045\'~dry~mass',
                             'C:N~ratio')))

# create depth profile figure of data with estimated mean and 95% CIs ----
ggplot() +
  facet_grid(basin ~ label, scales = 'free_x', labeller = label_parsed,
             switch = 'x') +
  geom_vline(xintercept = 0, color = 'dodgerblue', lwd = 1.5, alpha = 0.3) +
  geom_ribbon(aes(depth_m, ymin = lwr, ymax = upr), preds_g, alpha = 0.3,
              fill = 'darkorange') +
  geom_ribbon(aes(depth_m, ymin = lwr, ymax = upr), preds, alpha = 0.3) +
  geom_point(aes(depth_m, value), isotopes_long, alpha = 0.5) +
  geom_line(aes(depth_m, mu), preds_g, alpha = 0.5, color = 'darkorange') +
  geom_line(aes(depth_m, mu), preds, alpha = 0.5) +
  scale_x_reverse(lim = c(NA, 0)) + # have zero at the top
  labs(x = 'Water depth (m)', y = NULL) +
  coord_flip()

# create figure of data with estimated mean and 95% CIs ----
ggplot() +
  facet_grid(label ~ basin, scales = 'free_y', labeller = label_parsed,
             switch = 'y') +
  geom_ribbon(aes(depth_m, ymin = lwr, ymax = upr), preds_g, alpha = 0.3,
              fill = 'darkorange') +
  geom_ribbon(aes(depth_m, ymin = lwr, ymax = upr), preds, alpha = 0.3) +
  geom_point(aes(depth_m, value), isotopes_long, alpha = 0.5) +
  geom_point(aes(x, y), color = 'transparent',
             tibble(x = 5, y = 0,
                    label = factor(c('C:N~ratio', 'C~\'\045\'~dry~mass'),
                                   levels(isotopes_long$label)))) +
  geom_line(aes(depth_m, mu), preds_g, alpha = 0.5, color = 'darkorange') +
  geom_line(aes(depth_m, mu), preds, alpha = 0.5) +
  labs(x = 'Water depth (m)', y = NULL)

# save figure
ggsave('figures/isotope-models.png', scale = 1.5, width = 8, height = 4)
