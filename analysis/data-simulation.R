library('dplyr')   # for data wrangling; e.g., %>%, mutate(), select()
library('mgcv')    # for fitting GAMs
library('tidyr')   # for data wrangling; e.g., pivot_*()
library('assertr') # for ensuring simulated data is sensible
library('gratia')  # for simulation from a fitted GAM
library('ggplot2') # for fancy plots
source('functions/sim_data.R') # for simulating data
# NREPLICATES <- 4 # number of sample replicates
# NBASINS <- 4 # number of basins in Gall Lake

# change ggplot theme
theme_set(theme_bw() +
            theme(legend.position = 'top', # place the legend above plots
                  strip.background = element_blank(), # no box for the facets
                  strip.placement = 'outside', # place the facet text outside
                  panel.grid = element_blank()))

isotopes <-
  readxl::read_xlsx('gushulak-et-al/Gall Lake all data.xlsx') %>%
  select(`Depth(m)`:`%C`) %>% # subset to only isotope data
  rename(depth_m = `Depth(m)`, # change colnames to avoid needing backticks
         d15n = d15NAIR,
         d13c = d13CVPDB,
         perc_n = `%N`,
         perc_c = `%C`) %>%
  filter(!is.na(depth_m),  # remove rows with NAs
         depth_m != 0) %>% # no data at depth_m == 0
  mutate(c_n_ratio = perc_c / perc_n, # create C:N ratio
         frac_n = perc_n / 100, # change to fractions to use Beta distribution,
         frac_c = perc_c / 100) # which is appropriate for values bound [0, 1]
summary(isotopes)

# modeling ----
# NOTE: Gushulak et al. use incorrect distributions (families) for some models

# d15N and d13C can be positive or negative, so assuming normality is ok
m_d15n <- gam(d15n ~ s(depth_m, k = 15), # accounts for effect of depth only
              family = gaussian(link = 'identity'),
              data = isotopes,
              method = 'REML') # optimiziation method for the smoothness parameter

m_d13c <- gam(d13c ~ s(depth_m, k = 15),
              family = gaussian(link = 'identity'),
              data = isotopes,
              method = 'REML')

#' There is no family in `mgcv` for values that are strictly between 0 and 100,
#' so we need to convert %C and %N to the fraction of C and N in the sample.
#' We can then use the beta distribution for values in the interval [0, 1].
#' Otherwise, we violate assumptions of normality as above. We use the `logit`
#' link function to go from [0, 1] to [-Inf, Inf].
#' If we call the proportions `p` with `0 <= p <= 1`, we have:
#' `logit([0, 1)) = log(odds(0), odds(1)) = log(0, Inf) = (-Inf, Inf)`
m_frac_c <- gam(frac_c ~ s(depth_m, k = 15),
                family = betar(link = 'logit'),
                data = isotopes,
                method = 'REML')

m_frac_n <- gam(frac_n ~ s(depth_m, k = 15),
                family = betar(link = 'logit'),
                data = isotopes,
                method = 'REML')

# create simulated dataset ----
set.seed(2) # to avoid randomness and have consistent results

# new dataset
new_depths <-
  expand_grid(replicate = 1:16, # 4 replicates for 4 basins each
              depth_m = 1:16) %>% # a sample every meter
  mutate(basin = ceiling(replicate / 4), # add a column for basin
         basin = case_when(basin == 1 ~ 'North',
                           basin == 2 ~ 'West',
                           basin == 3 ~ 'South',
                           basin == 4 ~ 'East'),
         replicate = paste0(basin, '_', replicate)) # add basin to replicate col

new_data <-
  bind_rows(sim_data('m_d13c'),
            sim_data('m_d15n'),
            sim_data('m_frac_c'),
            sim_data('m_frac_n')) %>%
  pivot_wider(values_from = 'sim', names_from = 'parameter') %>%
  # ensure fractions are between 0 and 1
  chain_start() %>%
  assert(within_bounds(0, 1), c(frac_n, frac_c)) %>%
  chain_end() %>%
  mutate(c_n_ratio = frac_c / frac_n) # bound [0, Inf)

# check that the simulated data is similar to the original ----
# pivot the two datasets to long format and bind them together
datasets <-
  bind_rows(mutate(new_data,
                   perc_c = frac_c * 100, # change fractions to percentages
                   perc_n = frac_n * 100) %>%
              select(-frac_c, -frac_n) %>%
              pivot_longer(-c('depth_m', 'basin', 'replicate'),
                           names_to = 'parameter') %>%
              mutate(dataset = 'Simulated~data'),
            pivot_longer(isotopes,
                         -'depth_m', names_to = 'parameter') %>%
              mutate(dataset = 'Original~data')) %>%
  filter(! (parameter == 'frac_c' | parameter == 'frac_n')) %>%
  # add a column of better labels for ggplot
  mutate(label = case_when(parameter == 'd15n' ~ 'delta^{15}~N~\'\211\'',
                           parameter == 'perc_n' ~ 'N~\'\045\'~dry~mass',
                           parameter == 'd13c' ~ 'delta^{13}~C~\'\211\'',
                           parameter == 'perc_c' ~ 'C~\'\045\'~dry~mass',
                           parameter == 'c_n_ratio' ~ 'C:N~ratio'),
         # convert to factor to sort variables non-alphabetically
         label = factor(label,
                        levels = c('delta^{13}~C~\'\211\'',
                                   'delta^{15}~N~\'\211\'',
                                   'C~\'\045\'~dry~mass',
                                   'N~\'\045\'~dry~mass',
                                   'C:N~ratio')))

# check the similarity
ggplot(datasets, aes(depth_m, value)) +
  facet_grid(label ~ dataset, scales = 'free_y', labeller = label_parsed,
             switch = 'y') +
  geom_point(alpha = 0.5) +
  labs(x = 'Water depth (m)', y = NULL)

# save a plot of the smoothed values with simulations
smoothed <-
  mutate(datasets, dataset = gsub('~', ' ', dataset)) %>% # change '~' to ' '
  ggplot(aes(depth_m, value, color = dataset, fill = dataset)) +
  facet_wrap(. ~ label, scales = 'free_y', labeller = label_parsed,
             strip.position = 'left') +
  geom_line(aes(depth_m, value, group = replicate), inherit.aes = FALSE,
            filter(datasets, dataset == 'Simulated~data'), alpha = 0.15) +
  geom_smooth(method = 'gam', formula = y ~ s(x), se = TRUE) +
  scale_color_brewer(NULL, type = 'qual', palette = 6) +
  scale_fill_brewer(NULL, type = 'qual', palette = 6) +
  labs(x = 'Water depth (m)', y = NULL) +
  xlim(c(0, NA)) +
  theme(legend.position = c(0.85, 0.2)); smoothed

ggsave(filename = 'figures/smoothed-datasets.png', plot = smoothed, width = 4,
       height = 2, scale = 2)

# write the new dataset as a csv ----
new_data %>%
  select(-c(c_n_ratio)) %>% # remove column to make more realistic
  relocate(basin, replicate, depth_m) %>%
  readr::write_csv(file = 'data/simulated-data.csv')
