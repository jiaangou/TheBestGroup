library('dplyr')   # for data wrangling; e.g., %>%, mutate(), select()
library('mgcv')    # for fitting GAMs
library('tidyr')   # for data wrangling; e.g., pivot_*()
library('assertr') # for ensuring simulated data is sensible
library('gratia')  # for simulation from a fitted GAM
library('ggplot2') # for fancy plots
source('functions/sim_data.R') # for simulating data
NREPLICATES <- 4 # number of sample replicates

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

#' NOTE: Gushulak et al. use incorrect distributions (families) for some models,
#'       but this doesn't matter when generating random data, as long as the
#'       values are sensible.

isotopes_long <-
  isotopes %>%
  select(depth_m:perc_c) %>%
  pivot_longer(-'depth_m', names_to = 'parameter') %>%
  mutate(parameter = factor(parameter)) # needs to be a factor for gam()

# separate intercepts and smooths for each parameter
m <- gam(value ~ parameter + s(depth_m, by = parameter, k = 15),
         data = isotopes_long)
draw(m, scales = 'free')

# create simulated dataset ----
set.seed(7) # to avoid randomness and have consistent results

# new dataset
new_depths <- expand_grid(parameter = unique(isotopes_long$parameter),
                          depth_m = 1:16) # a sample every meter

new_data <-
  sim_data() %>%
  pivot_wider(values_from = 'sim', names_from = 'parameter') %>%
  mutate(perc_n = if_else(perc_n < 0, 4, perc_n), # remove the negative value
         c_n_ratio = perc_c / perc_n) %>% # bound [0, Inf)
  ungroup() %>%
  chain_start() %>%
  # ensure fractions are between 0 and 1
  assert(within_bounds(0, 100), c(perc_n, perc_c)) %>%
  # ensure ratio is positive
  assert(within_bounds(0, Inf), c_n_ratio) %>%
  chain_end()

# check that the simulated data is similar to the original ----
# pivot the two datasets to long format and bind them together
datasets <-
  bind_rows(pivot_longer(new_data,
                         -c('depth_m', 'replicate'), names_to = 'parameter') %>%
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
  facet_grid(label ~ ., scales = 'free_y', labeller = label_parsed,
             switch = 'y') +
  geom_smooth(method = 'gam', formula = y ~ s(x), se = TRUE) +
  geom_line(aes(depth_m, value, group = replicate), inherit.aes = FALSE,
            filter(datasets, dataset == 'Simulated~data'), alpha = 0.3) +
  scale_x_reverse() +
  scale_color_brewer(NULL, type = 'qual', palette = 6) +
  scale_fill_brewer(NULL, type = 'qual', palette = 6) +
  labs(x = 'Water depth (m)', y = NULL); smoothed

ggsave(filename = 'figures/smoothed-datasets.png', plot = smoothed, width = 4,
       height = 4, scale = 1.5)

# write the new dataset as a csv ----
new_data %>%
  select(-c(c_n_ratio)) %>% # remove column to make more realistic
  readr::write_csv(file = 'data/simulated-data.csv')
