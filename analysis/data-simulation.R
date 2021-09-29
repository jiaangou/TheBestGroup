library('dplyr')   # for data wrangling; e.g., %>%, mutate(), select()
library('mgcv')    # for fitting GAMs
library('tidyr')   # for data wrangling; e.g., pivot_*()
library('gratia')  # for simulation from a fitted GAM
library('ggplot2') # for fancy plots
theme_set(theme_bw() +
            theme(legend.position = 'top', # place the legend above plots
                  strip.background = element_blank(), # no box for the facets
                  strip.placement = 'outside')) # place the facet text outside

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

# NOTE: Gushulak et al. use incorrect distributions (families) for some models

# d15N and d13C can be positive or negative, so assuming normality is ok
m_n15 <- gam(d15n ~ s(depth_m), # predictor of depth
             family = gaussian(link = 'identity'),
             data = isotopes,
             method = 'REML') # optimiziation method for the smoothness parameter

m_c13 <- gam(d13c ~ s(depth_m),
             family = gaussian(link = 'identity'),
             data = isotopes,
             method = 'REML')

# not necessary for simulating new data ----------------------------------------
#' C:N is strictly > 0, so it violates assumptions of normality.
#' A Gamma  distribution is more appropriate. The `log` link function is
#' necessary for satisfying the assumptions of normality of the residuals. 
#' Note: `log(0) = -Inf` and `log(Inf) = Inf` => `log([0, Inf)) = (-Inf, Inf)`
m_c_n <- gam(c_n_ratio ~ s(depth_m),
             family = Gamma(link = 'log'),
             data = isotopes,
             method = 'REML')
# ------------------------------------------------------------------------------

#' There is no family in `mgcv` for values that are strictly between 0 and 100,
#' so we need to convert %C and %N to the fraction of C and N in the sample.
#' We can then use the beta distribution for values in the interval [0, 1].
#' Otherwise, we violate assumptions of normality as above. We use the `logit`
#' link function to go from [0, 1] to [-Inf, Inf].
#' If we call the proportions `p` with `0 <= p <= 1`, we have:
#' `logit([0, 1)) = log(0/(0 - 1), 1/(1 - 1)) = log(0, Inf) = (-Inf, Inf)`
m_frac_c <- gam(frac_c ~ s(depth_m),
                family = betar(link = 'logit'),
                data = isotopes,
                method = 'REML')

m_frac_n <- gam(frac_n ~ s(depth_m),
                family = betar(link = 'logit'),
                data = isotopes,
                method = 'REML')

# create simulated dataset ----
set.seed(4) # to avoid randomness and have consistent results

# new dataset
new_data <- tibble(depth_m = seq(1, 16, by = 0.25))

# create function for Brownian motion
noise <- function(x) {
  # vector of normal noise with SD a 5th of the vector's SD
  e <- rnorm(n = length(x), mean = 0, sd = sd(x) / 5)
  
  # sum the values sequentially, then reverse the order
  cumsum(e) %>% rev()
}

# function for quick predictions
sim_data <- function(model) {
  sim <-
    simulate(object = model,
             nsim = 1, # only return a single simulation
             seed = 3, # to have consistent results
             newdata = new_data) %>% # data for simulations
  as.numeric() # convert to a vector
  
  sim + noise(sim)
}

new_data <- mutate(new_data,
                   d15n = sim_data(m_n15),
                   d13c = sim_data(m_c13),
                   frac_n = sim_data(m_frac_n),
                   frac_c = sim_data(m_frac_c),
                   perc_n = frac_n * 100,
                   perc_c = frac_c * 100,
                   c_n_ratio = perc_c / perc_n)

# check that the simulated data is similar to the original ----
# pivot the two datasets to long format and bind them together
datasets <-
  bind_rows(pivot_longer(new_data, -'depth_m', names_to = 'parameter') %>%
              mutate(dataset = 'Simulated~data'),
            pivot_longer(isotopes, -'depth_m', names_to = 'parameter') %>%
              mutate(dataset = 'Original~data')) %>%
  filter(! (parameter == 'frac_c' | parameter == 'frac_n')) %>%
  # add a column of better labels for ggplot
  mutate(label = case_when(parameter == 'c.n.ratio' ~ 'C:N~ratio',
                           parameter == 'd15n' ~ 'delta^{15}~N~\'\211\'',
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
  scale_x_reverse() +
  labs(x = 'Water depth (m)', y = NULL)

# smoothed values only
smoothed <-
  mutate(datasets, dataset = gsub('~', ' ', dataset)) %>% # change '~' to ' '
  ggplot(aes(depth_m, value, color = dataset, fill = dataset)) +
  facet_grid(label ~ ., scales = 'free_y', labeller = label_parsed,
             switch = 'y') +
  geom_smooth(method = 'gam', formula = y ~ s(x), se = TRUE) +
  scale_x_reverse() +
  scale_color_brewer(NULL, type = 'qual', palette = 6) +
  scale_fill_brewer(NULL, type = 'qual', palette = 6) +
  labs(x = 'Water depth (m)', y = NULL); smoothed

# save a plot of the smooths
ggsave(filename = 'figures/smoothed-datasets.png', plot = smoothed, width = 4,
       height = 4, scale = 1.5)

# write the new dataset as a csv ----
# write.csv(new_data, file = 'data/simulated-data.csv')
