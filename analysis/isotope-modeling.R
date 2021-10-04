library('readr')   # for read_csv() and cols()
library('ggplot2') # for fancy plots
library('mgcv')    # for modeling
library('dplyr')   # for data wrangling, e.g., %>%, mutate(), select()
library('tidyr')   # for data wrangling, e.g., pivot_*()
library('gratia')  # for model diagnostics

# change ggplot theme
theme_set(theme_bw() +
            theme(legend.position = 'top', # place the legend above plots
                  strip.background = element_blank(), # no box for the facets
                  strip.placement = 'outside', # place the facet text outside
                  panel.grid = element_blank()))

# import data and add necessary columns ----
isotopes <-
  read_csv('data/simulated-data.csv',
           col_types = cols(.default = 'd', replicate = 'f', basin = 'f')) %>%
  mutate(c_n_ratio = frac_c / frac_n) # create C:N ratio

# convert to "long" format for easier plotting
isotopes_long <-
  isotopes %>%
  select(-c(replicate, basin)) %>%
  pivot_longer(-depth_m, names_to = 'parameter') %>%
  # fancy labels for ggplots
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

# plot the data
ggplot(isotopes_long, aes(value, depth_m)) +
  facet_grid(. ~ label, scales = 'free_x', labeller = label_parsed,
             switch = 'x') +
  geom_hline(yintercept = 0, color = 'dodgerblue') +
  geom_point(alpha = 0.5) +
  scale_y_reverse(lim = c(NA, 0)) + # have zero at the top
  labs(x = NULL, y = 'Water depth (m)')

# modeling ----
# NOTE: Gushulak et al. use incorrect distributions (families) for some models

# d15N and d13C can be positive or negative, so assuming normality is ok
m_d15n <- gam(d15n ~
                # average effect of depth across all basins and replicates
                s(depth_m, k = 15) +
                # effect of basin (keep complexity < that of global smooth)
                s(depth_m, basin, k = 10, bs = 'fs') +
                # effect of replicates (keep complexity < that of basin)
                s(depth_m, replicate, k = 5, bs = 'fs'),
              family = gaussian(link = 'identity'),
              data = isotopes,
              method = 'REML') # optimization method for the smoothness param

m_d13c <- gam(d13c ~
                # average effect of depth across all basins and replicates
                s(depth_m, k = 15) +
                # effect of basin
                s(depth_m, basin, k = 10, bs = 'fs') +
                # effect of replicates
                s(depth_m, replicate, k = 5, bs = 'fs'),
              family = gaussian(link = 'identity'),
              data = isotopes,
              method = 'REML')

#' C:N is strictly > 0, so it violates assumptions of normality.
#' A Gamma  distribution is more appropriate. The `log` link function is
#' necessary for satisfying the assumptions of normality of the residuals. 
#' Note: `log(0) = -Inf` and `log(Inf) = Inf`, so `log([0, Inf)) = (-Inf, Inf)`
m_c_n <- gam(c_n_ratio ~
               # average effect of depth across all basins and replicates
               s(depth_m, k = 15) +
               # effect of basin
               s(depth_m, basin, k = 10, bs = 'fs') +
               # effect of replicates
               s(depth_m, replicate, k = 5, bs = 'fs'),
             family = Gamma(link = 'log'),
             data = isotopes,
             method = 'REML')

#' There is no family in `mgcv` for values that are strictly between 0 and 100,
#' so we need to convert %C and %N to the fraction of C and N in the sample.
#' We can then use the beta distribution for values in the interval [0, 1].
#' Otherwise, we violate assumptions of normality as above. We use the `logit`
#' link function to go from [0, 1] to [-Inf, Inf].
#' If we call the proportions `p` with `0 <= p <= 1`, we have:
#' `logit([0, 1)) = log(odds(0), odds(1)) = log(0, Inf) = (-Inf, Inf)`
m_frac_c <- gam(frac_c ~
                  # average effect of depth across all basins and replicates
                  s(depth_m, k = 15) +
                  # effect of basin
                  s(depth_m, basin, k = 10, bs = 'fs') +
                  # effect of replicates
                  s(depth_m, replicate, k = 5, bs = 'fs'),
                family = betar(link = 'logit'),
                data = isotopes,
                method = 'REML')

m_frac_n <- gam(frac_n ~
                  # average effect of depth across all basins and replicates
                  s(depth_m, k = 15) +
                  # effect of basin
                  s(depth_m, basin, k = 10, bs = 'fs') +
                  # effect of replicates
                  s(depth_m, replicate, k = 5, bs = 'fs'),
                family = betar(link = 'logit'),
                data = isotopes,
                method = 'REML')

# model diagnostics ----
B <- round(log2(nrow(isotopes))) + 1 # number of histogram bins

appraise(m_d15n, n_bins = B) # no issues
appraise(m_d13c, n_bins = B) # no issues
appraise(m_c_n, n_bins = B) # no issues
appraise(m_frac_c, n_bins = B) # no issues, but low data density at low values
appraise(m_frac_n, n_bins = B) # no issues, but lower data density at low values

draw(m_d15n, residuals = TRUE, scales = 'fixed')
draw(m_d13c, residuals = TRUE, scales = 'fixed')
draw(m_c_n, residuals = TRUE, scales = 'fixed')
draw(m_frac_c, residuals = TRUE, scales = 'fixed')
draw(m_frac_n, residuals = TRUE, scales = 'fixed')

# make predictions from the fitted models ----
new_data <- expand_grid(depth_m = seq(1, 16, length.out = 400),
                        basin = unique(isotopes$basin),
                        replicate = 'ignore')

make_preds <- function(model) {
  link_inv <- model$family$linkinv # get inverse link function
  param <- names(model$model[1]) # name of the response variable
  
  # change "frac" to "perc"
  param <- if_else(grepl('frac', param), gsub('frac', 'perc', param), param)
  
  bind_cols(new_data,
            predict.gam(model, newdata = new_data, type = 'link',
                        se.fit = TRUE, exclude = c('s(depth_m,replicate)'))) %>%
    # convert fractions to percentages, and keep the rest the same
    mutate(parameter = param,
           mu = link_inv(fit), # move mean to the response scale
           lwr = link_inv(fit - 1.96 * se.fit), # create 95% CIs, then move back
           upr = link_inv(fit + 1.96 * se.fit), # to response scale
           # change fractions back to percentages
           lwr = lwr * if_else(grepl('perc', param), true = 100, false = 1),
           mu = mu * if_else(grepl('perc', param), true = 100, false = 1),
           upr = upr * if_else(grepl('perc', param), true = 100, false = 1),
           # fancy labels for ggplots
           label = case_when(parameter == 'd15n' ~ 'delta^{15}~N~\'\211\'',
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
}

preds <- bind_rows(make_preds(m_d15n),
                   make_preds(m_d13c),
                   make_preds(m_c_n),
                   make_preds(m_frac_c),
                   make_preds(m_frac_n))

ggplot() +
  facet_grid(label ~ basin, scales = 'free_y', labeller = label_parsed,
             switch = 'y') +
  geom_ribbon(aes(depth_m, ymin = lwr, ymax = upr), preds, alpha = 0.3) +
  geom_line(aes(depth_m, mu), preds, alpha = 0.5)
  
saveRDS(preds, 'analysis/predictions.rds')
