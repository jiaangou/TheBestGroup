library('dplyr')
library('mgcv')
library('ggplot2')

# change ggplot theme
theme_set(theme_bw() +
            theme(legend.position = 'top', # place the legend above plots
                  strip.background = element_blank(), # no box for the facets
                  strip.placement = 'outside', # place the facet text outside
                  panel.grid = element_blank()))

IS <- readxl::read_xlsx('gushulak-et-al/Gall Lake all data.xlsx') %>%
  select(`Depth(m)`:`%C`) %>%
  rename(N = `%N`,
         C = `%C`,
         depth_m = `Depth(m)`) %>%
  mutate(CN = C/N)

# plot the data ----
isotopes_long <- tidyr::pivot_longer(IS, -'depth_m', names_to = 'parameter') %>%
  mutate(label = case_when(parameter == 'CN' ~ 'C:N~ratio',
                           parameter == 'd15NAIR' ~ 'delta^{15}~N~\'\211\'',
                           parameter == 'N' ~ 'N~\'\045\'~dry~mass',
                           parameter == 'd13CVPDB' ~ 'delta^{13}~C~\'\211\'',
                           parameter == 'C' ~ 'C~\'\045\'~dry~mass') %>%
           # convert to factor to sort variables non-alphabetically
           factor(levels = c('delta^{13}~C~\'\211\'',
                             'delta^{15}~N~\'\211\'',
                             'C~\'\045\'~dry~mass',
                             'N~\'\045\'~dry~mass',
                             'C:N~ratio')))

ggplot() +
  facet_grid(label ~ ., scales = 'free_y', labeller = label_parsed,
             switch = 'y') +
  geom_point(aes(depth_m, value), isotopes_long, alpha = 0.5) +
  scale_x_reverse(lim = c(NA, 0)) + # have zero at the top
  labs(x = 'Water depth (m)', y = NULL)

#' the original code uses Date instead of Depth_m, not sure why
#' bs = 'tp' is the default basis, so it is omitted here for simplicity
#' note that the CN model should be Gamma, not Gaussian, and the models for
#' percentage of carbon and nitrogen should be Beta for proportions of C and N
N15 <- gam(d15NAIR ~ s(depth_m), data = IS, method = 'REML')
C13 <- gam(d13CVPDB ~ s(depth_m), data = IS, method = 'REML')
CN <- gam(CN ~ s(depth_m), data = IS, method = 'REML')
C <- gam(C ~ s(depth_m), data = IS, family = Gamma(link='log'), method = 'REML')
N <- gam(N ~ s(depth_m), data = IS, family = Gamma(link='log'), method = 'REML')

# predict from the models ----
new_data <- tibble(depth_m = seq(1, 16, length.out = 400))

make_preds <- function(model) {
  link_inv <- model$family$linkinv # find inverse link function
  param <- names(model$model[1]) # name of the response variable
  
  bind_cols(new_data,
            predict.gam(model, newdata = new_data, type = 'link',
                        se.fit = TRUE, terms = 's(depth_m)')) %>%
    # convert fractions to percentages, and keep the rest the same
    mutate(parameter = param,
           mu = link_inv(fit), # move mean to the response scale
           lwr = link_inv(fit - 1.96 * se.fit), # create 95% CIs, then move back
           upr = link_inv(fit + 1.96 * se.fit), # to response scale
           # fancy labels for ggplots (**fractions are actually percentages!**)
           label = case_when(parameter == 'CN' ~ 'C:N~ratio',
                             parameter == 'd15NAIR' ~ 'delta^{15}~N~\'\211\'',
                             parameter == 'N' ~ 'N~\'\045\'~dry~mass',
                             parameter == 'd13CVPDB' ~ 'delta^{13}~C~\'\211\'',
                             parameter == 'C' ~ 'C~\'\045\'~dry~mass'),
           # convert to factor to sort variables non-alphabetically
           label = factor(label,
                          levels = c('delta^{13}~C~\'\211\'',
                                     'delta^{15}~N~\'\211\'',
                                     'C~\'\045\'~dry~mass',
                                     'N~\'\045\'~dry~mass',
                                     'C:N~ratio')))
}

preds <- bind_rows(make_preds(N15),
                   make_preds(C13),
                   make_preds(CN),
                   make_preds(C),
                   make_preds(N))

ggplot() +
  facet_grid(. ~ label, scales = 'free_x', labeller = label_parsed,
             switch = 'x') +
  geom_vline(xintercept = 0, color = 'dodgerblue', lwd = 1.5, alpha = 0.3) +
  geom_ribbon(aes(depth_m, ymin = lwr, ymax = upr), preds, alpha = 0.3) +
  geom_point(aes(depth_m, value), isotopes_long, alpha = 0.5) +
  geom_line(aes(depth_m, mu), preds, alpha = 0.5) +
  scale_x_reverse(lim = c(NA, 0)) + # have zero at the top
  labs(x = 'Water depth (m)', y = NULL) +
  coord_flip()

saveRDS(preds, 'analysis/predictions-gushulak.rds')
