source('functions/noise.R') # for adding noise to simulated data

# function for quick predictions
sim_data <- function(model = m) {
  # bind new data and all simulations
  bind_cols(simulate(object = model,
                     nsim = NREPLICATES, # number of simulations or replicates
                     seed = 1, # to have consistent results
                     newdata = new_depths) %>% # data for simulations
              as_tibble(.name_repair = 'unique'),
            new_depths) %>%
    suppressMessages() %>% # prevent messages on change of column names
    # change to long format
    pivot_longer(-c('depth_m', 'parameter'),
                 names_to = 'replicate', values_to = 'sim') %>%
    # change replicate column from '...1'-'...4' to '1'-'4'
    mutate(replicate = gsub('...', '', replicate) %>% factor()) %>%
    noise()
}
