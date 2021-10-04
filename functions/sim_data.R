source('functions/noise.R') # for adding noise to simulated data

# function for quick predictions
sim_data <- function(model.name) {
  model <- get(model.name)
  
  # bind new data and all simulations
  bind_cols(new_depths, # add data used for simulations
            sim = simulate(object = model, # sims returned on response scale
                           nsim = 1, # replicates accounted for in `new_depths`
                           seed = 1, # to have consistent results
                           newdata = new_depths) %>% # data for simulations
              as.numeric()) %>% # convert from matrix to vector
    suppressMessages() %>% # prevent messages on change of column names
    mutate(replicate = factor(replicate), # factors are needed for models
           basin = factor(basin),
           # remove "m_" from the model name and pass it as the parameter
           parameter = substr(model.name, start = 3, nchar(model.name)),
           # move sims to link function to apply normal noise
           sim = model$family$linkfun(sim)) %>%
    noise() %>%
    # move sims back to response scale
    mutate(sim = model$family$linkinv(sim))
}
