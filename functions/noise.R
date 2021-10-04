# create function for Brownian motion by basin and by replicate
noise <- function(x) {
  x %>%
    group_by(basin) %>%
    mutate(
      # add normal noise to each basin with SD a fraction of the vector's SD
      e.basin = rnorm(n = n(), mean = 0, sd = sd(sim) / 4),
      # sum the values sequentially in each basin, then reverse the order
      e.basin = cumsum(e.basin) %>% rev(),
      # add the noise to the simulated data
      sim = sim + e.basin) %>%
    ungroup() %>%
    # do the same for each replicate within each basin
    group_by(replicate) %>%
    mutate(e.replicate = rnorm(n = n(),
                               mean = 0, sd = sd(sim) / 10),
           e.replicate = cumsum(e.replicate) %>% rev(),
           sim = sim + e.replicate) %>%
    ungroup() %>%
    select(-e.basin, -e.replicate)
}
