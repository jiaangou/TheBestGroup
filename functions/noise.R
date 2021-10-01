# create function for Brownian motion
noise <- function(x) {
  x %>%
    # add a vector of normal noise with SD a fraction of the vector's SD
    mutate(e = rnorm(n = nrow(x), mean = 0, sd = sd(sim) / 50)) %>%
    # sum the values sequentially in each replicate set, then reverse the order
    group_by(replicate) %>%
    mutate(e = cumsum(e) %>% rev(),
           sim = sim + e) %>%
    ungroup() %>%
    select(-e)
}
