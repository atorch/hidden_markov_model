library(ngspatial)  # For adjacency.matrix
library(plot.matrix)

simulate_ising <- function(n_pixels=100, beta=1, n_iter=1000) {
    ## Following http://statweb.stanford.edu/~jtaylo/courses/stats352/notes/ising.pdf
    values <- c(-1, 1)
    n_pixels_per_side <- sqrt(n_pixels)
    adjacency <- adjacency.matrix(m=n_pixels_per_side, n=n_pixels_per_side)  # TODO Is this slow?  Cache it
    z <- sample(values, size=n_pixels, replace=TRUE)
    for(iter in seq_len(n_iter)) {
        for(index in seq_len(n_pixels)) {
            neighbors <- which(adjacency[index, ] > 0)
            neighbor_sum <- sum(z[neighbors])
            odds <- exp(2 * beta * neighbor_sum)
            p <- odds / (1 + odds)
            z[index] <- sample(values, size=1, prob=c(1-p, p))
        }
    }
    return(z)
}

if(FALSE) {
    ## Example
    n_pixels_per_side <- 50
    z <- simulate_ising(n_pixels=n_pixels_per_side^2, beta=0.4, n_iter=100)
    plot(matrix(z, n_pixels_per_side, n_pixels_per_side))
}
