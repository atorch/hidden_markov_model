library(ngspatial)  # For adjacency.matrix
library(plot.matrix)

simulate_ising <- function(n_pixels, adjacency, beta, n_iter=1000) {

    stopifnot(dim(adjacency) == c(n_pixels, n_pixels))

    values <- c(-1, 1)
    z <- sample(values, size=n_pixels, replace=TRUE)

    ## Following http://statweb.stanford.edu/~jtaylo/courses/stats352/notes/ising.pdf
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
