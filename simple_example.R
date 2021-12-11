library(ggplot2)

get_transition_probs <- function(M_S) {
    return(t(M_S) / rowSums(t(M_S)))
}

expected_deforestation_rate_freq <- function(deforestation_rate=0.04, pr_y_diagonals=c(0.9, 0.8)) {

    ## Transition probabilities for S
    ## Convention: hidden state 1 is forest, hidden state 2 is something else (not forest)
    P <- rbind(c(1 - deforestation_rate, deforestation_rate),
               c(0.02, 0.98))
    
    ## Pr[Y | S] (transpose of the way it's defined in the paper)
    pr_y <- rbind(c(pr_y_diagonals[1], 1 - pr_y_diagonals[1]),
                  c(1 - pr_y_diagonals[2], pr_y_diagonals[2]))
    
    ## Initial distribution over S
    mu <- c(0.7, 0.3)
    
    ## Joint distribution of S_t+1, S_t
    M_S <- t(P * matrix(mu, length(mu), length(mu)))
    
    ## Joint distribution of Y_t+1, Y_t
    M_Y <- t(pr_y) %*% M_S %*% pr_y

    stopifnot(isTRUE(all.equal(sum(M_S), 1.0)))
    stopifnot(isTRUE(all.equal(sum(M_Y), 1.0)))
        
    ## Sanity check: this should be true (we correctly recover P from M_S)
    stopifnot(isTRUE(all.equal(P, get_transition_probs(M_S))))

    P_freq <- get_transition_probs(M_Y)
    
    deforestation_rate_freq <- P_freq[1, 2]
    return(deforestation_rate_freq)
}

true_deforestation_rates <- seq(0.0, 0.601, 0.01)
dfs <- list()
for(pr_y_diagonals in list(c(0.99, 0.99), c(0.95, 0.92), c(0.90, 0.85), c(0.85, 0.85))) {
    expected_deforestation_rates_freq <- sapply(true_deforestation_rates, expected_deforestation_rate_freq, pr_y_diagonals=pr_y_diagonals)
    df <- data.frame(true_deforestation_rate=true_deforestation_rates,
                     expected_deforestation_rate_freq=expected_deforestation_rates_freq)
    df$label <- sprintf("Pr[Y | S] diagonals: %s, %s", pr_y_diagonals[1], pr_y_diagonals[2])
    dfs[[length(dfs) + 1]] <- df
}

df <- do.call(rbind, dfs)

p <- (ggplot(df, aes(x=true_deforestation_rate, y=expected_deforestation_rate_freq, color=label)) +
      geom_point() +
      theme_bw() +
      scale_color_manual("", c()) +
      geom_abline(slope=1, lty=2, alpha=0.5))
ggsave(p, filename="example_expected_deforestation_rate_freq_given_true_deforestation_rate.png", width=6, height=4, units="in")
