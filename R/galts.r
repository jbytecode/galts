library(genalg)
library(DEoptim)
library(dfoptim)

ga.lts <- function(formula,
                   data,
                   h = NULL,
                   iters = 2,
                   popsize = 50,
                   lower,
                   upper,
                   csteps = 2,
                   method = "ga",
                   verbose = FALSE) {
    df_data <- data.frame(data)
    x <- model.matrix(formula, df_data)
    y <- model.frame(formula, df_data)[, 1]
    n <- length(y)
    p <- dim(x)[2]
    if (is.null(h)) {
        h <- floor(n / 2) + floor((p + 1) / 2)
    }
    if (is.null(method) ||
        (method != "de" && method != "ga") && method != "hj") {
        cat("Please select a method:
        'de' for differential evolution
        or 'ga' for genetic algorithm\n")
        return(NULL)
    }

    cstep <- function(candidates, csteps) {
        cmybetas <- candidates
        indices <- order(abs(y - x %*% cmybetas))[1:p]
        for (i in 1:csteps) {
            # ols <- lm(y[indices] ~ x[indices, ] - 1)
            ols <- lm(formula, df_data[indices, ])
            mybetas <- ols$coefficients
            # res <- y - x %*% mybetas
            # res2 <- abs(res)
            # o <- order(res2)
            o <- order(abs(y - x %*% mybetas))
            indices <- sort(o[1:h])
        }
        return(mybetas)
    }

    cost <- function(candidates) {
        newbetas <- cstep(candidates, csteps)
        res <- y - x %*% newbetas
        fitn <- sum(sort(res^2)[1:h])
        return(fitn)
    }

    best <- rep(0, p)
    if (method == "ga") {
        ga <- rbga(
            stringMin = rep(lower, p),
            stringMax = rep(upper, p),
            evalFunc = cost,
            iters = iters,
            popSize = popsize,
            verbose = verbose
        )
        best <- ga$population[1, ]
    } else if (method == "de") {
        de <- DEoptim(
            fn = cost,
            lower = rep(lower, p),
            upper = rep(upper, p),
            control = DEoptim.control(
                itermax = iters,
                NP = popsize,
                trace = verbose
            )
        )
        best <- de$optim$bestmem
    } else if (method == "hj") {
        hj_result <- hjkb(
            fn = cost,
            control = list(
                maximize = FALSE
            ),
            lower = rep(lower, p),
            upper = rep(upper, p),
            par = runif(p, min = lower, max = upper)
        )
        print("besting")
        best <- hj_result$par
    }
    newbetas <- cstep(best, csteps = csteps)
    res <- y - x %*% newbetas
    crit <- sum(sort(res^2)[1:h])
    result <- list(
        coefficients = as.double(newbetas),
        crit = crit,
        method = method
    )
    return(result)
}