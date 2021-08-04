library("testthat")
library("MASS")
library("dfoptim")

test_that("galts - Random regression for p = 2", {
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  e <- rnorm(n)
  y <- 5 + 5 * x1 + 5 * x2 + e
  mydata <- data.frame(x1, x2, y)
  result <- ga.lts(y ~ x1 + x2, mydata, lower = -100, upper = 100)
  expect_true(all(result$coefficients < c(6, 6, 6)))
  expect_true(all(result$coefficients > c(4, 4, 4)))
})


test_that("galts - phone data with ga", {
  # > ltsreg(calls ~ year, phones)
  #  Call:
  # lqs.formula(formula = calls ~ year, data = phones, method = "lts")
  # Coefficients:
  # (Intercept)         year
  #    -56.162        1.159
  #  Scale estimates 1.249 1.131
  result <- ga.lts(calls ~ year,
    phones,
    lower = -100,
    upper = 100,
    csteps = 10,
    popsize = 50,
    method = "ga"
  )
  # print(result)
  expect_true(all(result$coefficients < c(-40, 2)))
  expect_true(all(result$coefficients > c(-56, 1)))
})