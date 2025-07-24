# ReER
 Renewable Estimation for Expectile Regression
 
## Installation
You can install the **development** version from
[Github](https://github.com/Weiccao/ReER)

```s
# install.packages("remotes")
remotes::install_github("Weiccao/ReER")
```

## Usage

```s
library(expectreg)
library(MASS)
library(dplyr)

# generate data
N <- 1000; k <- 10
x1 <- runif(N); x2 <- runif(N)
Xmat <- cbind(x1, x2)
X <- cbind(1, Xmat); Xc <- cbind(1, scale(Xmat, scale = FALSE))
u <- rnorm(N)
beta <- c(2,1,2); gamma <- c(1, 0, 0)
Y <- X %*% beta + X %*% gamma * u
Yc <- Xc %*% beta + Xc %*% gamma * u

# generate stream datatset
indices <- split(1:N, cut(1:N, breaks = k, labels = FALSE))
X_split <- lapply(indices, function(idx) X[idx, , drop = FALSE])
Y_split <- lapply(indices, function(idx) Y[idx])
Xc_split <- lapply(indices, function(idx) Xc[idx, , drop = FALSE])
Yc_split <- lapply(indices, function(idx) Yc[idx])

# estimate
result <- online_alg(X = X_split, Y = Y_split,
                     Xc = Xc_split, Yc = Yc_split,
                     K = k, tau = 0.25,
                     regFormula = Y ~ x1 + x2,
                     dataName = c("Y", "x0", "x1", "x2"))

Result.df <- as.data.frame(rbind(result$DCER, result$PAER, result$ReER))
rownames(Result.df) <- c("DCER", "PAER", "ReER"); colnames(Result.df) <- c("beta0", "beta1", "beta2")
Result.df
```

## License

This package is free and open source software, licensed under GPL-3.
