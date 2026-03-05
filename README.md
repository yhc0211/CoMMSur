# CoMMSur


**CoMMSur** implements **CoMM-Surv**, a Bayesian survival mediation framework for identifying mediators between an exposure and a time-to-event outcome.

The method jointly models exposure–mediator and mediator–outcome pathways using a **configuration-based prior**. By placing a structured prior on the pair of coefficients corresponding to the exposure–mediator and mediator–outcome effects, CoMMSur estimates the posterior probability that each mediator is actively involved in the mediation pathway.


The model is implemented in **Stan** .

---

## Installation

```r
install.packages("remotes")
remotes::install_github("yhc0211/CoMMSur")

## Simulation example

The following example simulates a small dataset (N = 50, K = 3) and fits the model.

```r
set.seed(123)

N <- 50
K <- 3

# Exposure
X <- rnorm(N)

# Exposure-mediator effects 
alpha <- c(0.8, 0.0, -0.6)

# Mediators:
M <- sapply(1:K, function(k) alpha[k] * X + rnorm(N, 0, 1))
M <- as.matrix(M)
colnames(M) <- c("gene1", "gene2", "gene3")

# Mediator-outcome effects and direct effect (gamma)
beta  <- c(0.6, 0.0, -0.5)
gamma <- 0.2

# Linear predictor
eta <- as.vector(M %*% beta + gamma * X)

# Weibull AFT data generation (shape=rho, scale=lambda)
rho    <- 1.5
lambda <- 1.0

# Event times
time_true <- (-log(runif(N)) / (lambda * exp(eta)))^(1 / rho)

# Censoring
censor <- rexp(N, rate = 0.1)

event <- as.integer(time_true <= censor)
time  <- pmin(time_true, censor)

# Fit model
library(CoMMSur)

fit <- fit_CoMMSur(
  X = X, M = M, time = time, event = event,
  chains = 2, iter = 1000, warmup = 500, seed = 1
)

fit$results