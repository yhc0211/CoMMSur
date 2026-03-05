# CoMMSur

CoMM-Surv (Configuration-based Mediation Model for Survival Outcomes) is a Bayesian framework for mediation analysis with right-censored survival data and multiple mediators.
The model jointly characterizes the exposure–mediator and mediator–outcome pathways within a unified probabilistic framework. For each mediator, the corresponding pathway effects are assigned to one of three states—negative, null, or positive—resulting in nine latent configurations that define a structured prior on the pair of regression coefficients. This configuration-based formulation enables flexible modeling of heterogeneous mediation patterns while borrowing information across mediators.
Posterior probabilities derived from the model quantify the evidence that a mediator belongs to one of the non-null configurations corresponding to an active mediation pathway. These probabilities provide a principled basis for identifying mediators that transmit the exposure effect to survival outcomes.

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