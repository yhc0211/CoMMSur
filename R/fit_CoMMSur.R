#' Fit CoMM-Surv and return mediator-level posterior alternative probabilities
#'
#' @param X Numeric vector of length N (exposure).
#' @param M Numeric matrix of dimension N x K (mediators).
#' @param time Numeric vector of length N (observed time).
#' @param event Integer/logical vector of length N (1 = event, 0 = censored).
#' @param mediator_names Optional character vector of length K. If NULL, colnames(M) are used,
#'   otherwise defaults to "M1",...,"MK".
#' @param chains,iter,warmup,seed Stan sampling controls passed to rstan::sampling().
#' @param ... Additional arguments forwarded to rstan::sampling().
#'
#' @return A list with:
#'   \itemize{
#'     \item results: data.frame with columns mediator, p_null (= 1 - p_alt), and posterior mean
#'       probabilities for each of the 9 joint components (neg/null/pos x neg/null/pos)
#'     \item fit: rstan fit object
#'   }
#' @export
fit_CoMMSur <- function(X, M, time, event,
                        mediator_names = NULL,
                        chains = 4, iter = 2000, warmup = floor(iter/2),
                        seed = 1, ...) {

  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Package 'rstan' is required. Please install it first.")
  }

  # Basic input checks
  if (!is.numeric(X)) stop("X must be numeric.")
  if (!is.matrix(M)) stop("M must be a matrix (N x K).")
  if (!is.numeric(time)) stop("time must be numeric.")
  if (!(is.numeric(event) || is.logical(event) || is.integer(event))) stop("event must be 0/1 (numeric/integer/logical).")

  N <- length(X)
  if (nrow(M) != N) stop("nrow(M) must equal length(X).")
  if (length(time) != N) stop("length(time) must equal length(X).")
  if (length(event) != N) stop("length(event) must equal length(X).")

  K <- ncol(M)

  if (is.null(mediator_names)) {
    mediator_names <- colnames(M)
    if (is.null(mediator_names) || any(mediator_names == "")) {
      mediator_names <- paste0("M", seq_len(K))
    }
  } else {
    if (length(mediator_names) != K) stop("mediator_names must have length K = ncol(M).")
  }

  # Coerce event to integer 0/1
  event_int <- as.integer(event)
  if (any(!event_int %in% c(0L, 1L))) stop("event must be 0/1.")

  stan_data <- list(
    N = N,
    K = K,
    X = as.vector(X),
    M = M,
    time = as.vector(time),
    event = event_int
  )

  stan_file <- system.file("stan", "comm_surv_lv3.stan", package = "CoMMSur")
  if (stan_file == "") stop("Stan model file not found inside the package installation.")

  # Compile + sample
  sm <- rstan::stan_model(file = stan_file)

  fit <- rstan::sampling(
    object = sm,
    data = stan_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    seed = seed,
    ...
  )

  # Extract generated quantities
  ext <- rstan::extract(fit, pars = c("p_alt", "joint_comp_probs"), permuted = TRUE)
  p_alt_draws <- ext$p_alt                 # iterations x K
  joint_draws <- ext$joint_comp_probs      # iterations x K x 9

  # Ensure expected shapes
  if (is.null(dim(p_alt_draws))) p_alt_draws <- matrix(p_alt_draws, ncol = K)
  if (length(dim(joint_draws)) != 3L) stop("Unexpected shape for joint_comp_probs draws.")

  p_alt_mean <- colMeans(p_alt_draws)
  p_null_mean <- 1 - p_alt_mean

  # Posterior mean for each (mediator, component)
  # Result: K x 9
  joint_mean <- apply(joint_draws, c(2, 3), mean)

  comp_names <- c(
    "p_neg_neg", "p_neg_null", "p_neg_pos",
    "p_null_neg", "p_null_null", "p_null_pos",
    "p_pos_neg", "p_pos_null", "p_pos_pos"
  )
  colnames(joint_mean) <- comp_names

  out_df <- data.frame(
    mediator = mediator_names,
    p_null = as.numeric(p_null_mean),
    joint_mean,
    stringsAsFactors = FALSE
  )

  list(results = out_df, fit = fit)
}
