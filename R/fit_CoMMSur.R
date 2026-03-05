#' Fit CoMM-Surv and return mediator-level posterior probabilities
#'
#' Fits the CoMM-Surv Level-3 model in Stan via \code{rstan} and returns
#' mediator-level posterior summaries, along with an overall convergence
#' diagnostics summary (Rhat/ESS percentiles and divergent transition rate).
#'
#' @param X Numeric vector of length N (exposure).
#' @param M Numeric matrix of dimension N x K (mediators).
#' @param time Numeric vector of length N (observed time).
#' @param event Integer/logical vector of length N (1 = event, 0 = censored).
#' @param mediator_names Optional character vector of length K. If NULL, \code{colnames(M)} are used,
#'   otherwise defaults to \code{"M1"},...,\code{"MK"}.
#' @param chains,iter,warmup,seed Stan sampling controls passed to \code{rstan::sampling()}.
#' @param ... Additional arguments forwarded to \code{rstan::sampling()}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{results}: data.frame with columns \code{mediator}, \code{p_null (= 1 - p_alt)},
#'         and posterior mean probabilities for each of the 9 joint components
#'         (neg/null/pos \eqn{\times} neg/null/pos).
#'   \item \code{overall_stats}: tibble with overall sampling diagnostics, including
#'         Rhat percentiles (median, 90/95/99%), ESS percentiles (1/5/10/50/90%),
#'         and the number/rate of divergent transitions (post-warmup).
#'   \item \code{fit}: \code{rstan} fitted model object.
#' }
#'
#' @export
fit_CoMMSur <- function(X, M, time, event,
                        mediator_names = NULL,
                        chains = 4, iter = 8000, warmup = floor(iter/2),
                        seed = 1, ...) {
  
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Package 'rstan' is required. Please install it first.")
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required. Please install it first.")
  }
  
  # -------------------------
  # Basic input checks
  # -------------------------
  if (!is.numeric(X)) stop("X must be numeric.")
  if (!is.matrix(M)) stop("M must be a matrix (N x K).")
  if (!is.numeric(time)) stop("time must be numeric.")
  if (!(is.numeric(event) || is.logical(event) || is.integer(event))) {
    stop("event must be 0/1 (numeric/integer/logical).")
  }
  
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
  
  # -------------------------
  # Compile + sample
  # -------------------------
  sm <- rstan::stan_model(file = stan_file)
  
  # Quiet run: suppress console messages/warnings and suppress progress output
  fit <- suppressMessages(
    suppressWarnings(
      rstan::sampling(
        object = sm,
        data = stan_data,
        chains = chains,
        iter = iter,
        warmup = warmup,
        seed = seed,
        control = list(adapt_delta = 0.99, max_treedepth = 18),
        refresh = 0,
        ...
      )
    )
  )
  
  # -------------------------
  # overall_stats (Rhat/ESS percentiles + divergences)
  # -------------------------
  sum_fit <- rstan::summary(fit)$summary
  rhat <- sum_fit[, "Rhat"]
  ess  <- sum_fit[, "n_eff"]
  ok <- is.finite(rhat) & is.finite(ess)
  
  rhat_q <- stats::quantile(rhat[ok], probs = c(.5, .9, .95, .99), na.rm = TRUE)
  ess_q  <- stats::quantile(ess[ok],  probs = c(.01, .05, .1, .5, .9), na.rm = TRUE)
  
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  n_post <- sum(vapply(sp, nrow, integer(1)))
  divergent_n <- sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1)))
  divergent_rate <- if (n_post > 0) divergent_n / n_post else NA_real_
  
  overall_stats <- tibble::tibble(
    metric = c(
      "rhat_median","rhat_p90","rhat_p95","rhat_p99",
      "ess_p01","ess_p05","ess_p10","ess_median","ess_p90",
      "divergent_n","divergent_rate",
      "n_post_warmup_draws_total"
    ),
    value = c(
      unname(rhat_q[1]), unname(rhat_q[2]), unname(rhat_q[3]), unname(rhat_q[4]),
      unname(ess_q[1]),  unname(ess_q[2]),  unname(ess_q[3]),  unname(ess_q[4]),  unname(ess_q[5]),
      divergent_n, divergent_rate,
      n_post
    )
  )
  
  # -------------------------
  # Extract generated quantities
  # -------------------------
  ext <- rstan::extract(fit, pars = c("p_alt", "joint_comp_probs"), permuted = TRUE)
  p_alt_draws <- ext$p_alt                 # iterations x K
  joint_draws <- ext$joint_comp_probs      # iterations x K x 9
  
  if (is.null(dim(p_alt_draws))) p_alt_draws <- matrix(p_alt_draws, ncol = K)
  if (length(dim(joint_draws)) != 3L) stop("Unexpected shape for joint_comp_probs draws.")
  
  p_alt_mean <- colMeans(p_alt_draws)
  p_null_mean <- 1 - p_alt_mean
  
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
  
  list(
    results = out_df,
    overall_stats = overall_stats,
    fit = fit
  )
}