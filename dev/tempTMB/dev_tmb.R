#' Develop the other TMB models here!

load("data/mw.rda")
mw <- sf::st_as_sf(mw)

constant_tmb(mw, its = 1000)
iid_tmb(mw, its = 1000)
besag_tmb(mw, its = 1000)
bym2_tmb(mw, its = 1000)
fck_tmb(mw, its = 1000)
fik_tmb(mw, its = 1000)
system.time({ck_tmb(mw, its = 1000)})
system.time({ik_tmb(mw, its = 1000)})

 #' Test that the cross-validation enabling code works
cv_test <- iid_tmb(mw, its = 1000, ii = NULL)
summary(cv_test)

#' Test the IK model crashing for some geometries
#' Get the geometries it crashes for
#' Work out where it might be crashing
ik_test <- ik_tmb(mw, its = 1000)

sf <- mw
L <- 10
type <- "random" # Crashes for hexagonal but not for random: something strange about the covariance created by hexagonal?
ii <- NULL

n <- nrow(sf)
samples <- sf::st_sample(sf, type = type, exact = TRUE, size = rep(L, n))
S <- sf::st_distance(samples, samples)

# Parameters of the length-scale prior
param <- arealutils:::invgamma_prior(lb = 0.1, ub = max(as.vector(S)), plb = 0.01, pub = 0.01)

# Data structure for unequal number of points in each area
sample_index <- sf::st_intersects(sf, samples)
sample_lengths <- lengths(sample_index)
start_index <- sapply(sample_index, function(x) x[1])

dat <- list(n = nrow(sf),
            y = sf$y,
            m = sf$n_obs,
            left_out = !is.null(ii),
            ii = if(!is.null(ii)) ii else 0,
            a = param$a,
            b = param$b,
            sample_lengths = sample_lengths,
            total_samples = sum(sample_lengths),
            start_index = start_index,
            S = S)

ggplot(sf) +
  geom_sf(fill = "lightgrey") +
  geom_sf(data = samples, size = 0.5) +
  labs(x = "", y = "") +
  theme_minimal() +
  labs(fill = "") +
  theme_void()

param <- list(beta_0 = 0,
              u = rep(0, dat$n),
              log_sigma_u = 0,
              log_l = 0)

obj <- TMB::MakeADFun(
  data = c(model = "integrated", dat),
  parameters = param,
  random = c("beta_0", "u"),
  DLL = "arealutils_TMBExports"
)

quad <- aghq::marginal_laplace_tmb(ff = obj, k = k, startingvalue = obj$par)
