#' Can an SPDE version of centroid and integrated kernel models be developed?

#' Specify the spatial structure (e.g., Matérn)
spatial_model <- inla.spde2.matern(
  param = c(1, 1), # Range and smoothness hyperparameters
  control.family = list(theta = TRUE) # Include theta hyperparameter
)

# Define the lengthscale hyperparameter as a random effect (latent field)
hyperparam_formula <- formula(log_lengthscale ~ 1)  # Formula specifying the random effect
hyperparam_prior <- list(prec = list(initial = 1))  # Prior for the lengthscale

# Specify the likelihood for the observed data
likelihood <- list(
  link = "identity", # Identity link function (for Gaussian)
  precision = list(prior = "gamma") # Prior for the precision parameter
)  

# Combine the components to create the complete latent Gaussian model
fit <- inla(
  formula = y ~ 1,
  data = dat,
  family = "gaussian",
  control.predictor = list(compute = TRUE, link = 1),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  control.inla = list(
    h = 0.01,  # Resolution of the spatial mesh
    strategy = "adaptive",
    int.strategy = "eb" # Spatial integration strategy
  )
)  
