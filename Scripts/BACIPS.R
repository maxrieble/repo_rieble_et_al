################################################################################
# BACIPS.R (Model Definition)
#
# DESCRIPTION:
# Defines and runs the hierarchical Bayesian Progressive-Change BACIPS model
# using JAGS, based on Thiault et al. (2017). Called by '01_Run_BACIPS_Model.R'.
#
# USAGE / DEPENDENCIES:
# Before running (via source()), the following objects must exist:
#   - `TABLE_DATA_DELTA`: Data frame of movement deltas loaded by the calling script.
#   - `species_num`: Integer (1=Perch, 2=Carp) set by the calling script.
#
# OUTPUT:
# Creates two objects in the calling environment:
#   - `out`: The raw output from the R2jags::jags() function.
#   - `All`: A summary data frame of the model's posterior distributions
#            (i.e., out$BUGSoutput$summary).
#
# REFERENCE:
# Thiault, L. et al. (2017) Progressive-Change BACIPS. Methods Ecol Evol.
# DOI: 10.1111/2041-210X.12655
#
################################################################################

# --- 1. Load Libraries and Prepare Data ---
library(R2jags)

# Prepare time vectors
length.time.before <- length(unique(TABLE_DATA_DELTA$TIME[TABLE_DATA_DELTA$TIME < 0]))
length.time.after <- length(unique(TABLE_DATA_DELTA$TIME[TABLE_DATA_DELTA$TIME > 0]))
time.after <- seq(1:length.time.after)

# Format data list for JAGS
TABLE_DATA_DELTA <- as.data.frame(TABLE_DATA_DELTA)
DATA.rep <- list(
  Delta.before = NULL,
  ref.before = NULL,
  Delta.after = NULL,
  time.after = NULL,
  ID.fish = NULL
)

ind <- length(unique(TABLE_DATA_DELTA$id)) # Number of capture events

for (i in 1:ind) {
  temp_before <- which(TABLE_DATA_DELTA$id == i & TABLE_DATA_DELTA$TIME < 0)
  DATA.rep$Delta.before <- c(DATA.rep$Delta.before, mean(TABLE_DATA_DELTA[temp_before, "TDELTA2"], na.rm = TRUE))
  DATA.rep$ref.before <- c(DATA.rep$ref.before, TABLE_DATA_DELTA[temp_before, "TDELTA2"])
  
  temp_after <- which(TABLE_DATA_DELTA$id == i & TABLE_DATA_DELTA$TIME > 0)
  DATA.rep$Delta.after <- c(DATA.rep$Delta.after, TABLE_DATA_DELTA[temp_after, "TDELTA2"])
  
  DATA.rep$ID.fish <- c(DATA.rep$ID.fish, rep(i, length(temp_after)))
  DATA.rep$time.after <- c(DATA.rep$time.after, time.after)
}

DATA.rep$n.fish <- ind
DATA.rep$n.obs <- length(DATA.rep$Delta.after)

# Calculate range for priors based on reference period
Krange.i <- data.frame(Min = NA, Max = NA)
for (i in 1:ind) {
  ref_indices <- (length.time.before * (i - 1) + 1):(length.time.before * i)
  range_val <- tryCatch(
    range((DATA.rep$ref.before[ref_indices] + DATA.rep$Delta.before[i]) / 2, na.rm = TRUE),
    warning = function(w) c(NA, NA)
  )
  Krange.i[i, 1:2] <- ifelse(is.finite(range_val), range_val, c(NA, NA))
}
DATA.rep$Min <- ifelse(is.na(Krange.i$Min), 1000, Krange.i$Min) # Fallback default
DATA.rep$Max <- ifelse(is.na(Krange.i$Max), 14000, Krange.i$Max) # Fallback default

# --- 2. Define JAGS Model ---

model_string <- "
model {
  # Likelihood
  for( i in 1:n.obs) {
    Delta.after[i] ~ dnorm(Delta.after.hat[i], tol.resid)
    Delta.after.hat[i] <- ifelse(
      time.after[i] == 1,
      # Eq for t=1
      kone.i[id[i]] + rr.i[id[i]] * kone.i[id[i]] * (1 - (kone.i[id[i]] / ktwo.i[id[i]])) - (alpha.i[id[i]] / (1 + beta.i[id[i]] * time.after[i])) * kone.i[id[i]],
      # Eq for t>1
      Delta.after.hat[i-1] + rr.i[id[i]] * Delta.after.hat[i-1] * (1 - (Delta.after.hat[i-1] / ktwo.i[id[i]])) - (alpha.i[id[i]] / (1 + beta.i[id[i]] * time.after[i])) * Delta.after.hat[i-1]
    )
  }

  # Population-level Priors
  alpha ~ dunif(-2, 1)
  beta ~ dunif(0.000001, 1.2)
  rr ~ dunif(0.5, 1)
  tol.resid ~ dgamma(0.001, 0.001)
  tol.beta ~ dgamma(0.001, 0.001)

  # Individual-level Priors
  for(i in 1:n.fish) {
    kone.i.hat[i] <- Delta.before[i]
    kone.i[i] ~ dnorm(kone.i.hat[i], kone.i.sd[i]) T(1000, 14000)
    kone.i.sd[i] ~ dgamma(0.001, 0.001) T(0.0001, 10000)

    ktwo.i.hat[i] ~ dunif(Min[i], Max[i])
    ktwo.i[i] ~ dnorm(ktwo.i.hat[i] * tol.ktwo[i], log(tol.ktwo[i]))

    alpha.i[i] ~ dnorm(alpha * tol.alpha.i[i], tol.alpha.i[i]) T(-2.99, 1.249)
    beta.i[i] ~ dgamma(beta, tol.beta)
    rr.i[i] ~ dnorm(rr * tol.rr.i[i], tol.rr.i[i]) T(0, 2)

    # Species-specific priors inserted here
    %s
  }
}
"

# Inject species-specific priors
if (species_num == 1) { # Perch
  species_priors <- "
    tol.ktwo[i] ~ dunif(0.5, 2)
    tol.alpha.i[i] ~ dunif(0.01, 1.5)
    tol.rr.i[i] ~ dunif(0.001, 1.99)
  "
} else { # Carp
  species_priors <- "
    tol.ktwo[i] ~ dunif(0.5, 2)
    tol.alpha.i[i] ~ dunif(0.01, 1.5)
    tol.rr.i[i] ~ dunif(0.001, 1.99)
  "
}

# Write model to temporary file
cat(sprintf(model_string, species_priors), file = "model.txt")

# --- 3. Run JAGS Model ---

# Prepare data list for JAGS
data.jags <- list(
  n.obs = DATA.rep$n.obs,
  Delta.after = DATA.rep$Delta.after,
  Delta.before = DATA.rep$Delta.before,
  time.after = DATA.rep$time.after,
  n.fish = DATA.rep$n.fish,
  id = DATA.rep$ID.fish,
  Min = DATA.rep$Min,
  Max = DATA.rep$Max
)

# Parameters to monitor
params <- c("rr", "alpha", "beta",
            "kone.i", "ktwo.i", "alpha.i", "beta.i", "rr.i",
            "tol.resid")

# MCMC settings
ni <- 110000  # Total iterations
nt <- 10      # Thinning
nb <- 10000   # Burn-in
nc <- 5       # Chains

# Execute JAGS
out <- jags(
  data = data.jags,
  inits = NULL,
  parameters.to.save = params,
  model.file = "model.txt",
  n.chains = nc,
  n.thin = nt,
  n.iter = ni,
  n.burnin = nb,
  jags.seed = 123
)

# --- 4. Finalize Output ---

# Create the raw summary table
All <- out$BUGSoutput$summary

# Clean up temporary model file
if (file.exists("model.txt")) {
  file.remove("model.txt")
}