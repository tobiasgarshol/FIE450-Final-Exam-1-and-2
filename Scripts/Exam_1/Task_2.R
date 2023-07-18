rm(list=ls())
options(scipen = 999)
################################################################################
# ASSIGNMENT 1 - TASK 2
################################################################################

# Load data
euro_stoxx_daily <- read.csv("euro_stoxx_daily_1.csv")       # Source: Refinitiv Eikon
euro_stoxx_monthly <- read.csv("euro_stoxx_monthly_1.csv")   # Source: Refinitiv Eikon


#1. 
################################ VOLATILITY ###################################

###################### EXPLANATION OF VOLATILITY MEASURE ######################

## The product's maturity date is 2026-09-30, so we get a forecast period of 
## approx. 3,5 years. To get the most accurate estimate of the volatility of 
## the underlying asset, we should consider options with an expiration date that 
## is close to 3,5 years from now. Since there are no options traded on the 
## underlying asset with a maturity date in 3,5 years, we have decided to use 
## the GARCH model instead of option implied volatility to measure the
## volatility of the underlying asset. 

############################# GARCH VOLATILITY ################################

# Compute log returns
euro_stoxx_daily$r = c(NA, diff(log(euro_stoxx_daily$Close)))

# Remove NA's, and storing log returns in a vector
r = na.omit(euro_stoxx_daily$r)

# Function to compute GARCH variance
garch.var <- function(r, omega, alpha, beta) {
  sigma2 <- r[1]^2
  for (i in 2:length(r)) {
    sigma2 <- c(sigma2, omega + alpha*r[i]^2 + beta*sigma2[i - 1])
  }
  return(sigma2)
}

# Log-likelihood function for GARCH volatility
garch.ll.fun <- function(par, r) {
  omega <- par[1]
  alpha <- par[2]
  beta <- par[3]
  sigma2 <- garch.var(r, omega, alpha, beta)
  r <- r[-1] # Dropping the first observation because of missing forecast
  sigma2 <- sigma2[-length(sigma2)] # Dropping last observation 
  ll <- sum(-log(sigma2) - r^2/sigma2) # Log likelihood function
  return(-ll)
}

# Optimizing GARCH parameters
res = nlminb(c(0.001, 0.3, 0.3), garch.ll.fun, lower = 1e-06,
             upper = 1 - 1e-06, r = r)

# Saving the weights of the GARCH parameters
omega = res$par[1]
alpha = res$par[2]
beta = res$par[3]

# Compute gamma (alpha + beta + gamma = 1)
gamma = 1 - alpha - beta

# Compute long term variance
VL = omega/gamma

# Compute annualized GARCH volatility
garch.vol = sqrt(VL*250)
garch.vol

#2.
############################# EXPLANATION OF RF ###############################
# Ideally, it would be more appropriate to use a 3.5-year EURIBOR interest rate 
# swap as a benchmark, given that it aligns with the maturity of the option. 
# However, we encountered difficulty in obtaining this instrument and thus had 
# to resort to using a 3 year EURIBOR interest rate swap as a proxy. We used the 
# average of the bid ask prices on 31.01.2023 (0.03109 + 0.03119)/2 = 0.03114. 
# Source: Refinitiv Eikon; Ticker: EURAB6E3Y. 

###############################################################################


#3. 
################################ VALUATION ####################################

############################# PAYOFF FUNCTION #################################

K_u = 3449.997                 # Upper strike
K_l = 2638.233                 # Lower strike
p0 = 4058.82                   # Initial reference level
S = seq(1500, 4500, by = 1)    # underlying price
p_n = 100                      # Notional amount
max_p = 124                    # Max payoff
min_p = 100                    # Low cash payment

# Payoff function
payoff <- ifelse(S >= K_u, max_p, ifelse(S >= K_l, min_p, pmin(min_p,  p_n * (S / p0))))

# Payoff diagram
plot(S, payoff, type = "l", xlab = "EURO STOXX 50", ylab = "Payoff", 
     main = "Payoff Plot of Option", col = "red", lwd=3)


################################ MONTE CARLO ##################################

set.seed(1)
# Variables
S0 =  euro_stoxx_daily$Close[which(euro_stoxx_daily$Exchange.Date == "2023-01-31")] 
rf =  0.03114 # Yield on EURIBOR 3Y IRS as of 31.01.2023 (Source: Refinitiv Eikon; Ticker: EURAB6E3Y)
sigma = garch.vol #iv$estimate 
dt = 2/250
T = as.numeric(as.Date('2026-09-30') - as.Date('2023-01-31'))/365
n = 5000
f = 1/100


# Function to simulate prices:
simulate.paths.fast = function(S0, rf, sigma, dt, T, n) {
  t = seq(dt, T, by = dt)
  m = length(t)
  e = exp((rf - 0.5*sigma^2)*dt + sigma * sqrt(dt) * rnorm(m*n))
  e = matrix(e, nrow = m, ncol = n)
  S = S0 * apply(e, 2, cumprod)
  S = rbind(S0, S)
  return(S)
}

S = simulate.paths.fast(S0 = S0, 
                        rf = rf, 
                        sigma = sigma, 
                        dt = dt,
                        T = T, 
                        n = n)

# Plot of S (n number of simulations):
matplot(as.data.frame(S), type = 'l')

# Find dates:
b1_d <- ceiling((as.numeric(as.Date('2023-09-27') - as.Date('2023-01-31'))/365)/dt)
b2_d <- ceiling((as.numeric(as.Date('2024-09-25') - as.Date('2023-01-31'))/365)/dt)
b3_d <- ceiling((as.numeric(as.Date('2025-09-30') - as.Date('2023-01-31'))/365)/dt)

# Autocall barrier levels:
b1 <- 4058.82
b2 <- 3855.879
b3 <- 3652.938

# Autocall barrier payments:
auto1 <- 109.60
auto2 <- 114.40
auto3 <- 119.20

# Below we caluclate the payoff if the underlying hits any of the options (three) 
# remaining autocall barrier levels. The final payoff (if the option reaches maturity) is
# sim4. Sim1 contains the autocall payment that occours if the underlying hits autocall barrier
# level 1 (i.e., price of 4058.82 at 27. sept. 2023)

S = data.frame(S)

# Simulations hitting barrier 1
sim1 <- S[b1_d, which(S[b1_d, ] > b1)]
sim1[ ,1:ncol(sim1)] = auto1
sim1_0 <- exp(-rf*b1_d/365) *sim1

# Simulations hitting barrier 2
sim2 <- S[b2_d, setdiff(which(S[b2_d, ] > b2), which(S[b1_d, ] > b1))]
sim2[ ,1:ncol(sim2)] = auto2
sim2_0 <- exp(-rf*b2_d/365) *sim2

# Simulations hitting barrier 3
sim3 <- S[b3_d, setdiff(setdiff(which(S[b3_d, ] > b3), which(S[b1_d, ] > b1)), which(S[b2_d, ] > b2))]
sim3[ ,1:ncol(sim3)] = auto3
sim3_0 <- exp(-rf*b3_d/365) *sim3

# Simulations not hitting a barrier
sim4 <- S[nrow(S), setdiff(seq_len(ncol(S)), c(which(S[b1_d, ] > b1), which(S[b2_d, ] > b2), which(S[b3_d, ] > b3)))]
sim4 = ifelse(sim4[nrow(sim4), ] >= K_u, max_p, 
           ifelse(sim4[nrow(sim4), ] >= K_l, min_p, pmin(min_p,  p_n * sim4[nrow(sim4), ] / p0)))
sim4_0 <- exp(-rf*T)*sim4

# Compute the value of the option at time = 0
X0 <- cbind(sim1_0, sim2_0, sim3_0, data.frame(sim4_0))
X0 <- as.numeric(X0)

# Histogram: 
hist(X0, breaks = 100)

# Expected price: 
mean(X0)

# Calculate the Monte Carlo Error
SE = sd(X0)/sqrt(length(X0)) # Standard error of the simulation

# Confidence interval: 
alpha = 0.99                                  # Confidence level
z = -qnorm((1-alpha)/2)                       # Z-score
c(mean(X0) - z*SE, mean(X0), mean(X0) + z*SE) # Lower limit, mean, and upper limit

conf_int = matrix(c(mean(X0) - z*SE, mean(X0), mean(X0) + z*SE), nrow = 3, ncol = 1)
rownames(conf_int) <- c("lower limit", "mean", "upper limit")
colnames(conf_int) <- c("without variance reduction")
conf_int <- t(conf_int)
conf_int

uncertainty = ((mean(X0) + z*SE) - (mean(X0) - z*SE)) / (mean(X0))
uncertainty # We can decrease the uncertainty by increasing n

########################### WITH VARIANCE REDUCTION ###########################
# Function for MC w/ antithetic sampling:
simulate.paths.fast.as = function(S0, rf, sigma, dt, T, n) {
  
  n = n # Cut the random sampling in two
  
  t = seq(dt, T, by = dt)
  
  m = length(t)
  
  z = rnorm(n*m) # 2500 random numbers
  
  z.as = -z # antithetic numbers
  
  Z = matrix(c(z, z.as), nrow = m, ncol = n*2)
  
  e = exp((rf- 0.5*sigma^2)*dt + sigma*sqrt(dt)*Z)
  
  S = apply(e, # matrix
            2, # column 
            cumprod) # function
  S = rbind(S0, S0*S)
  
  return(S)
}

# Calculate option price
set.seed(1)
n = 2500
S.as = simulate.paths.fast.as(S0, rf, sigma, dt, T, n = n)

S.as = data.frame(S.as)

# Calculate the option value at maturity for each simulation: 
# Simulations hitting barrier 1
sim1 <- S.as[b1_d, which(S.as[b1_d, ] > b1)]
sim1[ ,1:ncol(sim1)] = auto1
sim1_0 <- exp(-rf*b1_d/365) *sim1

# Simulations hitting barrier 2
sim2 <- S.as[b2_d, setdiff(which(S.as[b2_d, ] > b2), which(S.as[b1_d, ] > b1))]
sim2[ ,1:ncol(sim2)] = auto2
sim2_0 <- exp(-rf*b2_d/365) *sim2

# Simulations hitting barrier 3
sim3 <- S.as[b3_d, setdiff(setdiff(which(S.as[b3_d, ] > b3), which(S.as[b1_d, ] > b1)), which(S.as[b2_d, ] > b2))]
sim3[ ,1:ncol(sim3)] = auto3
sim3_0 <- exp(-rf*b3_d/365) *sim3

# Simulations not hitting a barrier
sim4 <- S.as[nrow(S.as), setdiff(seq_len(ncol(S.as)), c(which(S.as[b1_d, ] > b1), which(S.as[b2_d, ] > b2), which(S.as[b3_d, ] > b3)))]
sim4 = ifelse(sim4[nrow(sim4), ] >= K_u, max_p, 
              ifelse(sim4[nrow(sim4), ] >= K_l, min_p, pmin(min_p,  p_n * sim4[nrow(sim4), ] / p0)))
sim4_0 <- exp(-rf*T)*sim4

# Compute the value of the option at time = 0
X0.as <- cbind(sim1_0, sim2_0, sim3_0, data.frame(sim4_0))
X0.as <- as.numeric(X0)

# Distribution:
hist(X0.as, breaks = 100)

# Estimated option price:
mean(X0.as)

# Store result for later:
X0_MC_AS <- mean(X0.as)

# Now we cant calculate the SE as above, as the antithetic values are dependent and not IID. 
# We can therefore do it as such where we remove the antithetic values and onlye use the "true"
# original values: 
X.pairs = (X0.as[1:n] + X0.as[(n + 1):(2*n)])/2 # Remove the antithetic values so that we are left with the originals.

SE.as = sd(X.pairs)/sqrt(n)
SE.as

# Confidence interval:
alpha <- 0.99
z <- -qnorm((1 - alpha)/2)
conf_int.as = matrix(c(mean(X0.as) - z * SE, mean(X0.as) ,mean(X0.as) + z * SE), nrow = 3, ncol = 1)
rownames(conf_int.as) <- c("lower limit", "mean", "upper limit")
colnames(conf_int.as) <- c("with antithetic sampling")
conf_int.as <- t(conf_int.as)

# Calculate the uncertainty in the estimate:
uncertainty.as = ((mean(X0.as) + z*SE.as) - (mean(X0.as) - z*SE.as)) / (mean(X0.as))
uncertainty.as # We can decrease the uncertainty by increasing n

# Compare confidence interval with and without antithetic sampling: 
comp <- rbind(conf_int, conf_int.as)
comp <- cbind(comp, c(uncertainty, uncertainty.as))
colnames(comp) <- c("lower limit", "mean", "upper limit", "uncertainty")
comp






