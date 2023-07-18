rm(list=ls())
options(scipen = 999)
################################################################################
# ASSIGNMENT 1 - TASK 1
################################################################################

euro_stoxx_daily <- read.csv("euro_stoxx_daily_1.csv")     # Source: Refinitiv Eikon
euro_stoxx_monthly <- read.csv("euro_stoxx_monthly_1.csv") # Source: Refinitiv Eikon

#1. 
######################### PAYOFF PLOT OF UNDERLYING ###########################
K_c = 3924.80                  # Call strike 
S = seq(0, 5500, by = 1)       # underlying price
p_n = 100                      # Notional amount
max_p = 124                    # Max payment
b = 3139.84                  # Put strike
min_p = 100                    # Cash payment

# Calculate the payoffs
payoff <- pmax(pmin(p_n*(S/K_c), max_p), pmin(min_p, p_n*(S/b)))

# Plot the payoffs
plot(S, payoff, type = "l", xlab = "EURO STOXX 50", ylab = "Payoff", 
     main = "Payoff Plot of Option", col = "red", lwd=3)

#2. 
########################### REPLICATING PORTFOLIO #############################
# REPLICATING PORTFOLIO: 
## 1. Long zero-strike (LEPO) call
## 2. Short a call with K = 3139.84
## 3. Long a call with K = 3924.80
## 4. Short call with K = cap = 4866.75 -- see calculations below
cap = (max_p*K_c)/p_n

# Variables for Black Scholes: 
s0 =  euro_stoxx_daily$Close[which(euro_stoxx_daily$Exchange.Date == "2023-01-31")] 
airbag = b
sigma = 0.13 # See below for estimation of sigma 
k = 3924.80
rf = 0.03363 # Yield on EURIBOR 18M IRS as of 31.01.2023 (Source: Refinitiv Eikon; Ticker: EURAB6E18M)
t = as.numeric(as.Date('2024-05-06') - as.Date('2023-01-31'))/365 # Time to maturity

# Black Scholes function:
call <- function(s0, sigma, k, rf, t) {
  
  d1 = (log(s0/k) + (rf + sigma^2/2)*t)/sigma/sqrt(t)
  d2 = d1 - sigma * sqrt(t)
  c = s0 * pnorm(d1) - exp(-rf*T)*k*pnorm(d2)
  return(c)
  
}

# Calculate the value of the portfolio:
set.seed(1)
call1 <- call(s0, sigma, k = 0, rf, t)
call2 <- call(s0, sigma, k = 3139.84, rf, t)
call3 <- call(s0, sigma, k = 3924.80, rf, t)
call4 <- call(s0, sigma, k = 4866.75, rf, t)

# We know that the strike of call1 is the same as the barrier (3139.84) and that the 
# parachute factor is 1.25, i.e. 3139.84*1.25 = 3924.80.
# Further, the slope of call1 and call2 is the barrier divided by 100 (one option controls 100 of...
# ... the underlying). And since the airbag has a parachute factor (1.25) the slope of call3 and call4
# is the barrier divided by 100 multiplied with 1.25. We illustrate this in the plot below:
plot(S, payoff, type = "l", xlab = "EURO STOXX 50", ylab = "Payoff", 
     main = "Payoff Plot of Option w/ Slopes", col = "red", lwd=3)
lines(S, S/(b/100), type = "l", col = "blue", lty=2, lwd = 3)
lines(S, S/(K_c/100), type = "l", col = "darkgreen", lty=2, lwd = 3)

# Calculate the ratio / slope of the option vs. underlying:
ratio1 = (b/100)
ratio2 = (K_c/100)
parachute_factor = 1/(ratio1/ratio2) 

# As we can see, this consistent with the product details which can be found here:
# https://www.xmarkets.db.com/LU/Product_Detail/DE000DB9U3W1?fbclid=IwAR0Q52ghQB2MCNq9mQTjkZ2JjnWxfT1qdx6Ocqt7UkfjZrnxiRidYH_Pdlo 

# Finally we calculate the price of the option portfolio (i.e., the price for one Airbag Certificate):
X0 = (call1 - call2)/ratio1 + (call3 - call4)/ratio2

# Estimated price: 
X0

# Store result for later
X0_BSM <- X0

#3.
################################ VOLATILITY ################################### 
#3.1 
############################ IMPLIED VOLATILITY ############################### 
# DEFINE PARAMETERS: 
s0 <- euro_stoxx_daily$Close[which(euro_stoxx_daily$Exchange.Date == "2023-01-31")] # price
k <- 4200 # strike 
t <- as.numeric(as.Date('2024-06-01') - as.Date('2023-01-31'))/365 # Time to maturity
rf = 0.03363 # Yield on EURIBOR 18M IRS as of 31.01.2023 (Source: Refinitiv Eikon; Ticker: EURAB6E18M)

# Define function to estimate call:
call <- function(s0, sigma, k, rf, t) {
  
  d1 = (log(s0/k) + (rf + sigma^2/2)*t)/sigma/sqrt(t)
  d2 = d1 - sigma * sqrt(t)
  c = s0 * pnorm(d1) - exp(-rf*T)*k*pnorm(d2)
  return(c)
  
}

# Market price of the call
c.market = 286.799987792969

# Argument: 
obj.fun <- function(sigma, c.market, s0, k, rf, t) {
  
  c.model <- call(s0, sigma, k, rf, t)
  eps <- (c.model - c.market)^2
  return(eps)
  
}

# Run minimize problem:
iv = nlm(obj.fun, 
         p = 0.2, 
         c.market = c.market, 
         s0 = s0,
         k = k, 
         rf = rf, 
         t = t)

# Implied volatility:
iv$estimate

# The option-implied volatility approach is considered to be the most 
# appropriate measure of future volatility of the underlying asset, due to its 
# forward-looking nature. This approach involves using options traded on the 
# underlying to estimate the volatility of the asset. We used an option with a 
# strike price that is close to the price of the underlying asset, as these 
# options typically are more liquid and therefore more accurately priced. 
# Additionally, we used an option with a maturity aligned with the maturity of 
# the option being valued, in order to accurately reflect the market's 
# expectations regarding the volatility of the underlying asset. 
# Source: Bloomberg; Ticker: SX5E 6/24 C4200.


#3.2
############################# HISTORICAL RETURNS ##############################
euro_stoxx_daily = euro_stoxx_daily[order(euro_stoxx_daily$Exchange.Date), ]
euro_stoxx_daily$ret = c(NA, diff(log(euro_stoxx_daily$Close)))

std_hist <- sd(euro_stoxx_daily$ret, na.rm = T)*sqrt(250)
std_hist

#4. 
############################# EXPLANATION OF RF ###############################
# Banks do not offer fixed deposit rates for 18 months, which is the maturity of
# the option. Thus, we use the 18M EURIBOR interest rate swap as the risk-free 
# interest rate. This rate reflect the maturity of the option better opposed to 
# using a 1y or 2y deposit rate. We used the average of the bid ask price: 
# (0.0336 + 0.03366)/2 = 0.03363. Source: Refinitiv Eikon; Ticker: EURAB6E18M.
# 

#5. 
################################ MONTE CARLO ##################################
# Variables
S0 = s0
rf = rf
sigma = iv$estimate
dt = 2/250
T = as.numeric(as.Date('2024-05-06') - as.Date('2023-01-04'))/365 # Divide by 365 since we are going to use all numbers in year fractions
n = 5000

#5.1
########################### WITHOUT VARIANCE REDUCTION ########################
# Function to simulate prices:
simulate.paths.fast = function(S0, rf, sigma, dt, T, n) {
  t = seq(dt, T, by = dt)
  m = length(t)
  e = exp((rf - 0.5*sigma^2)*dt + sigma * sqrt(dt) * rnorm(m*n))
  e = matrix(e, nrow = m, ncol = n)
  S = S0 * apply(e, 2, cumprod) # 2 = apply function columnwise, 1 = apply function rowwise
  S = rbind(S0, S)
  return(S)
}

set.seed(1)
S = simulate.paths.fast(S0 = S0, 
                    rf = rf, 
                    sigma = sigma, 
                    dt = dt,
                    T = T, 
                    n = n)

# Plot of S (n number of simulations):
matplot(as.data.frame(S), type = 'l')

# Calculate the option value at maturity for each simulation: 
X = pmax(pmin(p_n*S[nrow(S), ]/K_c, max_p), pmin(min_p, p_n*S[nrow(S), ]/b))
X0 = exp(-rf*T)*X   

# Distribution:
hist(X0, breaks = 10)

# Estimated option price: 
mean(X0)

# Store result for later
X0_MC <- mean(X0)

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

# Calculate the uncertainty in the estimate:
uncertainty = ((mean(X0) + z*SE) - (mean(X0) - z*SE)) / (mean(X0))
uncertainty # We can decrease the uncertainty by increasing n

#5.2
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

# Calculate the option value at maturity for each simulation: 
X.as = pmax(pmin(p_n*S.as[nrow(S.as), ]/K_c, max_p), pmin(min_p, p_n*S.as[nrow(S.as), ]/b))
X0.as = exp(-rf*T)*X.as

# Distribution:
hist(X0.as, breaks = 10)

# Estimated option price:
mean(X0.as)

class(X0.as)
X0.as

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

#6. 
############################## Assumptions ####################################

# 1.
# We assume that volatility is constant. According to the Black-Scholes model
# there is just one volatility that prices all options traded 
# on a given underlying asset.


# 2.
# We make a simplification by assuming that the price only changes once per day.
# The Black-Scholes formula as a stochastic differential equation (SDE) 
# describes continuously price changes (infinitesimal small changes).

# 3. 
# We make the assumption that we are in a risk-neutral world because we 
# can replicate our option with the underlying and a risk-free investment.
# This leads to two simplifications:
# 3.1 The risk-free rate is the expected return for any investment
# 3.2 The rate used for discounting the expected payoff will be the risk-free 
#    rate for any investment


# Further, we know that the actual market price is somewhat lower (around 102 EUR as of 31.01.2023).
# This is likely due to transaction costs, which we dont take into account.

################################### SUMMARY ###################################
summary <- matrix(c(X0_BSM, X0_MC, X0_MC_AS), nrow = 3, ncol = 1)
colnames(summary) <- "Estimated prices"
rownames(summary) <- c("Black Scholes w/ Replicating portfolio", "Monte Carlo", "Monte Carlo w/ Antithetic sampling")
summary
