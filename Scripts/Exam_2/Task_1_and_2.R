rm(list = ls())
options(scipen = 999)

# Loading packages
library(quadprog)
library(Matrix)

# Import data
df <- read.csv("Assignment-2-Data.csv")

#-------------------------------------------------------------------------------
# Clean data
#-------------------------------------------------------------------------------
# Set first row as column: 
colnames(df) <- df[1, ]

# Remove first row:
df <- df[-1, ]

# Clean dataframe:
summary(df)

df$Date <- as.Date(df$Date, format = "%d/%m/%Y")              # Fix date format
df[, 2:length(df)] <- lapply(df[, 2:length(df)], as.numeric)  # Fix numeric
df <- na.omit(df)                                             # Omit NA

summary(df)

# Convert USD figures to EUR
df[ ,c(2,3,5,7,8)] <- df[ ,c(2,3,5,7,8)]/df$EUR

#-------------------------------------------------------------------------------
# Task 1.1 | Plot mean-variance frontier
#-------------------------------------------------------------------------------
# Order by dates:
df <- df[order(df$Date), ]

# Estimate Returns: 
df[, c(-1, -10)] <- data.frame(apply(df[, c(-1, -10)],              # Dont apply to date column 
                             2,                                     # 2 = apply to columns
                             function(x) c(NA, diff(x) / lag(x))))  # Estimate the difference between x value and lag x value

# As EURIBOR is already annualized we want to convert it to percent: 
df[, "1M Euribor"] <- df[, "1M Euribor"]/100
df[, "1M Euribor"] <- df[, "1M Euribor"]/12                         # Calculate monthly NIBOR Rate

# Calculate expected risk free rate:
rf = mean(df$`1M Euribor`)

# Calculate excess returns
df[, 2:8] <- df[, 2:8] - (df[, 10])
# The reason we use excess returns and not total returns is that we are working
# with a global portfolio, and using total returns introduces a bias that will cause
# optimal portfolios to overweight domestic cash securities at the expense
# of higher-risk instruments and instruments abroad.
# Source: https://www.sciencedirect.com/science/article/abs/pii/S1058330002000472

# Remove first row as it is na: 
df <- na.omit(df)

# Plot the mean-variance frontier: 
R = df[, 2:8]
mu <- apply(R, 2, mean) * 12 
Rho <- cor(R, use = "pairwise.complete.obs")
Sigma <- cov(R, use = "pairwise.complete.obs") * 12

A <- t(rbind(1, mu))
d <- rep(0, length(mu))

mu.p.vec <- seq(0, 0.2, length = 100)
sigma.p.vec <- c()

Sigma2 <- nearPD(Sigma)$mat
Sigma2 <- as.matrix(Sigma2)

for(i in 1:length(mu.p.vec)) {
  mu.star <- mu.p.vec[i]
  b0 <- c(1, mu.star)
  res <- solve.QP(Dmat = Sigma2, dvec = d, Amat = A, bvec = b0, meq = 2)
  omega <- res$solution
  sigma.p.vec <- c(sigma.p.vec, sqrt(t(omega) %*% Sigma2 %*% as.matrix(omega)))
}

plot(sigma.p.vec, mu.p.vec, type = "l", xlim = c(0, max(sigma.p.vec)),
     xlab = "Volatility", ylab = "Expected return")

#-------------------------------------------------------------------------------
# Task 1.2 | Portfolio with target return of 7% p.a.
#-------------------------------------------------------------------------------
# NOTE: WE ASSUME THAT WE DONT HAVE LEVERAGE CONSTRAINT!
# I.E., PORTFOLIO WEIGHTS CAN BE >100%

mu.star <- 0.07
d <- rep(0, length(mu))

A <- t(rbind(1, mu))
b0 <- c(1, mu.star)
res <- solve.QP(Dmat = Sigma, dvec = d, Amat = A, bvec = b0, meq = 2)

# Weights:
omega <- res$solution
sum(omega)

# Calculate expected return: 
mu.p = t(omega) %*% as.matrix(mu)

# Calcualte volatility
sigma.p = sqrt(t(omega) %*% Sigma %*% as.matrix(omega))

# Calculate Sharpe: 
sr = mu.p / sigma.p
sr

# Plot the frontier:
mu.p.vec <- seq(0, 0.2, length = 100)
sigma.p.vec <- c()

for(i in 1:length(mu.p.vec)) {
  mu.star <- mu.p.vec[i]
  b0 <- c(1, mu.star)
  res <- solve.QP(Dmat = Sigma, dvec = d, Amat = A, bvec = b0, meq = 2)
  omega <- res$solution
  sigma.p.vec <- c(sigma.p.vec, sqrt(t(omega) %*% Sigma %*% as.matrix(omega)))
}

plot(sigma.p.vec, mu.p.vec, type = "l", xlim = c(0, max(sigma.p.vec)),
     xlab = "Volatility", ylab = "Expected return", lty=1, lwd = 1, col = "black")
points(sigma.p, mu.p, pch = 4, lwd = 4, col = "red")
text(sigma.p-0.025, mu.p+0.002, "Task 2", col = "red")


#-------------------------------------------------------------------------------
# Task 1.3 | Portfolio with target volatility of 15% p.a.
#-------------------------------------------------------------------------------
sigma.star <- 0.15                                           # Target volatility
mu.star <- sample(seq(0.1, 0.2, by = 0.01), 1)               # Random variable

# Create objective function: 
obj.fn <- function(mu.star, sigma.star) {
  
  A <- cbind(1, mu)
  b0 <- c(1, mu.star)
  
  res <- solve.QP(Dmat = Sigma, dvec = d, Amat = A, bvec = b0, meq = 2)
  
  omega <- res$solution
  
  sigma.model <- sqrt(t(omega) %*% Sigma %*% as.matrix(omega))
  
  eps <- (sigma.model - sigma.star)^2
  
  return(eps)

}

obj.fn(mu.star, sigma.star)
res <- nlminb(0.15, obj.fn, lower = 1e-06, upper = 1 - 1e-06, sigma.star = sigma.star)
res

b0 <- c(1,res$par)
A <- cbind(1, mu)

omega.1 <- solve.QP(Dmat = Sigma, dvec = d, Amat = A, bvec = b0, meq = 2)$solution

# weights:
omega.1
sum(omega.1)

# Calculate standard deviation:
sigma.p <- sqrt(t(omega.1) %*% Sigma %*% omega.1)
sigma.p

# Calculate expected return
mu.p = t(omega.1) %*% as.matrix(mu)
mu.p

# Calculate Sharpe: 
sr = mu.p / sigma.p
sr

# Adding this to our plot
points(sigma.p, mu.p, pch = 4, lwd = 4, col = "blue")
text(sigma.p+0.025, mu.p-0.002, "Task 3", col = "blue")


#-------------------------------------------------------------------------------
# Task 1.4 | Compute the minimum-variance portfolio
#-------------------------------------------------------------------------------
d <- rep(0, length(mu))
A <- cbind(1, sqrt(diag(Sigma)))
b0 <- c(1, 0) # We want to minimize the variance

res <- solve.QP(Dmat = Sigma, dvec = d, Amat = A, bvec = b0, meq = 1)

# weights:
omega <- res$solution
omega
sum(omega)

# Calculate standard deviation:
sigma.p <- sqrt(t(omega) %*% Sigma %*% omega)
sigma.p

# Calculate expected return
mu.p = t(omega) %*% as.matrix(mu)
mu.p

# Calculate Sharpe: 
sr = mu.p / sigma.p
sr

# Plot the minimum variance portfolio to our plot
points(sigma.p, mu.p, pch = 4, lwd = 4, col = "darkgreen")
text(sigma.p+0.023, mu.p, "Task 4", col = "darkgreen")


#-------------------------------------------------------------------------------
# Task 1.5 | Compute the maximum sharpe ratio portfolio
#-------------------------------------------------------------------------------
# Defining constraints 
mu.star = 0.1          # Random target return
A = t(rbind(mu))
b0 = c(mu.star)
d = rep(0, length(mu))

# Optimization
res <- solve.QP(Dmat = Sigma, dvec = d, Amat = A, bvec = b0, meq = 1)

# Save weights and re-scale so it sums to 1
w <- res$solution
omega <- w/sum(w)
sum(omega)

# Expected return of max sharpe portfolio
mu.max = t(omega) %*% as.matrix(mu)
mu.max

# Sigma of max sharpe portfolio
sigma.max = sqrt(t(omega) %*% Sigma %*% as.matrix(omega))
sigma.max

# Maximum sharpe ratio
sr = mu.max /sigma.max
sr


# Plot
points(sigma.max, mu.max, pch = 4, lwd = 4, col = "orange")
text(sigma.max+0.025, mu.max, "Task 5", col = "orange")


#-------------------------------------------------------------------------------
# Task 1.6 | Portfolio with target return of 7% w/ short constraints
#-------------------------------------------------------------------------------
mu.star = 0.07
d <- rep(0, length(mu))
A <- t(rbind(1, mu, diag(1, length(mu))))
b0 <- c(1, mu.star, rep(0, length(mu)))

res <- solve.QP(Dmat = Sigma, dvec = d, Amat = A, bvec = b0, meq = 2)

# weights:
omega <- res$solution
omega
sum(omega)

# Calculate standard deviation:
sigma.p <- sqrt(t(omega) %*% Sigma %*% omega)
sigma.p

# Calculate expected return
mu.p = t(omega) %*% as.matrix(mu)
mu.p

# Calculate Sharpe: 
sr = mu.p / sigma.p
sr

# Plot the minimum variance portfolio to our plot
points(sigma.p, mu.p, pch = 4, lwd = 4, col = "pink")
text(sigma.p, mu.p-0.01, "Task 6", col = "pink")

#-------------------------------------------------------------------------------
# Task 1.7 | Compute the minimum-variance 
#-------------------------------------------------------------------------------
mu.star = 0.07
d <- rep(0, length(mu))

# The below answer does not have a solution if we include a leverage constraint (!)
# i.e. weights must equal 100%. We therefore remove these constraint to be 
# able to solve the optimization problem -- and therefore assume that leverage is possible.
# This should mean that we get a higher Sharpe Ratio in task 1.7 than task 1.6 as we 
# allow for leverage in task 1.7. 

A <- cbind(mu, diag(length(mu)), -diag(length(mu)))
b0 <- c(mu.star, rep(0, length(mu)), -rep(0.3, length(mu)))

# Meq = 1, means that there is only one equality constraint (return = 7%) and that the weights must be <= 30%.
res <- solve.QP(Dmat = Sigma, dvec = d, Amat = A, bvec = b0, meq = 1)

# weights:
omega <- res$solution
omega
sum(omega)

# Calculate standard deviation:
sigma.p <- sqrt(t(omega) %*% Sigma %*% omega)
sigma.p

# Calculate expected return
mu.p = t(omega) %*% as.matrix(mu)
mu.p

# Calculate Sharpe: 
sr = mu.p / sigma.p
sr

# Plot the minimum variance portfolio to our plot
points(sigma.p, mu.p, pch = 4, lwd = 4, col = "brown")
text(sigma.p+0.01, mu.p+0.012, "Task 7", col = "brown")

#-------------------------------------------------------------------------------
# Task 2.1 | Backtest maximizing Sharpe Strategy 
#-------------------------------------------------------------------------------
w_mat <- matrix(NA, nrow = (nrow(df)-59), ncol = ncol(R))

# Calculate monthly weighting:
for(i in 1:(nrow(df) - 59)) {
  
  Sigma <- cov(df[i:(i+59), c(-1, -9,-10)], use = "pairwise.complete.obs") * 12
  mu <- apply(df[i:(i+59), c(-1, -9,-10)], 2, mean) * 12
  
  # Defining constraints 
  mu.star = 0.1          # Random target return
  A = t(rbind(mu))
  b0 = c(mu.star)
  d = rep(0, length(mu))
  
  # Optimization
  w <- solve.QP(Dmat = Sigma, dvec = d, Amat = A, bvec = b0, meq = 1)$solution
  
  omega <- w/sum(w)
  
  w_mat[i, ] <- omega 

}

# Check that all rows sum to 1:
rowSums(w_mat)

# Calculate monthly returns with weights: 
w = as.data.frame(w_mat)
rownames(w) <- as.Date(df[,1][60:nrow(df)])

head(w)
tail(w)

R = df[60:nrow(R), c(-9,-10)]
rownames(R) <- R[,1]
R <- R[,-1]

# Calculate position size in portfolio: 
port = w*R

# Calculate portfolio returns:
port$port <- rowSums(port)

# Add dates:
port$date <- df$Date[60:nrow(df)]

# Plot cumulative returns:
port$cumport <- cumprod(1 + port$port) - 1
plot(as.Date(port$date), port$cumport, type = "l")

# Expected return of portfolio: 
mu.p = mean(port$port) * 12
mu.max.sr = mu.p

# Expected sd of portfolio: 
sigma.p = sd(port$port) * sqrt(12)
sigma.max.sr = sigma.p

# SR of port: 
sr.p = mu.p/sigma.p
sr.max.sr = sr.p

#-------------------------------------------------------------------------------
# Task 2.2 | Backtest minimizing variance
#-------------------------------------------------------------------------------
w_mat <- matrix(NA, nrow = (nrow(df)-59), ncol = ncol(R))

# Calculate monthly weighting:
for(i in 1:(nrow(df) - 59)) {
 
  Sigma <- cov(df[i:(i+59), c(-1, -9,-10)], use = "pairwise.complete.obs") * 12
  mu <- apply(df[i:(i+59), c(-1, -9,-10)], 2, mean) * 12
  
  d <- rep(0, length(mu))
  A <- cbind(1, sqrt(diag(Sigma)))
  b0 <- c(1, 0) # We want to minimize the variance
  
  w <- solve.QP(Dmat = Sigma, dvec = d, Amat = A, bvec = b0, meq = 1)$solution
  
  w_mat[i, ] <- w 

}

# Check that all rows sum to 1:
rowSums(w_mat)

# Calculate monthly returns with weights: 
w = as.data.frame(w_mat)

# Calculate position size in portfolio: 
port = w*R

# Calculate portfolio returns:
port$port <- rowSums(port)

# Add dates:
port$date <- df$Date[60:nrow(df)]

# Plot cumulative returns:
port$cumport <- cumprod(1 + port$port) - 1
plot(as.Date(port$date), port$cumport, type = "l")

# Expected return of portfolio: 
mu.p = mean(port$port) * 12
mu.min.var = mu.p

# Expected sd of portfolio: 
sigma.p = sd(port$port) * sqrt(12)
sigma.min.var = sigma.p

# SR of port: 
sr.p = mu.p/sigma.p
sr.min.var = sr.p

#-------------------------------------------------------------------------------
# Task 2.3 | Backtest equal-weighting strategy
#-------------------------------------------------------------------------------
# Weightings:
w_mat <- matrix(1/ncol(R), nrow = (nrow(df)-59), ncol = ncol(R))

# Check that all rows sum to 1:
rowSums(w_mat)

# Calculate monthly returns with weights: 
w = as.data.frame(w_mat)

# Calculate position size in portfolio: 
port = w*R

# Calculate portfolio returns:
port$port <- rowSums(port)

# Add dates:
port$date <- df$Date[60:nrow(df)]

# Plot cumulative returns:
port$cumport <- cumprod(1 + port$port) - 1
plot(as.Date(port$date), port$cumport, type = "l")

# Expected return of portfolio: 
mu.p = mean(port$port) * 12
mu.ew = mu.p

# Expected sd of portfolio: 
sigma.p = sd(port$port) * sqrt(12)
sigma.ew = sigma.p

# SR of port: 
sr.p = mu.p/sigma.p
sr.ew = sr.p
#-------------------------------------------------------------------------------
# Task 2.4 | Which strategy dominates and why?
#-------------------------------------------------------------------------------
# Summary: 
summary <- matrix(c(mu.max.sr, mu.min.var, mu.ew,
                    sigma.max.sr, sigma.min.var, sigma.ew,
                    sr.max.sr, sr.min.var, sr.ew),
                  nrow = 3, ncol = 3)
colnames(summary) <- c("Mu", "Sigma", "SR")
rownames(summary) <- c("Max SR strategy", 
                       "Min Var. strategy", 
                       "Equal Weight strategy")
summary <- t(summary)

summary

# WHY THE EQUALLY WEIGHTED PORTFOLIO PERFORMS BEST:
# A potential explanation for the superior performance of the equally weighted
# portfolio can be attributed to its ability to function as a short volatility
# strategy. This approach involves the buying of securities that have
# demonstrated strong performance, while concurrently selling those that have
# displayed weaker performance. As a result of its short volatility strategy,
# the equally weighted portfolio generates a volatility risk premium, which can
# be a key contributor to its success.

# Conversely, the Maximum Sharpe Ratio strategy may be sub-optimal in this
# context, as it tends to allocate an excessive amount of assets to
# well-performing securities, while neglecting those with poorer performance.
# The potential drawbacks of this strategy can be observed in the plot from
# Task 2.2, which indicates that it is particularly susceptible to failure in
# market crashes.

# Although the Minimum Variance Portfolio strategy exhibits strong performance,
# it does not match the results of the equally weighted portfolio. Therefore,
# it can be inferred that the superior performance of the equally weighted
# portfolio may be primarily due to its short volatility strategy.
