### Single Growth Models Comparison: Exponential, Gompertz, Logistic

# Load required libraries
library(minpack.lm)
library(readxl)


# Read the data from Excel file
data <- read_excel("dataset.xlsx", sheet = "mortality")  # You may change the path accordingly

# Compute cumulative mortality
cumulative_mortality <- cumsum(data$daily_mortality[1:358])

# Prepare data frame for modeling
days <- 1:358
mortality_data <- data.frame(day = days, cumulative = cumulative_mortality)

# Setup two plots side by side
par(mfrow = c(1, 2))

# Plot 1: Cumulative mortality
plot(mortality_data$day, mortality_data$cumulative, type = "l", lwd = 2,
     xlab = "day", ylab = "cumulative number of deaths")

# Plot 2: Daily mortality 
plot(data$day[1:358], data$daily_mortality[1:358], type = "l", lwd = 2,
     xlab = "day", ylab = "number of deaths")

# Optional: reset graphic device
dev.off()

# Define growth model functions
exp_model <- function(x, a1, a2) {
  a1 * exp(a2 * x)
}
gompertz_model <- function(x, a1, a2, a3) {
  a1 * exp(-a2 * exp(-a3 * x))
}
logistic_model <- function(x, a1, a2, a3) {
  a1 * (1 + a2 * exp(-a3 * x))^(-1)
}

# Fit exponential model
fit_exp <- nlsLM(cumulative ~ exp_model(day, a1, a2), data = mortality_data, start = c(a1 = 10, a2 = 0.01))
coef_exp <- coef(fit_exp)
predict_exp <- function(x) { coef_exp["a1"] * exp(coef_exp["a2"] * x) }

# Fit Gompertz model
fit_gomp <- nlsLM(cumulative ~ gompertz_model(day, a1, a2, a3), data = mortality_data, start = c(a1 = 100000, a2 = 10, a3 = 0.01), control=nls.control(maxiter = 100))
coef_gomp <- coef(fit_gomp)
predict_gomp <- function(x) { coef_gomp["a1"] * exp(-coef_gomp["a2"] * exp(-coef_gomp["a3"] * x)) }

# Fit Logistic model
fit_log <- nlsLM(cumulative ~ logistic_model(day, a1, a2, a3), data = mortality_data, start = c(a1 = 10000, a2 = 1, a3 = 0.001))
coef_log <- coef(fit_log)
predict_log <- function(x) { coef_log["a1"] * (1 + coef_log["a2"] * exp(-coef_log["a3"] * x))^(-1) }

# Print model evaluation metrics
cat("Exponential Model: ", AIC(fit_exp), BIC(fit_exp), "\n")
cat("Gompertz Model:    ", AIC(fit_gomp), BIC(fit_gomp), "\n")
cat("Logistic Model:    ", AIC(fit_log), BIC(fit_log), "\n")

# Plot observed and estimated curves
plot(mortality_data$cumulative, type = "l", col = "black", lwd = 2,
     xlab = "day", ylab = "cumulative number of deaths", ylim = c(0, 32000))
lines(predict_exp(days), col = "blue", lwd = 2)
lines(predict_gomp(days), col = "orange", lwd = 2)
lines(predict_log(days), col = "green", lwd = 2)

legend("bottomright",
       legend = c("observed", "exponential", "gompertz", "logistic"),
       col = c("black", "blue", "orange", "green"),
       lwd = 2, bty = "n")

# Predict next seven days using the best-fitting model
predict_log <- function(x) { coef_log["a1"] * (1 + coef_log["a2"] * exp(-coef_log["a3"] * x))^(-1) }
pred_days <- 359:365
round(predict_log(pred_days), 0)





#### Change Point Determination and Piecewise Growth Model

# Read the data from Excel file
data <- read_excel("dataset.xlsx", sheet = "mortality2")  # You may change the path accordingly

# Compute cumulative mortality
cumulative_mortality <- cumsum(data$daily_mortality)

# Prepare data frame for modeling
days <- 1:594
mortality_data <- data.frame(day = days, cumulative = cumulative_mortality)

# Generalized growth function
f.grwth <- function(x, r, p, A) {
  (r * (1 - p) * x + A)^(1 / (1 - p))
}

# Initialize an empty data frame to store the results
results <- data.frame(
  k = integer(),
  AIC = numeric(),
  BIC = numeric(),
  r = numeric(),
  p = numeric(),
  A = numeric()
)

# Iterative fitting
for (k in 3:594) {
  tryCatch({
    # Fit the model using data up to time point k
    grwth1 <- nlsLM(
      cumulative ~ f.grwth(day[1:k], r, p, A),
      data = mortality_data[1:k, ],
      start = c(r = 1.2, p = 0.1, A = 10),
      control = nls.control(maxiter = 1000)
    )
    
    # Coefficients
    coefs <- as.numeric(coef(grwth1))
    
    # Append the results for this iteration to the results data frame
    results <- rbind(results, data.frame(
      k = k,
      AIC = AIC(grwth1),
      BIC = BIC(grwth1),
      r = coefs[1],
      p = coefs[2],
      A = coefs[3]
    ))
    
  }, error = function(e) {
    message(paste("Iteration", k, "failed:", e$message))
  })
}

# Results
write.csv(results, "growth_model_results.csv", row.names = FALSE)

# In this analysis, the first change point (d1) was identified.
# To detect the subsequent change points, the same procedure was repeated. 
# For instance, since d1 was defined as day 58, the next search for a change point was conducted using data after this day. 
# Following this approach, the subsequent change points were identified as d2 = 166, d3 = 246, d4 = 299, and d5 = 358.


# Define change points
i1 <- 58
i2 <- 166

# Split data based on change points
D1 <- mortality_data[1:i1, ]
D2 <- mortality_data[(i1 + 1):i2, ]

### Exponential model for D1
exp1 <- nlsLM(cumulative ~ exp_model(day, a1, a2), data = D1, start = c(a1 = 10, a2 = 0.01))
alpha.exp1 <- as.numeric(coef(exp1))
exp.est1 <- function(d) { alpha.exp1[1] * exp(alpha.exp1[2] * d) }
print(c(AIC(exp1), BIC(exp1)))

### Gompertz model for D1
gomp1 <- nlsLM(cumulative ~ gompertz_model(day, a1, a2, a3), data = D1, start = c(a1 = 100000, a2 = 10, a3 = 0.01))
alpha.gomp1 <- as.numeric(coef(gomp1))
gomp.est1 <- function(d) { alpha.gomp1[1] * exp(-alpha.gomp1[2] * exp(-alpha.gomp1[3] * d)) }
print(c(AIC(gomp1), BIC(gomp1)))

### Logistic model for D1
log1 <- nlsLM(cumulative ~ logistic_model(day, a1, a2, a3), data = D1, start = c(a1 = 10000, a2 = 1, a3 = 0.001))
alpha.log1 <- as.numeric(coef(log1))
log.est1 <- function(d) { alpha.log1[1] * (1 + alpha.log1[2] * exp(-alpha.log1[3] * d))^(-1) }
print(c(AIC(log1), BIC(log1)))

# Plot the results for D1 - Gompertz model (you can replace with other models if needed)
plot(mortality_data[1:58,1], mortality_data[1:58,2], pch = 19, col = "grey60", xlab = "day", ylab = "cumulative number of deaths")
points(1:58 , gomp.est1(1:58), type="l", col = "darkblue", lwd = 2)
legend("bottomright",
       legend = c("observed", "predicted"),
       lty = c(NA, 1), lwd = c(NA, 2),
       pch = c(19, NA), col = c("grey60", "darkblue"), bty = "n")



### Exponential model for D2
exp2 <- nlsLM(cumulative ~ exp_model(day, a1, a2), data = D2, start = c(a1 = 10, a2 = 0.01))
alpha.exp2 <- as.numeric(coef(exp2))
exp.est2 <- function(d) { alpha.exp2[1] * exp(alpha.exp2[2] * d) }
print(c(AIC(exp2), BIC(exp2)))

### Gompertz model for D2
gomp2 <- nlsLM(cumulative ~ gompertz_model(day, a1, a2, a3), data = D2, start = c(a1 = 100000, a2 = 10, a3 = 0.01))
alpha.gomp2 <- as.numeric(coef(gomp2))
gomp.est2 <- function(d) { alpha.gomp2[1] * exp(-alpha.gomp2[2] * exp(-alpha.gomp2[3] * d)) }
print(c(AIC(gomp2), BIC(gomp2)))

### Logistic model for D2
log2 <- nlsLM(cumulative ~ logistic_model(day, a1, a2, a3), data = D2, start = c(a1 = 10000, a2 = 1, a3 = 0.001))
alpha.log2 <- as.numeric(coef(log2))
log.est2 <- function(d) { alpha.log2[1] * (1 + alpha.log2[2] * exp(-alpha.log2[3] * d))^(-1) }
print(c(AIC(log2), BIC(log2)))

# Plot the results for D2 - Gompertz model
plot(mortality_data[1:166,1], mortality_data[1:166,2], pch = 19, col = "grey60", xlab = "day", ylab = "cumulative number of deaths")
points(gomp.est1(1:58) , type="l", col = "darkblue", lwd = 2)
points(59:166, gomp.est2(59:166), type="l", col = "green", lwd = 2)
legend("bottomright", legend = c("observed", "predicted for d1", "predicted for d2"),
       lty = c(NA, 1, 1), lwd = c(NA, 2, 2),
       pch = c(19, NA, NA), col = c("grey60", "darkblue", "green"), bty = "n")


# The suitability of growth models is evaluated for each change point, and the best-fitting model is selected accordingly. 
# For illustrative purposes, code is provided only for two out of five change points. 
# The same analysis should be repeated for the remaining change points to obtain complete results.