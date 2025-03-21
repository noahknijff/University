library(ggplot2)
library(gridExtra)
library(ggpubr)
library(MASS)
library(lmtest)
library(latex2exp)


options(echo = FALSE)

# Close the current plot device if one is open
if (dev.cur() > 1) {
   dev.off()
}

stormer <- get("stormer", envir = asNamespace("MASS"))


# EDA
print("EDA")


print(head(stormer))
print(summary(stormer))

# print(colnames(stormer))


cat("\n \n \n")



plot_residual_diagnostics <- function(model, exercise_name = "Exercise") {
  # Extract residuals and fitted values
  if ("nls" %in% class(model)) {
    fitted_vals <- fitted(model)
    residuals_vals <- residuals(model)
  } else {
    fitted_vals <- model$fitted.values
    residuals_vals <- residuals(model)
  }
  
  # Residuals vs Fitted (ggplot)
  p1 <- ggplot(data.frame(Fitted = fitted_vals, Residuals = residuals_vals), 
               aes(x = Fitted, y = Residuals)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste("(", exercise_name, ") Residuals vs Fitted"), x = "Fitted Values", y = "Residuals") +
    theme_minimal()
  
  # Normal Q-Q plot (ggplot)
  p2 <- ggplot(data.frame(Residuals = residuals_vals), aes(sample = Residuals)) +
    stat_qq() +
    stat_qq_line(color = "red") +
    labs(title = paste("(", exercise_name, ") Normal Q-Q")) +
    theme_minimal()
  
  # Arrange the plots side by side
  combined_plot <- ggarrange(p1, p2, ncol = 2)
  
  # Save with ggsave
  file_path <- paste0("Figures/residuals_diagnostics_", exercise_name, ".png")
  ggsave(file_path, plot = combined_plot, width = 12, height = 6, dpi = 300, bg = "white")
  
  return(combined_plot)
}




# Exercise A
print("Exercise A")


print(ggplot(stormer, aes(x = Viscosity, y = Time, color = Wt)) +
        geom_point(size = 3) +
        labs(title = "Stormer Viscometer Data",
             x = "Viscosity (v)",
             y = "Time (T)",
             color = "Weight") +
        theme_minimal() + 
        guides(color = guide_colorbar(barwidth = 1.5, barheight = 25)))
ggsave("Figures/time_viscosity_scatter.png", plot = last_plot(),
       width = 8, height = 6, dpi = 300, bg = "white")


print(ggplot(stormer, aes(x = Wt, y = Time, color = Viscosity)) +
  geom_point(size = 3) +
  labs(title = "Stormer Viscometer Data",
       x = "Weight (w)",
       y = "Time (T)",
       color = "Viscosity") +
  theme_minimal() + 
    guides(color = guide_colorbar(barwidth = 1.5, barheight = 25)))
ggsave("Figures/time_weight_scatter.png", plot = last_plot(), 
       width = 8, height = 6, dpi = 300, bg = "white")


# # Fit the nonlinear regression model
# nls_model <- nls(Time ~ (theta1 * Viscosity) / (Wt - theta2),
#                  data = stormer,
#                  start = list(theta1 = 50, theta2 = 5))  # Adjust initial values if needed
# 
# # Print estimated parameters
# print(summary(nls_model))
# 
# theta_estimates <- coef(nls_model)
# theta1_hat = theta_estimates["theta1"]
# theta2_hat = theta_estimates["theta2"]
# 
# sigma_hat <- summary(nls_model)$sigma
# se_theta1 <- summary(nls_model)$parameters["theta1", "Std. Error"]
# se_theta2 <- summary(nls_model)$parameters["theta2", "Std. Error"]
# 
# print(paste("Theta1 estimate:", theta1_hat))
# print(paste("Theta2 estimate:", theta2_hat))
# print(paste("Sigma estimate:", sigma_hat))
# print(paste("Sigma-sq estimate:", sigma_hat**2))



# Fit linear regression model to estimate theta1, theta2
linear_model <- lm(wT ~ Viscosity + Time, data = stormer)
print(summary(linear_model))

# Extract linear estimates as initial values
theta1_init <- coef(linear_model)["Viscosity"]
theta2_init <- coef(linear_model)["Time"]

print(paste("Theta1 init:", theta1_init))
print(paste("Theta2 init:", theta2_init))



# Fit nonlinear regression model using the linear estimates
nls_model <- nls(Time ~ (theta1 * Viscosity) / (Wt - theta2),
                 data = stormer,
                 start = list(theta1 = theta1_init, theta2 = theta2_init))  

# Print estimated parameters from NLS
print(summary(nls_model))

theta_estimates <- coef(nls_model)
theta1_hat = theta_estimates["theta1"]
theta2_hat = theta_estimates["theta2"]

sigma_hat <- summary(nls_model)$sigma
se_theta1 <- summary(nls_model)$parameters["theta1", "Std. Error"]
se_theta2 <- summary(nls_model)$parameters["theta2", "Std. Error"]

print(paste("Theta1 estimate:", theta1_hat))
print(paste("Theta2 estimate:", theta2_hat))
print(paste("Sigma estimate:", sigma_hat))
print(paste("Sigma-sq estimate:", sigma_hat**2))
print(paste("SE(theta1):", se_theta1))
print(paste("SE(theta2):", se_theta2))


print(paste("theta1 misestimation error :", 100*abs(theta1_hat-theta1_init)/theta1_init, "%"))
print(paste("theta2 misestimation error :", 100*abs(theta2_hat-theta2_init)/theta2_init, "%"))
            

# Compute wT (Weight * Time)
stormer$wT <- stormer$Wt * stormer$Time

linear_model <- lm(wT ~ Viscosity + Time, data = stormer)
print(summary(linear_model))


plot <- ggplot(stormer, aes(x = Viscosity, y = wT)) +
  geom_point(color = "blue", size = 3) +  # Scatter points
  geom_smooth(method = "lm", formula = y ~ x, color = "red") +  # Linear fit
  labs(title = "Scatter Plot of w•T vs. Viscosity with Fitted Model",
       x = TeX("$v$"),  
       y = TeX("$w \\cdot T$")) +
  theme_minimal()
ggsave("Figures/Wt_viscosity_scatter.png", plot = last_plot(), 
       width = 8, height = 6, dpi = 300, bg = "white")

print(plot)

diagnostics_A <- plot_residual_diagnostics(linear_model, "A")
print(diagnostics_A)

diagnostics_A2 <- plot_residual_diagnostics(nls_model, "A2")
print(diagnostics_A2)

cat("\n \n \n")


# Exercise B
print("Exercise B")

nls_model_restricted <- nls(wT ~ 25 * (1 + Viscosity / theta2), 
                            data = stormer, 
                            start = list(theta2 = 5))

print(waldtest(nls_model, nls_model_restricted))


diagnostics_B <- plot_residual_diagnostics(nls_model_restricted, "B")
print(diagnostics_B)




# Exercise C
print("Exercise C")

# Compute 92% confidence intervals
z_value <- qnorm(0.96)  # 92% confidence level

theta1_CI <- c(
  theta1_hat - z_value * se_theta1,
  theta1_hat + z_value * se_theta1
)

theta2_CI <- c(
  theta2_hat - z_value * se_theta2,
  theta2_hat + z_value * se_theta2
)

# Print results
print(paste("92% CI for Theta1:", round(theta1_CI[1], 4), "to", round(theta1_CI[2], 4)))
print(paste("92% CI for Theta2:", round(theta2_CI[1], 4), "to", round(theta2_CI[2], 4)))

cat("\n \n \n")



# Exercise D
print("Exercise D")

# 
# model <- lm(Time ~ Wt + Viscosity, data = stormer)
# print(summary(model))


viscosity_grid <- seq(10, 300, length.out = 100)

new_data <- data.frame(Wt = 50, Viscosity = viscosity_grid)

predicted_values <- predict(nls_model, newdata = new_data, se.fit = TRUE)

# Compute 94% confidence intervals
z_crit <- qnorm(0.97)
pred_df <- data.frame(
  Viscosity = viscosity_grid,
  Expected_Time = predicted_values,
  Lower_CI = predicted_values - z_crit * SE_T,
  Upper_CI = predicted_values + z_crit * SE_T
)

# Compute confidence interval width
pred_df$CI_Width <- pred_df$Upper_CI - pred_df$Lower_CI


print(ggplot(pred_df, aes(x = Viscosity, y = Expected_Time)) +
  geom_line(color = "blue", size = 1) +  # Expected Time line
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.3, fill = "blue") +  # Confidence band
  labs(title = "94% Confidence Intervals for Expected Time",
       x = "Viscosity",
       y = "Expected Time") +
  theme_minimal())
ggsave("Figures/expected_T_confidence.png", width = 6, height = 5, dpi = 300, bg = "white")


print(ggplot(pred_df, aes(x = Viscosity, y = CI_Width)) +
  geom_smooth(method = "loess", se = FALSE, color = "red") +  
  labs(x = "Viscosity", y = "Confidence Interval Width", 
       title = "Confidence Interval Width vs. Viscosity") +
  theme_minimal())
ggsave("Figures/viscosity_vs_CI.png", width = 6, height = 5, dpi = 300, bg = "white")


stormer$residuals <- resid(model)
stormer$fitted <- fitted(model)

diagnostics_D <- plot_residual_diagnostics(model, "D")
print(diagnostics_D)

cat("\n \n \n")



# Exercise E
print("Exercise E")


print(summary(nls_model_restricted))

diagnostics_E <- plot_residual_diagnostics(nls_model_restricted, "E")
print(diagnostics_E)


# Compute 95% confidence intervals
z_value <- qnorm(0.975)

theta1_CI <- c(
  theta1_hat - z_value * se_theta1,
  theta1_hat + z_value * se_theta1
)

theta2_CI <- c(
  theta2_hat - z_value * se_theta2,
  theta2_hat + z_value * se_theta2
)

# Print results
print(paste("95% CI for Theta1:", round(theta1_CI[1], 4), "to", round(theta1_CI[2], 4)))
print(paste("95% CI for Theta2:", round(theta2_CI[1], 4), "to", round(theta2_CI[2], 4)))

cat("\n \n \n")



# Extract unique weights (w) from the stormer dataset
unique_weights <- unique(stormer$Wt)

# Prepare a new dataframe for plotting
pred_data <- data.frame(Viscosity = rep(stormer$Viscosity, 3),  # Repeat viscosity values
                        Wt = rep(unique_weights, each = nrow(stormer)))  # Replicate weights for each viscosity value

# Compute the predicted T values using the formula T = theta1 * v / (w - theta2)
pred_data$Predicted_T <- (theta1_hat * pred_data$Viscosity) / (pred_data$Wt - theta2_hat)

# Add the original data points to the plot
print(ggplot(stormer, aes(x = Viscosity, y = Time, color = as.factor(Wt))) +
  geom_point(size = 3) +  # Scatter plot of original data
  geom_line(data = pred_data, aes(x = Viscosity, y = Predicted_T, color = as.factor(Wt)), size = 1, linetype = "dashed") +  # Fitted lines
  labs(title = "Viscosity vs Time with Fitted Curves for Different Weights",
       x = TeX("$v$"),
       y = TeX("$T$"),
       color = "Weight (w)") +
  scale_color_manual(values = c("red", "blue", "green")) +  # Customize colors for each weight
  theme_minimal() +
  theme(legend.position = "bottom"))  # Position the legend at the bottom
ggsave("Figures/Time_viscosity_scatter_fitted.png", plot = last_plot(), 
       width = 8, height = 6, dpi = 300, bg = "white")


# Function to calculate R-squared
calculate_r_squared <- function(actual, predicted) {
  ss_residuals <- sum((actual - predicted)^2)
  ss_total <- sum((actual - mean(actual))^2)
  r_squared <- 1 - (ss_residuals / ss_total)
  return(r_squared)
}

# Calculate R-squared for each unique weight
for (weight in unique_weights) {
  # Filter the data for the current weight
  subset_data <- stormer[stormer$Wt == weight, ]
  
  # Compute the predicted values based on the model
  predicted_values <- theta1_hat * subset_data$Viscosity / (weight - theta2_hat)
  
  # Calculate R-squared
  r_squared <- calculate_r_squared(subset_data$Time, predicted_values)
  
  print(paste("R-squared for weight", weight, ": ", round(r_squared, 3)))
}
