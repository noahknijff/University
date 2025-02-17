library(ggplot2)  
library(reshape2) 
library(gridExtra)

# Close the current plot device if one is open
if (dev.cur() > 1) {
  dev.off()
}

cholesterol = read.table("data/cholesterol.txt", header=TRUE)
before_data <- cholesterol$Before
after_data <- cholesterol$After8weeks


# Exercise A

# Convert data to long format for ggplot
cholesterol_long <- reshape2::melt(cholesterol, id.vars = NULL, variable.name = "Time", value.name = "Value")

print(ggplot(cholesterol_long, aes(x = Value, fill = Time)) +
  geom_histogram(alpha = 0.6, bins = 10, color = "black") +
  labs(title = "Cholesterol Levels Before and After 8 Weeks",
       x = "Cholesterol Level",
       y = "Count") +
  theme_minimal() +
  scale_fill_manual(values = c("blue", "red")) +
  facet_wrap(~ Time, ncol = 2) +  # Arrange plots in 2 columns
  theme(panel.spacing = unit(2, "lines")))  # Adjust spacing between panels
ggsave("Figures/cholesterol_histograms.png", plot = last_plot(),
       width = 10, height = 6, dpi = 300, bg = "white") 


# Check for normality

# QQ plot for "Before"
qq_before <- ggplot(subset(cholesterol_long, Time == "Before"), aes(sample = Value)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "QQ Plot for Cholesterol Levels Before")

# QQ plot for "After8weeks"
qq_after <- ggplot(subset(cholesterol_long, Time == "After8weeks"), aes(sample = Value)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "QQ Plot for Cholesterol Levels After 8 Weeks")

combined_qq_plot <- grid.arrange(qq_before, qq_after, ncol = 2)
ggsave("Figures/qq_plots_combined.png", plot = combined_qq_plot, 
       width = 12, height = 6, dpi = 300, bg = "white")


# Run Pearson correlation test
cor_test <- cor.test(before_data, after_data, method = "pearson")
print(cor_test)



# Exercise B

# Two paired t-test
t_test_result <- t.test(before_data, after_data, paired = TRUE)
print(t_test_result)

# Wilcoxon Signed-Rank Test
wilcox_test_result <- wilcox.test(before_data, after_data, paired = TRUE)
print(wilcox_test_result)

# Permutation Test
observed_diff <- mean(after_data - before_data)

# Combine the two groups and permute
combined_data <- c(before_data, after_data)
n <- length(before_data)

perm_diffs <- replicate(10000, {
  perm <- sample(combined_data)
  mean(perm[1:n] - perm[(n+1):(2*n)])
})

# Compute p-value
p_value <- mean(abs(perm_diffs) >= abs(observed_diff))
print(sprintf("Permutation test, p-value: %f", p_value))




# Exercise C

data <- after_data  # or after_data if you're analyzing that one
n <- length(data)
mean_data <- mean(data)
sd_data <- sd(data)
alpha <- 0.03  # For 97% CI
df <- n - 1

# t-distribution critical value
t_value <- qt(1 - alpha / 2, df)

# Confidence Interval
CI_parametric <- c(mean_data - t_value * sd_data / sqrt(n), mean_data + t_value * sd_data / sqrt(n))

# 2. Bootstrap CI for mu
bootstrap_samples <- 10000
bootstrap_means <- numeric(bootstrap_samples)

# Resampling with replacement
for (i in 1:bootstrap_samples) {
  bootstrap_sample <- sample(data, size = n, replace = TRUE)
  bootstrap_means[i] <- mean(bootstrap_sample)
}

# Calculating the 97% CI from the bootstrap distribution
CI_bootstrap <- quantile(bootstrap_means, c(0.015, 0.985))

# Output the results
print(paste("Parametric CI for mu (97%):", round(CI_parametric, 2)))
print(paste("Bootstrap CI for mu (97%):", round(CI_bootstrap, 2)))



# Exercise D


n <- length(after_data)  # sample size
theta_values <- seq(3, 12, by = 0.1)  # Range of theta to test
num_bootstrap <- 10000  # Number of bootstrap samples
alpha <- 0.05  # Significance level

# Observed test statistic (maximum value)
T_obs <- max(after_data)

# Function to perform bootstrap test
bootstrap_test <- function(data, theta, num_bootstrap) {
  bootstrap_max <- numeric(num_bootstrap)
  
  # Generate bootstrap samples and calculate maximum value for each
  for (i in 1:num_bootstrap) {
    bootstrap_sample <- runif(length(data), min = 3, max = theta)
    bootstrap_max[i] <- max(bootstrap_sample)
  }
  
  # p-value calculation: proportion of bootstrap max values >= observed max
  p_value <- mean(bootstrap_max >= T_obs)
  
  return(p_value)
}


# Function to perform KS test
ks_test <- function(data, theta) {
  # Generate a uniform distribution from 3 to theta
  uniform_data <- runif(length(data), min = 3, max = theta)
  
  # Perform the Kolmogorov-Smirnov test
  ks_result <- ks.test(data, uniform_data)
  
  # Return the p-value of the KS test
  return(ks_result$p.value)
}



# Run the bootstrap test and KS test for each theta
results <- data.frame(theta = theta_values, bootstrap_p_value = NA, ks_p_value = NA)

for (i in 1:length(theta_values)) {
  # Apply bootstrap test
  bootstrap_p_value <- bootstrap_test(after_data, theta_values[i], num_bootstrap)
  
  # Apply KS test
  ks_p_value <- ks_test(after_data, theta_values[i])
  
  # Store p-values
  results$bootstrap_p_value[i] <- bootstrap_p_value
  results$ks_p_value[i] <- ks_p_value
}

# Find the first non-rejected theta for both tests (where p-value > 0.05)
first_non_rejected_bootstrap <- results[which(results$bootstrap_p_value > 0.05)[1], ]
first_non_rejected_ks <- results[which(results$ks_p_value > 0.05)[1], ]

# Print the first non-rejected theta and p-values
print("First non-rejected theta and p-value (Bootstrap test):")
print(first_non_rejected_bootstrap)

print("First non-rejected theta and p-value (KS test):")
print(first_non_rejected_ks)

# Plot the p-values against theta values for both tests
library(ggplot2)

# Bootstrap test plot
p1 <- ggplot(results, aes(x = theta, y = bootstrap_p_value)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = expression(paste("Bootstrap test p-values vs. ", theta)),
       x = expression(theta),
       y = "p-value") +
  theme_minimal()

# KS test plot
p2 <- ggplot(results, aes(x = theta, y = ks_p_value)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = expression(paste("KS test p-values vs. ", theta)),
       x = expression(theta),
       y = "p-value") +
  theme_minimal()

# Save the plots as images
ggsave("Figures/bootstrap_pvalues.png", plot = p1, dpi = 300, bg = "white")
ggsave("Figures/ks_pvalues.png", plot = p2, dpi = 300, bg = "white")

# Print the plots
print(p1)
print(p2)



# Exercise E
# Perform Wilcoxon signed-rank test
wilcox_test_result <- wilcox.test(after_data, mu = 6, alternative = "less")
print(wilcox_test_result)


# Perform binomial test
count_below_4_5 <- sum(after_data < 4.5)
binom_test_result <- binom.test(count_below_4_5, length(after_data), p = 0.25, alternative = "greater")
print(binom_test_result)
