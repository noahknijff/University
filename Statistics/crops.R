crops = read.table("data/crops.txt", header=TRUE)
print(colnames(crops))

# Convert categorical variables to factors
crops$County <- as.factor(crops$County)
crops$Related <- as.factor(crops$Related)

# Fit the Two-Way ANOVA model (without Size)
anova_model <- aov(Crops ~ County * Related, data = crops)

# Display the ANOVA table
print(summary(anova_model))