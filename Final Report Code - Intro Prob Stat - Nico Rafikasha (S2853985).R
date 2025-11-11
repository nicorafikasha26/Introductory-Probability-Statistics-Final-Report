library(vcd) # Goodfit
library(dgof) #ks.test
library(fitdistrplus) #MLE
library(ggplot2)
library(dplyr)
library(EnvStats)
library(caret)
library(nnet)
library(pROC)
library(GGally)
library(VGAM)
library(readxl)

covariates <- read_excel("covariates.xlsx")
head(covariates)
hist(covariates$`Vas-12months`)

biomarkers <- read_excel("biomarkers.xlsx")
head(biomarkers)
hist(biomarkers$`VEGF-A`)

#First, set Variable factor for gender and smoker's status
covariates$`Sex (1=male, 2=female)`=as.factor(covariates$`Sex (1=male, 2=female)`)
covariates$`Smoker (1=yes, 2=no)`=as.factor(covariates$`Smoker (1=yes, 2=no)`)

#Statistic Descriptive
summary(covariates)

#Re-level Data Factor
covariates = covariates %>% mutate(`Sex (1=male, 2=female)` = relevel(`Sex (1=male, 2=female)`, ref = "1"))
covariates = covariates %>% mutate(`Smoker (1=yes, 2=no)` = relevel(`Smoker (1=yes, 2=no)`, ref = "2"))

#Statistic Descriptive
summary(biomarkers)

#Next, for biomarkers data, we want to separate "biomarker" column to patient id and incubation observation
library(tidyr)
biomarkers <- biomarkers %>% separate (`Biomarker`, into = c("PatientID", "Observation"), sep = "-", remove = FALSE)
biomarkers <- biomarkers %>% select(-`Biomarker`)
biomarkers

#Then, we want each observation under patient ID to be assigned manually to Patient ID's info in covariates
#First, make sure the columns name for two separate dataset are same for PatientID
colnames(covariates)[colnames(covariates) == "PatientID"] <- "PatientID"

#Pivot -> change the shape of biomarkers' table
value_cols <- setdiff(colnames(biomarkers), c("PatientID", "Observation"))
biomarkers_wide <- biomarkers %>% pivot_wider(names_from = Observation, values_from = all_of(value_cols)) %>% arrange(as.numeric(`PatientID`))

biomarkers_wide <- biomarkers_wide %>% mutate(`PatientID` = as.numeric(`PatientID`))

covariates_final <- covariates %>% left_join(biomarkers_wide, by="PatientID")
covariates_final


# Select baseline biomarkers' column
baseline_cols <- grep("0weeks", names(covariates_final), value = TRUE)

# Prepare result's dataframe
results_ttest <- data.frame(
  Protein = baseline_cols,
  Test_Used = "Welch's t-test",
  p_value = NA,
  Adj_p_value = NA
)

# Loop for each biomarker
for (i in seq_along(baseline_cols)) {
  prot <- baseline_cols[i]
  
  male_vals <- covariates_final[covariates_final$`Sex (1=male, 2=female)` == 1, prot, drop = TRUE]
  female_vals <- covariates_final[covariates_final$`Sex (1=male, 2=female)` == 2, prot, drop = TRUE]
  
  male_vals <- male_vals[!is.na(male_vals)]
  female_vals <- female_vals[!is.na(female_vals)]
  
  # Two-side Hypothesis Testing (Welch t-test, unequal variances)
  test <- t.test(male_vals, female_vals, var.equal = FALSE)
  results_ttest$p_value[i] <- test$p.value
}

# Bonferroni correction
results_ttest$Adj_p_value <- p.adjust(results_ttest$p_value, method = "bonferroni")

# Compute Type I Error Probability
alpha <- 0.05
m <- nrow(results_ttest)
type1_prob_uncorrected <- 1 - (1 - alpha)^m
type1_prob_bonferroni <- 1 - (1 - alpha/m)^m

cat("1. Probability(at least one Type I error, uncorrected):",
    round(type1_prob_uncorrected, 4), "\n")
cat("2. After Bonferroni correction:",
    round(type1_prob_bonferroni, 4), "\n\n")

print(results_ttest)


# Creating a Simple Linear Regression Model
y_var <- "Vas-12months"

#Select variable at inclusion only
covariates_final1 <- covariates_final %>% select(-matches("6weeks|12months"))

#Create automatic formula: y ~ all column except y and PatientID
predictors <- setdiff(names(covariates_final1), c("PatientID",y_var))
predictors_bt <- paste0("`", predictors, "`")
formula_auto <- as.formula(paste0("`", y_var, "` ~", paste(predictors_bt, collapse = " + ")))
formula_auto


# Run Regression Model
n <- nrow(covariates_final)
train_index <- sample(1:n, size = 0.8 * n, replace = FALSE)

train_data <- covariates_final[train_index, ]
test_data <- covariates_final[-train_index, ]

model1a = lm(formula_auto, data = train_data)
summary(model1a)


#Predicting
# Predict Testing Data

test_data$predicted <- predict(model1a, newdata = test_data)

# Evaluate Prediction Accuracy

correlation <- cor(test_data$`Vas-12months`, test_data$predicted)
mse <- mean((test_data$`Vas-12months` - test_data$predicted)^2)

cat("Correlation between actual and predicted:", round(correlation, 3), "\n")
cat("Mean Squared Error (MSE):", round(mse, 3), "\n")

# Plot predicted vs actual value

plot(test_data$`Vas-12months`, test_data$predicted,
     xlab = "Actual VAS-12months", ylab = "Predicted", main = "Simple Linear Regression Prediction",
     pch = 19, col = "blue")
abline(0, 1, col = "red", lwd = 2)



