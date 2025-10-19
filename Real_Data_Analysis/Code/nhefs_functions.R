
# remove(list = ls())
# # Load required libraries
library(dplyr)
# # Load the dataset
# nhefs <- read.csv("posterior_for_mu/Real_Data/NHEFS/nhefs.csv")
# # Basic data exploration
# str(nhefs)
# summary(nhefs)
# Function to properly identify and convert factor variables
prepare_nhefs_data <- function(data) {
  
  # Identify categorical variables that should be factors
  factor_vars <- c(
    # Demographics
    "sex",           # 0=male, 1=female
    "race",          # race categories
    "education",     # education levels  
    "marital",       # marital status categories
    "income",        # income brackets
    "birthplace",    # birthplace codes
    
    # Health conditions (binary)
    "asthma", "bronch", "tb", "hf", "hbp", "pepticulcer", 
    "colitis", "hepatitis", "chroniccough", "hayfever", 
    "diabetes", "polio", "tumor", "nervousbreak",
    
    # Lifestyle factors
    "alcoholtype",   # type of alcohol
    "alcoholfreq",   # frequency categories
    "exercise",      # exercise levels
    "active",        # activity levels
    "birthcontrol",  # birth control use
    
    # Symptoms/conditions
    "pica", "headache", "otherpain", "weakheart", 
    "allergies", "nerves", "lackpep", "hbpmed", 
    "boweltrouble", "wtloss", "infection"
  )
  
  # Convert specified variables to factors
  data[factor_vars] <- lapply(data[factor_vars], function(x) {
    if(!is.null(x)) factor(x) else x
  })
  
  # Keep qsmk as binary (0/1) - DO NOT convert to factor
  data$qsmk <- as.numeric(data$qsmk)
  
  return(data)
}

# Create model matrices
# Corrected function that explicitly excludes seqn from covariates
create_model_matrices <- function(data, covariate_set, treatment_var, outcome_var) {
  
  # Ensure seqn is not included in covariate_set
  if("seqn" %in% covariate_set) {
    warning("'seqn' detected in covariate_set. Removing it as it's an ID variable.")
    covariate_set <- covariate_set[covariate_set != "seqn"]
  }
  
  # Select analysis variables - seqn kept only for tracking
  analysis_vars <- c("seqn", treatment_var, outcome_var, covariate_set)
  analysis_data <- data[analysis_vars]
  
  # Remove rows with missing data for outcome and treatment
  analysis_data <- analysis_data[!is.na(analysis_data[[treatment_var]]) & 
                                   !is.na(analysis_data[[outcome_var]]), ]
  
  # Further remove rows with missing covariates
  analysis_data <- analysis_data[complete.cases(analysis_data), ]
  
  # Explicitly separate variables
  seqn <- analysis_data$seqn                    # ID variable (kept separate)
  A <- analysis_data[[treatment_var]]           # Treatment variable
  Y <- analysis_data[[outcome_var]]             # Outcome variable
  X_data <- analysis_data[covariate_set]        # ONLY covariates (no seqn, no treatment, no outcome)
  
  # Verify that seqn is NOT in X_data
  if("seqn" %in% names(X_data)) {
    stop("Error: seqn incorrectly included in covariate matrix!")
  }
  
  # Create design matrix with proper factor handling (seqn excluded)
  X_matrix <- model.matrix(~ ., data = X_data)
  X_matrix_no_intercept <- model.matrix(~ . - 1, data = X_data)
  
  # Get variable information for interpretation
  var_info <- data.frame(
    variable = colnames(X_matrix),
    type = sapply(colnames(X_matrix), function(x) {
      if(grepl("\\(Intercept\\)", x)) return("intercept")
      # Extract base variable name
      base_var <- gsub("\\d+$", "", x)
      base_var <- gsub("[^A-Za-z_]", "", base_var)
      if(base_var %in% names(X_data)) {
        if(is.factor(X_data[[base_var]])) return("factor_level")
        else return("continuous")
      }
      return("derived")
    })
  )
  
  # Print summary to verify correct structure
  cat("Data structure verification:\n")
  cat("- ID variable (seqn): length =", length(seqn), "\n")
  cat("- Treatment variable (", treatment_var, "): length =", length(A), "\n")
  cat("- Outcome variable (", outcome_var, "): length =", length(Y), "\n")
  cat("- Covariate matrix: dimensions =", dim(X_matrix), "\n")
  cat("- Covariates included:", paste(names(X_data), collapse = ", "), "\n")
  cat("- seqn in covariates?", "seqn" %in% names(X_data), "\n\n")
  
  return(list(
    X = X_matrix,
    X_no_int = X_matrix_no_intercept, 
    A = A,
    Y = Y,
    seqn = seqn,
    covariate_data = X_data,           # Raw covariate data (no seqn)
    analysis_data = analysis_data,     # Full analysis data (includes seqn for reference)
    var_info = var_info,
    n_obs = nrow(analysis_data),
    factor_levels = lapply(X_data[sapply(X_data, is.factor)], levels),
    covariates_used = names(X_data)    # Explicit list of covariates in model
  ))
}

# Simple function for quick positivity check
quick_positivity_check <- function(data, treatment_var = "qsmk") {
  
  # Get factor variables
  factor_vars <- names(data)[sapply(data, is.factor)]
  factor_vars <- factor_vars[factor_vars != treatment_var]
  
  issues <- list()
  
  for(var in factor_vars) {
    # Cross-tabulation
    tab <- table(data[[var]], data[[treatment_var]], useNA = "ifany")
    
    # Check each level
    levels_with_no_treated <- which(tab[, "1"] == 0)
    levels_with_no_control <- which(tab[, "0"] == 0)
    
    if(length(levels_with_no_treated) > 0 | length(levels_with_no_control) > 0) {
      issues[[var]] <- list(
        no_treated = names(levels_with_no_treated),
        no_control = names(levels_with_no_control),
        table = tab
      )
    }
  }
  
  if(length(issues) == 0) {
    cat("✓ All factor levels have both treated and control units!\n")
  } else {
    cat("⚠ Positivity violations found in:\n")
    for(var in names(issues)) {
      cat("  -", var, "\n")
      if(length(issues[[var]]$no_treated) > 0) {
        cat("    Levels with no treated units:", paste(issues[[var]]$no_treated, collapse = ", "), "\n")
      }
      if(length(issues[[var]]$no_control) > 0) {
        cat("    Levels with no control units:", paste(issues[[var]]$no_control, collapse = ", "), "\n")
      }
    }
  }
  
  return(issues)
}


# # Apply proper factor conversion first
# nhefs_clean <- prepare_nhefs_data(nhefs)
# # Check the structure
# str(nhefs_clean[c("sex", "race", "education", "asthma", "exercise")])
# #colnames of cleaned data
# colnames(nhefs_clean)
# 
# 
# # Create different analysis datasets
# # 1. Minimal set
# minimal_cov = c("sex", "age", "exercise", "education", "smokeintensity", "smokeyrs", "wt71")
# # Create model matrices for minimal dataset
# minimal_data <- create_model_matrices(nhefs_clean, minimal_cov)
# head(minimal_data$X_no_int)
# colSums(is.na(cbind(minimal_data$X_no_int, minimal_data$Y, minimal_data$Y)))
# 
# 
# 
# # 2. Extended set
# extended_data <- create_analysis_data(
#   nhefs,
#   c("sex", "age", "race", "education", "income", "marital",
#     "smokeintensity", "smokeyrs", "wt71", "ht", "sbp", "dbp",
#     "asthma", "bronch", "hf", "hbp", "diabetes",
#     "alcoholpy", "alcoholfreq", "exercise", "active")
# )
# 
# head(extended_data)
# colSums(is.na(extended_data))
# nrow(extended_data)
# ncol(extended_data)
# # Create model matrices for extended dataset
# model_data <- create_model_matrices(extended_data)
# 

