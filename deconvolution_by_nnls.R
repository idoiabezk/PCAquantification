---
  title: "R Notebook"
output: html_notebook
editor_options: 
  chunk_output_type: console
---

    
library(nnls)
library(tidyverse)
library(tidymodels)
library(readxl)
library(openxlsx)



## Setup automatic nnls deconvolution process from excel sheet


#### TODO ####
#filter(rsquared > 0.7) #cannot do this as some are still false positive and need to replace formula later with 0 for nnls 

# NEED TO ADD using case_when or elseif..
# if (sum(!is.na(filtered_dataA$Area)) == 0 || sum(!is.na(filtered_dataA$`Analyte Concentration`)) == 0) {
#       cat("No valid cases for fitting the model for molecule:", molecule_name, "\n")
# OR add a filter to remove those rsquared that are very low (perhaps below 0.7?)  

##############

# reading data from excel
Skyline_output <- read_excel("F:/LINKOPING/CP analysis/Samples_From_Orebro/Skyline/ResultsFromSkiline.xlsx") |>
  mutate(`Analyte Concentration` = as.numeric(`Analyte Concentration`)) |> 
  mutate(Area = as.numeric(Area)) |> 
  mutate(Area = replace_na(Area, 0)) |> # Replace missing values in the Response_factor column with 0
  mutate(Area = replace_na(Area, 0)) |>
  mutate(RatioQuanToQual = as.numeric(RatioQuanToQual)) |> #--< convert to numeric #-->
  mutate(RatioQualToQuan = as.numeric(RatioQualToQuan)) #--< convert to numeric #-->



#This is currently filtered for C10-C13 only to compare with the Perkons script. Will remove later or add as an argument in function
CPs_standards <- Skyline_output |> 
  filter(`Sample Type` == "Standard",
         Molecule != "IS",
         Molecule != "RS",
         `Isotope Label Type` == "Quan",
         Note != "NA") |>
  group_by(Note, Molecule) |>
  mutate(rel_int = Area/sum(Area)) |> #why ius it needed? Maybe can be removed
  nest() |> 
  # Remove groups where all 'Analyte Concentration' values are NA or there are no non-NA cases
  filter(map_lgl(data, ~sum(!is.na(.x$`Analyte Concentration`)) > 0)) |> 
  mutate(models = map(data, ~lm(Area ~ `Analyte Concentration`, data = .x))) |> 
  mutate(coef = map(models, coef)) |> 
  mutate(Response_factor = map(coef, pluck("`Analyte Concentration`"))) |> 
  mutate(intercept = map(coef, pluck("(Intercept)"))) |> 
  mutate(rsquared = map(models, summary)) |> 
  mutate(rsquared = map(rsquared, pluck("r.squared"))) |>
  select(-coef) |>  # remove coef variable since it has already been plucked
  unnest(c(Response_factor, intercept, rsquared)) |>  #removing the list type for these variables
  mutate(Response_factor = if_else(Response_factor < 0, 0, Response_factor)) |> # replace negative RF with 0
  mutate(rsquared = ifelse(is.nan(rsquared), 0, rsquared)) |> 
  mutate(Chain_length = paste0("C", str_extract(Molecule, "(?<=C)[^H]+"))) |> 
  #filter(Chain_length == "C18" | Chain_length == "C19" | Chain_length == "C20" | Chain_length == "C21" | Chain_length == "C22" | Chain_length == "C23" | Chain_length == "C24"| Chain_length == "C25"|Chain_length == "C26" | Chain_length == "C27"| Chain_length == "C28"| Chain_length == "C29" | Chain_length == "C30")|> #this will be remove later or added as arg in fn
  #filter(Note == "G" | Note == "H"| Note == "I" | Note == "J") |> #this will be remove later or added as arg in fn
  ungroup() |> 
  group_by(Note, Chain_length) |> 
  mutate(Sum_response_factor_chainlength = sum(Response_factor, na.rm = TRUE)) |> 
  ungroup()

#For SCCPs
CPs_standardsS<-CPs_standards |> 
filter(str_detect(Note, "S-")) |> 
mutate(Response_factor = if_else(Chain_length %in% c("C14", "C15", "C16", "C17", "C18", "C19", "C20", "C21", "C22", "C23", "C24", "C25", "C26", "C27", "C28", "C29", "C30"), 0, Response_factor))  
#For MCCPs
CPs_standardsM<-CPs_standards |> 
  filter(str_detect(Note, "M-")) |> 
  mutate(Response_factor = if_else(Chain_length %in% c("C10", "C11", "C12", "C13", "C18", "C19", "C20", "C21", "C22", "C23", "C24", "C25", "C26", "C27", "C28", "C29", "C30"), 0, Response_factor))  
#For MCCPs
CPs_standardsL<-CPs_standards |> 
  filter(str_detect(Note, "L-")) |> 
  mutate(Response_factor = if_else(Chain_length %in% c("C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17"), 0, Response_factor))  
#Together
CPs_standards<- rbind(CPs_standardsS, CPs_standardsM, CPs_standardsL)


CPs_samples <- Skyline_output |> 
  filter(`Sample Type` %in% c("Unknown", "Blank"),
         Molecule != "IS",
         Molecule != "RS",
         `Isotope Label Type` == "Quan") |> 
  mutate(Chain_length = paste0("C", str_extract(Molecule, "(?<=C)[^H]+"))) |> 
  #filter(Chain_length == "C18" | Chain_length == "C19" | Chain_length == "C20" | Chain_length == "C21" | Chain_length == "C22" | Chain_length == "C23" | Chain_length == "C24"| Chain_length == "C25"|Chain_length == "C26" | Chain_length == "C27"| Chain_length == "C28"| Chain_length == "C29" | Chain_length == "C30") |> #this will be remove later or added as arg in fn
  group_by(`Replicate Name`) |> 
  mutate(Relative_distribution = Area / sum(Area, na.rm = TRUE)) |> 
  select(Molecule, Area, Relative_distribution) |> 
  mutate(across(Relative_distribution, ~replace(., is.nan(.), 0))) |>    #replace NaN with zero
  #nest() |> 
  ungroup() 

CPs_standards_input <- CPs_standards |> 
  select(Molecule, Note, Response_factor) |> 
  pivot_wider(names_from = "Note", values_from = "Response_factor")

CPs_samples_input <- CPs_samples |> 
  select(Molecule, `Replicate Name`, Relative_distribution) |> 
  pivot_wider(names_from = "Replicate Name", values_from = "Relative_distribution")

# This step ensures that all values are corresponding to the same molecule for std and sample        
combined <- CPs_samples_input |> 
  right_join(CPs_standards_input, by = "Molecule")


############################################################################### DECONVOLUTION #############################################################################

# Ensure combined_matrix is correctly defined as a matrix prior to the deconvolution
combined_matrix <- CPs_standards_input  |> 
  select(-Molecule) |> 
  as.matrix()

# Ensure combined_sample is correctly defined with nested data frames prior to the deconvolution
combined_sample <- CPs_samples  |> 
  group_by(`Replicate Name`) |> 
  select(-Molecule, -Area) |> 
  nest() |> 
  ungroup()

# Function to perform deconvolution on a single data frame: I have many errors so I included more things but not sure if it is ok, or if I made it more complicated
perform_deconvolution <- function(df, combined_matrix) {
  df_matrix <- as.matrix(df)
  
  print(paste("df_matrix dimensions:", dim(df_matrix)))
  print(paste("combined_matrix dimensions:", dim(combined_matrix)))
  
  if (nrow(combined_matrix) != nrow(df_matrix)) {
    stop("Dimensions of combined_matrix and df are incompatible.")
  }
  
  # Reshape df_matrix if it has only one column 
  if (ncol(df_matrix) == 1) { 
    df_vector <- as.vector(df_matrix)
  } else {
    df_vector <- df_matrix
  }
  
  # Perform nnls
  deconv <- nnls(combined_matrix, df_vector)
  
  # Extract deconvolution results
  deconv_coef <- deconv$x
  deconv_resolved <- deconv$fitted.values
  deconv_reconst <- rowSums(combined_matrix %*% deconv_coef)
  
  # Ensure that values are positive for chi-square test
  if (any(deconv_resolved <= 0) || any(df_vector <= 0)) {
    warning("Non-positive values found, skipping chi-square test")
    chisq_result <- NULL
  } else {
    chisq_result <- chisq.test(deconv_resolved, p = df_vector, rescale.p = TRUE)
  }
  
  return(list(
    deconv_coef = deconv_coef,
    deconv_resolved = deconv_resolved,
    deconv_reconst = deconv_reconst,
    chisq_result = chisq_result
  ))
}

# Apply the perform_deconvolution function to each nested data frame
Deconvolution <- combined_sample  |> 
  mutate(result = map(data, ~ perform_deconvolution(.x, combined_matrix)))

# Extract deconv_coef from results and create a new data frame
deconv_coef_df <- Deconvolution  |> 
  mutate(deconv_coef = map(result, "deconv_coef"))  |> 
  select(`Replicate Name`, deconv_coef)  |> 
  unnest_wider(deconv_coef, names_sep = "_")

# View the result
print(deconv_coef_df)


########################################################## Calculate the concentration in ng/uL ###############################################################

#Remove the replicate name to generate vectors:
deconv_coef_df_matrix<- deconv_coef_df |> 
  select(-`Replicate Name`)

# Initialize an empty list to store results
result_list <- list()

# Iterate through each row of deconv_coef_df_matrix
for (i in 1:nrow(deconv_coef_df_matrix)) {
  
  # Extract row vector from deconv_coef_df_matrix
  deconv_coef_vector <- as.numeric(deconv_coef_df_matrix[i, ])
  
  # I had this one before, but to make sure
  combined_matrix <- CPs_standards_input |> 
    select(-Molecule)  |> 
    mutate(across(everything(), as.numeric)) |> 
    as.matrix()
  
  # Perform element-wise multiplication
  result <- sweep(combined_matrix, 2, deconv_coef_vector, `*`)
  
  # Create data frame with column names from CPs_standards_input
  result_df <- as.data.frame(result)
  colnames(result_df) <- colnames(CPs_standards_input)[-which(names(CPs_standards_input) == "Molecule")]
  
  # Assign name to the result_df from deconv_coef_df
  replicate_name <- deconv_coef_df$`Replicate Name`[i]
  
  # Store result in result_list with the corresponding name
  result_list[[replicate_name]] <- result_df
}

# Print the names of each data frame and their first few rows
for (name in names(result_list)) {
  cat("DataFrame Name:", name, "\n")
  print(head(result_list[[name]]))
  cat("\n")
}

# Combine all data frames into a single data frame with 'Replicate Name'
final_df <- do.call(rbind, Map(function(df, name) {
  df$`Replicate Name` <- name
  df <- df[, c("Replicate Name", setdiff(names(df), "Replicate Name"))]
  df
}, result_list, names(result_list)))

# Add CPs_standards_input$Molecule column to final_df
final_df$Molecule <- CPs_standards_input$Molecule

# Print the final combined data frame
print(final_df)


#Organize the data
final_df_tidy<-final_df|> 
  group_by(`Replicate Name`) |> 
  nest() |> 
  ungroup()


#Total sum the values for each replicate

#Remove the molecule
final_df_matrix<- final_df |> 
  select(-Molecule) |> 
  group_by(`Replicate Name`) |> 
  nest()

# Initialize an empty data frame to store results
total_sums_df <- data.frame(
  `Replicate Name` = character(),
  `Total Sum` = numeric(),
  stringsAsFactors = FALSE
)

# Iterate through each row of final_df_grouped
for (i in 1:nrow(final_df_matrix)) {
  # Extract nested data frame
  nested_df <- final_df_matrix$data[[i]]
  
  # Calculate total sum for the current `Replicate Name`
  `Replicate Name` <- final_df_matrix$`Replicate Name`[[i]]
  total_sum <- sum(colSums(nested_df[, -1]))  # Exclude the grouping column
  
  # Append results to total_sums_df
  total_sums_df <- rbind(total_sums_df, data.frame(
    `Replicate Name` = `Replicate Name`,
    `Total Sum` = total_sum
  ))
}

# Print the resulting data frame
print(total_sums_df)


################################################### FINAL RESULTS ####################################################################
CPs_samples<-CPs_samples |> 
rename(`Replicate.Name` = `Replicate Name`)

# Merge total_sums_df into CPs_samples based on Replicate Name
Concentration <- CPs_samples  |> 
  left_join(total_sums_df, by = "Replicate.Name")  |> 
  mutate(Concentration = `Relative_distribution` * `Total.Sum`)

print(Concentration)
Concentration<-Concentration |> 
  group_by(Replicate.Name) |> 
  distinct( `Molecule`, Concentration) |> 
  nest()

# Perform operations to reorganize data
reorganized_data <- Concentration  |> 
  unnest() |>  
  distinct(`Replicate.Name`, `Molecule`, .keep_all = TRUE)  |> 
  pivot_wider(names_from = `Molecule`, values_from = `Concentration`)
reorganized_data <- t(reorganized_data) #transpose

#Make the first row (replicate names) the column names
colnames(reorganized_data) <- reorganized_data[1, ]
Samples_Concentration <- reorganized_data[-1, ]
# Convert the result back to a data frame
Samples_Concentration <- as.data.frame(Samples_Concentration)
Samples_Concentration<- Samples_Concentration |> 
  mutate(Molecule = CPs_samples_input$Molecule)|> 
  relocate(Molecule, .before = everything()) 

######################################################### SAVE RESULTS ###################################################################

# Specify the file path where you want to save the Excel file
excel_file <- "F:/LINKOPING/CP analysis/Samples_From_Orebro/Samples_ConcentrationFromScript2.xlsx"

# Write 'Samples_Concentration' to Excel
write.xlsx(Samples_Concentration, excel_file, rowNames = FALSE)

# Confirm message that the file has been saved
cat("Excel file saved:", excel_file, "\n")

















################################### CHANGE SO IT TAKES ALL THE SAMPLES###################################
#CPs_samples_individual <- CPs_samples |> 
#filter(`Replicate Name` == "NIST_R1") |> 
#select(Molecule, Area, Relative_distribution)
#############################################################################################

#combined_matrix <- CPs_standards_input |> 
#select(-Molecule) |> 
#as.matrix()

#combined_sample <- CPs_samples |> 
# group_by(`Replicate Name`) |> 
#select(-Molecule, -Area) |> 
#nest() |> 
#ungroup() 

#combined_sample_list <- combined_sample |> 
#mutate(data = map(data, as.data.frame)) |> 
#pull(data)

# Using the non-negative least squares, nnls, package to deconvolute the distribution


#perform_deconvolution<- function(combined_matrix, combined_sample){
#deconv <- nnls(combined_matrix, combined_sample)

#deconv

#deconv_coef <- deconv$x

#deconv_resolved <- deconv[["fitted"]]

#deconv_reconst <- rowSums(combined_matrix %*% deconv_coef)

#chisq.test(deconv_resolved, p = combined_sample, rescale.p = TRUE)
#}

#Deconvolution <- combined_sample |> 
#mutate(combined_sample = map(data, ~ perform_deconvolution(combined_matrix, combined_sample)))



#Plot the concentration together with the patterns to see if it is ok
#par(mfrow = c(3,1))
#barplot(combined_sample, main = "Sample")
#barplot(deconv_reconst, main = "Reconstructed")
#barplot(CPs_samples_individual$Concentration, main = "Sample concentration ng/uL")#Plot the concentration


########### CONCENTRATION IN ng/uL

#CPs_standards <- CPs_standards|>
#mutate(Standard_response = CPs_standards$Response_factor * deconv_coef) #Calculate the concentration in the standards

#Standard_response <- sum(CPs_standards$Response_factor) #Sum the concentration of the homologues
#SumResponse <- sum(CPs_samples_individual$Area) #Sum the signal of the homologues in the sample
#SumConcentration <- SumResponse/Standard_response #Calculate the sum concentration in the sample

#CPs_samples_individual <- CPs_samples_individual |> 
# mutate(Concentration = SumConcentration* Relative_distribution) #Calculate the concentration of each homologue in the sample



################################################################################################################



#----JUST SOME TEST SCRIPTS-----
  
  
  ## test using Bogdal paper Figure 1


s <- c(0.1, 0.4, 0.4, 0.1)
y1 <- c(0.1, 0.4, 0.3, 0.2)
y2 <- c(0.5, 0.1, 0.2, 0.2)
y3 <- c(0.2, 0.4, 0.4, 0)
Y <- cbind(y1, y2, y3)
deconv <- nnls(Y, s)
deconv
coefficients <- deconv$x
resolved <- c(y1*deconv$x[1]+y2*deconv$x[2]+y3*deconv$x[3])
resolved2 <- deconv[["fitted"]]

chisq.test(resolved, p = s, rescale.p = TRUE)
qchisq(.5, df = 3)

reconstructed_pattern <- rowSums(Y%*%coefficients)

par(mfrow = c(2,1))
barplot(s)
barplot(reconstructed_pattern)




## Test using data from Gao et al std 57, 61, 63 and sample 74

y1 <- c(0,0,0.00890544,0.021319085,0.019430052,0.023855786,0,0.003022453,0.064173143,0.143404577,0.087543178,0.028875216,0,0.003562176,0.048197323,0.146265112,0.13800734,0.04015544,0,0.007016408,0.022398532,0.074805699,0.096556563,0.022506477)
y2<- c(0,	0.000420628,	0.014890216,0.020947253,0.019432994,0.016488601,0,0.015394969,0.091192059,0.125094641,0.070665433,0.022798015,0,0.014385463,0.083789013,0.143097501,0.103053756,0.028182048,0.001009506,0.012450576,0.041389754,0.081854126,0.07596534,0.017498107)
y3 <- c(0,0,0.009719222,0.020347846,0.018983744,0.024155962,0,0.002671365,0.05632602,0.137944754,0.08821189,0.029726043,0,0.001648289,0.039047403,0.140729794,0.139365693,0.045413209,0.007332045,0,0.017051267,0.071331136,0.106172559,0.043821757)
A <- as.matrix(cbind(y1,y2, y3))

s <- c(0,0.002484031,0.038549451,0.128585146,0.133542771,0.066473928,0,0.00241619,0.024986432,0.101563478,0.147058824,0.072579635,0.000741035,0.005437732,0.01244103,0.027548741,0.041403999,0.023180812,0.003846074,0.011256419,0.031332192,0.043131341,0.047916754,0.033523984)

a <- nnls(A,s)
a
coefficients <- a$x
resolved <- c(y1*a[["x"]][1] + y2*a[["x"]][2] + y3*a[["x"]][3])
chisq.test(resolved, p = s, rescale.p = TRUE)

sresolved <- cbind(s, resolved)
chisq.test(sresolved)




## test tidymodels 

library(RcppML)
library(Matrix)

data(biomass, package = "modeldata")

rec <- recipe(HHV ~ ., data = biomass) %>%
  update_role(sample, new_role = "id var") %>%
  update_role(dataset, new_role = "split variable") %>%
  step_nnmf_sparse(
    all_numeric_predictors(),
    num_comp = 2,
    seed = 473,
    penalty = 0.01
  ) %>%
  prep(training = biomass)

bake(rec, new_data = NULL)


bake(rec, new_data = NULL) %>%
  ggplot(aes(x = NNMF2, y = NNMF1, col = HHV)) +
  geom_point()





glmnet


library(glmnet)
mod2 <- glmnet(x, y, lambda = 0, lower.limits = 0, intercept = FALSE)
coef(mod2)




#Neural network


library(tidymodels)
library(brulee)

# Example data (replace with your own)
# Create synthetic data for illustration
set.seed(123)
n_obs <- 100
n_features <- 5

# Generate random features
features <- matrix(rnorm(n_obs * n_features), nrow = n_obs)

# Create response variable (classification example)
response <- sample(c("A", "B", "C"), size = n_obs, replace = TRUE)

# Combine features and response into data frames
cls_train <- data.frame(features, class = response)
cls_val <- data.frame(features, class = sample(c("A", "B", "C"), size = n_obs, replace = TRUE))

# Check the structure of cls_train and cls_val
str(cls_train)
str(cls_val)

# Create a recipe
biv_rec <- recipe(class ~ ., data = cls_train) %>%
  step_normalize(all_predictors())

# Specify the neural network model
nnet_spec <- mlp(epochs = 1000, hidden_units = 10, penalty = 0.01, learn_rate = 0.1) %>%
  set_engine("brulee", validation = 0) %>%
  set_mode("classification")

# Create a workflow
nnet_wflow <- biv_rec %>% workflow(nnet_spec)

# Fit the model
set.seed(987)
nnet_fit <- fit(nnet_wflow, cls_train)








# Install and load the nnls package (if not already installed)
# install.packages("nnls")
library(nnls)

# Example vectors
a <- c(1, 2, 3)
b <- c(3, 2, 1)
target_vector <- c(4, 5, 6)  # The vector you want to express

# Create a matrix from vectors a and b
ab_matrix <- matrix(c(a, b), nrow = length(a))

# Solve for coefficients
coefficients <- nnls(A = ab_matrix, b = target_vector)$x

# Print the coefficients
print(coefficients)

reconstr <- rowSums(coefficients*ab_matrix)


