
##############################################################################################################
################################ LOAD LIBRARIES, FILES AND CLEAN THE ENVRIONMENT #############################################
#############################################################################################################

#Clean the environment
rm(list = ls()) 

#Load the libraries
library(dplyr)
library(tidyr)
library(tidyverse)
library(readxl)
library(patchwork)


##############################################################################################################
################################ CALCULATE THE RECOVERY #####################################################
#############################################################################################################

RECOVERY <- read_excel("F:/LINKOPING/Manuscripts/Skyline/Skyline/OrbitrapDust.xlsx") |> 
  filter(`Isotope Label Type` == "Quan") |> 
  pivot_wider(id_cols = c(`Replicate Name`, `Sample Type`),
              names_from = Molecule, # name of the new column
              values_from = Area) #name of the values for the new columns


##########Filter and prepare the standards that we want to compare to, (I named them Std)
Std<-RECOVERY|>  
  select(`Replicate Name`, IS, RS)|> #Select the columns that we need
  filter(str_detect(`Replicate Name`, "Std", negate = FALSE))|> #Select the standards
  mutate(RatioStd = IS / RS) |> #Calculate the IS/RS `Normalized Area` for the standards 
  summarize(AverageRatio = mean(RatioStd, na.rm = TRUE)) #Calculate the average of the three standards

##############Calculate the recovery
RECOVERY<-RECOVERY|>
  filter( `Sample Type` != "Standard") |>  #Filter to remove the standards out of the data frame|>#Exclude the calibration standards
  mutate(RatioSample= IS / RS) |> #Divide IS/RS for all samples
  mutate(Recovery=RatioSample/as.numeric(Std))|> #Divide each sample ratio of IS/RS by the average ratio of the standards
  mutate(RecoveryPerc=Recovery*100) #The recovery in percentage


#################Plotting the recovery

# Define colors based on conditions: less than 50% red and above green
RECOVERY$color <- ifelse(RECOVERY$RecoveryPerc > 50, "#CD3333", "#9ACD32")

#Reorder the samples
RECOVERY$`Replicate Name` <- factor(RECOVERY$`Replicate Name`, levels = c("1", "1B", "1C", "2", "3", "4A", "4B", "5A", "5B", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "18", "19", "20", "D1", "D2", "D3", "21 Na2SO4", "23 24 FB", "25 LB", "26 LB", "27 LB", "28 LB", "29 MLB"))

# Bar plot
ggplot(RECOVERY, aes(x = `Replicate Name`, y = RecoveryPerc, fill = color)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#CD3333", "#9ACD32"), guide = FALSE) +
  labs(x = "Replicate Name", y = "RECOVERY", title = "RECOVERY") +
  theme_minimal()

#################################################################################################
##################################### SCCPs##################################################
##################################################################################################

##############################################################################################################
################################ PREPARE THE STANDARDS FOR CALIBRATION #####################################################
#############################################################################################################

########################### STANDARDS A #############################################################
# Load the file
TESTING <- read_excel("F:/LINKOPING/Manuscripts/Skyline/Skyline/OrbitrapDust.xlsx") |>
  mutate(`Analyte Concentration` = as.numeric(`Analyte Concentration`)) 

# Replace missing values in the Response_factor column with 0
TESTING <- TESTING |> 
         mutate(`Normalized Area` = as.numeric(`Normalized Area`),  # Convert to numeric
                `Normalized Area` = ifelse(is.na(`Normalized Area`), 0, `Normalized Area`))  # Replace NAs with 0

# Function to create molecule names
create_molecule_name <- function(i, j) {
  molecule_name <- paste0("C", i, "H", (2 * i) + 2 - j, "Cl", j)
  return(molecule_name)
}

# Create an empty list to store plots and calibration results
plotsA <- list()
calibration_resultsA <- data.frame(STD_code = character(),
                                   Reference_standard = character(),
                                   Chain_length = character(),
                                   Type = character(),
                                   Homologue = character(),
                                   Response_factor = numeric(),
                                   Intercept = numeric(),
                                   R_squared = numeric(),
                                   stringsAsFactors = FALSE)

# Define the range of i and j
i_values <- 10:13
j_values <- 3:12
# Nested loops to iterate over i and j
for (i in i_values) {
  for (j in j_values) {
    # Create molecule name
    molecule_name <- create_molecule_name(i, j)
    
    # Filter data for the current molecule
    filtered_dataA <- TESTING |>
      filter(`Isotope Label Type` == "Quan", Molecule == molecule_name, Note == "A")|> 
      filter(str_detect(Molecule, "C9", negate = TRUE))|>#Exclude the calibration vSCCPs
      filter(str_detect(Molecule, "C14", negate = TRUE))|>#Exclude the calibration MCCPs
      filter(str_detect(Molecule, "C15", negate = TRUE))|>#Exclude the calibration MCCPs
      filter(str_detect(Molecule, "C16", negate = TRUE))|>#Exclude the calibration MCCPs
      filter(str_detect(Molecule, "C17", negate = TRUE)) |>#Exclude the calibration MCCPs
      filter(str_detect(Molecule, "IS", negate = TRUE)) |> #Exclude the calibration IS
      filter(str_detect(Molecule, "IS", negate = TRUE)) #Exclude the calibration RS 
    
    # Check if there are any non-NA cases in the filtered data
    if (sum(!is.na(filtered_dataA$`Normalized Area`)) == 0 || sum(!is.na(filtered_dataA$`Analyte Concentration`)) == 0) {
      cat("No valid cases for fitting the model for molecule:", molecule_name, "\n")
    } else {
      # Fit linear model to the data
      lm_modelA <- lm(`Normalized Area` ~ `Analyte Concentration`, data = filtered_dataA)
      
      # Extract slope and intercept from the model
      slopeA <- coef(lm_modelA)[2]
      interceptA <- coef(lm_modelA)[1]
      
      # Calculate R-squared
      R_squaredA <- summary(lm_modelA)$r.squared
      
      # Determine the type based on the value of i
      if (i >= 10 && i <= 13) {
        type <- "SCCPs"
      } else if (i >= 14 && i <= 17) {
        type <- "MCCPs"
      } else if (i >= 18 && i <= 30) {
        type <- "LCCPs"
      } else {
        type <- "Unknown"
      }
      
      # Store the calibration results in the data frame, create the columns and values that we need to fit with the next script (the one for the deconvolution paper published by I don't remember)
      calibration_resultsA <- rbind(calibration_resultsA, 
                                    data.frame(STD_code = paste("A-", type),
                                               Reference_standard = paste(type, "C", i, "Cl%", "52"),
                                               Chain_length = paste("C",i),
                                               Type = type,
                                               Homologue = molecule_name,
                                               Response_factor = round(slopeA),
                                               Intercept = interceptA,
                                               R_squared = R_squaredA))
      
      
      # Create plot for the current molecule
      plotA <- ggplot(filtered_dataA, aes(x = `Analyte Concentration`, y = `Normalized Area`)) +
        geom_point() +
        geom_smooth(method = "lm", 
                    se = FALSE, 
                    colour = "orange",
                    size = 1) +
        theme_classic() +
        ylab("`Normalized Area`") +
        ggtitle(paste("Molecule:", molecule_name))
      
      # Store the plot in the list
      plotsA[[paste("Molecule", molecule_name)]] <- plotA
    }
  }
}
# Print the data frame with calibration results
print(calibration_resultsA)


# Arrange plots in a grid
calibration_curves_gridA <- plotsA[[1]]
for (i in 2:length(plotsA)) {
  calibration_curves_gridA <- calibration_curves_gridA + plotsA[[i]]
}

# Print the grid of calibration curves
calibration_curves_gridA


########################### STANDARDS B #####################################################################################



# Replace missing values in the Response_factor column with 0
TESTING <- TESTING |> 
  mutate(`Normalized Area` = replace_na(`Normalized Area`, 0)) 

# Function to create molecule names
create_molecule_name <- function(i, j) {
  molecule_name <- paste0("C", i, "H", (2 * i) + 2 - j, "Cl", j)
  return(molecule_name)
}

# Create an empty list to store plots and calibration results
plotsB <- list()
calibration_resultsB <- data.frame(STD_code = character(),
                                   Reference_standard = character(),
                                   Chain_length = character(),
                                   Type = character(),
                                   Homologue = character(),
                                   Response_factor = numeric(),
                                   Intercept = numeric(),
                                   R_squared = numeric(),
                                   stringsAsFactors = FALSE)

# Define the range of i and j
i_values <- 10:13
j_values <- 3:12
# Nested loops to iterate over i and j
for (i in i_values) {
  for (j in j_values) {
    # Create molecule name
    molecule_name <- create_molecule_name(i, j)
    
    # Filter data for the current molecule
    filtered_dataB <- TESTING |>
      filter(`Isotope Label Type` == "Quan", Molecule == molecule_name, Note == "B") |> 
      filter(str_detect(Molecule, "C9", negate = TRUE))|>#Exclude the calibration vSCCPs
      filter(str_detect(Molecule, "C14", negate = TRUE))|>#Exclude the calibration MCCPs
      filter(str_detect(Molecule, "C15", negate = TRUE))|>#Exclude the calibration MCCPs
      filter(str_detect(Molecule, "C16", negate = TRUE))|>#Exclude the calibration MCCPs
      filter(str_detect(Molecule, "C17", negate = TRUE))|>#Exclude the calibration MCCPs
      filter(str_detect(Molecule, "IS", negate = TRUE)) |> #Exclude the calibration IS
      filter(str_detect(Molecule, "IS", negate = TRUE)) #Exclude the calibration RS
    
    # Check if there are any non-NA cases in the filtered data
    if (sum(!is.na(filtered_dataB$`Normalized Area`)) == 0 || sum(!is.na(filtered_dataB$`Analyte Concentration`)) == 0) {
      cat("No valid cases for fitting the model for molecule:", molecule_name, "\n")
    } else {
      # Fit linear model to the data
      lm_modelB <- lm(`Normalized Area` ~ `Analyte Concentration`, data = filtered_dataB)
      
      # Extract slope and intercept from the model
      slopeB <- coef(lm_modelB)[2]
      interceptB <- coef(lm_modelB)[1]
      
      # Calculate R-squared
      R_squaredB <- summary(lm_modelB)$r.squared
      
      # Determine the type based on the value of i
      if (i >= 10 && i <= 13) {
        type <- "SCCPs"
      } else if (i >= 14 && i <= 17) {
        type <- "MCCPs"
      } else if (i >= 18 && i <= 30) {
        type <- "LCCPs"
      } else {
        type <- "Unknown"
      }
      
      # Store the calibration results in the data frame, create the columns and values that we need to fit with the next script (the one for the deconvolution paper published by I don't remember)
      calibration_resultsB <- rbind(calibration_resultsB, 
                                    data.frame(STD_code = paste("B-", type),
                                               Reference_standard = paste(type, "C", i, "Cl%", "52"),
                                               Chain_length = paste("C",i),
                                               Type = type,
                                               Homologue = molecule_name,
                                               Response_factor = round(slopeB),
                                               Intercept = interceptB,
                                               R_squared = R_squaredB))
      
      
      # Create plot for the current molecule
      plotB <- ggplot(filtered_dataB, aes(x = `Analyte Concentration`, y = `Normalized Area`)) +
        geom_point() +
        geom_smooth(method = "lm", 
                    se = FALSE, 
                    colour = "purple",
                    size = 1) +
        theme_classic() +
        ylab("`Normalized Area`") +
        ggtitle(paste("Molecule:", molecule_name))
      
      # Store the plot in the list
      plotsB[[paste("Molecule", molecule_name)]] <- plotB
    }
  }
}
# Print the data frame with calibration results
print(calibration_resultsB)


# Arrange plots in a grid
calibration_curves_gridB <- plotsB[[1]]
for (i in 2:length(plotsB)) {
  calibration_curves_gridB <- calibration_curves_gridB + plotsB[[i]]
}

# Print the grid of calibration curves
calibration_curves_gridB


########################### COMBINE STNADARDS A AND B ######################################################################

combined_df <- rbind(calibration_resultsA, calibration_resultsB)


##############################################################################################################
################################ PREPARE THE SAMPLES FOR CALIBRATION #####################################################
#############################################################################################################

# Load the data
TESTINGB <- TESTING |> 
  filter(`Isotope Label Type` == "Quan") |> 
  mutate(
    `Normalized Area` = as.numeric(`Normalized Area`),
    Chain_length = paste("C", sub(".*-C(9|10|11|12|13|14|15|16|17)$", "\\1", `Molecule List`)), #Create the columns that contains the chain length as for example C10 or C11 or C12 etc.
    type = ifelse(
      Chain_length >= 10 & Chain_length <= 13, "SCCPs",
      ifelse(
        Chain_length >= 14 & Chain_length <= 17, "MCCPs",
        ifelse(Chain_length >= 18 & Chain_length <= 30, "LCCPs", "Unknown") #Create the column type which should contain SCCPs, MCCPs and LCCPs depending on the chain length
      )
    )
  ) |> 
  rename(Homologue = Molecule) |> 
  filter( `Sample Type` != "Standard" #Filter to remove the standards out of the data frame
  ) |> 
  filter(str_detect(Homologue, "C9", negate = TRUE))|>#Exclude the calibration vSCCPs
  filter(str_detect(Homologue, "C14", negate = TRUE))|>#Exclude the calibration MCCPs
  filter(str_detect(Homologue, "C15", negate = TRUE))|>#Exclude the calibration MCCPs
  filter(str_detect(Homologue, "C16", negate = TRUE))|>#Exclude the calibration MCCPs
  filter(str_detect(Homologue, "C17", negate = TRUE))|>#Exclude the calibration MCCPs
  filter(str_detect(Homologue, "IS", negate = TRUE)) |> #Exclude the calibration IS
  filter(str_detect(Homologue, "RS", negate = TRUE)) #Exclude the calibration RS

# Group samples for creating one data frame for each: the CALIBRATION script requires to have one data frame for each sample, so out of the current data frame we should create a new of for each of the samples
list_of_samples <- split(TESTINGB, TESTINGB$`Replicate Name`)



##############################################################################################################
########################################## CALIBRATION #####################################################
#############################################################################################################

#####GROUP STANDARD MIXTURES USED#####

#####Set working directory#####
working.directory <- "F:/LINKOPING/Manuscripts/Skyline/Skyline/"

##########################################PREPARE DATASET FOR PATTERN RECONSTRUCTION#################
{
  # Rename the file of the standards that we have created previous for calibration so it will fit the script
  input <- combined_df
  input$STD_code <- as.factor(input$STD_code)
  
  # Create all possible binary combinations between given calibration sets
  Combinations <- combn(x = levels(input$STD_code), m = 2, FUN = NULL, simplify = TRUE)
  
  # Store sum RFs for each group CP standard
  input <- input |> 
    group_by(Reference_standard) |> 
    mutate(Sum_response_factor = sum(Response_factor, na.rm = TRUE))
  input[1:5] <- lapply(input[1:5], as.factor)
  input$Response_factor[is.na(input$Response_factor)] <- 0
}

############################################SECTION FOR SAMPLE PROCESSING##########################
## The script was made for quantifying one sample at the time, so I changed it to iterate for all the data frames of the samples

# Initialize an empty list to store results for all samples
all_results <- list()
all_plots <- list()

# Extract unique sample names, it will give the to each "data frame"
unique_sample_names <- unique(TESTINGB$`Replicate Name`)

# Iterate over each unique sample name
for (sample_name in unique_sample_names) {
  # Get the corresponding sample data frame, now I will use sample_df for further procesing
  sample_df <- list_of_samples[[sample_name]]
  
  # Set sample name, to see which samples it will calibrate
  print(paste("Processing sample:", sample_name))
  
  
  
  ####################################RUN PATTERN RECONSTRUCTION FOR SELECTED (LOADED) SAMPLE####################
  #This section is mostly how it was given, I only included the info to fill in the empty list that I created at the begging so I will get the results for each sample. I also combined all results into one data frame so I would not have one data frame for ach sample. 
  
  sample_df <- sample_df |> 
    mutate(
      Chain_length = as.factor(Chain_length),
      `Normalized Area` = as.numeric(`Normalized Area`),
      Relative_distribution = `Normalized Area` / sum(`Normalized Area`, na.rm = TRUE)
    )
  
  # Calculate relative `Normalized Area` distribution within each homologue group
  sample_df$Relative_distribution <- NA
  sample_df$`Normalized Area`[is.na(sample_df$`Normalized Area`)] <- 0
  sample_df <- sample_df |>  mutate(Relative_distribution = `Normalized Area` / sum(`Normalized Area`, na.rm = TRUE))
  results <- sample_df
  results[c("Comp_1", "Comp_2", "Fraction_Comp_1", "Simulated_pattern")] <- NA
  
  # Deconvolution of homologue patterns (adjust your loop accordingly)
  REF <- sample_df$Relative_distribution
  Distance <- 100
  
  for (z in 1:length(Combinations[1,])) {
    C_1 <- subset(input, subset = (STD_code == Combinations[1, z]))
    C_2 <- subset(input, subset = (STD_code == Combinations[2, z]))
    
    for (j in 1:100) {
      Combo <- (C_1$Response_factor * j + C_2$Response_factor * (100 - j)) / sum((C_1$Response_factor * j + C_2$Response_factor * (100 - j)), na.rm = TRUE)
      
      # Check for NA values in REF or Combo before comparison
      if (!any(is.na(REF)) && !any(is.na(Combo)) && Distance > (sum(sqrt((REF - Combo)^2)))) {
        results$Comp_1 <- as.character(C_1$STD_code)
        results$Comp_2 <- as.character(C_2$STD_code)
        results$Fraction_Comp_1 <- j
        results$Simulated_pattern <- Combo
        Distance <- sum(sqrt((REF - Combo)^2))
      }
    }
  }
  
  
  # Calculate concentrations (ng per microliter)
  results$RF_1st <- NA
  results$RF_2nd <- NA
  
  for (i in 1:nrow(results)) {
    results$RF_1st[i] <- input$Sum_response_factor[input$STD_code == results$Comp_1[i]]
    results$RF_2nd[i] <- input$Sum_response_factor[input$STD_code == results$Comp_2[i]]
  }
  
  results <- results |> 
    mutate(
      RF_1st = as.numeric(RF_1st),
      RF_2nd = as.numeric(RF_2nd),
      Concentration = sum(`Normalized Area`) / (RF_1st * (Fraction_Comp_1 / 100) + RF_2nd * ((100 - Fraction_Comp_1) / 100))
    )
  
  # Store the results for the current sample in the list
  all_results[[sample_name]] <- results
  
  # Visualization of results
  plot_table <- data.frame(
    Distribution = c(results$Relative_distribution, results$Simulated_pattern),
    Homologue = results$Homologue,
    Chain_length = results$Chain_length,
    Origin = rep(as.factor(c("Measured", "Simulated")), each = nrow(results))
  )
  
  plot_table$Homologue <- factor(plot_table$Homologue, levels = unique(plot_table$Homologue))
  
  plot <- ggplot(plot_table, aes(x = Homologue, y = Distribution * 100, fill = Origin, colour = Origin)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8, size = .4) +
    theme(panel.background = element_blank()) +
    scale_fill_manual(values = c("darkolivegreen3", "darkslategray4")) +
    scale_color_manual(values = c("darkolivegreen4", "darkslategray")) +
    ggtitle(label = paste(sample_name, " - Distribution of CP homologues")) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) +
    xlab("") + ylab("Relative `Normalized Area` distribution, %") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(
      legend.key.size = unit(0.15, "in"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(fill = NA, colour = "grey20"),
      panel.grid.major.y = element_line(colour = "grey50"),
      panel.grid.minor.y = element_line(colour = "grey80"),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey20"),
      strip.text = element_text(colour = "white", face = "bold")
    ) +
    facet_wrap(. ~ Chain_length, scales = "free", nrow = 4, ncol = 4)
  
  # Store the plot for the current sample in the list
  all_plots[[sample_name]] <- plot
  
  results_output_SCCPs <- results |> 
    summarise(median(Concentration)) |> 
    mutate(
      Type = results$Type[1],
      Sample = sample_name,
      Comment = paste("The best match:", results$Fraction_Comp_1[1], "% of ", results$Comp_1[1], " and ", 100 - results$Fraction_Comp_1[1], "% of ", results$Comp_2[1])
    ) |> 
    rename("Total concentration, ng/µL" = "median(Concentration)")
  
  # Combine results for all samples into a single dataframe
  all_results_df_SCCPs <- bind_rows(all_results, .id = "Sample")
  
  # Print or further process the combined results dataframe
  print(all_results_df_SCCPs)
  # Print results output
  print(results_output_SCCPs)
  
}

#########################################VIEW RESULTS############################################

#View Overview
print(all_results_df_SCCPs)

#View graph, REPLACE THE NAME OF THE SAMPLE YOU WANT TO SEE, another option is to open the list from the Environment and from there open each plot
all_plots[["NIST_R1"]]
all_plots[["NIST_R2"]]
all_plots[["NIST_R3"]]



#################################################################################################
##################################### MCCPs##################################################
##################################################################################################

##############################################################################################################
################################ PREPARE THE STANDARDS FOR CALIBRATION #####################################################
#############################################################################################################

########################### STANDARDS A #############################################################
# Load the file
TESTING <- read_excel("F:/LINKOPING/Manuscripts/Skyline/Skyline/OrbitrapDust.xlsx") |>
  mutate(`Analyte Concentration` = as.numeric(`Analyte Concentration`)) 

# Replace missing values in the Response_factor column with 0
TESTING <- TESTING %>%
  mutate(`Normalized Area` = ifelse(is.na(`Normalized Area`) | is.nan(`Normalized Area`), 0, `Normalized Area`)) #Replace with 0
   
         

# Function to create molecule names
create_molecule_name <- function(i, j) {
  molecule_name <- paste0("C", i, "H", (2 * i) + 2 - j, "Cl", j)
  return(molecule_name)
}

# Create an empty list to store plots and calibration results
plotsA <- list()
calibration_resultsA <- data.frame(STD_code = character(),
                                   Reference_standard = character(),
                                   Chain_length = character(),
                                   Type = character(),
                                   Homologue = character(),
                                   Response_factor = numeric(),
                                   Intercept = numeric(),
                                   R_squared = numeric(),
                                   stringsAsFactors = FALSE)

# Define the range of i and j
i_values <- 14:17
j_values <- 3:12
# Nested loops to iterate over i and j
for (i in i_values) {
  for (j in j_values) {
    # Create molecule name
    molecule_name <- create_molecule_name(i, j)
    
    # Filter data for the current molecule
    filtered_dataA <- TESTING |>
      filter(`Isotope Label Type` == "Quan", Molecule == molecule_name, Note == "D")|> 
      filter(str_detect(Molecule, "C9", negate = TRUE))|>#Exclude the calibration vSCCPs
      filter(str_detect(Molecule, "C10", negate = TRUE))|>#Exclude the calibration SCCPs
      filter(str_detect(Molecule, "C11", negate = TRUE))|>#Exclude the calibration SCCPs
      filter(str_detect(Molecule, "C12", negate = TRUE))|>#Exclude the calibration SCCPs
      filter(str_detect(Molecule, "C13", negate = TRUE)) |>#Exclude the calibration SCCPs
      filter(str_detect(Molecule, "IS", negate = TRUE)) |> #Exclude the calibration IS
      filter(str_detect(Molecule, "RS", negate = TRUE)) #Exclude the calibration RS
    
    # Check if there are any non-NA cases in the filtered data
    if (sum(!is.na(filtered_dataA$`Normalized Area`)) == 0 || sum(!is.na(filtered_dataA$`Analyte Concentration`)) == 0) {
      cat("No valid cases for fitting the model for molecule:", molecule_name, "\n")
    } else {
      # Fit linear model to the data
      lm_modelA <- lm(`Normalized Area` ~ `Analyte Concentration`, data = filtered_dataA)
      
      # Extract slope and intercept from the model
      slopeA <- coef(lm_modelA)[2]
      interceptA <- coef(lm_modelA)[1]
      
      # Calculate R-squared
      R_squaredA <- summary(lm_modelA)$r.squared
      
      # Determine the type based on the value of i
      if (i >= 10 && i <= 13) {
        type <- "SCCPs"
      } else if (i >= 14 && i <= 17) {
        type <- "MCCPs"
      } else if (i >= 18 && i <= 30) {
        type <- "LCCPs"
      } else {
        type <- "Unknown"
      }
      
      # Store the calibration results in the data frame, create the columns and values that we need to fit with the next script (the one for the deconvolution paper published by I don't remember)
      calibration_resultsA <- rbind(calibration_resultsA, 
                                    data.frame(STD_code = paste("A-", type),
                                               Reference_standard = paste(type, "C", i, "Cl%", "52"),
                                               Chain_length = paste("C",i),
                                               Type = type,
                                               Homologue = molecule_name,
                                               Response_factor = round(slopeA),
                                               Intercept = interceptA,
                                               R_squared = R_squaredA))
      
      
      # Create plot for the current molecule
      plotA <- ggplot(filtered_dataA, aes(x = `Analyte Concentration`, y = `Normalized Area`)) +
        geom_point() +
        geom_smooth(method = "lm", 
                    se = FALSE, 
                    colour = "orange",
                    size = 1) +
        theme_classic() +
        ylab("`Normalized Area`") +
        ggtitle(paste("Molecule:", molecule_name))
      
      # Store the plot in the list
      plotsA[[paste("Molecule", molecule_name)]] <- plotA
    }
  }
}
# Print the data frame with calibration results
print(calibration_resultsA)


# Arrange plots in a grid
calibration_curves_gridA <- plotsA[[1]]
for (i in 2:length(plotsA)) {
  calibration_curves_gridA <- calibration_curves_gridA + plotsA[[i]]
}

# Print the grid of calibration curves
calibration_curves_gridA


########################### STANDARDS B #####################################################################################


# Replace missing values in the Response_factor column with 0
TESTING <- TESTING |> 
  mutate(`Normalized Area` = replace_na(`Normalized Area`, 0))

# Function to create molecule names
create_molecule_name <- function(i, j) {
  molecule_name <- paste0("C", i, "H", (2 * i) + 2 - j, "Cl", j)
  return(molecule_name)
}

# Create an empty list to store plots and calibration results
plotsB <- list()
calibration_resultsB <- data.frame(STD_code = character(),
                                   Reference_standard = character(),
                                   Chain_length = character(),
                                   Type = character(),
                                   Homologue = character(),
                                   Response_factor = numeric(),
                                   Intercept = numeric(),
                                   R_squared = numeric(),
                                   stringsAsFactors = FALSE)

# Define the range of i and j
i_values <- 14:17
j_values <- 3:12
# Nested loops to iterate over i and j
for (i in i_values) {
  for (j in j_values) {
    # Create molecule name
    molecule_name <- create_molecule_name(i, j)
    
    # Filter data for the current molecule
    filtered_dataB <- TESTING |>
      filter(`Isotope Label Type` == "Quan", Molecule == molecule_name, Note == "E") |> 
      filter(str_detect(Molecule, "C9", negate = TRUE))|>#Exclude the calibration vSCCPs
      filter(str_detect(Molecule, "C10", negate = TRUE))|>#Exclude the calibration SCCPs
      filter(str_detect(Molecule, "C11", negate = TRUE))|>#Exclude the calibration SCCPs
      filter(str_detect(Molecule, "C12", negate = TRUE))|>#Exclude the calibration SCCPs
      filter(str_detect(Molecule, "C13", negate = TRUE)) |>#Exclude the calibration SCCPs
      filter(str_detect(Molecule, "IS", negate = TRUE)) |> #Exclude the calibration IS
      filter(str_detect(Molecule, "RS", negate = TRUE)) #Exclude the calibration RS 
    
    # Check if there are any non-NA cases in the filtered data
    if (sum(!is.na(filtered_dataB$`Normalized Area`)) == 0 || sum(!is.na(filtered_dataB$`Analyte Concentration`)) == 0) {
      cat("No valid cases for fitting the model for molecule:", molecule_name, "\n")
    } else {
      # Fit linear model to the data
      lm_modelB <- lm(`Normalized Area` ~ `Analyte Concentration`, data = filtered_dataB)
      
      # Extract slope and intercept from the model
      slopeB <- coef(lm_modelB)[2]
      interceptB <- coef(lm_modelB)[1]
      
      # Calculate R-squared
      R_squaredB <- summary(lm_modelB)$r.squared
      
      # Determine the type based on the value of i
      if (i >= 10 && i <= 13) {
        type <- "SCCPs"
      } else if (i >= 14 && i <= 17) {
        type <- "MCCPs"
      } else if (i >= 18 && i <= 30) {
        type <- "LCCPs"
      } else {
        type <- "Unknown"
      }
      
      # Store the calibration results in the data frame, create the columns and values that we need to fit with the next script (the one for the deconvolution paper published by I don't remember)
      calibration_resultsB <- rbind(calibration_resultsB, 
                                    data.frame(STD_code = paste("B-", type),
                                               Reference_standard = paste(type, "C", i, "Cl%", "52"),
                                               Chain_length = paste("C",i),
                                               Type = type,
                                               Homologue = molecule_name,
                                               Response_factor = round(slopeB),
                                               Intercept = interceptB,
                                               R_squared = R_squaredB))
      
      
      # Create plot for the current molecule
      plotB <- ggplot(filtered_dataB, aes(x = `Analyte Concentration`, y = `Normalized Area`)) +
        geom_point() +
        geom_smooth(method = "lm", 
                    se = FALSE, 
                    colour = "purple",
                    size = 1) +
        theme_classic() +
        ylab("`Normalized Area`") +
        ggtitle(paste("Molecule:", molecule_name))
      
      # Store the plot in the list
      plotsB[[paste("Molecule", molecule_name)]] <- plotB
    }
  }
}
# Print the data frame with calibration results
print(calibration_resultsB)


# Arrange plots in a grid
calibration_curves_gridB <- plotsB[[1]]
for (i in 2:length(plotsB)) {
  calibration_curves_gridB <- calibration_curves_gridB + plotsB[[i]]
}

# Print the grid of calibration curves
calibration_curves_gridB


########################### COMBINE STNADARDS A AND B ######################################################################

combined_df <- rbind(calibration_resultsA, calibration_resultsB)


##############################################################################################################
################################ PREPARE THE SAMPLES FOR CALIBRATION #####################################################
#############################################################################################################

# Load the data
TESTINGB <- TESTING |> 
  filter(`Isotope Label Type` == "Quan") |> 
  mutate(
    `Normalized Area` = as.numeric(`Normalized Area`),
    Chain_length = paste("C", sub(".*-C(9|10|11|12|13|14|15|16|17)$", "\\1", `Molecule List`)), #Create the columns that contains the chain length as for example C10 or C11 or C12 etc.
    type = ifelse(
      Chain_length >= 10 & Chain_length <= 13, "SCCPs",
      ifelse(
        Chain_length >= 14 & Chain_length <= 17, "MCCPs",
        ifelse(Chain_length >= 18 & Chain_length <= 30, "LCCPs", "Unknown") #Create the column type which should contain SCCPs, MCCPs and LCCPs depending on the chain length
      )
    )
  ) |> 
  rename(Homologue = Molecule) |> 
  filter( `Sample Type` != "Standard" #Filter to remove the standards out of the data frame
  ) |> 
  filter(str_detect(Homologue, "C9", negate = TRUE))|>#Exclude the calibration vSCCPs
  filter(str_detect(Homologue, "C10", negate = TRUE))|>#Exclude the calibration SCCPs
  filter(str_detect(Homologue, "C11", negate = TRUE))|>#Exclude the calibration SCCPs
  filter(str_detect(Homologue, "C12", negate = TRUE))|>#Exclude the calibration SCCPs
  filter(str_detect(Homologue, "C13", negate = TRUE)) |>#Exclude the calibration SCCPs
  filter(str_detect(Homologue, "IS", negate = TRUE)) |> #Exclude the calibration IS
  filter(str_detect(Homologue, "RS", negate = TRUE)) #Exclude the calibration RS

# Group samples for creating one data frame for each: the CALIBRATION script requires to have one data frame for each sample, so out of the current data frame we should create a new of for each of the samples
list_of_samples <- split(TESTINGB, TESTINGB$`Replicate Name`)



##############################################################################################################
########################################## CALIBRATION #####################################################
#############################################################################################################

#####GROUP STANDARD MIXTURES USED#####

#####Set working directory#####
working.directory <- "F:/LINKOPING/Manuscripts/Skyline/Skyline/"

##########################################PREPARE DATASET FOR PATTERN RECONSTRUCTION#################
{
  # Rename the file of the standards that we have created previous for calibration so it will fit the script
  input <- combined_df
  input$STD_code <- as.factor(input$STD_code)
  
  # Create all possible binary combinations between given calibration sets
  Combinations <- combn(x = levels(input$STD_code), m = 2, FUN = NULL, simplify = TRUE)
  
  # Store sum RFs for each group CP standard
  input <- input |> 
    group_by(Reference_standard) |> 
    mutate(Sum_response_factor = sum(Response_factor, na.rm = TRUE))
  input[1:5] <- lapply(input[1:5], as.factor)
  input$Response_factor[is.na(input$Response_factor)] <- 0
}

############################################SECTION FOR SAMPLE PROCESSING##########################
## The script was made for quantifying one sample at the time, so I changed it to iterate for all the data frames of the samples

# Initialize an empty list to store results for all samples
all_results <- list()
all_plots <- list()

# Extract unique sample names, it will give the to each "data frame"
unique_sample_names <- unique(TESTINGB$`Replicate Name`)

# Iterate over each unique sample name
for (sample_name in unique_sample_names) {
  # Get the corresponding sample data frame, now I will use sample_df for further procesing
  sample_df <- list_of_samples[[sample_name]]
  
  # Set sample name, to see which samples it will calibrate
  print(paste("Processing sample:", sample_name))
  
  
  
  ####################################RUN PATTERN RECONSTRUCTION FOR SELECTED (LOADED) SAMPLE####################
  #This section is mostly how it was given, I only included the info to fill in the empty list that I created at the begging so I will get the results for each sample. I also combined all results into one data frame so I would not have one data frame for ach sample. 
  
  sample_df <- sample_df |> 
    mutate(
      Chain_length = as.factor(Chain_length),
      `Normalized Area` = as.numeric(`Normalized Area`),
      Relative_distribution = `Normalized Area` / sum(`Normalized Area`, na.rm = TRUE)
    )
  
  # Calculate relative `Normalized Area` distribution within each homologue group
  sample_df$Relative_distribution <- NA
  sample_df$`Normalized Area`[is.na(sample_df$`Normalized Area`)] <- 0
  sample_df <- sample_df |>  mutate(Relative_distribution = `Normalized Area` / sum(`Normalized Area`, na.rm = TRUE))
  results <- sample_df
  results[c("Comp_1", "Comp_2", "Fraction_Comp_1", "Simulated_pattern")] <- NA
  
  # Deconvolution of homologue patterns (adjust your loop accordingly)
  REF <- sample_df$Relative_distribution
  Distance <- 100
  
  for (z in 1:length(Combinations[1,])) {
    C_1 <- subset(input, subset = (STD_code == Combinations[1, z]))
    C_2 <- subset(input, subset = (STD_code == Combinations[2, z]))
    
    for (j in 1:100) {
      Combo <- (C_1$Response_factor * j + C_2$Response_factor * (100 - j)) / sum((C_1$Response_factor * j + C_2$Response_factor * (100 - j)), na.rm = TRUE)
      
      # Check for NA values in REF or Combo before comparison
      if (!any(is.na(REF)) && !any(is.na(Combo)) && Distance > (sum(sqrt((REF - Combo)^2)))) {
        results$Comp_1 <- as.character(C_1$STD_code)
        results$Comp_2 <- as.character(C_2$STD_code)
        results$Fraction_Comp_1 <- j
        results$Simulated_pattern <- Combo
        Distance <- sum(sqrt((REF - Combo)^2))
      }
    }
  }
  
  
  # Calculate concentrations (ng per microliter)
  results$RF_1st <- NA
  results$RF_2nd <- NA
  
  for (i in 1:nrow(results)) {
    results$RF_1st[i] <- input$Sum_response_factor[input$STD_code == results$Comp_1[i]]
    results$RF_2nd[i] <- input$Sum_response_factor[input$STD_code == results$Comp_2[i]]
  }
  
  results <- results |> 
    mutate(
      RF_1st = as.numeric(RF_1st),
      RF_2nd = as.numeric(RF_2nd),
      Concentration = sum(`Normalized Area`) / (RF_1st * (Fraction_Comp_1 / 100) + RF_2nd * ((100 - Fraction_Comp_1) / 100))
    )
  
  # Store the results for the current sample in the list
  all_results[[sample_name]] <- results
  
  # Visualization of results
  plot_table <- data.frame(
    Distribution = c(results$Relative_distribution, results$Simulated_pattern),
    Homologue = results$Homologue,
    Chain_length = results$Chain_length,
    Origin = rep(as.factor(c("Measured", "Simulated")), each = nrow(results))
  )
  
  plot_table$Homologue <- factor(plot_table$Homologue, levels = unique(plot_table$Homologue))
  
  plot <- ggplot(plot_table, aes(x = Homologue, y = Distribution * 100, fill = Origin, colour = Origin)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8, size = .4) +
    theme(panel.background = element_blank()) +
    scale_fill_manual(values = c("darkolivegreen3", "darkslategray4")) +
    scale_color_manual(values = c("darkolivegreen4", "darkslategray")) +
    ggtitle(label = paste(sample_name, " - Distribution of CP homologues")) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) +
    xlab("") + ylab("Relative `Normalized Area` distribution, %") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(
      legend.key.size = unit(0.15, "in"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(fill = NA, colour = "grey20"),
      panel.grid.major.y = element_line(colour = "grey50"),
      panel.grid.minor.y = element_line(colour = "grey80"),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey20"),
      strip.text = element_text(colour = "white", face = "bold")
    ) +
    facet_wrap(. ~ Chain_length, scales = "free", nrow = 4, ncol = 4)
  
  # Store the plot for the current sample in the list
  all_plots[[sample_name]] <- plot
  
  results_output_MCCPs <- results |> 
    summarise(median(Concentration)) |> 
    mutate(
      Type = results$Type[1],
      Sample = sample_name,
      Comment = paste("The best match:", results$Fraction_Comp_1[1], "% of ", results$Comp_1[1], " and ", 100 - results$Fraction_Comp_1[1], "% of ", results$Comp_2[1])
    ) |> 
    rename("Total concentration, ng/µL" = "median(Concentration)")
  
  # Combine results for all samples into a single dataframe
  all_results_df_MCCPs <- bind_rows(all_results, .id = "Sample")
  
  # Print or further process the combined results dataframe
  print(all_results_df_MCCPs)
  # Print results output
  print(results_output_MCCPs)
  
}

#########################################VIEW RESULTS############################################

#View Overview
print(all_results_df_MCCPs)

#View graph, REPLACE THE NAME OF THE SAMPLE YOU WANT TO SEE, another option is to open the list from the Environment and from there open each plot
all_plots[["NIST_R1"]]
all_plots[["NIST_R2"]]
all_plots[["NIST_R3"]]


########################### COMBINE SCCPs AND MCCPs ######################################################################

combined_RESULTS <- rbind(all_results_df_SCCPs, all_results_df_MCCPs)








##############################################################################################################
################################ NORMALIZE THE `Normalized Area` #####################################################
#############################################################################################################

###I think it would be better if we normalize the `Normalized Area` in Skyline, so it is easier to decide which labeled standard to use for each homologue etc. 

#Load the file
#GC_qToF_CPsFoodResults<-read_excel("F:/OREBRO/CP analysis/Food project/Result/Skyline/GC_qToF_CPsFoodResults.xlsx")

#Wide the data so each sample is in one row
#GC_qToF_CPsFoodResultsB <- GC_qToF_CPsFoodResults |>  
#filter(Quantitative == "YES") |> 
#pivot_wider(id_cols = `Replicate Name`,
#names_from = Molecule, # name of the new column
#values_from = `Normalized Area` #name of the values for the new columns
#) #data type       

#Make the data in the data frame numeric
#pattern <- "^C\\d+" # This pattern matches column names starting with "C" followed by numbers
#GC_qToF_CPsFoodResultsB <- mutate_at(GC_qToF_CPsFoodResultsB, vars(matches(pattern)), as.numeric) # Apply as.numeric() to columns matching the pattern
#GC_qToF_CPsFoodResultsB$IS = as.numeric(GC_qToF_CPsFoodResultsB$IS)# Define the IS as numeric 

#Divide the `Normalized Area` of the homologues by the IS
#Corrected`Normalized Area` <- mutate_at(GC_qToF_CPsFoodResultsB, vars(matches(pattern)), 
#~ #. / GC_qToF_CPsFoodResultsB$IS) #It will now do it for all the columns that have the pattern (start the name with "C")




