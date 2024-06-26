
#Load the libraries
library(dplyr)
library(tidyr)
library(tidyverse)
library(readxl)
library(patchwork)
library(stringr)

#Clean the environment
rm(list = ls()) 

# Load the file
TESTING <- read_excel("F:/LINKOPING/Manuscripts/Skyline/Skyline/OrbitrapDust.xlsx") |>
  mutate(`Analyte Concentration` = as.numeric(`Analyte Concentration`)) 

# Replace missing values in the Response_factor column with 0
TESTING <- TESTING |> 
  mutate(`Normalized Area` = ifelse(is.na(`Normalized Area`), 0, `Normalized Area`))  # Replace NAs with 0


##############################################################################################################
################################ PREPARE THE STANDARDS FOR CALIBRATION #####################################################
#############################################################################################################

#Filter data for 
filtered_dataA <- TESTING |>
  filter(`Sample Type` == "Standard") |> #Remove the samples etc. and only take the calibration std
  filter(`Isotope Label Type` == "Quan") |>  #Remove the qual ions
  filter((str_detect(`Replicate Name`, "^SCCP"))) |> #Take only the SCCPs
  filter(str_detect(Molecule, "C9", negate = TRUE))|>#Exclude the calibration vSCCPs
  filter(str_detect(Molecule, "C14", negate = TRUE))|>#Exclude the calibration MCCPs
  filter(str_detect(Molecule, "C15", negate = TRUE))|>#Exclude the calibration MCCPs
  filter(str_detect(Molecule, "C16", negate = TRUE))|>#Exclude the calibration MCCPs
  filter(str_detect(Molecule, "C17", negate = TRUE))|>#Exclude the calibration MCCPs
  filter(str_detect(Molecule, "IS", negate = TRUE)) |> #Exclude the calibration IS
  filter(str_detect(Molecule, "RS", negate = TRUE)) #Exclude the calibration RS


filtered_dataB <- TESTING |>
  filter(`Sample Type` == "Standard") |> #Remove the samples etc. and only take the calibration std
  filter(`Isotope Label Type` == "Quan") |>  #Remove the qual ions
  filter((str_detect(`Replicate Name`, "^MCCP"))) |> #Take only the MCCPs
  filter(str_detect(Molecule, "C9", negate = TRUE))|>#Exclude the calibration vSCCPs
  filter(str_detect(Molecule, "C10", negate = TRUE))|>#Exclude the calibration MCCPs
  filter(str_detect(Molecule, "C11", negate = TRUE))|>#Exclude the calibration MCCPs
  filter(str_detect(Molecule, "C12", negate = TRUE))|>#Exclude the calibration MCCPs
  filter(str_detect(Molecule, "C13", negate = TRUE))|>#Exclude the calibration MCCPs
  filter(str_detect(Molecule, "IS", negate = TRUE)) |> #Exclude the calibration IS
  filter(str_detect(Molecule, "RS", negate = TRUE)) #Exclude the calibration RS


# Combine filtered data from A and B
filtered_data <- rbind(filtered_dataA, filtered_dataB)

# Initialize an empty list to store results for all standards
all_results <- list()

# Group the standards
list_of_std <- split(filtered_data, filtered_data$Note)
unique_std_names <- unique(filtered_data$Note)

for (std_name in unique_std_names) {
  std_df <- list_of_std[[std_name]]
  
  # Group the Molecules
  list_of_molecule <- split(std_df, std_df$Molecule)
  unique_molecule_names <- unique(std_df$Molecule)
  
  for (molecule_name in unique_molecule_names) {
    molecule_df <- list_of_molecule[[molecule_name]]
    
    # Build calibration curves
    cal <- lm(`Normalized Area` ~ `Analyte Concentration`, data = molecule_df)
    
    # Identify the carbon length and the chlorine content
    # Use a regular expression to extract the number of carbon atoms
    matches <- regmatches(molecule_name, regexec("C(\\d+)", molecule_name))
    carbon_count <- as.numeric(matches[[1]][2])
    
    # Assign the carbon count to the variable i
    i <- carbon_count
    
    # Use a regular expression to extract the number of chlorine atoms
    matches_cl <- regmatches(molecule_name, regexec("Cl(\\d+)", molecule_name))
    chlorine_count <- as.numeric(matches_cl[[1]][2])
    
    # Assign the chlorine count to the variable j
    j <- chlorine_count
    
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
    
    # Construct the Reference_standard
    replicate_name <- unique(molecule_df$`Replicate Name`)
    replicate_suffix <- sapply(strsplit(as.character(replicate_name), "_"), "[", 2)
    reference_standard <- paste0(type, "C", i, "Cl%", replicate_suffix, sep=" ")
    
    # Save results
    result <- data.frame(
      STD_code = paste0(std_name, "-", type),
      Reference_standard = reference_standard,
      Chain_length = paste0("C",i),
      Type = type,
      Homologue = molecule_name,
      Response_factor = coef(cal)[2],
      Intercept = coef(cal)[1],
      R_squared = summary(cal)$r.squared,
      stringsAsFactors = FALSE
    )
    
    # Append the result to the list
    all_results <- append(all_results, list(result))
  }
}

# Combine all results into a single data frame
cal_results <- do.call(rbind, all_results)
cal_results <- distinct(cal_results)
cal_results

##############################################################################################################
################################ PREPARE THE SAMPLES FOR CALIBRATION #####################################################
#############################################################################################################

# Load the data
TESTINGB <- TESTING |> 
  filter(`Isotope Label Type` == "Quan") |> 
  mutate(
    `Normalized Area` = as.numeric(`Normalized Area`),
    Chain_length = str_extract(`Molecule List`, "C\\d+"),
    Carbon_number = as.numeric(str_remove(Chain_length, "C")),
    Type = case_when(
      Carbon_number < 10 ~"vSCCPs",
      Carbon_number >= 10 & Carbon_number <= 13 ~"SCCPs",
      Carbon_number >= 14 & Carbon_number <= 17 ~"MCCPs",
      Carbon_number >= 18 & Carbon_number <= 30 ~"LCCPs",
      Carbon_number > 30 ~"vLCCPs") #Create the column type which should contain SCCPs, MCCPs and LCCPs depending on the chain length
  )|> 
  rename(Homologue = Molecule) |> 
  filter( `Sample Type` != "Standard" #Filter to remove the standards out of the data frame
  ) |> 
  filter(str_detect(Type, "vSCCPs", negate = TRUE))|>#Exclude the calibration vSCCPs
  filter(str_detect(Homologue, "IS", negate = TRUE)) |> #Exclude the calibration IS
  filter(str_detect(Homologue, "RS", negate = TRUE)) |>  #Exclude the calibration RS
  filter(str_detect(`Replicate Name`, "Std", negate = TRUE)) #Exclude the calibration standards

# Group samples for creating one data frame for each: the CALIBRATION script requires to have one data frame for each sample, so out of the current data frame we should create a new of for each of the samples
list_of_samples <- split(TESTINGB, TESTINGB$`Replicate Name`)



##############################################################################################################
########################################## CALIBRATION #####################################################
#############################################################################################################

#####GROUP STANDARD MIXTURES USED#####


##########################################PREPARE DATASET FOR PATTERN RECONSTRUCTION#################
{
  # Rename the file of the standards that we have created previous for calibration so it will fit the script
  input <- cal_results
  input$Type <- as.factor(input$Type)
  
  # Initialize an empty list to store results for all standards
  SCCP_MCCP_combinations <- list()
  
  # Group the standards
  input_list <- split(input, input$Type)
  SCCPs_MCCPs <- unique(input$Type)
  
  for (i in SCCPs_MCCPs) {
    SCCPs_MCCPs_df <- input_list[[i]]
    SCCPs_MCCPs_df$STD_code <- as.factor(SCCPs_MCCPs_df$STD_code)
    Combinations <- combn(x = levels(SCCPs_MCCPs_df$STD_code), m = 2, FUN = NULL, simplify = TRUE)
    SCCP_MCCP_combinations <- append(SCCP_MCCP_combinations, list(Combinations))
    
  }
  
  # Combine all results into a single data frame
  Combinations <- do.call(cbind, SCCP_MCCP_combinations)
  
  
  # Store sum RFs for each group CP standard
  input <- input |> 
    group_by(Reference_standard) |> 
    mutate(Sum_response_factor = sum(Response_factor, na.rm = TRUE))
  input[c(1:5)] <- lapply(input[c(1:5)], as.factor) 
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
for (sample_name in 1:length(unique_sample_names)) {
  # Get the corresponding sample data frame, now I will use sample_df for further procesing
  sample_df <- list_of_samples[[sample_name]]
  
  # Set sample name, to see which samples it will calibrate
  print(paste("Processing sample:", sample_name))
  
  
  ####################################RUN PATTERN RECONSTRUCTION FOR SELECTED (LOADED) SAMPLE####################
  #This section follows the structure that was in the original script
  
  
  ######################################################If I prepare combinations in columns 1:3 for SCCPs and 4:6 for MCCPs#############################################
  ####################################################Separating different  depending on conditioning #########################################
  ## The script was made for quantifying one sample at the time, so I changed it to iterate for all the data frames of the samples
  
  # Initialize an empty list to store results for all samples
  all_results <- list()
  all_plots <- list()
  
  # Extract unique sample names, it will give the to each "data frame"
  unique_sample_names <- unique(TESTINGB$`Replicate Name`)
  
  # Iterate over each unique sample name
  for (sample_name in 1:length(unique_sample_names)) {
    # Get the corresponding sample data frame, now I will use sample_df for further procesing
    sample_df <- list_of_samples[[sample_name]]
    
    # Set sample name, to see which samples it will calibrate
    print(paste("Processing sample:", sample_name))
    
    # Convert Type to factor and Normalized Area to numeric
    
    sample_df <- sample_df |> 
      mutate(
        Type = as.factor(Type),
        Chain_length = as.factor(Chain_length),
        `Normalized Area` = as.numeric(`Normalized Area`),
        Relative_distribution = `Normalized Area` / sum(`Normalized Area`, na.rm = TRUE)
      )
    
    # Calculate relative 'Normalized Area' distribution within each homologue group
    sample_df$Relative_distribution <- NA
    sample_df$`Normalized Area`[is.na(sample_df$`Normalized Area`)] <- 0
    
    # Calculate relative 'Normalized Area' distribution within each homologue group in SCPPs and in MCCPs separately
    sample_df <- sample_df  |>  
      group_by(Type) |>  
      mutate(Relative_distribution = `Normalized Area` / sum(`Normalized Area`, na.rm = TRUE))
    
    results <- sample_df
    results[c("Comp_1", "Comp_2", "Fraction_Comp_1", "Simulated_pattern")] <- NA
    
    #Type as factor in input
    input <- input |> 
      mutate(Type = as.factor(Type))
    
    # Deconvolution of homologue patterns
    
    for (i in 1:length(sample_df$Type)) {
      REF <- sample_df$Relative_distribution[sample_df$Type == levels(sample_df$Type)[i]]
      Distance <- 100
      
      for (z in 1:length(Combinations[1, ])) { 
        
        C_1 <- subset(input, subset = (STD_code == Combinations[1, z] & Type == sample_df$Type[i]))
        C_2 <- subset(input, subset = (STD_code == Combinations[2, z] & Type == sample_df$Type[i]))
        
        for (j in 1:100) {
          Combo <- (C_1$Response_factor * j + C_2$Response_factor * (100 - j)) / sum((C_1$Response_factor * j + C_2$Response_factor * (100 - j)), na.rm = TRUE)
          
          if (Distance > sum(sqrt((REF - Combo)^2))) {
            results$Comp_1[results$Type == levels(sample_df$Type)[i]] <- as.character(C_1$STD_code)
            results$Comp_2[results$Type == levels(sample_df$Type)[i]] <- as.character(C_2$STD_code)
            results$Fraction_Comp_1[results$Type == levels(sample_df$Type)[i]] <- j
            results$Simulated_pattern[results$Type == levels(sample_df$Type)[i]] <- Combo
            Distance <- sum(sqrt((REF - Combo)^2))
          }
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
}
#########################################VIEW RESULTS############################################

#View Overview
print(all_results_df_MCCPs)

#View graph, REPLACE THE NAME OF THE SAMPLE YOU WANT TO SEE, another option is to open the list from the Environment and from there open each plot
all_plots[["NIST_R1"]]
all_plots[["NIST_R2"]]
all_plots[["NIST_R3"]]
