utils::globalVariables(c(
  "X", "species", "Variable", "Summary", "value", "PC1", "PC2",
  "LD1", "LD2", "treatment", "Letters"
))

#' Run Orangutan
#'
#' This function runs the full Orangutan analysis pipeline.
#'
#' @importFrom dplyr %>% select arrange filter all_of add_row group_by group_split pull bind_rows mutate
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom ggplot2 ggplot aes geom_boxplot labs scale_fill_manual theme_minimal theme element_text ggsave geom_point scale_color_manual element_blank element_rect geom_polygon geom_text
#' @importFrom RColorBrewer brewer.pal
#' @importFrom adegenet dapc
#' @importFrom stats aov sd lm residuals prcomp shapiro.test bartlett.test TukeyHSD as.formula kruskal.test setNames
#' @importFrom utils read.csv write.csv
#' @importFrom grDevices chull
#' @importFrom methods Summary
#' @export

run_orangutan <- function() {
  
  # ------------------- i. Data -------------------
  data_path <- readline(prompt = "Enter dataset (path/to/file_name.csv): ")
  data <- read.csv(data_path)
  output_dir <- dirname(data_path)
  if ("X" %in% colnames(data)) data <- data %>% select(-X)
  cleaned_data <- data %>% arrange(species)

  # ------------------- ii. Outlier removal -------------------
  remove_outliers <- readline("Do you want to remove outliers? [yes/no]: ")
  
  if (tolower(remove_outliers) == "yes") {
    vars_input <- readline("Enter variable names to clean outliers from (comma-separated): ")
    vars <- trimws(unlist(strsplit(vars_input, ",")))
    tail_pct <- as.numeric(readline("Enter tail proportion for outlier removal (e.g., 0.05 for 5%): "))
    for (v in vars) {
      if (v %in% colnames(cleaned_data)) {
        Q_low <- quantile(cleaned_data[[v]], tail_pct, na.rm = TRUE)
        Q_high <- quantile(cleaned_data[[v]], 1 - tail_pct, na.rm = TRUE)
        IQR <- Q_high - Q_low
        lower <- Q_low - 1.5 * IQR
        upper <- Q_high + 1.5 * IQR
        cleaned_data <- cleaned_data[cleaned_data[[v]] >= lower & cleaned_data[[v]] <= upper, ]
      } else {
        cat(sprintf("Variable '%s' not found in dataset.\n", v))
      }
    }
    # Save the new cleaned_data as a CSV
    write.csv(cleaned_data, file = file.path(output_dir, "cleaned_data_outliers_removed.csv"), row.names = FALSE)
  }
  
  # ------------------- iii. Automating colors -------------------
  if (!"species" %in% colnames(cleaned_data)) stop("The 'species' column is not present in the dataset.")
  palette_name <- "Paired"
  max_palette_colors <- 12
  num_colors <- length(unique(cleaned_data$species))
  colors <- brewer.pal(min(num_colors, max_palette_colors), palette_name)
  if (num_colors > max_palette_colors) {
    colors <- rep(colors, length.out = num_colors)
    message("Warning: More than ", max_palette_colors, " species detected. Color assignments will repeat.")
  }
  unique_species <- unique(cleaned_data$species)
  custom_colors <- setNames(colors, unique_species)
  print(custom_colors)
  
  # ------------------- 1. Summary Statistics -------------------
  variables <- colnames(cleaned_data)[sapply(cleaned_data, is.numeric)]
  summary_stats <- data.frame(Species = character(), Variable = character(), Summary = character(), stringsAsFactors = FALSE)
  for (var in variables) {
    for (sp in unique(cleaned_data$species)) {
      subset_data <- cleaned_data %>% filter(species == sp) %>% select(all_of(var))
      mean_val <- mean(subset_data[[var]], na.rm = TRUE)
      std_dev <- sd(subset_data[[var]], na.rm = TRUE)
      min_val <- min(subset_data[[var]], na.rm = TRUE)
      max_val <- max(subset_data[[var]], na.rm = TRUE)
      summary_string <- paste0(round(mean_val, 2), " \u00b1 ", round(std_dev, 2), 
                               " (", round(min_val, 2), "-", round(max_val, 2), ")")
      summary_stats <- summary_stats %>% add_row(Species = sp, Variable = var, Summary = summary_string)
    }
  }
  summary_stats_reshaped <- summary_stats %>% pivot_wider(names_from = Variable, values_from = Summary)
  write.csv(summary_stats_reshaped, file = file.path(output_dir, "summary_statistics.csv"), row.names = FALSE)
  
  # ------------------- 2. Identification of Non-Overlapping Variables and Plotting ------------------- 
  # Create a function to check overlap of min and max values across species pairs
  check_no_overlap <- function(species1, species2, variable) {
    # Extract values for the variable for each species
    sp1_values <- cleaned_data %>% filter(species == species1) %>% pull(variable)
    sp2_values <- cleaned_data %>% filter(species == species2) %>% pull(variable)
    
    # Get min and max values for each species
    sp1_min <- min(sp1_values, na.rm = TRUE)
    sp1_max <- max(sp1_values, na.rm = TRUE)
    sp2_min <- min(sp2_values, na.rm = TRUE)
    sp2_max <- max(sp2_values, na.rm = TRUE)
    
    # Check if the max of one species is less than the min of the other species
    no_overlap <- sp1_max < sp2_min || sp2_max < sp1_min
    return(no_overlap)
  }
  
  # Create a dataframe to store results
  no_overlap_results <- data.frame(species1 = character(0), species2 = character(0), variable = character(0), no_overlap = logical(0))
  
  # Iterate through all species pairs and variables to check no overlap
  for (var in variables) {
    species_list <- unique(cleaned_data$species)
    
    for (i in 1:(length(species_list) - 1)) {
      for (j in (i + 1):length(species_list)) {
        sp1 <- species_list[i]
        sp2 <- species_list[j]
        
        no_overlap <- check_no_overlap(sp1, sp2, var)
        
        # Store the result
        no_overlap_results <- rbind(no_overlap_results, data.frame(species1 = sp1, species2 = sp2, variable = var, no_overlap = no_overlap))
      }
    }
  }
  
  # Filter the species pairs that do not overlap
  non_overlapping_pairs <- no_overlap_results %>% filter(no_overlap == TRUE)
  
  # Check if there are no non-overlapping pairs
  if (nrow(non_overlapping_pairs) == 0) {
    cat("All variables overlap for the species pairs. No plots were produced.\n")
  } else {
    # Reshape the data for plotting
    cleaned_data_long <- cleaned_data %>%
      pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
    
    # Ensure species are ordered according to the custom_colors vector
    species_order <- unique(cleaned_data$species)
    
    # Create a named vector that maps species to colors from custom_colors
    species_colors <- setNames(custom_colors[1:length(species_order)], species_order)
    
    # Loop through each non-overlapping pair and plot
    for (i in seq_len(nrow(non_overlapping_pairs))) {
      pair <- non_overlapping_pairs[i, ]
      
      # Filter data for the current pair and variable
      pair_data <- cleaned_data_long %>%
        filter(
          species %in% c(pair$species1, pair$species2) &
            variable == pair$variable
        )
      
      # Create the plot
      plot <- ggplot(pair_data, aes(x = species, y = value, fill = species)) +
        geom_boxplot() +
        labs(
          title = paste("Non-Overlapping Pair:", pair$species1, "vs", pair$species2, "for", pair$variable),
          x = "Species",
          y = pair$variable
        ) +
        scale_fill_manual(values = species_colors) +  # Apply the correct species colors
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)
        )
      
      # Display the plot
      print(plot)
      
      # Save the plot
      ggsave(filename = file.path(output_dir, paste0("plot_non_overlap_", pair$species1, "_vs_", pair$species2, "_", pair$variable, ".tiff")), 
      plot = plot, width = 7.5, height = 6, dpi = 300, device = "tiff")
    }
  }
  
  # ------------------- iii. Allometric transformation -------------------
  # Prompt user to decide whether to apply allometric transformation
  apply_allometric <- tolower(readline("Do you want to apply allometric transformation? (yes/no): "))
  
  if (apply_allometric == "yes") {
    message("Applying allometric transformation...")
    
    # 1. Dynamically identify the main_length column
    main_length_column <- grep("main_length", colnames(cleaned_data), ignore.case = TRUE, value = TRUE)
    
    # If main_length column is not found, stop the process
    if (length(main_length_column) == 0) {
      stop("main_length column not found in the dataset!")
    }
    
    # Extract main_length data
    main_length <- cleaned_data[[main_length_column]]
    
    # 2. Dynamically identify mensural variables based on decimal numbers, excluding main_length
    mensural_columns <- sapply(names(cleaned_data), function(col_name) {
      # Check if the column is not main_length and contains numeric values with decimals
      if (!grepl("main_length", col_name, ignore.case = TRUE)) {
        col <- cleaned_data[[col_name]]  # Access the column data
        
        # Check if the column is numeric and contains decimal numbers
        is_numeric <- is.numeric(col)
        has_decimal <- any(grepl("\\.", col))  # Check for numbers with decimals
        
        # Return TRUE if the column is numeric, has decimals, and is not main_length
        return(is_numeric && has_decimal)
      }
      return(FALSE)  # Return FALSE if it's main_length
    })
    
    # Extract mensural variables (excluding main_length)
    mensural_vars <- cleaned_data[, mensural_columns]
    
    # 3. Apply allometric adjustment
    log_main_length <- log(main_length)  # Log-transformed main_length
    log_mensural_vars <- log(mensural_vars)  # Log-transformed mensural variables
    
    # Perform allometric adjustment using regression on each mensural variable
    adjusted_data <- apply(log_mensural_vars, 2, function(x) {
      model <- lm(x ~ log_main_length)  # Linear regression with log-transformed main_length as the predictor
      residuals(model)  # Extract residuals from the model
    })
    
    # Convert the result into a data frame
    adjusted_data <- as.data.frame(adjusted_data)
    
    # 4. Combine the adjusted data with the species classification
    adjusted_data_with_species <- cbind(cleaned_data[, 1], adjusted_data)  # Assuming the first column is 'Species'
    colnames(adjusted_data_with_species) <- c("species", colnames(adjusted_data))
    
    # Add the non-mensural variables (FALSE in mensural_columns) to adjusted_data_with_species
    non_mensural_vars <- cleaned_data[, !mensural_columns]  # Select columns that are NOT mensural variables
    non_mensural_vars <- non_mensural_vars[, !grepl("species", colnames(non_mensural_vars), ignore.case = TRUE)]  # Exclude species
    
    # Concatenate the non-mensural variables to the adjusted data (excluding the species column)
    cleaned_allometric_data <- cbind(adjusted_data_with_species, non_mensural_vars)
    
    message("Now cleaned_allometric_data contains species, mensural (adjusted), and meristic variables.")
  } else {
    message("Skipping allometric transformation. Using original cleaned_data.")
    cleaned_allometric_data <- cleaned_data
  }
  
  # ------------------- 3.1 Multivariate Analysis: PCA -------------------
  # Keep only numeric columns for PCA
  variables_matrix <- as.matrix(cleaned_allometric_data[, sapply(cleaned_allometric_data, is.numeric)])
  
  # Perform PCA
  pca_variables <- prcomp(variables_matrix, center = TRUE, scale. = TRUE)
  
  # Create a dataframe with PCA results
  pca_df <- as.data.frame(pca_variables$x)
  
  # Add species column to the PCA dataframe
  pca_df$species <- cleaned_allometric_data$species
  
  # Calculate explained variance
  explained_variance <- round(pca_variables$sdev^2 / sum(pca_variables$sdev^2) * 100, 2)
  
  # Dynamically get top contributing variables to PC1
  loading_scores <- pca_variables$rotation[, 1]
  absolute_values <- abs(loading_scores)
  variable_score_ranked <- sort(absolute_values, decreasing = TRUE)
  top_n <- min(10, length(variable_score_ranked))
  top_vars <- names(variable_score_ranked[1:top_n])
  
  # Display loading scores of top variables
  cat("Top contributing variables to PC1:\n")
  print(pca_variables$rotation[top_vars, 1])
  
  # Prompt user for species to encircle
  species_input <- readline(prompt = "Enter species to encircle (comma-separated): ")
  # Split the input into a vector of species names
  species_to_encircle <- strsplit(species_input, ",")[[1]] %>% 
    trimws()  # Remove any leading or trailing whitespace from each species name
  
  # Print the selected species
  cat("Species to encircle:", species_to_encircle, "\n")
  
  # ------------------- 3.1 Multivariate Analysis: PCA -------------------
  # Prepare convex hulls for PCA if species_to_encircle is not NULL or empty
  if (!is.null(species_to_encircle) && length(species_to_encircle) > 0) {
    pca_df_polygon <- pca_df %>%
      filter(species %in% species_to_encircle) %>%
      group_by(species) %>%
      group_split() %>%
      lapply(function(df) {
        df[chull(df$PC1, df$PC2), ] %>%
          rbind(df[chull(df$PC1, df$PC2)[1], ])  # Close the polygon
      }) %>%
      bind_rows()
  } else {
    pca_df_polygon <- NULL
  }
  
  # PCA plot
  ggsave(filename = file.path(output_dir, "plot_multivar_pca.tiff"), 
  plot = (ggplot(pca_df, aes(x = PC1, y = PC2, color = species)) +
      {if (!is.null(pca_df_polygon)) 
        geom_polygon(data = pca_df_polygon, aes(x = PC1, y = PC2, color = species),
                     fill = "lightgrey", alpha = 0.2, linewidth = 0.5)} +
      geom_point(size = 3.5) +
      scale_color_manual(
        breaks = sort(unique(pca_df$species)),
        values = custom_colors
      ) +
      labs(
        x = paste("PC1 (", explained_variance[1], "%)", sep = ""),
        y = paste("PC2 (", explained_variance[2], "%)", sep = "")
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 1),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
  ), width = 7.5, height = 6, dpi = 300, device = "tiff", bg = "white")
  
  # ------------------- 3.2 Multivariate Analysis: DAPC -------------------
  # Convert species column to factor
  cleaned_allometric_data$species <- as.factor(cleaned_allometric_data$species)
  
  # Perform DAPC using the same variables_matrix from PCA
  dapc_result <- dapc(variables_matrix, 
                      n.pca = 10, 
                      grp = cleaned_allometric_data$species, 
                      n.da = 5, 
                      jackknife = TRUE)
  
  # Extract LD scores and add species info
  LD_scores <- as.data.frame(dapc_result$ind.coord)
  LD_scores$species <- cleaned_allometric_data$species
  
  # Extract eigenvalues and calculate variance explained by LD1 and LD2
  eigenvalues <- dapc_result$eig
  var_explained <- eigenvalues[1:2] / sum(eigenvalues)
  
  # Prepare convex hulls for DAPC if species_to_encircle is not NULL or empty
  if (!is.null(species_to_encircle) && length(species_to_encircle) > 0) {
    LD_scores_polygon <- LD_scores %>%
      filter(species %in% species_to_encircle) %>%
      group_by(species) %>%
      group_split() %>%
      lapply(function(df) {
        df[chull(df$LD1, df$LD2), ] %>%
          rbind(df[chull(df$LD1, df$LD2)[1], ])  # Close the polygon
      }) %>%
      bind_rows()
  } else {
    LD_scores_polygon <- NULL
  }
  
  # Sort species factor levels alphabetically
  LD_scores$species <- factor(LD_scores$species, levels = sort(unique(LD_scores$species)))
  if (!is.null(LD_scores_polygon) && nrow(LD_scores_polygon) > 0) {
    LD_scores_polygon$species <- factor(LD_scores_polygon$species, levels = sort(unique(LD_scores_polygon$species)))
  }
  
  # DAPC plot
  ggsave(filename = file.path(output_dir, "plot_multivar_dapc.tiff"),
   plot = (ggplot(LD_scores, aes(x = LD1, y = LD2, color = species)) +
      {if (!is.null(LD_scores_polygon)) 
        geom_polygon(data = LD_scores_polygon, aes(x = LD1, y = LD2, color = species),
                     fill = "lightgrey", alpha = 0.2, linewidth = 0.5)} +
      geom_point(size = 3.5) +
      scale_color_manual(
        breaks = sort(unique(LD_scores$species)),
        values = custom_colors
      ) +
      labs(
        x = paste("LD1 (", round(var_explained[1] * 100, 2), "%)"),
        y = paste("LD2 (", round(var_explained[2] * 100, 2), "%)")
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 1),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
  ), width = 7.5, height = 6, dpi = 300, device = "tiff", bg = "white")
  
  # ------------------- 3.3 Confusion table for DAPC -------------------
  predicted_groups <- as.factor(dapc_result$assign)
  conf_table <- table(cleaned_allometric_data$species, predicted_groups)
  write.csv(conf_table, file = file.path(output_dir, "confusion_table.csv"), row.names = TRUE)
  
  # ------------------- 4. Univariate Analyses (ANOVA and Tukey HSD, Kruskal-Wallis and Dunn's test) -------------------
  # Initialize an empty vector to store variables that failed ANOVA assumptions
  failed_anova_assumptions <- c()
  
  # Initialize an empty list to store results
  anova_results_table <- data.frame(
    Variable = character(),
    F_value = numeric(),
    DF_between = numeric(),
    DF_within = numeric(),
    P_value = numeric(),
    Shapiro_p = numeric(),
    Bartlett_p = numeric(),
    Assumptions_Met = character(),
    stringsAsFactors = FALSE
  )
  
  # Function to generate Tukey HSD labels for post-hoc analysis
  generate_label_df <- function(TUKEY, variable) {
    tukey_result <- TUKEY[[variable]][[1]]
    Tukey.labels <- data.frame(multcompView::multcompLetters(tukey_result[, "p adj"])['Letters'])
    Tukey.labels$treatment <- rownames(Tukey.labels)
    Tukey.labels <- Tukey.labels[order(Tukey.labels$treatment), ]
    return(Tukey.labels)
  }
  
  # Create an empty list to store Tukey HSD results
  TUKEY <- list()
  
  # Function to conduct ANOVA and Tukey HSD test, and store variables failing assumptions
  check_anova_assumptions <- function(data, variable) {
    anova_result <- aov(data[[variable]] ~ data$species)
    shapiro_test <- shapiro.test(residuals(anova_result))
    bartlett_test <- bartlett.test(data[[variable]], data$species)
    anova_summary <- summary(anova_result)
    
    f_statistic <- anova_summary[[1]][["F value"]][[1]]
    df_between <- anova_summary[[1]][["Df"]][[1]]
    df_within <- anova_summary[[1]][["Df"]][[2]]
    p_value <- anova_summary[[1]][["Pr(>F)"]][[1]]
    
    cat("Variable:", variable, "\n")
    cat("F:", f_statistic, "\n")
    cat("Degrees of Freedom (Between):", df_between, "\n")
    cat("Degrees of Freedom (Within):", df_within, "\n")
    cat("P-value:", p_value, "\n")
    cat("Shapiro-Wilk normality test p-value:", shapiro_test$p.value, "\n")
    cat("Bartlett test p-value:", bartlett_test$p.value, "\n")
    
    assumptions_met <- if (shapiro_test$p.value > 0.05 & bartlett_test$p.value > 0.05) {
      cat("ANOVA assumptions met for", variable, "\n\n")
      "Yes"
    } else {
      cat("ANOVA assumptions NOT met for", variable, "\n\n")
      failed_anova_assumptions <<- c(failed_anova_assumptions, variable)
      "No"
    }
    
    # Store results into the table
    anova_results_table <<- rbind(anova_results_table, data.frame(
      Variable = variable,
      F_value = f_statistic,
      DF_between = df_between,
      DF_within = df_within,
      P_value = p_value,
      Shapiro_p = shapiro_test$p.value,
      Bartlett_p = bartlett_test$p.value,
      Assumptions_Met = assumptions_met,
      stringsAsFactors = FALSE
    ))
    
    # Post-hoc if significant
    if (p_value < 0.05) {
      posthoc_result <- TukeyHSD(anova_result)
      TUKEY[[variable]] <<- posthoc_result
    } else {
      cat("No significant difference found for", variable, "\n\n")
    }
  }
  
  # Automatically select numeric columns and apply function
  numeric_columns <- sapply(cleaned_allometric_data, is.numeric)
  variable_names <- names(cleaned_allometric_data)[numeric_columns]
  
  for (variable in variable_names) {
    check_anova_assumptions(cleaned_allometric_data, variable)}
  
  # Save results to a CSV file
  write.csv(anova_results_table, file = file.path(output_dir, "summary_anova.csv"), row.names = FALSE)
  
  # Plotting ANOVA Results with Tukey HSD Labels
  # Loop over variables that passed ANOVA assumptions
  anova_variables <- setdiff(variable_names, failed_anova_assumptions)
  
  for (variable in anova_variables) {
    # Check if a Tukey HSD result exists for the variable
    if (is.null(TUKEY[[variable]])) {
      cat("Skipping Tukey labels for", variable, "because no TukeyHSD result was stored.\n")
      next
    }
    
    # Generate Tukey HSD labels
    LABELS <- generate_label_df(TUKEY, variable)
    
    # Create the base boxplot using ggpubr
    p <- ggpubr::ggboxplot(cleaned_allometric_data,
                           x = "species", y = variable,
                           color = "species",
                           legend = "none", 
                           add = "jitter") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10)
        ) +
        scale_color_manual(values = custom_colors)
    
    # Calculate y-position for the labels
    max_y <- max(cleaned_allometric_data[[variable]], na.rm = TRUE)
    over <- 0.1 * diff(range(cleaned_allometric_data[[variable]], na.rm = TRUE))
    
    # Add Tukey HSD labels to the plot
    p <- p + geom_text(data = LABELS,  
                       aes(x = as.numeric(factor(treatment, 
                                                 levels = unique(cleaned_allometric_data$species))), 
                           y = max_y + over, 
                           label = Letters),
                       size = 6, color = "black")
    
    print(p)
    ggsave(filename = file.path(output_dir, paste0("plot_anova_", variable, ".tiff")), 
           plot = p, width = 7.5, height = 6, dpi = 300, device = "tiff")
  }
  
  # Kruskal-Wallis and Dunn's Test for Variables Failing ANOVA Assumptions #######
  # Now use the dynamically stored variables that failed ANOVA assumptions for Kruskal-Wallis and Dunn's tests
  variables <- failed_anova_assumptions
  
  # Generating Dunn labels based on significant pairwise comparisons
  generate_dunn_label_df <- function(dunn_test, species) {
    # Create a matrix of pairwise p-values
    pairwise_p <- matrix(1, nrow = length(species), ncol = length(species))
    rownames(pairwise_p) <- colnames(pairwise_p) <- species
    
    # Populate the matrix with adjusted p-values from the Dunn test
    comparisons <- dunn_test$comparisons
    p_adjusted <- dunn_test$P.adjusted
    
    for (i in seq_along(comparisons)) {
      pair <- unlist(strsplit(comparisons[i], " - "))
      pairwise_p[pair[1], pair[2]] <- p_adjusted[i]
      pairwise_p[pair[2], pair[1]] <- p_adjusted[i]
    }
    
    # Generate labels using multcompLetters
    dunn_letters <- multcompView::multcompLetters(pairwise_p, threshold = 0.05)
    label_df <- data.frame(Letters = dunn_letters$Letters, species = names(dunn_letters$Letters))
    
    return(label_df)
  }
  
  # Initialize a data frame to store the summary of Dunn's test
  summary_table <- data.frame(
    Variable = character(),
    Kruskal_p_value = numeric(),
    Kruskal_Chi_Squared = numeric(),
    Significant_Pairs = character(),
    stringsAsFactors = FALSE
  )
  
  # Conduct Kruskal-Wallis and Dunn's post hoc tests on the dynamically selected variables
  dunn_labels_list <- list()  # To store Dunn's labels for each variable
  
  for (var in variables) {
    # Conduct Kruskal-Wallis test
    kruskal_test <- kruskal.test(as.formula(paste(var, "~ species")), data = cleaned_allometric_data)
    
    # Add Kruskal-Wallis results to the summary table
    summary_table <- summary_table %>%
      add_row(Variable = var, 
              Kruskal_p_value = kruskal_test$p.value,
              Kruskal_Chi_Squared = kruskal_test$statistic)
    
    # Perform Dunn's test if Kruskal-Wallis is significant
    if (kruskal_test$p.value < 0.05) {
      dunn_test <- dunn.test::dunn.test(cleaned_allometric_data[[var]], cleaned_allometric_data$species, method = "bonferroni")
      
      # Generate Dunn's test labels
      dunn_labels <- generate_dunn_label_df(dunn_test, unique(cleaned_allometric_data$species))
      dunn_labels_list[[var]] <- dunn_labels
      
      # Store significant pairs
      significant_pairs <- paste(dunn_test$comparisons[dunn_test$P.adjusted < 0.05], collapse = "; ")
      summary_table$Significant_Pairs[summary_table$Variable == var] <- significant_pairs
    } else {
      summary_table$Significant_Pairs[summary_table$Variable == var] <- "None"
    }
  }
  
  # Save the summary table
  write.csv(summary_table, file = file.path(output_dir, "summary_kruskalwallis.csv"), row.names = FALSE)
  
  # Plotting with Dunn Labels
  # Loop over variables to generate plots for Dunn's test labels
  for (var in variables) {
    if (!is.null(dunn_labels_list[[var]])) {
      dunn_labels <- dunn_labels_list[[var]]
      
      # Create boxplot with ggplot2
      p <- ggpubr::ggboxplot(cleaned_allometric_data, 
                             x = "species", y = var, 
                             color = "species", 
                             legend = "none", 
                             add = "jitter") + 
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10)
        ) + 
        scale_color_manual(values = custom_colors)
      
      # Calculate y-position for the labels
      max_y <- max(cleaned_allometric_data[[var]], na.rm = TRUE)
      over <- 0.1 * diff(range(cleaned_allometric_data[[var]], na.rm = TRUE))
      
      # Add Dunn's test labels to the plot
      p <- p + geom_text(data = dunn_labels, 
                         aes(x = as.numeric(factor(species, levels = unique(cleaned_allometric_data$species))), 
                             y = max_y + over, 
                             label = Letters), 
                         size = 6, color = "black")
      
      print(p)
      ggsave(filename = file.path(output_dir, paste0("plot_kruskalwallis_", var, ".tiff")), 
             plot = p, width = 7.5, height = 6, dpi = 300, device = "tiff")
    }
  }
  
message("\nOrangutan run completed.\nAll output files are saved in:\n", output_dir)
}
