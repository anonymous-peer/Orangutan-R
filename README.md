# Orangutan

**Orangutan** is an R package designed for analyzing and visualizing phenotypic data in the context of species descriptions or, more broadly, comparing groups of individuals. The package automates group comparisons by generating summary statistics, identifying non-overlapping variables, running multivariate and univariate analyses.

## Features

- **Analysis Pipeline**: Runs the complete pipeline in one function.
- **Outlier Removal**: To remove individuals with extreme values (optional).
- **Summary Statistics**: Includes mean, standard deviation, min, and max values for each variable.
- **Multivariate Analyses**:  
  - Principal Component Analysis (PCA)  
  - Discriminant Analysis of Principal Components (DAPC)  
  - Cross-classification table
- **Univariate Tests** (_the appropriate test is selected automatically_):
  - ANOVA with Tukey's HSD as post hoc Test
  - Kruskal-Wallis with Dunn's as post hoc Test
- **Interactive Console Prompts**: Dynamically handles user input for dataset paths, outlier removal, species encircling in plots, and allometric transformations.

---

## Installation

Install the package from GitHub using the `devtools` package:  
   ```R
   # Install devtools if not already installed
   install.packages("devtools")
   ```
   ```R
   # Install Orangutan
   devtools::install_github("anonymous-peer/Orangutan-R")
   ```

---

## Usage

### 1. Running the Orangutan Analysis Pipeline

The primary function of the package is `run_orangutan()`, which executes the full analysis pipeline. 
Once the function is called, it will guide the user through the steps interactively via console prompts.

```R
# Execute the analysis pipeline
run_orangutan()
```

### 2. Console Prompts

When the `run_orangutan()` function is executed, the user will be prompted to provide the following inputs:

1. **Dataset Path**:  
   - A CSV file containing the phenotypic data. Ensure the dataset includes a `species` column and numeric variables with no missing data.  
   - Example:  
     ```
     Enter dataset (path/to/file_name.csv): path/to/dataset.csv
     ```
2. **Outlier Removal**:  
   - The user will be asked whether to remove outliers from the dataset.  
   - Example:  
     ```
     Do you want to remove outliers? [yes/no]: yes
     ```
   - If you answer "yes", you will be prompted to further specify which variable(s) to process, and the lower and upper cutoff value for removal.
   - Example:  
     ```
     Enter variable names to clean outliers from (comma-separated): variable1, variable2
     ```
     ```
     Enter tail proportion for outlier removal (e.g., 0.05 for 5%): 0.05
     ```
3. **Allometric Transformation**:  
   - The user will be asked whether to apply allometric transformation (adjusts variables based on body size).
   - Example:  
     ```
     Do you want to apply allometric transformation? (yes/no): yes
     ```
4. **Species for Encircling**:  
   - For the PCA and DAPC plots, optionally specify species to highlight by encircling their data points.  
   - Make sure species names have no trailing spaces in the input dataset (e.g., "setosa" is not the same as "setosa ").  
   - Example:  
     ```
     Enter species to encircle (comma-separated): species1, species2
     ```
     or
     ```
     Enter species to encircle (comma-separated): species1,species2
     ```

### 3. Output Files

The analysis generates the following output files in the directory of the dataset:

- **Summary Statistics**: `summary_statistics.csv`
- **Outlier Removal** (if ran): `cleaned_data_outliers_removed.csv`
- **PCA Plot**: `plot_multivar_pca.tiff`
- **DAPC Plot**: `plot_multivar_dapc.tiff`
- **Confusion Table**: `confusion_table.csv`
- **ANOVA Summary**: `summary_anova.csv`
- **Kruskal-Wallis Summary**: `summary_kruskalwallis.csv`
- **Additional Plots**: plots for non-overlapping variables, ANOVA and/or Kruskal-Wallis results.

---

## Example Dataset Format

The dataset should be a CSV file with the following structure:

### If there is no need for allometric transformations:
| species    | variable1 | variable2 | variable3 | ... |
|------------|-----------|-----------|-----------|-----|
| species1   | value1    | value2    | value3    | ... |
| species2   | value4    | value5    | value6    | ... |

### If allometric transformations are needed:
| species    | main_length | variable2 | variable3 | ... |
|------------|-------------|-----------|-----------|-----|
| species1   | value1      | value2    | value3    | ... |
| species2   | value4      | value5    | value6    | ... |

Tested example datasets can be found in the [`example_datasets/`](example_datasets) directory of this repository.

---

## Notes and Warnings

- **Color Assignments**: If more than 12 species are present, color assignments will repeat due to palette limitations.
- **Data Requirements**:  
  - The `species` column is mandatory, and numeric variables are required for statistical analyses.  
  - The `main_length` column is mandatory if allometric adjustment is needed.
- **Error Handling**: The function will stop execution if required columns are missing or if invalid data is provided.

---

## Recommended Citation

###. (2025). Orangutan: an R package for analyzing and visualizing phenotypic data in the context of ecology and systematics. GitHub repository: [https://github.com/metalofis/Orangutan-R].
