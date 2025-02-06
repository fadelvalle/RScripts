# Key Steps:
# 1. Load necessary packages: ggplot2, dplyr, and tidyverse for data manipulation and visualization.
# 2. Manually select two CSV files ("new.csv" and "old.csv") using file.choose().
# 3. Read the CSV files into data frames.
# 4. Merge the data frames by the shared column "SYMBOL" and rename relevant columns for clarity.
# 5. Reshape the data into long format using pivot_longer() for easier plotting with ggplot.
# 6. Add a new column for adjusted p-values (adjP) depending on the dataset.
# 7. Remove rows with missing logFC values.
# 8. Add a new column for color coding based on significance (p < 0.05).
# 9. Plot the data using ggplot with customized aesthetics for color, shape, and axes.
# 10. Add a reference line to the plot for comparison (at y = -1.3).
# 11. Display the plot.

# Load necessary packages
library(ggplot2)   # for creating plots
library(dplyr)      # for data manipulation
library(tidyverse)  # for data manipulation and visualization

# Manually select files (these will open file dialogs to choose the files)
new_file <- file.choose() # Select the "new.csv" file
old_file <- file.choose() # Select the "old.csv" file

# Read the CSV files
new_data <- read.csv(new_file)   # Read data from "new.csv"
old_data <- read.csv(old_file)   # Read data from "old.csv"

# Combine the data by the "SYMBOL" column, renaming columns for clarity
combined_data <- full_join(new_data %>% rename(logFC_new = logFC, adjP_new = adj.P.Val),
                           old_data %>% rename(logFC_old = logFC, adjP_old = adj.P.Val),
                           by = "SYMBOL")  # Join the datasets by the SYMBOL column

# Convert data to long format for ggplot
plot_data <- combined_data %>%
  pivot_longer(cols = c(logFC_new, logFC_old), 
               names_to = "Dataset", 
               values_to = "logFC") %>%  # Reshape the data for plotting
  mutate(adjP = ifelse(Dataset == "logFC_new", adjP_new, adjP_old),  # Assign the correct p-value column based on the dataset
         Dataset = factor(Dataset, levels = c("logFC_new", "logFC_old"))) %>%  # Set the order of the 'Dataset' factor
  filter(!is.na(logFC))  # Remove rows where logFC is missing

# Add a color column based on the significance of the p-value
plot_data <- plot_data %>%
  mutate(Color = ifelse(adjP < 0.05, "p < 0.05", "p > 0.05"))  # Color by p-value significance

# Create a plot comparing the logFC values from both datasets
genes_plot <- ggplot(plot_data, aes(x = SYMBOL, y = logFC, color = Color, shape = Dataset)) +
  geom_point(size = 3, alpha = 0.8) +  # Plot points with custom size and transparency
  scale_color_manual(values = c("p < 0.05" = "green", "p > 0.05" = "red")) +  # Customize color scale for significance
  scale_shape_manual(values = c("logFC_new" = 15, "logFC_old" = 16)) +  # Customize shapes for new vs. old datasets (square and circle)
  labs(title = "Comparison of logFC between New and Old",
       x = "SYMBOL",  # x-axis label
       y = "logFC",   # y-axis label
       color = "Significance",  # legend title for color
       shape = "Dataset") +  # legend title for shape
  theme_minimal() +  # Apply minimal theme to the plot
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  # Rotate x-axis labels for better visibility
  geom_hline(yintercept = -1.3, linetype = "dashed", color = "blue", size = 1)  # Add a reference line at y = -1.3
genes_plot  # Display the plot
