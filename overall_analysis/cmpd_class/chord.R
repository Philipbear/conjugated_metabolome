# Install required packages if not already installed
if (!require("circlize")) install.packages("circlize")
if (!require("RColorBrewer")) install.packages("RColorBrewer")

# Load libraries
library(circlize)
library(RColorBrewer)

# Read the combined summary data
combined_data <- read.csv("/Users/shipei/Documents/projects/conjugated_metabolome/overall_analysis/cmpd_class/data/combined_summary.csv")

# Get unique classes in the specified order, excluding 'Others'
class_order <- c(
  'Fatty acids',
  'Shikimates and Phenylpropanoids',
  'Terpenoids',
  'Alkaloids',
  'Amino acids and Peptides',
  'Carbohydrates',
  'Polyketides'
)

# Filter out 'Others' class
filtered_data <- combined_data[combined_data$class_1 != "Others" & combined_data$class_2 != "Others", ]

# Create a matrix for the chord diagram
# Initialize an empty matrix with zeros
chord_matrix <- matrix(0, nrow = length(class_order), ncol = length(class_order))
rownames(chord_matrix) <- class_order
colnames(chord_matrix) <- class_order

# Fill the matrix with counts from the data
for (i in 1:nrow(filtered_data)) {
  # Check if both classes are in our filtered list
  if (filtered_data$class_1[i] %in% class_order && filtered_data$class_2[i] %in% class_order) {
    row_idx <- which(class_order == filtered_data$class_1[i])
    col_idx <- which(class_order == filtered_data$class_2[i])
    
    # Add the count to the matrix
    chord_matrix[row_idx, col_idx] <- filtered_data$count[i]
    
    # Make the matrix symmetric (if the same class appears on both sides)
    if (row_idx != col_idx) {
      chord_matrix[col_idx, row_idx] <- filtered_data$count[i]
    }
  }
}

# Set up colors for the plot using the provided custom colors
sector_colors <- c(
  "Fatty acids" = "#4e639e", 
  "Shikimates and Phenylpropanoids" = "#e54616", 
  "Terpenoids" = "#dba053",
  "Alkaloids" = "#ff997c", 
  "Amino acids and Peptides" = "#7fbfdd", 
  "Carbohydrates" = "#96a46b",
  "Polyketides" = "#760f00"
)

# Function to create the chord diagram
create_chord_diagram <- function() {
  # Reset the circos parameters
  circos.clear()
  
  # Create the chord diagram with customized parameters
  chordDiagram(
    chord_matrix,
    grid.col = sector_colors,
    transparency = 0.4,
    directional = FALSE,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.1),
    annotationTrackHeight = c(0.05, 0.1)
  )
}

# Create the SVG version with transparent background and smaller size
svg("chord_diagram_np_pathways.svg", width = 3.1, height = 3.1, bg = "transparent")
par(mfrow = c(1,1), mar = c(1, 1, 1, 1))
create_chord_diagram()
dev.off()

# Print confirmation message
cat("Chord diagram saved as chord_diagram_np_pathways.svg\n")