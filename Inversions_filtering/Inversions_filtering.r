library(webr)
library(scales)
library(dplyr)
library(tibble)
library(reshape2)
library(tidyr)
library(ggplot2)

# Read the table contained in the .txt for the genotypes info. 
class <- read.table("Genotyping_info.txt", header=FALSE, sep=",", stringsAsFactors=FALSE)

#Transform information in the third column to specific values. 
class <- class %>%
  mutate(V3 = case_when(
    V3 %in% c("0", "SD") ~ "SD-flanked",
    V3 %in% c("IR") ~ "IR-flanked",
    V3 %in% c("LINE", "MEI-flanked", "SINE", "SVA", "TEs", "LTR") ~ "MEI-flanked",
    TRUE ~ V3
  ))

#Count the number of SD-flanked inversions that are monomorphic, only 1 sample genotyped, or are in Chr Y.
SD_monomorphic <- sum(class$V3 == "SD-flanked" & class$V2 == "monomorphic", na.rm = TRUE)
SD_only_1 <- sum(class$V3 == "SD-flanked" & class$V2 == "only 1 sample genotyped", na.rm = TRUE)
SD_Y <- sum(class$V3 == "SD-flanked" & class$V2 == "Chr Y inversion", na.rm = TRUE)

#Count the number of non-SD-flanked inversions that are monomorphic, only 1 sample genotyped, or are in Chr Y.
non_SD_monomorphic <- sum(class$V3 == "non-SD-flanked" & class$V2 == "monomorphic", na.rm = TRUE)
non_SD_only_1 <- sum(class$V3 == "non-SD-flanked" & class$V2 == "only 1 sample genotyped", na.rm = TRUE)
non_SD_Y <- sum(class$V3 == "non-SD-flanked" & class$V2 == "Chr Y inversion", na.rm = TRUE)

#Count the number of MEI-flanked inversions that are monomorphic, only 1 sample genotyped, or are in Chr Y.
MEI_monomorphic <- sum(class$V3 == "MEI-flanked" & class$V2 == "monomorphic", na.rm = TRUE)
MEI_only_1 <- sum(class$V3 == "MEI-flanked" & class$V2 == "only 1 sample genotyped", na.rm = TRUE)
MEI_Y <- sum(class$V3 == "MEI-flanked" & class$V2 == "Chr Y inversion", na.rm = TRUE)

#Count the number of IR-flanked inversions that are monomorphic, only 1 sample genotyped, or are in Chr Y.
IR_monomorphic <- sum(class$V3 == "IR-flanked" & class$V2 == "monomorphic", na.rm = TRUE)
IR_only_1 <- sum(class$V3 == "IR-flanked" & class$V2 == "only 1 sample genotyped", na.rm = TRUE)
IR_Y <- sum(class$V3 == "IR-flanked" & class$V2 == "Chr Y inversion", na.rm = TRUE)

#Read the file containing the flanking info: original data.
flanking_info <- read.csv("data/Inversions_coordinates_porubsky_v3.csv", header = FALSE, sep = " ")

#Transform information in the seventh column to specific values. 
flanking_info <- flanking_info %>%
  mutate(V7 = case_when(
    V7 %in% c("0", "SD") ~ "SD-flanked",
    V7 %in% c("IR") ~ "IR-flanked",
    V7 %in% c("LINE", "MEI-flanked", "SINE", "SVA", "TEs", "LTR") ~ "MEI-flanked",
    TRUE ~ V7
  ))

#Count the number of inversions in each category:
SD_flanked_count <- sum(flanking_info$V7 == "SD-flanked", na.rm = TRUE)
non_SD_flanked_count <- sum(flanking_info$V7 == "non-SD-flanked", na.rm = TRUE)
IR_flanked_count <- sum(flanking_info$V7 == "IR-flanked", na.rm = TRUE)
MEI_flanked_count <- sum(flanking_info$V7 == "MEI-flanked", na.rm = TRUE)

#Merge all data in one single data frame. 
data_frame_merged <- data.frame(
  Category = c("SD-flanked", "non-SD-flanked", "MEI-flanked", "IR-flanked"),
  Original_counts = c(SD_flanked_count,non_SD_flanked_count,MEI_flanked_count,IR_flanked_count),
  Counts_polymorphic = c(SD_flanked_count-SD_monomorphic-SD_only_1-SD_Y,non_SD_flanked_count-non_SD_monomorphic-non_SD_only_1-non_SD_Y,MEI_flanked_count-MEI_monomorphic-MEI_only_1-MEI_Y,IR_flanked_count-IR_monomorphic-IR_only_1-IR_Y),
  Counts_monomorphic =  c(SD_monomorphic,non_SD_monomorphic,MEI_monomorphic,IR_monomorphic),
  Counts_only_1_sample = c(SD_only_1,non_SD_only_1,MEI_only_1,IR_only_1),
  Counts_chrY = c(SD_Y,non_SD_Y,MEI_Y,IR_Y)
)

#Define colours for proper visualization:
category_colors <- c(
  "SD-flanked" = "lightgreen",
  "non-SD-flanked" = "sienna",
  "MEI-flanked" = "dodgerblue3",
  "IR-flanked" = "goldenrod1"
)

# Convert to long format
data_long <- data_frame_merged %>%
  pivot_longer(cols = -c(Category, Original_counts), 
               names_to = "Count_Type", 
               values_to = "Count") %>%
  mutate(Proportion = Count / Original_counts,
         # Clean up count type names
         Count_Type = gsub("Counts_", "", Count_Type),
         # Create a combined factor for color mapping
         Color_Group = paste(Category, Count_Type, sep = "_"))

#Create a gradient of colors for the different categories:
gradient_colors <- unlist(lapply(names(category_colors), function(cat) {
  base_color <- category_colors[cat]
  colorRampPalette(c("white", base_color))(4)[2:4] # Generate 3 shades, skipping white
}))

# Assign names to the gradient colors:
names(gradient_colors) <- paste(rep(names(category_colors), each = 3),
                                c("polymorphic", "monomorphic", "only_1_sample", "chrY")[1:3],
                                sep = "_")

# Add chrY colors (using darker shade of each category)
gradient_colors <- c(gradient_colors, 
                     "SD-flanked_chrY" = "darkgreen",
                     "non-SD-flanked_chrY" = "sienna4",
                     "MEI-flanked_chrY" = "dodgerblue4",
                     "IR-flanked_chrY" = "goldenrod3")

# Create the plot with proper formatting
filtering <- ggplot(data_long, aes(x = Category, y = Proportion, fill = Color_Group)) +
  geom_bar(stat = "identity",color = "black",linewidth=0.3) +
  scale_fill_manual(
    values = gradient_colors,
    breaks = paste(rep(names(category_colors), each = 4),
                   rep(c("polymorphic", "monomorphic", "only_1_sample", "chrY"), 4),
                   sep = "_"),
    labels = rep(c("Polymorphic", "Monomorphic", "Only 1 sample", "chrY"), 4),
    name = "Count Type"
  ) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proportion of Different Count Types by Flanking Category",
    x = "Flanking Category",
    y = "Proportion"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )
ggsave(filename = "Barplot_filtering_no_low.png", plot = filtering, width = 10, height = 10, bg = "white")
