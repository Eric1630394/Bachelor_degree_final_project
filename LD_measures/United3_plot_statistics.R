# Load necessary library
library(ggplot2)
library(dplyr)

# Read the tab-separated file where information for LD values is stored. 
r2_values <- read.table("United3.ld", header=TRUE, sep="\t", stringsAsFactors=FALSE)

#Read the file containing the flanking information: original data.
flanking_info <- read.csv("data/Inversions_coordinates_porubsky_v3.csv", header = FALSE, sep = " ")

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

#Merge all data in one single data frame. 
data_frame_merged <- data.frame(
  Category = c("SD-flanked", "non-SD-flanked", "MEI-flanked", "IR-flanked"),
  Original_counts = c(SD_flanked_count,non_SD_flanked_count,MEI_flanked_count,IR_flanked_count),
  Counts_polymorphic = c(SD_flanked_count-SD_monomorphic-SD_only_1-SD_Y,non_SD_flanked_count-non_SD_monomorphic-non_SD_only_1-non_SD_Y,MEI_flanked_count-MEI_monomorphic-MEI_only_1-MEI_Y,IR_flanked_count-IR_monomorphic-IR_only_1-IR_Y),
  Counts_monomorphic =  c(SD_monomorphic,non_SD_monomorphic,MEI_monomorphic,IR_monomorphic),
  Counts_only_1_sample = c(SD_only_1,non_SD_only_1,MEI_only_1,IR_only_1),
  Counts_chrY = c(SD_Y,non_SD_Y,MEI_Y,IR_Y)
)

#Assign flanking information for every inversion in the LD dataset. 
r2_values$Type <- NA

for (i in 1:nrow(r2_values)) {
  matching_row <- flanking_info[flanking_info$V1 == r2_values$Inv[i], ]
  corresponding_name <- matching_row$V7
  r2_values$Type[i] <- corresponding_name
}

#Transform information in the Type column to specific values. 
r2_values <- r2_values %>%
  mutate(Type = case_when(
    Type %in% c("0", "SD") ~ "SD-flanked",
    Type %in% c("IR") ~ "IR-flanked",
    Type %in% c("LINE", "MEI-flanked", "SINE", "SVA", "TEs", "LTR") ~ "MEI-flanked",
    TRUE ~ Type
  ))

#Sort the data by Type.
sorted_r2_values <- r2_values %>%
  mutate(Inv = factor(Inv, levels = unique(Inv[order(Type)])))

#Create a dot plot to visualize the r2 values for inversions and SNPs.
plot_united <- ggplot(data = sorted_r2_values, mapping = aes(x = Inv, y = GLB, colour = Type)) +
  geom_point(size=3) + 
  labs(title = "r2 values for inversions and SNPs",
                      x = "Inversions",
                      y = "r2 values") +
  scale_color_manual(values = c("goldenrod1", "dodgerblue3","sienna","lightgreen")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 30,hjust = 0.5),
    axis.text.x = element_text(angle = 75, hjust = 1, size = 14, margin = margin(t = 10, b = 40)),  # Rotate and adjust x-axis labels
    axis.title.x = element_text(margin = margin(t = 10),size = 25),  # Add space below x-axis title
    axis.title.y = element_text(margin = margin(r = 10), size = 20),   # Add space to the left of y-axis title
    legend.text = element_text(size=15),
    legend.title = element_text(size=20),
)

ggsave(filename = "United3_ld_plot.png", plot = plot_united, width = 49, height = 20, bg = "white")

#Create a boxplot to know the distribution of r2 values for each inversion type.
boxplot <- ggplot(data = sorted_r2_values, mapping = aes(y = GLB, x = Type,fill=Type)) +
  geom_boxplot() + theme_bw() + scale_fill_manual(values = c("goldenrod1", "dodgerblue3","sienna","lightgreen"))

ggsave(filename = "Boxplot.png", plot = boxplot, width = 10, height = 10, bg = "white")

#Also, to know the distribution of r2 values for each inversion type, create a histogram.
histogram <- ggplot(data = sorted_r2_values, mapping = aes(x = GLB, fill = Type)) +
  geom_histogram(color="black",binwidth = 0.015,linewidth=0.2) + scale_fill_manual(values = c("goldenrod1", "dodgerblue3","sienna","lightgreen")) +
  scale_x_continuous(name = "r2 values") + 
  ggtitle("Histogram on different r2 values",) + 
  theme_light() +
  theme(plot.title = element_text(size = 13,hjust = 0.5))

ggsave(filename = "Histogram.png", plot = histogram, width = 10, height = 10, bg = "white")

#Generate an output containing five random rows from the dataset, just to visualize how data is organized. 
random_rows <- sorted_r2_values %>% sample_n(5)

#Inside the dataset, create a data frame containing those inversions with at least one tagSNP.  
ss <- subset(sorted_r2_values, sorted_r2_values$GLB == 1)

#Count the number of SNPs with r2 value = 1 for each inversion. 
count_data <- ss %>%
  group_by(Inv) %>%
  summarise(count = n())

count_data$Type <- NA

#Classify SNPs depending on the category of the inversion they are being assigned. 
for (i in 1:nrow(count_data)) {
  matching_row <- sorted_r2_values[sorted_r2_values$Inv == count_data$Inv[i], ]
  corresponding_name <- matching_row$Type
  count_data$Type[i] <- corresponding_name
}

#Create a barplot to visualize the results:
r2_values_perfect_LD <- ggplot(count_data, aes(x = Inv, y = count, fill = Type)) +
  geom_bar(stat = "identity",color="black",linewidth=0.2) +
  theme_minimal() +
  labs(title = "Inversions with r2 value = 1", x = "Inversions", y = "Total number") +
  scale_fill_manual(values = c("dodgerblue3","sienna","lightgreen")) +
  theme(
    plot.title = element_text(size = 10,hjust = 0.5),
    axis.text.x = element_text(angle = 75, hjust = 1, size = 5, margin = margin(t = 10, b = 40)),
)

ggsave(filename = "r2_values_perfect_LD.png", plot = r2_values_perfect_LD, width = 8, height = 8, bg = "white")

#It is also interesting to study the representation of each inversion category among original data and inversions with at least 1 tagSNP.

#First, count the number of inversions for each category inside those with at least 1 tagSNP. 
ss <- subset(sorted_r2_values, sorted_r2_values$GLB == 1)
r2_1_SD_flanked <- c()
r2_1_non_SD_flanked <- c()
r2_1_MEI_flanked <- c()
r2_1_IR_flanked <- c()
for (i in 1:nrow(ss)) {
  type <- ss$Type[i]
  if (type == "SD-flanked"){
    r2_1_SD_flanked <- c(r2_1_SD_flanked,ss$Inv[i])
  } else if (type == "non-SD-flanked"){
    r2_1_non_SD_flanked <- c(r2_1_non_SD_flanked,ss$Inv[i])
  } else if (type == "MEI-flanked") {
    r2_1_MEI_flanked <- c(r2_1_MEI_flanked,ss$Inv[i])
  } else {
    r2_1_IR_flanked <- c(r2_1_IR_flanked,ss$Inv[i])
  }
}

SD_flanked_r2_1_count <- length(unique(r2_1_SD_flanked))
non_SD_flanked_r2_1_count <- length(unique(r2_1_non_SD_flanked))
MEI_flanked_r2_1_count <- length(unique(r2_1_MEI_flanked))
IR_flanked_r2_1_count <- length(unique(r2_1_IR_flanked))

#Then two vectors are defined depending on if original data or inversions with tagSNPs is considered.  
original_count <- c(SD_flanked_count,non_SD_flanked_count,MEI_flanked_count,IR_flanked_count)
r2_1_count <- c(SD_flanked_r2_1_count,non_SD_flanked_r2_1_count,MEI_flanked_r2_1_count,IR_flanked_r2_1_count)

#Create a data frame with counts stored in the previous two vectors. 
original_perfect_ld <- data.frame(
  Category = c("SD-flanked", "non-SD-flanked", "MEI-flanked", "IR-flanked"),
  Original_count = original_count,
  Perfect_LD_count = r2_1_count
)

#Since it is interesting to plot proportions and not absoulte counts we then create:
original_proportions <- c()
perfect_ld_proportions <- c()

#Calculate proportion of inversions for each category among original data (after filtered) and inversions with at least 1 tagSNP.
for (i in 1:nrow(original_perfect_ld)) {
  original_proportions <- c(original_proportions, data_frame_merged$Counts_polymorphic[i]/sum(data_frame_merged$Counts_polymorphic))
  perfect_ld_proportions <- c(perfect_ld_proportions, original_perfect_ld$Perfect_LD_count[i]/sum(r2_1_count))
}

#Create the data frame with the corresponding proportions. 
proportions <- data.frame(
  Category = c("SD-flanked", "non-SD-flanked", "MEI-flanked", "IR-flanked"),
  Original_proportion = original_proportions,
  Perfect_LD_proportion = perfect_ld_proportions
)

#Convert the data frame for the proportion to long format. 
data_long <- proportions %>%
  tidyr::pivot_longer(
    cols = c(Original_proportion, Perfect_LD_proportion),
    names_to = "Count_Type",
    values_to = "Count"
  )

#Define colours for proper visualization:
category_colors <- c(
  "SD-flanked" = "lightgreen",
  "non-SD-flanked" = "sienna",
  "MEI-flanked" = "dodgerblue3",
  "IR-flanked" = "goldenrod1"
)

#Make the proper visualization in a barplot. 
barplot <- ggplot(data_long, aes(x = Category, y = Count, fill = Category, alpha = Count_Type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,color="black",linewidth=0.5) +
  scale_fill_manual(
    name = "Category",
    values = category_colors,
    guide = guide_legend(order = 1)  # Main legend first
  ) +
  scale_alpha_manual(
    name = "Count (number of inversions)",
    values = c(Original_proportion = 1, Perfect_LD_proportion = 0.5),
    labels = c("Original", "At least 1 tagSNP"),
    guide = guide_legend(order = 2)  # Secondary legend second
  ) +
  labs(
    title = "Number of original inversions vs inversions with at least 1 tagSNP",
    x = "Category",
    y = "Count (number of inversions)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.box = "horizontal",
    legend.spacing.x = unit(0.5, "cm"),
    plot.title = element_text(hjust=0.5)
  ) +
  guides(
    fill = guide_legend(title.position = "top"),
    alpha = guide_legend(title.position = "top")
  )

ggsave(filename = "Barplot_comparison.png", plot = barplot, width = 10, height = 10, bg = "white")

#Test statistical differences with a Fisher's exact test. 

# Initialize vector to store raw p-values
raw_p_values <- numeric(nrow(original_perfect_ld))

for (i in 1:nrow(original_perfect_ld)) {
  # In-category counts
  a <- original_perfect_ld$Perfect_LD_count[i]
  b <- data_frame_merged$Counts_polymorphic[i]
  
  # Out-of-category counts (the rest)
  rest_total <- sum(data_frame_merged$Counts_polymorphic)-b
  rest_perfect <- sum(original_perfect_ld$Perfect_LD_count)-a
  c <- rest_perfect
  d <- rest_total
  
  # Construct contingency table
  contingency <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                        dimnames = list(c("This_category", "Others"),
                                        c("Perfect_LD", "Not_Perfect")))
  
  # Run Fisher's Exact Test
  test_result <- fisher.test(contingency)
  raw_p_values[i] <- test_result$p.value
  cat("\nCategory:", original_perfect_ld$Category[i], "\n")
  print(contingency)
  print(fisher.test(contingency))
}

# Apply Benjamini-Hochberg correction
adjusted_p_values <- p.adjust(raw_p_values, method = "BH")

# Add adjusted p-values to your data frame
original_perfect_ld$BH_adjusted_p <- adjusted_p_values

# Comparison with Giner et al. for inversions originated by NH. 
Giner_invs_NH <- read.table("evolutionary_impact_NH_invs.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
vector_Giner_invs_NH <- c(Giner_invs_NH$V1)
for (inv in vector_Giner_invs_NH) {
  matching_rows <- sorted_r2_values[sorted_r2_values$Inv == inv, ]
  max_value <- max(matching_rows$GLB)
  print(max_value)
}

# Comparison with Giner et al. for inversions originated by NAHR.
Giner_invs_NAHR <- read.table("evolutionary_impact_NAHR_invs.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
vector_Giner_invs_NAHR <- c(Giner_invs_NAHR$V1)
for (inv in vector_Giner_invs_NAHR) {
  matching_rows <- sorted_r2_values[sorted_r2_values$Inv == inv, ]
  max_value <- max(matching_rows$GLB)
  print(max_value)
}
