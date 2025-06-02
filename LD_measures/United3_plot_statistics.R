# Load necessary library
library(ggplot2)
library(dplyr)
library(grid)

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

category_colors <- c(
  "SD-flanked" = "seagreen3",
  "non-SD-flanked" = "paleturquoise4",
  "MEI-flanked" = "turquoise3",
  "IR-flanked" = "sienna3"
)

#Create a dot plot to visualize the r2 values for inversions and SNPs.
plot_united <- ggplot(data = sorted_r2_values, mapping = aes(x = Inv, y = GLB, colour = Type)) +
  geom_point(size=3) + 
  labs(x = "Inversions",y = "r2 values") +
  scale_color_manual(name = "Category", values = category_colors,  guide = guide_legend(byrow = TRUE)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 30,hjust = 0.5),
    axis.text.x = element_text(angle = 75, hjust = 1, size = 17, margin = margin(t = 10, b = 40)),  # Rotate and adjust x-axis labels
    axis.title.x = element_text(margin = margin(t = 10),size = 30),  # Add space below x-axis title
    axis.title.y = element_text(margin = margin(r = 10), size = 30),   # Add space to the left of y-axis title
    legend.text = element_text(size=25),
    legend.title = element_text(size=27),
    legend.spacing.y = unit(20, "cm"),
    legend.box = "vertical", 
    axis.text.y = element_text(size=25)
)

plot_united <- ggplot(data = sorted_r2_values, mapping = aes(x = Inv, y = GLB, colour = Type)) +
  geom_point(size = 3) + 
  labs(x = "Inversions", y = "r2 values", title = NULL) +  # Remove title
  scale_color_manual(
    name = "Category",
    values = category_colors,
    guide = guide_legend(byrow = TRUE)
  ) +
  theme_minimal() +
  theme(
    plot.title = element_blank(),  # Ensure no title is shown
    axis.text.x = element_text(angle = 75, hjust = 1, size = 17, margin = margin(t = 10, b = 40)),
    axis.title.x = element_text(margin = margin(t = 10), size = 30),
    axis.title.y = element_text(margin = margin(r = 10), size = 30),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 27),
    legend.spacing.y = unit(20, "cm"),
    legend.box = "vertical", 
    axis.text.y = element_text(size = 25)
  )


ggsave(filename = "United3_ld_plot.png", plot = plot_united, width = 49, height = 20, bg = "white")

#Also, to know the distribution of r2 values for each inversion type, create a histogram.
histogram <- ggplot(data = sorted_r2_values, mapping = aes(x = GLB, fill = Type)) +
  geom_histogram(color="black",binwidth = 0.03,linewidth=0.2) + scale_fill_manual(name = "Category", values = category_colors) + 
  labs(x = "LD measures",
       y = "Absolute count of SNPs") + 
  theme_light() +
  theme(
    axis.text.x = element_text(size=20),
    axis.text.y = element_text(size=20),
    legend.text = element_text(size=15),
    axis.title = element_text(size=25),
    axis.title.x = element_text(size = 25),
    legend.title = element_text(size=20)
)

ggsave(filename = "Histogram.png", plot = histogram, width = 15, height = 15, bg = "white")

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
  labs(x = "Inversions", y = "Absolute count of SNPs in perfect LD") +
  scale_fill_manual(name = "Category", values = c("turquoise3","paleturquoise4","seagreen3")) +
  theme(
    plot.title = element_text(size = 10,hjust = 0.5),
    axis.text.x = element_text(angle=45,hjust = 1, size = 11.5, margin = margin(t = 10, b = 40)),
    axis.text.y = element_text(size=7),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    legend.title = element_text(size=15),
    legend.text=element_text(size=12)
)

ggsave(filename = "r2_values_perfect_LD.png", plot = r2_values_perfect_LD, width = 22, height = 7,bg="white")

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
original_count <- c(data_frame_merged$Counts_polymorphic[1],data_frame_merged$Counts_polymorphic[2],data_frame_merged$Counts_polymorphic[3],data_frame_merged$Counts_polymorphic[4])
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

#Make the proper visualization in a barplot. 
barplot <- ggplot(data_long, aes(x = Category, y = Count, fill = Category, alpha = Count_Type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.5) +
  scale_fill_manual(
    name = "Category",
    values = category_colors,
    guide = guide_legend(order = 1)
  ) +
  scale_alpha_manual(
    name = "Group of inversions",
    values = c(Original_proportion = 1, Perfect_LD_proportion = 0.5),
    labels = c("Original", "At least 1 tagSNP"),
    guide = guide_legend(order = 2)
  ) +
  labs(
    x = NULL,  # Removes x-axis title
    y = "Proportion among the total of inversions",
    title = NULL  # Removes plot title
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 15),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_blank(),  # Also ensures x-axis title is blank
    legend.position = "top",
    legend.box = "horizontal",
    legend.spacing.x = unit(1, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 13),
    plot.title = element_blank(),  # Ensures plot title is blank
    axis.title.y = element_text(size = 15)
  ) +
  guides(
    fill = guide_legend(title.position = "top"),
    alpha = guide_legend(title.position = "top")
  )


ggsave(filename = "Barplot_comparison.png", plot = barplot, width = 10, height = 10,bg="white")

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
