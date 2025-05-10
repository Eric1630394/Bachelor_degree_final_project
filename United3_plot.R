# Load necessary library
library(ggplot2)
library(dplyr)

# Read the tab-separated file where information for LD values is stored. 
r2_values <- read.table("United3.ld", header=TRUE, sep="\t", stringsAsFactors=FALSE)

#Read the file containing the flanking info: original data.
flanking_info <- read.csv("data/Inversions_coordinates_porubsky_v3.csv", header = FALSE, sep = " ")

#Assign flanking information for every inversion in the LD dataset. 
r2_values$Type <- NA

for (i in 1:nrow(r2_values)) {
  matching_row <- flanking_info[flanking_info$V1 == r2_values$Inv[i], ]
  corresponding_name <- matching_row$V7
  r2_values$Type[i] <- corresponding_name
}

#Transform information in the seventh column to specific values. 
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

#We can generate an output containing five random rows from the dataset, just to visualize how data is organized. 
random_rows <- sorted_r2_values %>% sample_n(5)

#Inside the dataset, it can be created a data frame containing those inversions with at least one tagSNP.  
ss <- subset(sorted_r2_values, sorted_r2_values$GLB == 1)

#Count the number of SNPs with r2 value = 1.
count_data <- ss %>%
  group_by(Inv) %>%
  summarise(count = n())

count_data$Type <- NA

#Classify SNPs depending on the category of the inversion they are being analysed. 
for (i in 1:nrow(count_data)) {
  matching_row <- sorted_r2_values[sorted_r2_values$Inv == count_data$Inv[i], ]
  corresponding_name <- matching_row$Type
  count_data$Type[i] <- corresponding_name
}

#Create a barplot to visalize results:
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
