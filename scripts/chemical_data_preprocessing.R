library(ggplot2)

#read in chemical data and reformat
chemicals <- read.delim("/rds/homes/a/acd446/Computational_Biology_Group7/water_chemicals.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
row.names(chemicals) <- chemicals$ChemName
chemicals_pca <- chemicals [,-1]
chemicals_pca <- chemicals_pca [,-1]

##### PCA to visualise variation within chemical data #####
#run pca of chemicals per sample
pca_result <- prcomp(chemicals_pca, scale. = TRUE)

# Plot PCA
pca_df <- as.data.frame(pca_result$x)
pca_df$Site <- rownames(chemicals_pca)  # Add site names

pca_chemicals <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Site)) +
  geom_point(size = 3, color = "blue") +
  geom_text(vjust = 1, hjust = 1) +
  labs(title = "PCA of Chemical Data", x = "PC1", y = "PC2") +
  theme_minimal()
ggsave("figures/preprocessing/pca_chemicals.png", plot =  pca_chemicals, width = 6, height = 4, dpi = 300)

#run pca of chemical data per site
pca_result_site <- prcomp(t(chemicals_pca))

# Plot PCA
pca_df_site <- as.data.frame(pca_result_site$x)
pca_df_site$Site <- rownames(t(chemicals_pca))  # Add site names

pca_chemicals_sites <- ggplot(pca_df_site, aes(x = PC1, y = PC2, label = Site)) +
  geom_point(size = 3, color = "blue") +
  geom_text(vjust = 1, hjust = 1) +
  labs(title = "PCA of Chemical Data", x = "PC1", y = "PC2") +
  theme_minimal()

ggsave("figures/preprocessing/pca_chemicals_sites.png", plot = pca_chemicals_sites, width = 6, height = 4, dpi = 300)

#filter any chemicals that are 0 (undetected) across all samples
chemicals_filtered <- chemicals_pca[rowSums(chemicals_pca != 0) > 0, ]
chemicals_filtered_t <- t(chemicals_filtered)


control_row <- data.frame(matrix(0, nrow = 1, ncol = ncol(chemicals_filtered_t)))
rownames(control_row) <- "Control"
colnames(control_row) <- colnames(chemicals_filtered_t)  # Ensure column names match

# Add the new row to the dataframe
chemicals_filtered <- rbind(control_row, chemicals_filtered_t)

#chemicals_scaled <- scale(chemicals_filtered)
#chemicals_log <- log10(chemicals_scaled +1)

#set rownames as sample ID
chemicals_filtered$Sample_ID <- row.names(chemicals_filtered)

#### Generating synthetic replicates for chemical data ####
# Function to generate replicates while maintaining order
generate_replicates <- function(sample_row) {
  synthetic <- as.data.frame(
    apply(sample_row[,-77], 2, function(x) jitter(rep(x, 5), amount = 0.1))
  )
  synthetic$Sample_ID <- rep(sample_row$Sample_ID, 5)
  synthetic$Replicate <- paste0("Rep_", 1:5)
  return(synthetic)
}

# Generate synthetic data for all samples
synthetic_data <- do.call(rbind, lapply(1:nrow(chemicals_filtered), function(i) generate_replicates(chemicals_filtered[i, ])))
chemicals_filtered2 <- chemicals_filtered[!rownames(chemicals_filtered) %in% "Control", ]
synthetic_data_2 <- do.call(rbind, lapply(1:nrow(chemicals_filtered2), function(i) generate_replicates(chemicals_filtered2[i, ])))

# Combine original and synthetic data
chemicals_filtered$Replicate <- "Rep_0"
chemicals_filtered2$Replicate <- "Rep_0"

final_data <- rbind(chemicals_filtered, synthetic_data)
final_data2 <- rbind(chemicals_filtered, synthetic_data_2)

# Order to keep replicates grouped
final_data <- final_data[order(final_data$Sample_ID, final_data$Replicate), ]
final_data2 <- final_data2[order(final_data2$Sample_ID, final_data2$Replicate), ]

final_data_final <- final_data[,-78]
final_data_final <- final_data_final[,-78]
final_data_final2 <- final_data2[,-78]
final_data_final2 <- final_data_final2[,-78]

finalfinaldata <- rbind(final_data_final, final_data_final2)

finalfinaldata <- finalfinaldata[!rownames(finalfinaldata) %in% "Control2", ]

#load before sample data
before_samples<- read.csv("sample_sheet.csv", row.names = 1, header = TRUE)

row.names(finalfinaldata) <- row.names(before_samples)

#remove missing samples
missing_samples <- c("D06A2", "D07A3", "D06B3", "D08B1", "D10B4", "D01B4", "D05B4")

filtered_chemicals <- subset(finalfinaldata, !(row.names(finalfinaldata) %in% missing_samples))

filtered_chemicals <- filtered_chemicals[,-77]
str(filtered_chemicals)

chemical_dataa <- log1p(filtered_chemicals)
final_chemicals <- scale(chemical_dataa)
hist(t(final_chemicals))

# Perform PCA
pca_result <- prcomp(final_chemicals, scale. = TRUE)

# Plot PCA
pca_df <- as.data.frame(pca_result$x)
pca_df$Site <- rownames(final_chemicals)  # Add site names

pca_final <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Site)) +
  geom_point(size = 3, color = "blue") +
  geom_text(vjust = 1, hjust = 1) +
  labs(title = "PCA of Chemical Data", x = "PC1", y = "PC2") +
  theme_minimal()
ggsave("figures/preprocessing/pca_final_chemicals.png", plot = pca_final, width = 6, height = 4, dpi = 300)

pca_result_site <- prcomp(t(final_chemicals))

# Plot PCA
library(ggplot2)
pca_df <- as.data.frame(pca_result_site$x)
pca_df$Site <- rownames(t(chemicals_pca))  # Add site names

pca_final_site <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Site)) +
  geom_point(size = 3, color = "blue") +
  geom_text(vjust = 1, hjust = 1) +
  labs(title = "PCA of Chemical Data", x = "PC1", y = "PC2") +
  theme_minimal()

ggsave("figures/preprocessing/pca_final_chemicals_site.png", plot = pca_final_site, width = 6, height = 4, dpi = 300)

write.csv(final_chemicals, "chemicals_mofa.csv")

