library(ggplot2)
library(dplyr)

### RNA Seq Preprocessing ###
#read in rna seq counts file
raw_rna_seq <- read.csv("/rds/homes/a/acd446/Computational_Biology_Group7/rna_raw_counts.csv", header = TRUE, row.names =1 ,stringsAsFactors = FALSE)

#calculate total read depth per sample
read_depth <- colSums(raw_rna_seq)

#view read depth
print(read_depth)

#convert to dataframe for plotting
read_depth_df <- data.frame(Sample = names(read_depth), ReadDepth = read_depth)
summary(read_depth_df$ReadDepth)

#Plot read depth distribution
read_depth <- ggplot(read_depth_df, aes(x = Sample, y = ReadDepth)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "RNA-seq Library Size per Sample", x = "Sample", y = "Total Reads") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("figures/preprocessing/read_depth.png", plot = read_depth, width = 6, height = 4, dpi = 300)


#Define a threshold (e.g., 5 million reads)
threshold <- 5000000

#Remove samples with low read depth
rna_counts_filtered <- read_depth_df[read_depth_df$ReadDepth > threshold,]

#Check dimensions after filtering
dim(rna_counts_filtered)

### Metabolite Data Preprocessing###
#Read in polar neg and polar pos csv
polar_neg_pqn <- read.csv("polar_neg_pqn_imputed.csv", header = TRUE, row.names =1 ,stringsAsFactors = FALSE)
polar_pos_pqn <- read.csv("polar_pos_pqn_imputed.csv", header = TRUE, row.names =1 ,stringsAsFactors = FALSE)

#Join polar neg and polar pos 
merged_data <- full_join(polar_neg_pqn, polar_pos_pqn, by = "mz", suffix = c("_pos", "_neg"))

#Count occurrences of each m/z value
mz_counts <- table(merged_data$mz)

#Find m/z values that appear more than once
duplicated_mz <- names(mz_counts[mz_counts > 1])

#Print duplicated m/z values
print(duplicated_mz)

#save cleaned RNA_seq and metabolites files 
write.csv(rna_counts_filtered, "filtered_rna.csv", row.names = FALSE)
write.csv(merged_data, "filtered_metabolites.csv", row.names = FALSE)
```