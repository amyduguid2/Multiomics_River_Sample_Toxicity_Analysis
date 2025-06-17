library(iRF)
library(AUC)
library(ggplot2)

#sample metadata
sample_data <- read.csv("/rds/homes/a/acd446/Computational_Biology_Group7/sample_data_filtered.csv", row.names =1, header = TRUE, stringsAsFactors = FALSE)

#rna filtered
rna_vst <- read.csv("/rds/homes/a/acd446/Computational_Biology_Group7/filtered_rna.csv", row.names = 1, header = TRUE, stringsAsFactors = FALSE)

#metabolites
metabolites <- read.csv("/rds/homes/a/acd446/Computational_Biology_Group7/filtered_metabolites.csv", row.names = 1, header = TRUE, stringsAsFactors = FALSE)

#chemical info
chemicals <- read.delim("/rds/homes/a/acd446/Computational_Biology_Group7/water_chemicals.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


##############RNA Seq IRF##################

X <- as.data.frame(t(rna_vst))
Y <- sample_data$REF

# prepare inputs
X <- X[Y != 'Control',]
Y <- as.factor(Y[Y != 'Control'])
Y <- as.factor(as.numeric(relevel(Y, '1x'))-1)

n <- dim(X)[1]
p <- dim(X)[2]

# split training-testing set
train.id <- sample(1:n, size = 8*round(n/10))
test.id <- setdiff(1:n, train.id)

# fit iRF without interaction importance estimation
sel.prob <- rep(1/p, p)

fit <- iRF(x = X[train.id,], 
           y = Y[train.id], 
           xtest = X[test.id,], 
           ytest = Y[test.id],
           n.iter = 5, 
           iter.return = 1:5,
           n.core = 4
)

# plot ROC/AUC
rf <- fit$rf.list

png("figures/roc_curve_rna.png", width = 1000, height = 800)  # open device

plot(0:1, 0:1, type = 'l', lty = 2, xlab = 'FPR', ylab = 'TPR', main = 'IRF ROC Curve for RNA-Seq Data')
for (iter in 1:5){
  cat(paste('iter = ', iter, ':: '))
  roc.info <- roc(rf[[iter]]$test$votes[,2], Y[test.id])
  lines(roc.info$fpr, roc.info$tpr, type = 'l', col = iter, lwd = 2)
  cat(paste('AUROC: ', round(100*auc(roc.info), 2), '%\n', sep = ''))
}
legend('bottomright', legend=paste('iter:', 1:iter), col=1:iter, lwd=2, bty='n')

dev.off()  # close device

### Extract important features based on Gini value ###

# Extract Gini importance values
gini_values <- fit$rf.list[[1]]$importance[, "MeanDecreaseGini"]

# Compute mean and standard deviation
mean_gini <- mean(gini_values)
sd_gini <- sd(gini_values)

# Set significance threshold (mean + 1 SD)
significant_threshold <- mean_gini + sd_gini

# Filter significant features
significant_features_rna <- names(gini_values[gini_values > significant_threshold])

#top 50 features rnaseq
top_50_features_rna <- significant_features_rna[1:50]

# View significant features
length(significant_features_rna)

#filter rna seq data for significant features
filtered_rna_seq <- rna_vst[rownames(rna_vst) %in% significant_features_rna,]

#filter top 50 features
top_50_rna <- rna_vst[rownames(rna_vst) %in% top_50_features_rna,]

#check dimensions
dim(filtered_rna_seq)

#pca of filtered RNA seq data
pca_rna <- prcomp(t(filtered_rna_seq), center = TRUE, scale. = TRUE)

pca_rna_df <- data.frame(PC1 = pca_rna$x[,1], PC2 = pca_rna$x[,2], Group = sample_data$REF, Site = sample_data$Site, SampleID = sample_data$SampleID)

# Extract variance explained
variance_explained <- pca_rna$sdev^2 / sum(pca_rna$sdev^2) * 100  # Convert to %

# Format axis labels with variance explained
x_label_rna <- paste0("PC1 (", round(variance_explained[1], 2), "%)")
y_label_rna <- paste0("PC2 (", round(variance_explained[2], 2), "%)")

pca_rna <- ggplot(pca_rna_df, aes( x = PC1, y = PC2, color = Site, shape = Group))+
  geom_point(size = 4, alpha = 0.8)+
  theme_minimal()+
  labs(title = "PCA of RNA Seq VST data After IRF Feature Selection", x = x_label_rna, y = x_label_rna)

ggsave("figures/pca_rna.png", plot = pca_rna, width = 6, height = 4, dpi = 300)


############## Metabolomics IRF ##################
X2 <- as.data.frame(t(metabolites))
Y2 <- sample_data$REF

# prepare inputs
X2 <- X2[Y2 != 'Control',]
Y2 <- as.factor(Y2[Y2 != 'Control'])
Y2 <- as.factor(as.numeric(relevel(Y2, '1x'))-1)

n2 <- dim(X2)[1]
p2 <- dim(X2)[2]

# split training-testing set
train.id2 <- sample(1:n2, size = 8*round(n2/10))
test.id2 <- setdiff(1:n2, train.id2)


# fit iRF without interaction importance estimation
sel.prob2 <- rep(1/p2, p2)

fit2 <- iRF(x = X2[train.id2,], 
            y = Y2[train.id2], 
            xtest = X2[test.id2,], 
            ytest = Y2[test.id2],
            n.iter = 5, 
            iter.return = 1:5,
            n.core = 4
)

# plot ROC/AUC

png("figures/roc_curve_metabolites.png", width = 1000, height = 800)  # open device

rf <- fit2$rf.list
plot(0:1, 0:1, type = 'l', lty = 2, xlab = 'FPR', ylab = 'TPR', main = 'IRF ROC Curve for Metabolite Data')
for (iter in 1:5){
  cat(paste('iter = ', iter, ':: '))
  roc.info <- roc(rf[[iter]]$test$votes[,2], Y[test.id])
  lines(roc.info$fpr, roc.info$tpr, type = 'l', col = iter, lwd = 2)
  cat(paste('AUROC: ', round(100*auc(roc.info), 2), '%\n', sep = ''))
}
legend('bottomright', legend=paste('iter:', 1:iter), col=1:iter, lwd=2, bty='n')

dev.off()  # close device

# Extract Gini importance values
gini_values2 <- fit2$rf.list[[1]]$importance[, "MeanDecreaseGini"]

# Compute mean and standard deviation
mean_gini2 <- mean(gini_values2)
sd_gini2 <- sd(gini_values2)

# Set significance threshold (mean + 1 SD)
significant_threshold2 <- mean_gini2 + sd_gini2

# Filter significant features
significant_features2 <- names(gini_values2[gini_values2 > significant_threshold2])

top_50_features_metabolites <- significant_features2[1:50]

# View significant features
print(significant_features2)

filtered_metabolites <- t(metabolites)[, colnames(t(metabolites)) %in% significant_features2]

top_50_metabolites <-  t(metabolites)[, colnames(t(metabolites)) %in% top_50_features_metabolites]
dim(filtered_metabolites)

pca_metabolites <- prcomp(filtered_metabolites, center = TRUE, scale. = TRUE)

pca_metabolites_df <- data.frame(PC1 = pca_metabolites$x[,1], PC2 = pca_metabolites$x[,2], Group = sample_data$REF, Site = sample_data$Site)

# Extract variance explained
variance_explained <- pca_metabolites$sdev^2 / sum(pca_metabolites$sdev^2) * 100  # Convert to %

# Format axis labels with variance explained
x_label_metabolites <- paste0("PC1 (", round(variance_explained[1], 2), "%)")
y_label_metabolites <- paste0("PC2 (", round(variance_explained[2], 2), "%)")

pca_metabolites <- ggplot(pca_metabolites_df, aes( x = PC1, y = PC2, color = Site, shape = Group))+
  geom_point(size = 4, alpha = 0.8)+
  theme_minimal()+
  labs(title = "PCA of Metabolite Data After IRF Feature Selection", x = x_label_metabolites, y = y_label_metabolites)

ggsave("figures/pca_metabolites.png", plot = pca_metabolites, width = 6, height = 4, dpi = 300)

top_50_metabolites <- as.data.frame(top_50_metabolites)
top_50_rna <- as.data.frame(t(top_50_rna))

# Ensure both datasets are numeric
#filtered_metabolites[] <- lapply(filtered_metabolites, as.numeric)
#filtered_rna_seq[] <- lapply(filtered_rna_seq, as.numeric)

combined_data <- cbind(top_50_metabolites, top_50_rna)
dim(combined_data)
str(combined_data)
# Ensure combined_data is a dataframe
combined_data_df <- as.data.frame(combined_data)

# Add group labels
group <- as.factor(sample_data$REF)


# Fit MANOVA model
manova_model <- manova(as.matrix(combined_data) ~ group)

#View summary with Wilks' Lambda
summary(manova_model, test = "Wilks")

#write significant genes and metabolites to a csv file
significant_genes_irf <- write.csv(filtered_rna_seq, "significant_genes_irf.csv", row.names = TRUE)
significant_metabolites_irf <- write.csv(filtered_metabolites, "significant_metabolites_irf.csv", row.names = TRUE)


