library(data.table)
library(MOFA2)
library(ggplot2)
library(dplyr)

#load rna data
metabolites <- read.csv("/rds/homes/a/acd446/Computational_Biology_Group7/significant_metabolites_irf.csv", row.names = 1, header = TRUE, stringsAsFactors = FALSE)

#load metabolites
rna <- read.csv("/rds/homes/a/acd446/Computational_Biology_Group7/significant_genes_irf.csv", row.names = 1, header = TRUE, stringsAsFactors = FALSE)

#load sample data
sample_data <- read.csv("/rds/homes/a/acd446/Computational_Biology_Group7/sample_data_filtered.csv", row.names = 1, header = TRUE, stringsAsFactors = FALSE)

#chemical data
chemicals <- read.csv("chemicals_mofa.csv", row.names =1, header = TRUE, stringsAsFactors = FALSE)

#create empty list to store data
combined_data <- list()

#combne all data types into a list of lists
combined_data[["RNA_Seq"]] <- as.matrix(rna)
combined_data[["Metabolites"]] <- as.matrix(t(metabolites))
combined_data[["Chemicals"]] <- as.matrix(t(chemicals))

#create metadata dataframe
metadata <- data.frame(
  sample = sample_data$SampleID,
  REF = sample_data$REF,
  site = sample_data$Site,
  group = sample_data$REF
)

#create MOFA object
MOFAobject <- create_mofa(combined_data, groups = sample_data$REF)

data_opts <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
train_opts <- get_default_training_options(MOFAobject)


#set metadata in MOFA object
samples_metadata(MOFAobject) <- metadata

#plot overview of MOFA object
png("figures/mofa/overview_mofa.png", width = 800, height = 600)  # open device
plot_data_overview(MOFAobject)
dev.off()  # close device

#initialise mofa object
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

#train MOFA model and save to output file
outfile = file.path(getwd(),"output/mofa_model.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile)

#plot overview of trained model
png("figures/mofa/overview_mofa_trained.png", width = 800, height = 600)  # open device
plot_data_overview(MOFAobject.trained)
dev.off()  # close device

#plot variance explained 
png("figures/mofa/variance_explained.png", width = 800, height = 600)  # open device
plot_variance_explained(MOFAobject.trained, x="view", y="factor")
dev.off()  # close device

#plot factors
png("figures/mofa/factor_plot.png", width = 800, height = 600)  # open device
plot_factor(MOFAobject.trained, 
            factor = c(1,2,3,4),
            color_by = "site",
            shape_by = "REF"
)
dev.off()  # close device

#plot factors ggplot
library(ggplot2)
p <- plot_factor(MOFAobject.trained, 
                 factors = c(1,2,3,4),
                 color_by = "REF",
                 dot_size = 3,        # change dot size
                 dodge = T,           # dodge points with different colors
                 legend = F,          # remove legend
                 add_violin = T,      # add violin plots,
                 violin_alpha = 0.25  # transparency of violin plots
)

# The output of plot_factor is a ggplot2 object that we can edit
p <- p + 
  scale_color_manual(values=c("Control"="black", "1x"="red", "10x" = "blue")) +
  scale_fill_manual(values=c("Control"="black", "1x"="red", "10x" = "blue"))

ggsave("figures/mofa/factor_plot2.png", plot = p, width = 6, height = 4, dpi = 300)

#plot factors again
png("figures/mofa/factor_plot3.png", width = 800, height = 600)  # open device
plot_factors(MOFAobject.trained, 
             factors = c(1,2,3),
             color_by = "site"
)
dev.off()  # close device

#plot factor weights
png("figures/mofa/chemicals_weights_plot.png", width = 800, height = 600)  # open device
plot_weights(MOFAobject.trained,
             view = "Chemicals",
             factor = 3,
             nfeatures = 20,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)
dev.off()  # close device

#plot top weights for chemicals
png("figures/mofa/chemicals_top_weights_plot.png", width = 800, height = 600)  # open device
plot_top_weights(MOFAobject.trained,
                 view = "Chemicals",
                 factor = 4,
                 nfeatures = 15
)
dev.off()  # close device

#plot top weights for metabolites
png("figures/mofa/metabolites_top_weights_plot.png", width = 800, height = 600)  # open device
plot_top_weights(MOFAobject.trained,
                 view = "Metabolites",
                 factor = 3,
                 nfeatures = 15
)
dev.off()  # close device

#plot top weights for RNA
png("figures/mofa/RNA_top_weights_plot.png", width = 800, height = 600)  # open device
plot_top_weights(MOFAobject.trained,
                 view = "RNA_Seq",
                 factor = 3,
                 nfeatures = 15
)
dev.off()  # close device

#plot chemicals heatmap
png("figures/mofa/chemicals_heatmap.png", width = 800, height = 600)  # open device
plot_data_heatmap(MOFAobject.trained,
                  view = "Chemicals",         # view of interest
                  factor = 2,             # factor of interest
                  features = 20,          # number of features to plot (they are selected by weight)
                  
                  # extra arguments that are passed to the `pheatmap` function
                  cluster_rows = TRUE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE
)
dev.off()  # close device

#plot rna heatmap
png("figures/mofa/rna_heatmap.png", width = 800, height = 600)  # open device
plot_data_heatmap(MOFAobject.trained,
                  view = "RNA_Seq",         # view of interest
                  factor = 2,             # factor of interest
                  features = 20,          # number of features to plot (they are selected by weight)
                  
                  # extra arguments that are passed to the `pheatmap` function
                  cluster_rows = TRUE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = TRUE,
)
dev.off()  # close device

#plot metabolites heatmap
png("figures/mofa/metabolite_heatmap.png", width = 800, height = 600)  # open device
plot_data_heatmap(MOFAobject.trained,
                  view = "Metabolites",         # view of interest
                  factor = 3,             # factor of interest
                  features = 20,          # number of features to plot (they are selected by weight)
                  
                  # extra arguments that are passed to the `pheatmap` function
                  cluster_rows = TRUE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
)
dev.off()  # close device

#run tsne to visualise clustering following mofa
set.seed(42)
#model <- run_umap(model)
model <- run_tsne(MOFAobject.trained)

#plot tsne
png("figures/mofa/tsne.png", width = 800, height = 600)  # open device
plot_dimred(model,
            method = "TSNE",  # method can be either "TSNE" or "UMAP"
            color_by = "site",
            shape_by = "REF",
            title = "tSNE coloured by site"
)
dev.off()  # close device

#run umap to visualise clustering following mofa
set.seed(42)
model <- run_umap(MOFAobject.trained)
#model <- run_tsne(MOFAobject.trained)

#plot umap
png("figures/mofa/umap.png", width = 800, height = 600)  # open device
plot_dimred(model,
            method = "UMAP",  # method can be either "TSNE" or "UMAP"
            color_by = "site",
            shape_by = "REF",
            title = "tSNE coloured by site"
)
dev.off()  # close device

#extract factors from mofa object (latent variables that show the sources of variation across samples)
factors <- get_factors(MOFAobject.trained, factors = "all")
lapply(factors,dim)

#extract weights from mofa object (how much each featurecontributes to each factor for each data view)
weights <- get_weights(MOFAobject.trained, views = "all", factors = "all")
lapply(weights,dim)

#extract weights for each data view
chemical_weights <- as.data.frame(weights$"Chemicals")
rna_weights <- as.data.frame(weights$"RNA_Seq")
metabolites_weights <- as.data.frame(weights$"Metabolites")

#arrange to create ranked list
rna_arranged <- rna_weights %>% arrange(desc(rna_weights$Factor3))
metabolites_arranged <- metabolites_weights %>% arrange(desc(metabolites_weights$Factor3))
metabolites_arranged2 <- metabolites_weights %>% arrange(desc(metabolites_weights$Factor1))

#saved ranked list
write.csv(rna_arranged, "rna_ranked_list.csv")
write.csv(metabolites_arranged2, "metabolites_ranked_list2.csv")