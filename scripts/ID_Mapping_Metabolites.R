library(dplyr)
library(tidyr)
library(stringr)
# read target IDs
# df <- read.csv("../data/rna_norm_counts.csv", header = T, row.names = 1, check.names = F)
df <- read.csv('metabolites_ranked_list.csv', header = T, row.names = 1, check.names = F)

#select top 50 metabolites
df <- df[1:50,]

#select factor3 weights
rownames(df) <- str_remove(rownames(df), "^X")
df <- df %>% select(Factor3)
ids <- row.names(df)
length(ids)

#read in hmdb annotations
id.maps <- read.csv("polar_pos_pkl_to_hmdb_annotations.tsv", header = F, row.names = 1, sep = '\t', stringsAsFactors = F)
id.maps2 <- read.csv("polar_neg_pkl_to_hmdb_annotations.tsv", header = F, row.names = 1, sep = '\t', stringsAsFactors = F)

#combine polar pos and polar neg annotations
id.maps <- rbind(id.maps, id.maps2)

#mapping
id.maps.matched <- subset(id.maps, row.names(id.maps) %in% ids)
#id.maps.matched$V2 <- as.character(id.maps.matched$V2)
#id.maps.matched <- id.maps.matched[id.maps.matched$V2 != '',]
id.maps.matched <- subset(id.maps.matched, V3 != "")

#rank by factor 3 scores
id.maps.matched.ranked <- id.maps.matched[rownames(df)[rownames(df) %in% rownames(id.maps.matched)], , drop = FALSE]

#keep only unique ids
id.ranked.unique <- id.maps.matched.ranked %>% distinct(V3, .keep_all = TRUE)
df_ranked <- df %>% filter(rownames(.) %in% rownames(id.maps.matched.ranked)) %>% mutate(Metabolite = id.maps.matched.ranked$V3 )%>% select(Metabolite, Factor3) %>% separate_rows(Metabolite, sep = ";")

#save to files for input to hmdb
write.table(df_ranked, file = "ranked_metabolites.rnk", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.csv(df_ranked, "ranked_metabolites.csv")
write.table(df_ranked, "metabolites_ranked.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(df_ranked$Metabolite, "metabolites_ranked_ID.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)