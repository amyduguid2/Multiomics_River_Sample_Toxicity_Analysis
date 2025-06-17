library(dplyr)
library(tidyr)

# read target IDs
df <- read.csv('rna_ranked_list.csv', header = T, row.names = 1, check.names = F)

ids <- row.names(df)
length(ids)

# load mapping files
id.maps <- read.csv("rna_dma_to_hsa_gn_mappings.tsv", header = F, row.names = 1, sep = '\t', stringsAsFactors = F)
#id.maps <- read.csv("polar_pos_pkl_to_kegg_annotations.tsv", header = F, row.names = 1, sep = '\t', #stringsAsFactors = F)

# mapping
id.maps.matched <- subset(id.maps, row.names(id.maps) %in% ids)
#id.maps.matched$V2 <- as.character(id.maps.matched$V2)
#id.maps.matched <- id.maps.matched[id.maps.matched$V2 != '',]
id.maps.matched <- subset(id.maps.matched, V2 != "")

id.maps.matched.ranked <- id.maps.matched[rownames(df)[rownames(df) %in% rownames(id.maps.matched)], , drop = FALSE]

#keep only unique ids
id.ranked.unique <- id.maps.matched.ranked %>% distinct(V2, .keep_all = TRUE)
mapped.ids <- unique(unlist(strsplit(id.maps.matched.ranked$V2, split = ';')))
length(mapped.ids)

df_ranked <- df %>% filter(rownames(.) %in% rownames(id.maps.matched.ranked)) %>% mutate(Gene_ID = id.maps.matched.ranked$V2 )%>% select(Gene_ID, Factor3) %>% separate_rows(Gene_ID, sep = ";")

# save to file
write.table(df_ranked, file = "ranked_genes.rnk", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)