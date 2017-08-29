library(ape)
library(phytools)
library(picante)

# reading in the tree
shark_tree <- read.nexus('./data/611_taxa_tree topology for Emmaline.txt')
shark_tree
plot(shark_tree)
plot(shark_tree, type ='f', cex=.25, show.tip.label = F)

# creating species list from range data
species_names <- dir('./data/polygon')
species_names <- sub(pattern = '.json', "", species_names)

# removing unnecessary tips
shark_tips <- shark_tree$tip.label
shark_bin <- as.character(sapply(sapply(shark_tips, function(x) strsplit(x, '_', 
                         fixed = T)), function(y) paste(y[[1]], y[[2]])))
shark_bin_clean <- gsub("'", '', shark_bin, fixed = T)
shark_bin_clean <- gsub('.', '', shark_bin_clean, fixed = T)
shark_tree$tip.label <- shark_bin_clean
shark_tree
shark_bin_discard <- grep('[2-5]', shark_bin_clean, value = T)
cf_names <- grep('cf', shark_tree$tip.label, value = T)
shark_tree_clean <- drop.tip(shark_tree, shark_bin_discard)
shark_tree_clean$tip.label = sub('cf', '', shark_tree_clean$tip.label)
sum(species_names %in% shark_tree_clean$tip.label)
drop_names <- shark_tree_clean$tip.label[!(shark_tree_clean$tip.label %in% species_names)]
shark_tree_clean <- drop.tip(shark_tree_clean, drop_names)
shark_tree_clean$tip.label
genus_names <- sapply(strsplit(species_names, ' '), function(x) x[[1]])
tree_genera <- sapply(strsplit(shark_tree_clean$tip.label, ' '), function(x) x[[1]])
sum(unique(genus_names) %in% unique(tree_genera))
length(unique(genus_names))
genus_names[!(genus_names %in% tree_genera)]

# rooting the shark tree
plot(shark_tree_clean, show.tip.label = F)
identify.phylo(shark_tree_clean)
is.rooted(shark_tree_clean)
shark_tree_clean <- root(shark_tree_clean, node = 575, resolve.root = T)
is.rooted(shark_tree_clean)
is.binary(shark_tree_clean)

# potential method to make tree ultrmetric
is.binary.tree(shark_tree_clean)
attempt_tree <- multi2di(shark_tree)
is.binary.tree(attempt_tree)
is.rooted(attempt_tree)
attempt_tree <- chronopl(attempt_tree, lambda = 0)
plot(attempt_tree, show.tip.label = F)

# phylogenetic diversity
val_res_list <- vector("list", length = 7)
mat_list <- vector("list", length = 7)
val_list <- vector("list", length = length(sp_files))
for (j in seq_along(c(1, 2, 3, 4, 5, 6, 7))) {
for (i in seq_along(sp_files)) {
  val_list[[i]] <- raster_res_list[[j]][[i]]@data@values
}
  val_res_list[[j]] <- val_list
  pixels <- (1:length(val_res_list[[j]][[j]]))
  com_mat <- data.frame(val_res_list[[j]], pixels, row.names = pixels)
  mat_list[[j]] <- com_mat
}


