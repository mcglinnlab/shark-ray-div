library(ape)
shark_tree <- read.nexus('./data/611_taxa_tree topology for Emmaline.txt')
shark_tree
plot(shark_tree)
plot(shark_tree, type ='f', cex=.25, show.tip.label = F)

shark_tips <- shark_tree$tip.label
shark_bin <- as.character(sapply(sapply(shark_tips, function(x) strsplit(x, '_', 
                         fixed = T)), function(y) paste(y[[1]], y[[2]])))
shark_bin_clean <- gsub("'", '', shark_bin, fixed = T)
shark_bin_clean <- gsub('.', '', shark_bin_clean, fixed = T)
shark_tree$tip.label <- shark_bin_clean
shark_tree
shark_bin_discard <- grep('[2-5]', shark_bin_clean, value = T)
shark_tree_clean <- drop.tip(shark_tree, shark_bin_discard)
shark_tree_clean$tip.label
plot(shark_tree_clean, cex = .25)
