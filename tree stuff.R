library(ape)
library(phytools)

# reading in the tree
shark_tree <- read.nexus('./data/611_taxa_tree topology for Emmaline.txt')
shark_tree
plot(shark_tree)
plot(shark_tree, type ='f', cex=.25, show.tip.label = F)

# removing unnecessary tips
shark_tips <- shark_tree$tip.label
shark_bin <- as.character(sapply(sapply(shark_tips, function(x) strsplit(x, '_', 
                         fixed = T)), function(y) paste(y[[1]], y[[2]])))
shark_bin_clean <- gsub("'", '', shark_bin, fixed = T)
shark_bin_clean <- gsub('.', '', shark_bin_clean, fixed = T)
shark_tree$tip.label <- shark_bin_clean
shark_tree
shark_bin_discard <- grep('[2-5]', shark_bin_clean, value = T)
shark_tree_clean <- drop.tip(shark_tree, shark_bin_discard)

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
