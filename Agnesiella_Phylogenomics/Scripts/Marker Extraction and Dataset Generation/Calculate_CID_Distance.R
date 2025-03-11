#!/usr/bin/env Rscript
# Author: WJJ v2024.12.1

# Usage:
# Rscript Calculate_CID_Distance.R <working_directory> <tree_file>

library("TreeTools")
library("TreeDist")

# Retrieve command-line arguments
args <- commandArgs(TRUE)

# Set working directory
setwd(args[1])

# Read tree file
trees <- read.tree(args[2])


# Remove node labels and edge lengths
trees[] <- lapply(trees, "[[<-", "node.label", NULL)
trees[] <- lapply(trees, "[[<-", "edge.length", NULL)
trees <- TreeTools::Preorder(trees)

# Calculate clustering information distance
distance <- ClusteringInfoDistance(trees)

# Convert distance to matrix
distance_matrix <- as.matrix(distance)

# Write distance matrix to CSV
write.csv(distance_matrix, "3-CID_distance_matrix.csv", row.names = TRUE)