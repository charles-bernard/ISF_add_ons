library(data.table)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(stats)

##############################################################
# SETUP ######################################################
##############################################################
# PATH TO THE DIRECTORY THAT STORES THE COG ENRICHMENT OF ALL
# GENE FAMILIES
enrich_dir <- "/home/charles/clones/ISF_add_ons/foutoir/enrich";
# PATH TO THE OUTPUT DIRECTORY
output_dir <- "/home/charles/clones/ISF_add_ons/foutoir/consensus_clustering";

# Retrieve enrichment files
enrich_files <- Sys.glob(file.path(enrich_dir, "*_COG_enrichment.txt"));

# Retrieve name of each gene family base on enrichment filenames
basenames <- basename(enrich_files);
family_names <- unlist(strsplit(basenames, "*_COG_enrichment.txt"));

# Load statistics into a matrix
mat = matrix(data = 0, nrow = 26, ncol = length(enrich_files));
for(i in 1:length(enrich_files)) {
  curr_data <- fread(enrich_files[i], sep = "\t", header = TRUE);
  mat[, i] <- curr_data[, 2][[1]];
}
rownames(mat) <- paste(curr_data[, 1][[1]], curr_data[, 4][[1]]);
colnames(mat) <- family_names;

##############################################################
# CONSENSUS CLUSTERING #######################################
##############################################################
wd <- getwd();
setwd(output_dir);

# simulated matrix
# mat <- matrix(data = 0, ncol = 20, nrow = nrow(curr_data));
# rownames(mat) <- paste(curr_data[, 4][[1]], curr_data[, 1][[1]]);
# for(i in 1:ncol(mat)) {
#   x <- runif(nrow(curr_data));
#   mat[, i] <- x/sum(x)*100;
# }
# colnames(mat) <- paste("family", 1:20)

ccl <- ConsensusClusterPlus(mat, reps = 500, pItem = 0.8, maxK = 10, 
                            clusterAlg = 'km', distance = 'euclidean',
                            plot = 'pdf', writeTable = TRUE);
setwd(wd);

###############################################################
# Take a look at consensus.pdf in the output directory        #
# and set the optimal number of cluster 'k' below             #
###############################################################
k <- 2;

##############################################################
# AGGLOMERATIVE CLUSTERING ###################################
##############################################################
# Create the partition into k clusters with kmean
cl <- kmeans(t(mat), k, iter.max = 500);
partition <- cl$cluster;


###############################################################
# HEATMAP WITH HIERARCHICAL CLUSTERING OF COLUMNS #############
###############################################################
# WEITHER THE COG CATEGORY (APPEARING AS ROWS IN THE HEATMAP)
# MUST BE CLUSTERED ACCORDING TO THEIR FONCTIONAL SUPER
# CATEGORY OR NOT
cluster_rows <- TRUE;
group_COG <- TRUE;

supercat = curr_data[, 3][[1]];
supercat[which(supercat == "NO COG ASSOCIATED")] <- "NA\n";
supercat[which(supercat == "INFORMATION STORAGE AND PROCESSING")] <- "Information storage\nand processing";
supercat[which(supercat == "CELLULAR PROCESSES AND SIGNALING")] <- "Cellular processes\nand signaling";
supercat[which(supercat == "METABOLISM")] <- "Metabolism\n";
supercat[which(supercat == "POORLY CHARACTERIZED")] <- "Poorly\nChara.";

supercat_col = c("white", "lightblue", "salmon", "lightgreen", "lightgrey");
names(supercat_col) <- unique(supercat);

if(group_COG) {
  split <- supercat;
} else {
  split <- NULL;
}

cluster_colors <- rainbow(n=k); 
names(cluster_colors) <- sort(unique(partition));

column_ha <- HeatmapAnnotation('Family clusters' = as.character(partition),
                               col = list('Family clusters' = cluster_colors),
                               annotation_legend_param = 
                                 list(nrow = 1, egend_direction = "horizontal"));

pdf(file.path(output_dir, "heatmap.pdf"),
    width = 10,
    height = 9);
hm <- 
  Heatmap(matrix = mat,
          name = "% of enrichment",
          split = split,
          row_gap = unit(2, "mm"),
          border = TRUE,
          rect_gp = gpar(col = "black", lwd = 1),
          cluster_rows = cluster_rows,
          clustering_distance_rows = "pearson",
          clustering_method_rows = "ward.D2",
          show_row_dend = FALSE,
          row_title_side = "left",
          row_title_gp = gpar(fill = "white", fontsize = 8.5, alpha = 0),
          row_names_side = "right",
          row_names_max_width = unit(10, "cm"),
          row_names_gp = gpar(fontsize = 9),
          column_title = "COG category enrichment (in %) per family of genes",
          cluster_columns = TRUE,
          clustering_distance_columns = "pearson",
          clustering_method_columns = "ward.D2",
          column_names_gp = gpar(fontsize = 11),
          top_annotation = column_ha,
          heatmap_legend_param = list(legend_direction = "horizontal"))
hm <- 
  draw(hm, merge_legend = TRUE,
       heatmap_legend_side = "bottom", 
       annotation_legend_side = "bottom");
ro = row_order(hm);

for(i in 1:length(unique(supercat))) {
  decorate_row_title("% of enrichment", {
    grid.rect(gp = gpar(fill = supercat_col[names(ro)[i]]));
    grid.text(names(ro)[i], rot = 90, gp = gpar(fontsize = 8.5));
  }, slice = i)
}
junk <- dev.off();
