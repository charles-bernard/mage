# Wrapper of all the previous internal functions
# -----------------------------------------------------------------
#' @title get_classified_asso
#'
#' @description After having called the \code{explore_clusters} function, this function proposes
#' to retrieve ALL the associations corresponding to any assigned classification in the 
#' 'associations_classification_template.csv' file.
#'
#' @param trees_info the trees info returned by \code{explore_clusters}
#' @param classification_file the filled in classification template produced by \code{explore_clusters}
#' @param scores_tab the table of scores
#'
#' @return A list of tables, each of them corresponding to the associations of a specific classification 
#'
#' @importFrom data.table fread
get_classified_asso <- function(trees_info, classification_file, scores_tab) {
  classif_data <- fread(classification_file, header = T, sep = ",");
  classifs <- unique(clust_classif[`Classification` != ""]$Classification);
  classified_asso <- list();
  
  for(i in 1:length(classifs)) {
    curr_classif <- classifs[i];
    curr_classif_data <- classif_data[`Classification` == curr_classif, ];
    classif_clusts <- sort(unique(curr_classif_data$`Cluster`));
    classif_ix <- NULL;
    for(j in 1:length(classif_clusts)) {
      curr_nodes <- 
        curr_classif_data[`Cluster` == classif_clusts[j], `Node`];
      curr_tree_ix <- which(names(trees_info) == classif_clusts[j]);
      classif_ix <- c(classif_ix, which(trees_info[[curr_tree_ix]]$pt_nodes %in% curr_nodes));
    }
    classified_asso[[i]] <- data.table(scores_tab[classif_ix, ]);
    names(classified_asso)[i] <- curr_classif;
  }
  
  return(classified_asso);
}