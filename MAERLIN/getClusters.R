getClusters <- function(trcl, ...){
########################################################
# Extracts clusters of records that have been declared #
# Matches by the achieved solution of a Deduplication  #
# task, after a transitive closure.                    #
########################################################
UseMethod("getClusters")
}

getClusters.default <- function(trcl, ...){
####################################
# Default method: raises an error. #
####################################
stop("Input must be the transitive closure of Deduplication task solution!")
}

getClusters.trcl <- function(trcl, file, value = c("data.frame", "list"), ...){
#######################################
# Reconstruct a dataframe (or a list) #
# of clusters from a transitive       #
# closure.                            #
#######################################

# Conformity check: TO BE IMPLEMENTED!!!!!!!
# conformF(trcl, file)
value <- match.arg(value)
dim.clust <- sapply(trcl, length)
n <- length(trcl)
clust.id <- rep(1:n, times=dim.clust)
rec.id <- unlist(trcl)
df.clust <- data.frame(cluster.id = clust.id, record.id = rec.id, file[rec.id, , drop = FALSE])
switch(value, "data.frame" = df.clust,
              "list"       = split(df.clust, f = clust.id))
}