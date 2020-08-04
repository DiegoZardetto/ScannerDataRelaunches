trcl <- function(tab){
################################################################
# Finds the transitive closure of a set of match-pairs arising #
# from a N:M DEDUPLICATION task.                               #
# NOTE: Match-pairs are internally represented as a two        #
#       column table ('tab'): each row identifies a pair.      #
# NOTE: Since 'tab' arises from a Deduplication task, we know  #
#       that surely tab[i,1]!=tab[i,2] for all i.              #
# NOTE: The return value is a list whose components are        #
#       clusters of Matching elements appearing in the         #
#       transitive closure.                                    #
################################################################

n <- nrow(tab)
# initialize first cluster
clust <- list(tab[1, ])
# loop on pairs
for (i in 2:n){
     # verify if elements of the pair already appear in some cluster
	 p.r <- tab[i, 1]
	 cl.r <- which(sapply(clust, function(cl) p.r %in% cl))
	 p.c <- tab[i, 2]
	 cl.c <- which(sapply(clust, function(cl) p.c %in% cl))
     if ( length(cl.r) && !length(cl.c)) clust[[cl.r]] <- c(clust[[cl.r]], p.c)
     if (!length(cl.r) &&  length(cl.c)) clust[[cl.c]] <- c(clust[[cl.c]], p.r)
     if ( length(cl.r) &&  length(cl.c)) {
	     if (cl.r != cl.c){
	         clust[[cl.r]] <- c(clust[[cl.r]], clust[[cl.c]])
		     clust[[cl.c]] <- NULL
            }
        }
     if (!length(cl.r) && !length(cl.c)) clust[[length(clust)+1]] <- c(p.r, p.c)
    }
clust
}


clust.df <- function(clust, df.mails){
#######################################
# Reconstruct a dataframe of clusters #
# from a transitive closure.          #
#######################################
dim.clust <- sapply(clust, length)
n <- length(clust)
clust.id <- rep(1:n, times=dim.clust)
clust.list <- lapply(clust, function(e) df.mails[e,])
df.clust <- NULL
sapply(clust.list, function(e) df.clust <<- rbind(df.clust, e))
df.clust <- data.frame(cluster.id = clust.id, order=rownames(df.clust), df.clust)
df.clust
}