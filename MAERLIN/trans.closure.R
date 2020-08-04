trans.closure <- function(object, ...){
################################################################
# Finds the transitive closure of a set of match-pairs arising #
# from a DEDUPLICATION task.                                   #
# Methods must handle the N:M and 1:1 matching cases.          #
################################################################
UseMethod("trans.closure")
}

trans.closure.default <- function(object, ...){
####################################
# Default method: raises an error. #
####################################
stop("Transitive Closure only applies to solutions of Deduplication tasks!")
}

trans.closure.dedup.mix.par <- function(object, D.mat,
                                        merge = c("auto", "never", "always")){
################################################################
# Finds a transitive (quasi)closure of the set of match-pairs  #
# arising from a N:M DEDUPLICATION task.                       #
# NOTE: 'object' must be the solution to some Deduplication    #
#       task with N:M constraints.                             #
#       Match-pairs are internally represented as a two        #
#       column table ('tab'): each row identifies a pair.      #
# NOTE: Since 'tab' arises from a Deduplication task, we know  #
#       that surely tab[i,1]!=tab[i,2] for all i.              #
# NOTE: The return value is a list whose components are        #
#       clusters of Matching elements appearing in the         #
#       transitive closure.                                    #
# NOTE: Parameter 'merge' specifies how the quasi-closure      #
#       dynamic algorithm handles the choice of merging two    #
#       already generated clusters:                            #
#       - "auto"   -> clusters can be merged when a strong     #
#                     'link' exists (i.e. they are linked by   #
#                     a pair with high posterior probability)  #
#                     and the overall 'affinity' of the pairs  #
#                     in the merged cluster is not too low.    #
#       - "always" -> gives the ordinary transitive closure.   #
#       - "never"  -> clusters cannot be merged, they can at   #
#                     most grow by capturing new records.      #
################################################################

# Conformity check
conformD(object, D.mat)

# Extract solution in matrix form
tab <- getMatSol(object)

# Are there declared positives?
n <- nrow(tab)
if (n==0) {
     cat("No Matches have been declared!\n")
     return(invisible(NULL))
    }

# Retrieve selected merge method
merge <- match.arg(merge)

# Posterior Match probability of declared matches
ppost.tab <- object[["fm.d"]](D.mat[tab])

# Order solution by decreasing ppost
tab <- tab[order(ppost.tab, decreasing = TRUE), ]

# Order ppost by decreasing ppost: each pair is tied to its actual ppost
ppost.tab <- ppost.tab[order(ppost.tab, decreasing = TRUE)]

# initialize first cluster (the more likely to be a true one)
clust <- list(tab[1, ])

# THIS ONLY FOR DEBUG/TUNING!!!!!!!!!!!!!!
#t0 <- proc.time()[3]
#t.iter <- 0
#nclust.iter <- 1
#avgsize.iter <- 2
#maxsize.iter <- 2
#merge.iter <- NULL
#merge.repuls <- NULL
#merge.link <- NULL

# loop on pairs
for (i in 2:n){
     # verify if elements of the pair already appear in some cluster
     p.r <- tab[i, 1]
     cl.r <- which(sapply(clust, function(cl) p.r %in% cl))
     p.c <- tab[i, 2]
     cl.c <- which(sapply(clust, function(cl) p.c %in% cl))
     # If either already belongs to a cluster BUT the other does not, add the
     # the second to the preexisting cluster
     if ( length(cl.r) && !length(cl.c)) clust[[cl.r]] <- c(clust[[cl.r]], p.c)
     if (!length(cl.r) &&  length(cl.c)) clust[[cl.c]] <- c(clust[[cl.c]], p.r)
     # If both already belong to some cluster BUT the clusters differ,
     # verify if the clusters have to be merged or not and with what method
     if ( ( merge != "never" ) && ( length(cl.r) &&  length(cl.c) ) ) {
         if (cl.r != cl.c){
             # find what cluster comes first (the earlier, the more reliable) 
             cl.1 <- min(cl.r, cl.c)
             cl.2 <- max(cl.r, cl.c)
             # candidate cluster: the merged one
             clmerge <- c(clust[[cl.1]], clust[[cl.2]])
             if ( merge == "auto" ){
                 # retrieve the ppost of the pair that "pushes" to merge the clusters
                 plink <- ppost.tab[i]
                 # generate all the meaningful pairs in the candidate merged cluster
                 merge.pairs <- t(combn(clmerge, 2, sort, decreasing = TRUE))
                 # compute its overall affinity score and the corresponding
                 # repulsion between the original clusters:
                 affinity <- mean(object[["fm.d"]](D.mat[merge.pairs]))
                 repulsion <- 3/2 - affinity 

                 # THIS ONLY FOR DEBUG/TUNING!!!!!!!!!!!!!!
                 #merge.repuls <- c(merge.repuls, repulsion)
                 #merge.link <- c(merge.link, plink)
                 #merge.iter <- c(merge.iter, i)

                 # compare plink with repulsion and decide accordingly
                 if (plink > repulsion) {
                     clust[[cl.1]] <- clmerge
                     clust[[cl.2]] <- NULL
                    }
                }
             else {
                 clust[[cl.1]] <- clmerge
                 clust[[cl.2]] <- NULL
                }
            }
        }
     # If neither already belongs to a cluster, put the pair into a NEW cluster
     if (!length(cl.r) && !length(cl.c)) clust[[length(clust)+1]] <- c(p.r, p.c)

     # THIS ONLY FOR DEBUG/TUNING!!!!!!!!!!!!!!
     #t.iter <- c(t.iter, proc.time()[3] - t0)
     #nclust.iter <- c(nclust.iter, length(clust))
     #avgsize.iter <- c(avgsize.iter, mean(sapply(clust,length)))
     #maxsize.iter <- c(maxsize.iter, max(sapply(clust,length)))
	}
# THIS ONLY FOR DEBUG/TUNING!!!!!!!!!!!!!!
#attr(clust,"repuls") <- merge.repuls
#attr(clust,"link") <- merge.link
#attr(clust,"iter") <- merge.iter
#attr(clust,"t.iter") <- t.iter
#attr(clust,"nclust.iter") <- nclust.iter
#attr(clust,"avgsize.iter") <- avgsize.iter
#attr(clust,"maxsize.iter") <- maxsize.iter

class(clust) <- c("trcl",class(clust))
clust
}


trans.closure.dedup.ea.out <- function(object, D.mat, mix.par){
################################################################
# Finds a transitive (quasi)closure of the set of match-pairs  #
# arising from a 1:1 DEDUPLICATION task.                       #
# Note: if object has class dedup.ea.out, this means that 1:1  #
#       deduplication must hold. But, even from a 1:1 solution #
#       from a "matrix standpoint", clusters of size greater   #
#       than 2 may arise via transitive closure.               #
#                                                              #
#       To avoid such wrong result, the merge = "never"        #
#       option of trans.closure.dedup.mix.par is applied,      #
#       and in addition clusters are not free to grow by       #
#       capturing new records.                                 #
################################################################ 

# Conformity check
conformD(object, D.mat)
conformD(mix.par, D.mat)

# Extract solution in matrix form
tab <- getMatSol(object)

# Are there declared positives?
n <- nrow(tab)
if (n==0) {
     cat("No Matches have been declared!\n")
     return(invisible(NULL))
    }

# Posterior Match probability of declared matches
ppost.tab <- mix.par[["fm.d"]](D.mat[tab])

# Order solution by decreasing ppost
tab <- tab[order(ppost.tab, decreasing = TRUE), ]

# Order ppost by decreasing ppost: each pair is tied to its actual ppost
ppost.tab <- ppost.tab[order(ppost.tab, decreasing = TRUE)]

# initialize first cluster (the more likely to be a true one)
clust <- list(tab[1, ])

# loop on pairs
for (i in 2:n){
     # verify if elements of the pair already appear in some cluster
     p.r <- tab[i, 1]
     cl.r <- which(sapply(clust, function(cl) p.r %in% cl))
     p.c <- tab[i, 2]
     cl.c <- which(sapply(clust, function(cl) p.c %in% cl))
     # 1) If either already belongs to a cluster BUT the other does not,
     #    DO NOT add the second to the preexisting cluster
     # 2) Even if both already belong to some cluster BUT the clusters differ,
     #    clusters must not be merged!

     # If neither already belongs to a cluster, put the pair into a NEW cluster
     if (!length(cl.r) && !length(cl.c)) clust[[length(clust)+1]] <- c(p.r, p.c)
	}
class(clust) <- c("trcl",class(clust))
clust
}