getSolution <- function(object){
####################################################
# Extracts the data structure storing the achieved #
# solution for a RL job.                           #
# Methods must handle distinct objects, e.g.       #
# fitted parameters for N:M matching, or           #
# Clustering EA return objects.                    #
####################################################
UseMethod("getSolution")
}

getSolution.default <- function(object){
####################################
# Default method: raises an error. #
####################################
stop("No solution found inside this object!")
}

getSolution.ea.out <- function(object){
###########################################
# getSolution method for class ea.out.    #
# NOTE: Return value is the best genotype #
#       found by the Clustering EA, hence #
#       a vector whose length depends on  #
#       matching constraints.             #
###########################################
# Extract best genotype from EA output
best.fit <- max(object[["last"]][["fit"]])
best.ind <- which.max(object[["last"]][["fit"]])
best.sol <- object[["last"]][["pop"]][best.ind,]
attr(best.sol,"fitness") <- best.fit
best.sol 
}

getSolution.reclink.mix.par.sol.NtoM <- function(object){
##################################################
# getSolution method for class mix.par.sol.NtoM. #
# NOTE: Return value is the MAP matrix.          #
##################################################
# If RL job was N:M, then return MAP pairs table
object$MAP.l
}

getSolution.reclink.mix.par.sol.1to1 <- function(object){
##################################################
# getSolution method for class mix.par.sol.1to1. #
# NOTE: Return value is a vector.                #
##################################################
# If RL job was 1:1, extract best genotype from mix.par.sol
# directly.
object$solution
}

getSolution.reclink.mix.par.sol.1toN <- getSolution.reclink.mix.par.sol.Nto1 <- getSolution.reclink.mix.par.sol.1to1
####################################################
# getSolution methods for classes mix.par.sol.1toN #
# and mix.par.sol.Nto1 (identical to 1:1).         #
# NOTE: Return value is a vector.                  #
####################################################

getSolution.dedup.mix.par.sol.NtoM <- function(object){
##################################################
# getSolution method for class mix.par.sol.NtoM. #
# NOTE: Return value is the MAP matrix.          #
##################################################
# If RL job was N:M, then return MAP pairs table
object$MAP.l
}

getSolution.dedup.mix.par.sol.1to1 <- function(object){
##################################################
# getSolution method for class mix.par.sol.1to1. #
# NOTE: Return value is a vector.                #
##################################################
# If RL job was 1:1, extract best genotype from mix.par.sol
# directly.
object$solution
}


getMatSol <- function(object){
####################################################
# Extracts the data structure storing the achieved #
# solution for a RL job, in matrix form (identical #
# to the one of MAPs).                             #
####################################################
sol <- getSolution(object)
sol.to.mat(sol, attr(object, "D.mat.dim"))
}



sol.to.mat <- function(sol, D.mat.dim, check = FALSE){
#####################################################
# Converts the data structure storing the achieved  #
# solution for a RL job into matrix form (identical #
# to the one of MAPs).                              #
# NOTE: If 'sol' does not come from a fit.mix.par   #
#       or EA.cluster output (e.g. it is some       #
#       benchmark solution built elsewhere) a check #
#       is made:                                    #
#    i) on its type and dimensions                  #
#   ii) if a matrix, on possible duplicated rows.   #
#####################################################
# SHOULD write some checks on input types, dimensions...

# If sol is a matrix (e.g. it stores MAPs)
# return sol...
if (is.matrix(sol)) {
     if (check){
         if (any(duplicated(sol))) {
             sol <- unique(sol)
             warning("Input solution contained duplicated rows: they have been removed!", immediate. = TRUE)
            }
        }
     colnames(sol) <- c("row", "col")
     return(sol)
    }
else {
     #... else transform sol into a MAP-like
     # matrix:
     on <- sol!=0
     n.link <- sum(on)
     if (length(sol)==D.mat.dim[1]) {
         a.l <- which(on)
         b.l <- sol[on]
        }
     else {
         b.l <- which(on)
         a.l <- sol[on]
        }
     mat <- matrix(c(a.l, b.l), n.link, 2)
     colnames(mat) <- c("row", "col")
     return(mat)
    }
}