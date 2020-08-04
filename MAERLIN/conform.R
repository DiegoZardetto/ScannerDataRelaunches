conformD <- function(object, D.mat){
###################################################
# Generic function. Checks if an object returned  #
# by a function whose input was a distance matrix #
# conforms to the passed distance matrix.         #
# NOTE: object may be the output of fit.mix.par,  #
#       EA.cluster, ...                           #
###################################################
if (!is.matrix(D.mat))
     stop("Wrong input object: second argument must be a matrix")
UseMethod("conformD")
}

conformD.reclink.mix.par <- function(object, D.mat){
#####################################################
# Conformity check with distance matrix: method for #
# RECORD LINKAGE fitted parameters.                 #
#####################################################
if (inherits(D.mat, "dedup.Dist"))
    stop("Mismatch: distance matrix built for deduplication, parameters fitted for record linkage")

if (!is.null(attr(object, "file.a")) && !is.null(attr(D.mat, "file.a"))) {
     if (!identical(attr(object, "file.a"),attr(D.mat, "file.a")) || 
         !identical(attr(object, "file.b"),attr(D.mat, "file.b"))) {
         stop("Mismatch: parameters were fit on a different distance matrix")
        }
    }

if (!identical(attr(object,"D.mat.dim"), dim(D.mat)))
    stop("Mismatch: parameters were fit on a different distance matrix")
}

conformD.dedup.mix.par <- function(mix.par, D.mat){
#####################################################
# Conformity check with distance matrix: method for #
# DEDUPLICATION fitted parameters.                  #
#####################################################
if (inherits(D.mat, "reclink.Dist"))
    stop("Mismatch: distance matrix built for record linkage, parameters fitted for deduplication")

if (!is.null(attr(mix.par, "file.a")) && !is.null(attr(D.mat, "file.a"))) {
     if (!identical(attr(mix.par, "file.a"), attr(D.mat, "file.a"))) {
         stop("Mismatch: parameters were fit on a different distance matrix")
        }
    }

if (!identical(attr(mix.par,"D.mat.dim"), dim(D.mat)))
    stop("Mismatch: parameters were fit on a different distance matrix")
}

conformD.reclink.ea.out <- function(object, D.mat){
#####################################################
# Conformity check with distance matrix: method for #
# RECORD LINKAGE solution from the Clustering EA.   #
#####################################################
if (inherits(D.mat, "dedup.Dist"))
    stop("Mismatch: distance matrix built for deduplication, solution found for record linkage")

if (!is.null(attr(object, "file.a")) && !is.null(attr(D.mat, "file.a"))) {
     if (!identical(attr(object, "file.a"),attr(D.mat, "file.a")) ||
         !identical(attr(object, "file.b"),attr(D.mat, "file.b"))) {
         stop("Mismatch: solution built using a different distance matrix")
        }
    }

if (!identical(attr(object,"D.mat.dim"), dim(D.mat)))
    stop("Mismatch: solution built using a different distance matrix")
}

conformD.dedup.ea.out <- function(object, D.mat){
#####################################################
# Conformity check with distance matrix: method for #
# RECORD LINKAGE solution from the Clustering EA.   #
#####################################################
if (inherits(D.mat, "reclink.Dist"))
    stop("Mismatch: distance matrix built for record linkage, solution found for deduplication")

if (!is.null(attr(object, "file.a")) && !is.null(attr(D.mat, "file.a"))) {
     if (!identical(attr(object, "file.a"),attr(D.mat, "file.a")) ||
         !identical(attr(object, "file.b"),attr(D.mat, "file.b"))) {
         stop("Mismatch: solution built using a different distance matrix")
        }
    }

if (!identical(attr(object,"D.mat.dim"), dim(D.mat)))
    stop("Mismatch: solution built using a different distance matrix")
}


conformF <- function(object, file.a, file.b){
###################################################
# Generic function. Checks if an object returned  #
# by a function whose input was tied to 2 files   #
# conforms to the passed files.                   #
# NOTE: object may be the output of Distance,     #
#       fit.mix.par, EA.cluster, ...              #
###################################################
if (!is.data.frame(file.a) || !is.data.frame(file.b))
     stop("Wrong input object: 2nd and 3rd arguments must be data.frames")
UseMethod("conformF")
}

conformF.reclink.Dist <- function(object, file.a, file.b){
##################################
# Method for RL distance matrix. #
##################################
m <- nrow(file.a) 
n <- nrow(file.b)
if (m > n){
     if (m==ncol(object) && n==nrow(object)) {
         stop("File.a has more rows than file.b: please exchange them")
        }
    }
else {
     if (m!=nrow(object) || n!=ncol(object)) {
         stop("Distance matrix does not conform to files")
        }
    }
}

conformF.dedup.Dist <- function(object, file.a, file.b){
#####################################
# Method for DEDUP distance matrix. #
#####################################
if (!identical(file.a, file.b))
     stop("Distance matrix built for deduplication: input files must be equal")
if (nrow(file.a)!=nrow(object))
     stop("Distance matrix does not conform to files")
}

conformF.reclink.mix.par <- function(object, file.a, file.b){
####################################
# Method for RL fitted parameters. #
####################################
m <- nrow(file.a) 
n <- nrow(file.b)
if (m > n){
     if (m == attr(object, "D.mat.dim")[2] && n == attr(object, "D.mat.dim")[1]) {
         stop("File.a has more rows than file.b: please exchange them")
        }
    }
else {
     if (m != attr(object, "D.mat.dim")[1] || n != attr(object, "D.mat.dim")[2]) {
         stop("Fitted parameters do not conform to files")
        }
    }
}

conformF.dedup.mix.par <- function(object, file.a, file.b){
#######################################
# Method for DEDUP fitted parameters. #
#######################################
if (!identical(file.a, file.b))
     stop("Distance matrix built for deduplication: input files must be equal")
if (nrow(file.a)!=attr(object, "D.mat.dim")[1])
     stop("Fitted parameteres matrix does not conform to files")
}

conformF.ea.out <- conformF.reclink.mix.par
####################################
# Method for Clustering EA output. #
####################################


#conformF.trcl <- function(trcl, file){
#########################################
# Method for transitive closure output. #
# MUST BE WRITTEN!!!!!!!!!!!!!!!!!!!!!  #
#########################################
#}