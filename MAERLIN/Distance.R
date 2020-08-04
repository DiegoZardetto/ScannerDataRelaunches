Distance <- function(file.a, file.b, match.vars = "All",
                     dist.funs = "Levenshtein", weights = NULL,
                     print = TRUE, ...){
############################################################
# Given two files and a set of common character columns as #
# matching variables, computes a normalized weighted       #
# distance matrix.                                         #
# If the input files are distinct a subsequent RECORD      #
# LINKAGE task is assumed, whereas if they are identical   #
# a DEDUPLICATION task on a single source is expected.     #
# Accordingly, the class of the output matrix will be      #
# either "reclink.Dist" or "dedup.Dist" (while still       #
# inheriting from class "matrix").                         #
# NOTE: If matching variables are not supplied by default  #
#       'All' common character columns are used.           #
# NOTE: 'dist.funs' is a character vector specifying the   #
#       distance to be applied to each matching variable   #
#       (recycled if necessary).                           #
#        Allowed distances are currently stored inside     #
#        vector dist.funs.list.                            #
# NOTE: If given, weights must be as many as the matching  #
#       variables (otherwise they get recycled).           #
#       If weights are not supplied by the user, the       #
#       function by default assigns the same weight to all #
#       the matching variables (w_i = 1).                  #
# NOTE: The dot-dot-dot argument stands for further        #
#       potential arguments that could be passed to the    #
#       selected 'dist.funs' (e.g. 'qcut' for emailBody).  # 
# NOTE: The actual formula for the computed normalized     #
#       weighted distance is obtained as the sequence      #
#       of the following two transformations:              #
#                                                          #
#       1. - Weigthed sum/mean of standardized individual  #
#            distances:                                    #
#            y = sum_i[ (x_i - m_i) / s_i ] * w_i          #
#                                                          #
#       2. - Normalization of y in [0,1]:                  #
#            z = ( y - min(y) ) / ( max(y) - min(y))       #
#                                                          #
#       where m_i and s_i stand for the mean and standard  #
#       deviation of the observed distances on the i-th    #
#       matching variable.                                 #
#       Notice that if there exist:                        #
#       1) a pair such that x_i = 0 for all i              #
#       2) a pair such that x_i = 1 for all i              #
#       then z equals the weigthed mean of individual      #
#       distances x_i with weights given by ww_i=w_i/s_i:  #
#                                                          #
#       z = sum_i[x_i * ( w_i / s_i )]/sum_i( w_i / s_i )  #
#                                                          #
# NOTE: If a subsequent deduplication task is assumed,     #
#       only the lower triangle of the distance matrix is  #
#       actually meaningful and will be returned           #
#       (remaining elements are NA). All statistics m_i    #
#       and s_i are accordingly computed only on the lower #
#       triangle.                                          #
# NOTE: White spaces are silently removed and letters      #
#       converted to lower-case, before computing any      #
#       distance values.                                   #
# NOTE: The distance matrix eventually undergoes a slight  #
#       regularization treatment: this is actually         #
#       mandatory for the subsequent fitting phase (since  #
#       Beta domain is (0,1) not [0,1]).                   #
#       The effect is that 0 (1) values are mapped into    #
#       value epsilon (eta), namely at the midpoint        #
#       between 0 (1) and the minimum (maximum) distance   #
#       value greater (smaller) than 0 (1).                #
# NOTE: An explicit normalization statement (like the one  #
#       at point 2. above) is actually needed ONLY if any  #
#       of the distance functions acting on matching       #
#       variables happens to be NOT ALREADY NORMALIZED IN  #
#       [0,1]. This should be avoided, whenever possible.  #
############################################################
dist.funs.list <- c("Levenshtein", "JaroWinkler", "Equality", "Three.grams",
                    "Cos.3.TfIdf", "VSM.TfIdf", "emailBody", "Diffabs", "FormatSD", "InRange", "simCos.3.TfIdf", "simVSM.TfIdf")

if (!is.data.frame(file.a) || !is.data.frame(file.b))
     stop("Files to be matched must be passed as data.frames")
if (nrow(file.a) > nrow(file.b))
    stop("File.a has more rows than file.b: please exchange them")
vars.a <- names(file.a)
vars.b <- names(file.b)
chars.a <- sapply(vars.a, function(var) is.character(file.a[, var]))
chars.b <- sapply(vars.b, function(var) is.character(file.b[, var]))
common <- intersect(vars.a[chars.a],vars.b[chars.b])
if (missing(match.vars) || match.vars=="All"){
    if (length(common) == 0)
        stop("No common character variables in files")
    }
else {
      if (!any(match.vars %in% common))
          stop("Matching variables must be character and common to both files")
      common <- match.vars[match.vars %in% common]
    }
if (!is.null(weights)){
    if (!is.numeric(weights))
        stop("If given, weights must be numeric")
    # Avoid negative weights
    weights <- abs(weights)
    w <- rep(weights, length.out = length(common))
    if (length(common) == 1){
        warning("Only one matching variable: its weight is 1", immediate. = TRUE)
        w <- 1
        }
    }
if (is.null(weights)){
    w <- rep(1, length.out = length(common))
    }
names(w) <- common

dist.funs <- rep(match.arg(dist.funs, dist.funs.list, several.ok = TRUE), length.out = length(common))
names(dist.funs) <- common

# Check if the function has been called from the command line
# or by another function (if so I cannot retrieve file names
# by using substitute())
directly <- !(length(sys.calls()) > 1)

# Check if a Deduplication job is in order
dedup <- identical(file.a, file.b)

# Compute a list with the distance matrices as components 
D.list  <- lapply(common, function(var)
                          do.call(what = dist.funs[var],
                                  args = list(x = file.a[, var], y = file.b[, var], dedup = dedup, ...)))

# Retrieve standard deviations vector
stdev <- sapply(D.list, function(el) attr(el, "stdev"))

# Retrieve NA frequency vector
NA.freq <- sapply(D.list, function(el) attr(el, "NAf"))

# Convert distance matrices list into an array
D.array <- array(unlist(D.list),dim=c(nrow(file.a),nrow(file.b),length(common)))

# Asses variable's weights.
# 1. If some stdev happens to be zero (i.e. all records from both files
#    are identical w.r.t. the variable, hence all distances are zero),
#    then leave the corresponding variable weight unchanged; # 18/10/2019 FIX 
#    then put the corresponding variable weight to zero;     # OLD VERSION

# 2. As mapping NAs to a given value may concentrate too much the
#    distance distribution, must soften the importance of variables
#    according to their NA rate;
ww <- (1 - NA.freq) * w/stdev
# ww[!is.finite(ww)] <- 0                # OLD VERSION
ww[!is.finite(ww)] <- w[!is.finite(ww)]  # 18/10/2019 FIX FOR SCANNER DATA
if (length(ww) > 1) { nww <- ww/sum(ww) } else { nww <- 1 }

# If more than 1 matching variable selected, compute distance as the weighted
# average of individual distances (current version uses array algebra, this is
# ~ 2.7 times faster than apply + weighted.mean)
D.array <- D.array * array(rep(nww, each=prod(dim(D.array)[1:2])), dim=dim(D.array))
D.mat <- rowSums(D.array, dims = 2)

# First normalize D.mat in [0,1], then if needed regularize it in (0,1)
has.cutoff <- FALSE
D.mat.min <- min(D.mat)
D.mat.max <- max(D.mat)
  # If Deduplication is in order, then restrict to lower triangle of matrix D.mat
  if (dedup) {
      D.mat.min <- min(D.mat[row(D.mat) > col(D.mat)])
      D.mat.max <- max(D.mat[row(D.mat) > col(D.mat)])
      }
# If all distances are equal (i.e. all zeros) leave D.mat unchanged...
if (!isTRUE(all.equal(D.mat.min, D.mat.max))){
    # ... otherwise, if needed, mormalize D.mat in [0,1]
    if ( (D.mat.min < 0) || (D.mat.min > 1) || (D.mat.max < 0) || (D.mat.max > 1) ) {
         # THIS CHECK ADDED 8/11/2010: if distances are natively big, should not
         # decrease them by translation, unless necessary...
         D.mat <- ( D.mat - D.mat.min ) / ( D.mat.max - D.mat.min )
        }
    # Don't regularize (very rare) cases where D.mat
    # has only values 0 and 1
    D.trim <- D.mat[ ( D.mat > 0 ) & ( D.mat < 1) ]
    if (dedup) {
        D.trim <- D.mat[ ( D.mat > 0 ) & ( D.mat < 1) & !is.na(D.mat) ]
        }
    npairs.trim <- length(D.trim)
    if (npairs.trim > 0){
        # Regularization in (0,1): perturb the original distance matrix
        # D.mat by mapping 0 -> epsilon and 1 -> eta
        epsilon <- 0.5 * min(D.trim)
        eta <- 1 - 0.5 * ( 1 - max(D.trim) )
        cutoff <- c(epsilon = epsilon, eta = eta)
        has.cutoff <- TRUE
        D.mat[D.mat <= 0] <- epsilon
        D.mat[D.mat >= 1] <- eta
        }
    }
# Thus output matrix will have a cutoff attribute
# if (and only if) its trimmed version (i.e. the
# one obtained by excluding 0 and 1 values) is
# NON-EMPTY. 

if (directly){
     attr(D.mat,"file.a") <- substitute(file.a)
     if (!dedup) {
         attr(D.mat,"file.b") <- substitute(file.b)
        }
    }
attr(D.mat,"match.vars") <- common
attr(D.mat,"dist.funs") <- dist.funs
attr(D.mat,"in.weights") <- w
attr(D.mat,"out.weights") <- ww
attr(D.mat,"NA.freq") <- NA.freq
attr(D.mat,"cutoff") <- if (has.cutoff) cutoff
attr(D.mat,"Call") <- sys.call()
class(D.mat) <- c(if (!dedup) "reclink.Dist" else "dedup.Dist", class(D.mat))
if (isTRUE(print)) {
     print(D.mat)
    }
invisible(D.mat)
}

print.reclink.Dist <- function(x, ...){
#########################################
# Print method for class "reclink.Dist" #
#########################################
cat("Matrix of pairwise distances for Record Linkage\n")
match.vars <- attr(x, "match.vars")
dist.funs <- attr(x, "dist.funs")
in.weights <- round(attr(x, "in.weights"), 3)
out.weights <- round(attr(x, "out.weights"), 3)
out.weights <- round(attr(x, "out.weights"), 3)
NA.freq <- round(attr(x, "NA.freq"), 3)
if (length(match.vars) == 1) {
     cat("- Matching variable:\t", match.vars, "\n")
     cat("- Distance function:\t", dist.funs, "\n")
    }
else{
     cat("- Matching variables and Distance functions:\n\n")
     vars.df <- data.frame(`Variable` = match.vars,
                           `Distance` = dist.funs,
                           `NA.freq` = NA.freq,
                           `Input.Weight` = in.weights,
                           `Output.Weight` = out.weights, stringsAsFactors = FALSE)
     vars.df <- vars.df[order(vars.df[["Output.Weight"]], decreasing = TRUE), ]
	 rownames(vars.df) <- NULL
     print(vars.df)
     cat("\n")
    }
cat(paste("- Matrix dimension:\t [", paste(dim(x), collapse = ", "),"]\n", sep = ""))
cat("- Pairs:\t\t", prod(dim(x)), "\n")
if (length(match.vars) == 1) {
     cat("- NA frequency:\t\t", NA.freq, "\n")
    }
if (!is.null(attr(x, "file.a"))) {
     cat("- Match files:\t\t ", deparse(attr(x, "file.a")), ", ", deparse(attr(x, "file.b")),"\n", sep="")
    }
#cat("- Distance distribution quartiles:\n")
#print(signif(quantile(x), 4))
cat("\n")
}

print.dedup.Dist <- function(x, ...){
#######################################
# Print method for class "dedup.Dist" #
#######################################
cat("Matrix of pairwise distances for Deduplication\n")
match.vars <- attr(x, "match.vars")
dist.funs <- attr(x, "dist.funs")
in.weights <- round(attr(x, "in.weights"), 3)
out.weights <- round(attr(x, "out.weights"), 3)
NA.freq <- round(attr(x, "NA.freq"), 3)
if (length(match.vars) == 1) {
     cat("- Matching variable:\t", match.vars, "\n")
     cat("- Distance function:\t", dist.funs, "\n")
    }
else{
     cat("- Matching variables and Distance functions:\n\n")
     vars.df <- data.frame(`Variable` = match.vars,
                           `Distance` = dist.funs,
                           `NA.freq` = NA.freq,
                           `Input.Weight` = in.weights,
                           `Output.Weight` = out.weights, stringsAsFactors = FALSE)
     vars.df <- vars.df[order(vars.df[["Output.Weight"]], decreasing = TRUE), ]
	 rownames(vars.df) <- NULL
     print(vars.df)
     cat("\n")
    }
cat(paste("- Matrix dimension:\t [", paste(dim(x), collapse = ", "),"]\n", sep = ""))
cat("- Pairs:\t\t", prod(nrow(x)/2, nrow(x) - 1), "\n")
if (length(match.vars) == 1) {
     cat("- NA frequency:\t\t", NA.freq, "\n")
    }
if (!is.null(attr(x, "file.a"))) {
     cat("- File:\t\t\t", deparse(attr(x, "file.a")),"\n")
    }
#cat("- Distance distribution quartiles:\n")
#print(signif(quantile(x[!is.na(x)]), 4))
cat("\n")
}