performance <- function(result, true.sol, 
                        file.a = NULL, file.b = NULL, D.mat = NULL, mix.par = NULL,
                        print = TRUE, plot = complete, complete = FALSE,
                        pairs = c("data.frame", "list")){
####################################################
# Compute the quality of the achieved RL result by #
# comparing it to THE true solution (tipically a   #
# benchmark solution assumed as "gold-standard").  #
# NOTE: 'result' can be an output of either        #
#       fit.mix.par or EA.cluster.                 #
####################################################


# Get from input object the achieved solution (in matrix form)
# and the Distance Matrix dimension
sol <- getMatSol(result)
D.mat.dim <- attr(result,"D.mat.dim")

# Convert true.sol into matrix form
t.sol <- sol.to.mat(true.sol, D.mat.dim)

# Compute basic quality measures
Q <- perf(sol, t.sol, D.mat.dim)

# Compare fitness values for exact and best genetic solutions
# exact.fit <- newfit(true.sol, mix.par, D.mat)
# fitness.comp <- c(best.fit, exact.fit)
# names(fitness.comp) <- c("Achieved solution", "True solution")

if  (!complete){
#    Q <- c(Q, list(fitness = fitness.comp))
     print(Q)
    }
else{
     # Check conformity among mandatory arguments
     if (is.null(file.a) || is.null(file.b))
         stop("When 'complete' is TRUE must specify both file.a and file.b")
     has.mix.par <- FALSE
     has.D.mat <- FALSE
     conformF(result, file.a, file.b)
     if (inherits(result, "ea.out")){
         if ( xor(is.null(D.mat), is.null(mix.par)) )
              stop("For results of class 'ea.out' D.mat and mix.par must be passed both or be both NULL!")
        }
     if (!is.null(D.mat)) {
         conformD(result, D.mat)
         has.D.mat <- TRUE
        }
     if (!is.null(mix.par)) {
         conformF(mix.par, file.a, file.b)
         actual.mix.par <- mix.par
         has.mix.par <- TRUE
        }
     if (inherits(result, "reclink.mix.par") || inherits(result, "dedup.mix.par")) {
         actual.mix.par <- result
         has.mix.par <- TRUE
         if (!is.null(mix.par)) {
             if (!identical(actual.mix.par, mix.par))
                 warning("Passed mix.par does not agree with input object: ignored", immediate. = TRUE)
            }
        }

     pairs <- match.arg(pairs)
     common <- intersect(names(file.a), names(file.b))

     getPairs <- function(mat){
         #######################################
         # Extracts record pairs corresponding #
         # to rows of a MAP-like matrix.       #
         # NOTE: most inputs via lex. scoping. #
         #######################################
         n.links <- nrow(mat)
         a.l <- mat[, 1]
         fa.l <- file.a[a.l, common]
         rownames(fa.l) <- paste("a", rownames(fa.l), sep=".")
         b.l <- mat[, 2]
         fb.l <- file.b[b.l, common]
         rownames(fb.l) <- paste("b", rownames(fb.l), sep=".")
         links.df <- rbind(fa.l, fb.l)
         if (has.mix.par && has.D.mat) {
             fm.d <- sapply( seq(along.with = a.l),
                                 function (i) {
                                 actual.mix.par[["fm.d"]](D.mat[a.l[i], b.l[i]])
                                }
                            )
             links.df[["fm.d"]] <- rep(fm.d, 2)
            }
         links.df <- links.df[rep(1:n.links, each = 2) + c(0, n.links),]
         switch(pairs, "data.frame" = links.df,
                       "list"       = lapply(seq(1, nrow(links.df), by = 2),
					                         function(i) links.df[c(i, i+1), ]))
        }

     ################################################
     # Build a tp, fn and fp lists of record pairs. #
     ################################################
     # Add convenience a third column to sol and t.sol matrices
     sol   <- cbind(  sol, 1)
     t.sol <- cbind(t.sol, 2)
     mm <- merge(t.sol, sol, by=1:2, all=TRUE)
     # Compute lists, sizes and pairs 
       # tp list
       tp <- mm[!is.na(mm[,3]) & !is.na(mm[,4]), 1:2]
       n.tp <- nrow(tp)
       tp.list <- if (n.tp > 0) getPairs(tp)
       # fn list
       fn <- mm[is.na(mm[,4]), 1:2]
       n.fn <- nrow(fn)
       fn.list <- if (n.fn > 0) getPairs(fn)
       # fp list
       fp <- mm[is.na(mm[,3]), 1:2]
       n.fp <- nrow(fp)
       fp.list <- if (n.fp > 0) getPairs(fp)

     # Bind performance measures and output lists
     Q <- c(Q, list(tp.list = tp.list, fn.list = fn.list, fp.list = fp.list))
     if (print) print(Q[1:5])
    }
invisible(Q)
}

perf <- function(sol, true.sol, D.mat.dim){
#####################################################
# Compute the quality of a RL solution by comparing #
# it to THE true solution (tipically a benchmark    #
# solution assumed as "gold-standard").             #
# NOTE: 'sol' and 'true.sol' format is free, both   #
#       vector and matrix representations work.     #
# VALUE: a list with components:                    #
#        1) Confusion Matrix                        #
#        2) Precision                               #
#        3) Recall                                  #
#        4) False Positives Rate                    #
#        5) F-Measure                               #
#####################################################
# Transform true and achieved solution into matrix form
t.sol <- sol.to.mat(true.sol, D.mat.dim)
sol <- sol.to.mat(sol, D.mat.dim)

# Initialize confusion matrix to be filled
confusion <- matrix(NA, 3, 3, dimnames = list(c("p", "n", "total"), c("M", "U", "total")))
	
# Compute confusion matrix pieces
  # Pairs
    npairs <- confusion[3, 3] <- prod(D.mat.dim)
  # True Match
    M <- confusion[3, 1] <- nrow(t.sol)
  # True Unmatch
    U <- confusion[3, 2] <- npairs - M

# Compute remaining confusion matrix pieces
  # Merge true and achieved solution matrices
    tp.mat <- merge(t.sol, sol)
  # Positives
    p <- confusion[1, 3] <- nrow(sol)
  # True positives
    tp <- confusion[1, 1] <- nrow(tp.mat)
  # False positives
    fp <- confusion[1, 2] <- p - tp
  # Negatives
    n <- confusion[2, 3] <- npairs - p
  # True negatives
    tn <- confusion[2, 2] <- U - fp
  # False negatives
    fn <- confusion[2, 1] <- M - tp

# Compute performance measures
  # precision = tp/(tp+fp)
    prec <- tp/(tp+fp)
    names(prec) <- "precision = tp/(tp+fp)"
  # recall = tp rate = tp/M
    rec <- tp / M
    names(rec) <- "recall = tp rate = tp/M"
  # fp rate = fp/U
    fp.rate <- fp / U
    names(fp.rate) <- "fp rate = fp/U"
  # F-measure = 2/(1/precision+1/recall)
    F <- 2 / (1/prec + 1/rec)
    names(F) <- "F-measure = 2/(1/precision+1/recall)"
Q <- list("confusion matrix" = confusion, prec, rec, fp.rate, F)
Q
}