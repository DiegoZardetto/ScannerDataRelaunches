screen <- function (D.mat, ...){
#######################################################
# Generic function. Diagnostics on computed pairwise  #
# distances when analysed by using the true solution  #
# of the (RL or DD) problem.                          #
# An underlying RECORD LINKAGE or DEDUPLICATION task  #
# is assumed.                                         #
#######################################################
UseMethod("screen")
}




screen.reclink.Dist <- function(D.mat, true.sol, cut = NULL,
                                print = TRUE, plot = TRUE, breaks = seq(0, 1, 0.02), ...){
#################################################################
# Method for RECORD LINKAGE.                                    #
# If cut = NULL the function has 2 effects:                     #
# 1. If plot = TRUE, plots the histograms of the observed       #
#    distances between pairs classified by the true assignments #
#    to the M and U classes.                                    #
# 2. Returns a list containing the vectors with the observed    #
#    distances for M and U pairs.                               #
#                                                               #
# If cut = dist.th, were dist.th is a distance threshold        #
# in [0,1], the function:                                       #
# 1. Returns the performance (in terms of Precision, Recall     #
#    and F-measure) of the binary classifier wich declares to   #
#    be M all pairs with distance less than or equal to the     #
#    given threshold dist.th, and U the remaining pairs.        #
# 2. If plot = TRUE, plots the histograms of the observed       #
#    distances between pairs classified by the true assignments #
#    to the M and U classes + the cut value vertical line + the #
#    counts of the False Positives and False Negatives          #
#    determined by the above classifier.                        #
#################################################################

# If needed prepare the plot window
if (plot) par(mfcol=c(2,1))

D.mat.dim <- dim(D.mat)

# Transform true solution into matrix form
t.sol <- sol.to.mat(true.sol, D.mat.dim)

# Match pairs distances
n.M <- nrow(t.sol)
D.M <- as.numeric(D.mat[t.sol])
D.mat2 <- D.mat
D.mat2[t.sol] <- NA

# Unmatch pairs distances
n.U <- prod(D.mat.dim) - n.M
D.U <- as.numeric(D.mat[!is.na(D.mat2)])

# Put all into a list
out <- list(D.M=D.M, D.U=D.U)

if (is.null(cut)){
    if (plot) {
        # Histograms
        main.M <- paste("Matches\n[", n.M," pairs]",sep="")
        hist(D.M, col = "lavender", xlim = c(0,1), xlab = "distance", main = main.M, breaks = breaks, ...)
        main.U <- paste("Unmatches\n[", n.U," pairs]",sep="")
        hist(D.U, col =  "peachpuff", xlim = c(0,1), xlab = "distance", main = main.U, breaks = breaks, ...)
        # Density plots
        par(mfcol=c(1,1))
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
        U.dens <- density(D.U)
        M.dens <- density(D.M)
        ylim <- range(U.dens$y, M.dens$y) 
        plot(U.dens, xlim=0:1, ylim=ylim, xlab="Distance", ylab="Density",
             main="Density plot for Matches (blue) and Unmatches (red)",
             lty=2, col="red", lwd=2)
        lines(M.dens, lty=1, col="blue", lwd=2)
        }
    }
else {
      names(D.M) <- rep("M", n.M)
      names(D.U) <- rep("U", n.U)
      # Put together labelled distances 
      dist <- c(D.M, D.U)
      
      classifier <- function(cut, plot){
      #########################################################################
      # Evaluate the performance of a classifier that declares a positive for #
      # every pair with distance less than or equal to a given cut threshold  #
      # Precision, Recall and F-measure scores are returned.                  #
      #########################################################################
      
	  positives <- (dist <= cut)
      # Performance measures
      p <- sum(positives)
      tp <- sum(names(dist[positives])=="M")
      prec <- tp / p
      rec  <- tp / n.M
      F <- 2 * (tp / (p + n.M))
      n.fn <- sum(names(dist[!positives])=="M")
      n.fp <- sum(names(dist[positives])=="U")
      if (print){
          cat("Distance Threshold = ", cut, "\n")
          cat("Precision = ", prec, "\n")
          cat("Recall    = ", rec, "\n")
          cat("F-measure = ", F, "\n")
          cat("\n")
          }
      if (plot) {
          main.fn <- paste("False Negatives\n[", n.fn," pairs]",sep="")
          hist(D.M, col = "lavender", xlim = c(0,1), xlab = "distance", main = main.fn, breaks = breaks, ...)
          abline(v=cut, col="red", lwd=2, lty=2)
          main.fp <- paste("False Positives\n[", n.fp," pairs]",sep="")
          hist(D.U, col =  "peachpuff", xlim = c(0,1), xlab = "distance", main = main.fp, breaks = breaks, ...)
          abline(v=cut, col="red", lwd=2, lty=2)
          }
      if (print & ( (n.fn == 0) & (n.fp == 0) ) ) {
         cat("Perfect Decision!\n")
         }
      score <- c(prec, rec, F)
      names(score) <- c("Precision", "Recall", "F-measure")
	  score
      }

     # Now compute PRF scores for all thresholds of the cut vector
     score <- sapply(cut, function(dist.th) classifier(cut = dist.th, plot = plot))
     score <- t(score)
     }
# Return output object
if (is.null(cut)) invisible(out) else invisible(score)
}


screen.dedup.Dist <- function(D.mat, true.sol, cut = NULL,
                              print = TRUE, plot = TRUE, breaks = seq(0, 1, 0.02), ...){
#################################################################
# Method for DEDUPLICATION.                                     #
# If cut = NULL the function has 2 effects:                     #
# 1. If plot = TRUE, plots the histograms of the observed       #
#    distances between pairs classified by the true assignments #
#    to the M and U classes.                                    #
# 2. Returns a list containing the vectors with the observed    #
#    distances for M and U pairs.                               #
#                                                               #
# If cut = dist.th, were dist.th is a distance threshold        #
# in [0,1], the function:                                       #
# 1. Returns the performance (in terms of Precision, Recall     #
#    and F-measure) of the binary classifier wich declares to   #
#    be M all pairs with distance less than or equal to the     #
#    given threshold dist.th, and U the remaining pairs.        #
# 2. If plot = TRUE, plots the histograms of the observed       #
#    distances between pairs classified by the true assignments #
#    to the M and U classes + the cut value vertical line + the #
#    counts of the False Positives and False Negatives          #
#    determined by the above classifier.                        #
#################################################################

# If needed prepare the plot window
if (plot) par(mfcol=c(2,1))

D.mat.dim <- dim(D.mat)

# Transform true solution into matrix form
t.sol <- sol.to.mat(true.sol, D.mat.dim)

# Match pairs distances
n.M <- nrow(t.sol)
D.M <- as.numeric(D.mat[t.sol])
D.mat2 <- D.mat
D.mat2[t.sol] <- NA

# Unmatch pairs distances
n.U <- prod(D.mat.dim[1], D.mat.dim[1] - 1)/2 - n.M
D.U <- as.numeric(D.mat[!is.na(D.mat2) & (row(D.mat) > col(D.mat))])

# Put all into a list
out <- list(D.M=D.M, D.U=D.U)

if (is.null(cut)){
    if (plot) {
        # Histograms
        main.M <- paste("Matches\n[", n.M," pairs]",sep="")
        hist(D.M, col = "lavender", xlim = c(0,1), xlab = "distance", main = main.M, breaks = breaks, ...)
        main.U <- paste("Unmatches\n[", n.U," pairs]",sep="")
        hist(D.U, col =  "peachpuff", xlim = c(0,1), xlab = "distance", main = main.U, breaks = breaks, ...)
        # Density plots
        par(mfcol=c(1,1))
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
        U.dens <- density(D.U)
        M.dens <- density(D.M)
        ylim <- range(U.dens$y, M.dens$y) 
        plot(U.dens, xlim=0:1, ylim=ylim, xlab="Distance", ylab="Density",
             main="Density plot for Matches (blue) and Unmatches (red)",
             lty=2, col="red", lwd=2)
        lines(M.dens, lty=1, col="blue", lwd=2)
        }
    }
else {
      names(D.M) <- rep("M", n.M)
      names(D.U) <- rep("U", n.U)
      # Put together labelled distances 
      dist <- c(D.M, D.U)
      
      classifier <- function(cut, plot){
      #########################################################################
      # Evaluate the performance of a classifier that declares a positive for #
      # every pair with distance less than or equal to a given cut threshold  #
      # Precision, Recall and F-measure scores are returned.                  #
      #########################################################################
      
	  positives <- (dist <= cut)
      # Performance measures
      p <- sum(positives)
      tp <- sum(names(dist[positives])=="M")
      prec <- tp / p
      rec  <- tp / n.M
      F <- 2 * (tp / (p + n.M))
      n.fn <- sum(names(dist[!positives])=="M")
      n.fp <- sum(names(dist[positives])=="U")
      if (print){
          cat("Distance Threshold = ", cut, "\n")
          cat("Precision = ", prec, "\n")
          cat("Recall    = ", rec, "\n")
          cat("F-measure = ", F, "\n")
          cat("\n")
          }
      if (plot) {
          main.fn <- paste("False Negatives\n[", n.fn," pairs]",sep="")
          hist(D.M, col = "lavender", xlim = c(0,1), xlab = "distance", main = main.fn, breaks = breaks, ...)
          abline(v=cut, col="red", lwd=2, lty=2)
          main.fp <- paste("False Positives\n[", n.fp," pairs]",sep="")
          hist(D.U, col =  "peachpuff", xlim = c(0,1), xlab = "distance", main = main.fp, breaks = breaks, ...)
          abline(v=cut, col="red", lwd=2, lty=2)
          }
      if (print & ( (n.fn == 0) & (n.fp == 0) ) ) {
         cat("Perfect Decision!\n")
         }
      score <- c(prec, rec, F)
      names(score) <- c("Precision", "Recall", "F-measure")
	  score
      }

     # Now compute PRF scores for all thresholds of the cut vector
     score <- sapply(cut, function(dist.th) classifier(cut = dist.th, plot = plot))
     score <- t(score)
     }
# Return output object
if (is.null(cut)) invisible(out) else invisible(score)
}