OTFS <- function(D.mat, true.sol, step = 0.005, plot = TRUE, stream.plot = FALSE, ...){
######################################################################
# Runs the OTFS (OPTIMAL THRESHOLD FULLY-SUPERVISED) Classifier.     #
# Given the distance matrix and the true solution to the RL problem, #
# determines the BEST distance threshold for a binary classifier     #
# wich declares to be M all pairs with distance less than or equal   #
# to the given threshold, and U the remaining pairs.                 #
# The BEST threshold is defined as the one such that the classifier  #
# above shows it's maximum F-measure.                                #
#                                                                    #
# NOTE: The BEST threshold is found by letting a cut value slide     #
#       along the [0,1] interval with steps of 'step' length.        #
#       All Precision, Recall and F-measure values obtained during   #
#       the sliding search are recorded and returned in output.      #
#       Moreover an F-measure plot and a Precision-Recall graph is   #
#       returned (if plot = TRUE).                                   #
######################################################################
steps <- seq(0, 1, by = step)
PRF <- screen(D.mat, true.sol, cut = steps, print = FALSE, plot = stream.plot)
F.max <- max(PRF[,"F-measure"], na.rm = TRUE)
i.max <- which.max(PRF[,"F-measure"])
d.max <- steps[i.max]
P.max <- PRF[i.max,"Precision"]
R.max <- PRF[i.max,"Recall"]
if (stream.plot) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
   }

PRF.list <- list(
                 PRF   = PRF,
                 steps = steps,
                 F.max = F.max,
                 d.max = d.max,
                 P.max = P.max,
                 R.max = R.max
                )
class(PRF.list) <- "OTFS"
if (plot)  plot(PRF.list, D.mat, true.sol, ...)
print(PRF.list)
}

plot.OTFS <- function(PRF, D.mat, true.sol, ...) {
################################
# Plot method for OTFS class.  #
################################
par(mfcol=c(1,2))
PRF.mat <- PRF$PRF
steps   <- PRF$steps
F.max   <- PRF$F.max
d.max   <- PRF$d.max
P.max   <- PRF$P.max
R.max   <- PRF$R.max
plot(steps, PRF.mat[,"F-measure"], xlim=0:1, ylim=0:1, ty="o", col="blue",
     xlab="Distance Threshold", ylab="F-measure", main="F-measure\nfor varying Distance Threshold")
points(d.max, F.max, col="red", pch=19, cex=4/3)
plot(PRF.mat[,"Recall"], PRF.mat[,"Precision"], xlim=0:1, ylim=0:1, ty="o", col="red",
     xlab="Recall", ylab="Precision", main="Precision and Recall\nfor varying Distance Threshold")
points(R.max, P.max, col="blue", pch=19, cex=4/3)
oask <- devAskNewPage(TRUE)
on.exit(devAskNewPage(oask))
invisible(screen(D.mat, true.sol, cut = d.max, print = FALSE, ...))
}

print.OTFS <- function(x, ...){
################################
# Print method for OTFS class. #
################################
PRF.mat <- x$PRF
steps   <- x$steps
F.max   <- x$F.max
d.max   <- x$d.max
P.max   <- x$P.max
R.max   <- x$R.max
cat("Best Distance Threshold:\n")
cat("Distance  = ", d.max, "\n")
cat("Precision = ", P.max, "\n")
cat("Recall    = ", R.max, "\n")
cat("F-measure = ", F.max, "\n")
cat("\n")
invisible(x)
}