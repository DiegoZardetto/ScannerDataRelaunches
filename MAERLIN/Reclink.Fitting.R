fit.mix.par <- function (D.mat, ...){
#################################################
# Generic function. Fits the parameters of my   #
# inflated-beta mixture model on a specified    #
# distance matrix by exploiting my perturbative #
# ML optimization algorithm.                    #
# An underlying RECORD LINKAGE or DEDUPLICATION #
# task is assumed.                              #
#################################################
UseMethod("fit.mix.par")
}


fit.mix.par.reclink.Dist <-
function(D.mat, match.constr = c("1:1", "1:N", "N:1", "N:M"),
         guessM = c(c.am = 1/2, c.bm = 2), force = TRUE,
         verbose = FALSE, plot = TRUE, print = TRUE, ...){
###################################################################
# Given a normalized distance matrix for record pairs, uses my    #
# perturbative fitting technique to compute Maximum Likelihood    #
# estimates for the parameters modeling the following             #
# two components inflated-beta mixture (a RECORD LINKAGE task is  #
# assumed):                                                       #
#                                                                 #
# f(x) = pm * fm(x; am, bm, p0) + (1 - pm) * fu(x; au, bu, p1)    #
#                                                                 #
# where:                                                          #
#                                                                 #
# fm(x; am, bm, p0) = p0 * delta(x,0) + (1-p0) * dbeta(x; am, bm) #
# fu(x; au, bu, p1) = p1 * delta(x,1) + (1-p1) * dbeta(x; au, bu) #
#                                                                 #
# Thus the model to fit includes from 5 (pm, am, bm, au, bu) up   #
# to 7 (the former 5 plus p0 and p1) parameters, being P0 (P1)    #
# non-null only if zero (one) pairwise distances exist.           #
#                                                                 #
# NOTE: First order perturbation theory is used to split the      #
#       original constrained optimization problem into two        #
#       major sequential steps:                                   #
#       STEP 1: neglect Matches contributions and estimate the 2  #
#               parameters of Unmatches f.u distribution (au,bu); #
#       STEP 2: use estimates obtained in STEP 1 to estimate the  #
#               3 remaining f.m parameters (pm,am,bm) focusing on #
#               parameter-space regions where Matches are more    #
#               likely to be found.                               #
#       Lastly, if needed, the algorithm estimates P0 and P1.     #
# NOTE: 'match.constr' specifies matching constraints.            #
#       One-to-One, One-to-Many, Many-to-One, and Many-to-Many    #
#       matching restrictions can all be handled.                 #
#       These map to assertions on whether file.a and file.b have #
#       duplicates:                                               #
#       - 1:1  ->  ( (file.a !has.dup) & (file.b !has.dup) )      #
#       - 1:N  ->  ( (file.a !has.dup) & (file.b  has.dup) )      #
#       - N:1  ->  ( (file.a  has.dup) & (file.b !has.dup) )      #
#       - N:M  ->  ( (file.a  has.dup) & (file.b  has.dup) )      #
#       Each one of these constraints translates into specific    #
#       restrictions in parameter-space (see below).              #
# NOTE: 'guessM' selects starting values for STEP 2 parameters.   #
#       See below for the default value, which is expected to be  #
#       a good choice almost always, since results depend only    #
#       very weakly on it.                                        #
# NOTE: 'force' asks to return the best achieved approximation    #
#       when STEP 2 optimization execution of nlminb() reaches    #
#       either the iterations or the function evaluations limits. #
# NOTE: 'verbose' enables printing of optimizations return codes. #
# NOTE: 'plot' enables plotting useful graphics.                  #
# NOTE: 'eps' is a lower, positive bound for beta distributions   #
#       shape parameters (if a=0 or b=0 dbeta() returns NaN).     #
# NOTE: Imposed constraints ensure that:                          #
#       STEP 1: au in (eps, Inf)                                  #
#               bu in (eps, Inf)                                  #
#       STEP 2: pm in [0, pm.max]                                 #
#               am in (eps, ~bu)                                  #
#               bm in (~au, Inf)                                  #
#       where last two bounds of both steps guarantee desirable   #
#       behaviour of M and U distribution (i.e. M dominating U at #
#       dist ~ 0, U dominating M at dist ~ 1, and M and U having  #
#       a small overlap).                                         #
#       pm.max is the maximum match probability given files       #
#       (i.e. D.mat) dimensions and matching constraints.         #
#       Roughly speaking:                                         #
#       - pm.max = 1/ncol(D.mat) if match.constr = "1:1"          #
#       - pm.max = 1/nrow(D.mat) if match.constr = "1:N"          #
#       - pm.max = 1/ncol(D.mat) if match.constr = "N:1"          #
#       - pm.max = 1             if match.constr = "N:M"          #
# NOTE: By default uses method of moments estimates as starting   #
#       point for STEP 1 optimization:                            #
#       STEP 1: au=E.MOM(au), bu=E.MOM(bu)                        #
#       For STEP 2, if 'guessM' is not provided, default initial  #
#       guess is as follows:                                      #
#       STEP 2: pm=pm.max/2, am=c.am*bu=bu/2, bm=c.bm*au=2*au     #
# NOTE: '...' further parameters to be passed to nlminb().        #
###################################################################

# Check and retrieve files dimensions
if ( ( m <- nrow(D.mat) ) > ( n <- ncol(D.mat) ) )
    stop("Rows must be less than or equal to columns: please transpose the Distance Matrix")
npairs <- prod(m, n)

# Beta parameters a and b must be greater than 0  
eps <- 1e-4

# Retrieve requested matching constraints
match.constr <- match.arg(match.constr)

# Check for guessM type, length, and values
if (!missing(guessM)){
     if (!is.numeric(guessM) || (length(guessM) != 2) ||
         guessM[1] < eps || guessM[1] > 1 || guessM[2] < 1){
         stop("guessM parameters must be numeric and must satisfy: 0 < c.am <= 1 and c.bm >= 1")
        }
}

# Two possible alternatives:
#  (i) D.mat contains only values 0 and 1     (very rare cases)
# (ii) trimmed version of D.mat is non-empty  (almost always)
# Must check and manage properly both cases (i) and (ii).
cutoff <- attr(D.mat, "cutoff")
if (!is.null(cutoff)){
     # Standard cases (ii)
     zero <- epsilon <- cutoff["epsilon"]
     one  <- eta     <- cutoff["eta"]
    }
else {
     # Rare cases (i)
     zero <- 0; epsilon <- 1e-4
     one  <- 1; eta     <- 1 - epsilon
    }

# Check for possible pairs at distance 0 or 1:
# if any, trim the distance distribution
n.0 <- sum(D.mat <= zero)
n.1 <- sum(D.mat >= one)
D.trim <- D.mat[ ( D.mat > zero ) & ( D.mat < one) ]
npairs.trim <- length(D.trim)

# Compute the maximum allowed value for the expected number of Matches,
# given current matching constraints
warn.msg <- "Too many zero distance pairs: either too few match variables were used or files contain duplicates!"
n.M.max <- switch(match.constr,
                  "1:1" = if (n.0 <= m) m else {warning(warn.msg, immediate. = TRUE); (npairs - n.1)},
                  "1:N" = if (n.0 <= n) n else {warning(warn.msg, immediate. = TRUE); (npairs - n.1)},
                  "N:1" = if (n.0 <= m) m else {warning(warn.msg, immediate. = TRUE); (npairs - n.1)},
                  "N:M" = npairs - n.1)

# Deduce the upper bound for the mixing parameter
# w.r.t. the trimmed distance distribution
pm.max  <- ( n.M.max - n.0 ) / npairs.trim

if (npairs.trim > 0){
# -- Standard treatment (case (ii)) START
# Use MOM estimates for (au, bu) as starting point
# for STEP 1 optimisation

MOM <- function(d){
############################
# MOM estimates for beta   #
# shape parameters (a, b). #
############################
n <- length(d)
# Handle the case of just one trimmed distance point
if (n == 1) {
     d <- rep(d, 3)
     n <- 3
    }
scale  <- (n - 1) / n
# Handle the case of identical trimmed distances
if (!(var(d) > 0)) {
     d <- jitter(d, amount = eps)
    }
v <- scale * var(d)
me <- mean(d)
onemme <- 1 - me
estimate <- numeric(2)
estimate[1] <-     me * (me * onemme / v - 1)
estimate[2] <- onemme * (me * onemme / v - 1)
estimate
}

guessU <- MOM(as.vector(D.trim))

# Check for STEP 1 starting point sanity and soundness
if (guessU[1] < eps || guessU[2] < eps) 
    stop("Bad MOM estimates for U parameters: au = ", guessU[1], " , bu = ", guessU[2])
if (guessU[1] < guessU[2]) 
    warning("MOM estimates for U conditional distribution put more mass at dist ~ 0 than at dist ~ 1: au = ",
            guessU[1], " , bu = ", guessU[2], immediate. = TRUE)

# Some useful 'global' values which must not be re-computed in optimization loops:
log.dist <- log(D.trim)
sum.log.dist <- sum(log.dist)
log.1mdist <- log(1-D.trim)
sum.log.1mdist <- sum(log.1mdist)

# Define STEP 1 objective:
# turn STEP1 max logLikelihood problem into minimization
fit.u.par <- function(par){
au <- par[1]
bu <- par[2]
log.p <- dbeta(D.trim,au,bu,log=TRUE)
-sum(log.p)
}
# Define STEP 1 objective gradient
grad.fit.u.par <- function(par){
au <- par[1]
bu <- par[2]
do.dau <- - sum.log.dist   + npairs.trim*(digamma(au) - digamma(au+bu))
do.dbu <- - sum.log.1mdist + npairs.trim*(digamma(bu) - digamma(au+bu))
c(do.dau, do.dbu)
}

# STEP 1 optimization: estimate au and bu
  # Faesible region boundaries
  lower1 <- c(eps, eps)
step1.out <- nlminb(start = guessU, objective = fit.u.par, gradient = grad.fit.u.par,
                    lower = lower1, ...)
  # Check for converge
  if (step1.out$convergence!=0){
      print(step1.out)
      stop("STEP1 optimization failed: try different starting values")
      }
u.par.best <- step1.out$par
au <- u.par.best[1]
bu <- u.par.best[2]
  # Check that estimates are meaningful and sound
  if (au < eps || bu < eps)
      stop("Bad U parameters: au = ", au, " , bu = ", bu,
           " (eps = ", eps, ")")
  if (au < bu)
      warning("U conditional distribution puts more mass at dist ~ 0 than at dist ~ 1: au = ", au, " , bu = ", bu, immediate. = TRUE)
  # Add out gradient value to optimization log
  step1.out <- c(step1.out, list( out.gradient = grad.fit.u.par(c(au, bu)) ) )
  # Check for high (i.e. greater than 'little') out gradient values over the scale npairs.trim
  # when estimated parameters are far (i.e. more than 'near' away) from boundaries
  near   <- 1E-4
  little <- 1E-2
  if (
      all( abs( u.par.best - lower1 ) > near ) && 
      any( abs( step1.out[["out.gradient"]] )/npairs.trim > little )
     ) warning("Oddly high out gradient in STEP1 optimization: check optimization log!", immediate. = TRUE)
# STEP 1 optimization: END


# Useful 'global' values which must not be re-computed in STEP 2 optimization loop:
f.u  <- dbeta(D.trim,au,bu)

# Define STEP 2 objective:
# turn STEP 2 max logLikelihood problem into minimization
fit.m.par <- function(par){
pm <- par[1]
am <- par[2]
bm <- par[3]
p <- pm*dbeta(D.trim,am,bm)+(1-pm)*f.u
-sum(log(p))
}
# Define STEP 2 objective gradient
grad.fit.m.par <- function(par){
pm <- par[1]
am <- par[2]
bm <- par[3]
f.m <- dbeta(D.trim,am,bm)
f <- pm*f.m+(1-pm)*f.u
do.dpm <- - sum( ( f.m-f.u ) / f )
do.dam <- - sum( (pm*f.m / f) * ( log.dist   - (digamma(am) - digamma(am+bm)) ) )
do.dbm <- - sum( (pm*f.m / f) * ( log.1mdist - (digamma(bm) - digamma(am+bm)) ) )
c(do.dpm, do.dam, do.dbm)
}

# STEP 2 optimization: estimate pm, am and bm
  # Faesible region boundaries
  offset <- 0.1 * ((au-bu)/2)
  lower2  <- c(     0,          eps,  au - offset)
  upper2  <- c(pm.max,  bu + offset,          Inf)
  
  # If no starting guess is provided for STEP 2 optimisation
  # use default values (2 -> 10 in former settings...)
  guessM.in <- guessM
  guessM    <- numeric(3)
  guessM[1] <- min( pm.max / 2, 1 / sqrt(npairs.trim) )
  guessM[2] <- bu * guessM.in[1]
  guessM[3] <- au * guessM.in[2]
  names(guessM) <- c("pm","am","bm")

step2.out <-  nlminb(start = guessM, objective = fit.m.par, gradient = grad.fit.m.par, 
                     lower = lower2, upper = upper2, ...)
optim.out <- list(step1.out = step1.out, step2.out = step2.out)
  # Check for converge
  if (step2.out$convergence!=0){
      if ( force && (step2.out$message %in% c("iteration limit reached without convergence (9)",
                                              "function evaluation limit reached without convergence (9)")) ) {
          warning("STEP2 limit reached", immediate. = TRUE)
          }
      else {
           if ( !(step2.out$message %in% c("absolute function convergence (6)",
                                           "singular convergence (7)")) ){
               print(optim.out)
               stop("STEP2 optimization failed: try different starting values for guessM")
               }
           }
    }
m.par.best <- step2.out$par
pm <- m.par.best[1]
am <- m.par.best[2]
bm <- m.par.best[3]
  # Check that estimates are meaningful and sound
  if (pm < 0 || pm > pm.max || am < eps)
      stop("Bad M parameters: pm = ", pm, " , am = ", am, " , bm = ", bm,
              " (pm.max = ", pm.max, " , eps = ", eps, ")")
  if (bm < bu)
      warning("M conditional distribution dominates U at dist ~ 1: bm = ",
              bm, " , bu = ", bu, immediate. = TRUE)
  if (am > au)
      warning("U conditional distribution dominates M at dist ~ 0: am = ",
              am, " , au = ", au, immediate. = TRUE)
  if (am > bm)
      warning("M conditional distribution puts more mass at dist ~ 1 than at dist ~ 0: am = ",
              am, " , bm = ", bm, immediate. = TRUE)
  # Add out gradient value to optimization log
  step2.out <- c(step2.out, list(out.gradient = grad.fit.m.par(c(pm, am, bm)) ) )
  optim.out <- list(step1.out = step1.out, step2.out = step2.out)
  # Check for high (i.e. greater than 'little') out gradient values over the scale npairs.trim
  # when estimated parameters are far (i.e. more than 'near' away) from boundaries
  if (
      all( abs( m.par.best - lower2) > near ) &&
      all( abs( m.par.best - upper2) > near ) &&
      any( abs( step2.out[["out.gradient"]] )/npairs.trim > little )
     ) warning("Oddly high out gradient in STEP2 optimization: check optimization log!", immediate. = TRUE)
# STEP 2 optimization: END

if (isTRUE(verbose)) print(optim.out) 

# Build output pieces
  # Basic fitted parameters
    best.par <- c(m.par.best,u.par.best)
    names(best.par) <- c("pm","am","bm","au","bu")
  # Now compute parameters describing the "spikes" of M and U distributions
  # at distance 0 and 1 respectively  
    best.par["pm"] <- n.0 / npairs + pm * ( npairs.trim / npairs )
    pm <- best.par["pm"]
    p0 <- if (pm > 0) n.0 / ( pm * npairs ) else 0
    p1 <- n.1 / ( (1 - pm ) * npairs )
  # Update best parameters vector
    best.par <- c(best.par, p0, p1)
    names(best.par) <- c("pm","am","bm","au","bu","p0","p1")
# -- Standard treatment END
}
else {
# -- Rare treatment (case (i)) START
# (cases with npairs.trim=0)
  pm <- n.0/npairs
  # If all distances are zero then define pm slightly less then 1
  # (otherwise the logLikelihood has only values Inf and -Inf)
  if (pm >= 1) pm <- 1 - 1E-6
  # This are convenience values for zero overall mass betas
  # (never used in what follows, e.g. by the clustering EA gen())
  am <- 1
  bm <- 1
  au <- 1
  bu <- 1
  # As above: p0 and p1 cannot be exactly equal to 1
  # (otherwise the logLikelihood has only values Inf and -Inf)
  p0 <- 1 - 1E-6
  p1 <- p0
  best.par <- c(pm,am,bm,au,bu,p0,p1)
  names(best.par) <- c("pm","am","bm","au","bu","p0","p1")
  guessU <- c(au, bu) 
  guessM <- c(pm, am, bm)
  optim.out <- "Empty trimmed distance matrix: only 0,1 found"
# Rare treatment END
}

  # Fitted distributions and loglikelihood ratio
    f.m    <- function(d) ifelse(d > epsilon, (1 - p0) * dbeta(d, am, bm), p0)
    f.u    <- function(d) ifelse(d < eta, (1 - p1) * dbeta(d, au, bu), p1)
    f      <- function(d) pm*f.m(d)+(1-pm)*f.u(d)
    fm.d   <- function(d) pm*f.m(d)/f(d)
    fu.d   <- function(d) (1-pm)*f.u(d)/f(d)
    loglik <- function(d) log( ( pm * f.m(d) ) / ( (1-pm) * f.u(d) ) )
  # Maximum A-Posteriori (MAP) expected number of Matches given observed distances
    MAP.l <- which( fm.d(D.mat) >= fu.d(D.mat), arr.ind = TRUE )
    # Check if MAP pairs fulfill matching constraints: if it is the case
    # the RL solution has been FOUND!

      # Initialization: suppose there's no solution
      solution <- NULL
      # Search duplicated file.a records in MAPs
      MAP.a <- table(MAP.l[, 1])
      MAP.dup.a <- MAP.a[MAP.a > 1]
      a.has.dup <- length(MAP.dup.a) > 0
      if (a.has.dup) attr(MAP.l, "MAP.dup.a") <- MAP.dup.a
      # Search duplicated file.b records in MAPs
      MAP.b <- table(MAP.l[, 2])
      MAP.dup.b <- MAP.b[MAP.b > 1]
      b.has.dup <- length(MAP.dup.b) > 0
      if (b.has.dup) attr(MAP.l, "MAP.dup.b") <- MAP.dup.b
      # Now recall that duplicated a(b) records in MAP
      # imply that b(a) has duplicates!
# T T
      if (a.has.dup && b.has.dup) {
          if (match.constr=="N:M") {
              solution <- "MAPs"
              cls <- c("reclink.mix.par.sol.NtoM", "reclink.mix.par.sol")
            }
        }
# T F
      if (a.has.dup && !b.has.dup) {
          if (match.constr=="1:N") {
              solution <- rep(0, n)
              mapply( function(i, j) solution[i] <<- j, MAP.l[, 2], MAP.l[, 1] )
              cls <- c("reclink.mix.par.sol.1toN", "reclink.mix.par.sol")
            }
          if (match.constr=="N:M") {
              solution <- "MAPs"
              cls <- c("reclink.mix.par.sol.NtoM", "reclink.mix.par.sol")
            }
        }
# F T
      if (!a.has.dup && b.has.dup) {
          if (match.constr=="N:1") {
              solution <- rep(0, m)
              mapply( function(i, j) solution[i] <<- j, MAP.l[, 1], MAP.l[, 2] )
              cls <- c("reclink.mix.par.sol.Nto1", "reclink.mix.par.sol")
            }
          if (match.constr=="N:M") {
              solution <- "MAPs"
              cls <- c("reclink.mix.par.sol.NtoM", "reclink.mix.par.sol")
            }
        }
# F F
      if (!a.has.dup && !b.has.dup) {
          if (match.constr=="1:1") {
              solution <- rep(0, m)
              mapply( function(i, j) solution[i] <<- j, MAP.l[, 1], MAP.l[, 2] )
              cls <- c("reclink.mix.par.sol.1to1", "reclink.mix.par.sol")
            }
          if (match.constr=="N:1") {
              solution <- rep(0, m)
              mapply( function(i, j) solution[i] <<- j, MAP.l[, 1], MAP.l[, 2] )
              cls <- c("reclink.mix.par.sol.Nto1", "reclink.mix.par.sol")
            }
          if (match.constr=="1:N") {
              solution <- rep(0, n)
              mapply( function(i, j) solution[i] <<- j, MAP.l[, 2], MAP.l[, 1] )
              cls <- c("reclink.mix.par.sol.1toN", "reclink.mix.par.sol")
            }
          if (match.constr=="N:M") {
              solution <- "MAPs"
              cls <- c("reclink.mix.par.sol.NtoM", "reclink.mix.par.sol")
            }
        }

    EC.M <- nrow(MAP.l)
  # Unconditional expected number of Matches given files (i.e. D.mat) dimensions
    E.M <- ceiling(pm*npairs)
# Build output list
out <- list(
     best.par = best.par,
     f      = f,
     f.m    = f.m,
     f.u    = f.u,
     fm.d   = fm.d,
     fu.d   = fu.d,
     loglik = loglik,
     MAP.l  = MAP.l,
     solution = solution,
     EC.M   = EC.M,
     E.M    = E.M
    )

if (!is.null(attr(D.mat, "file.a"))) {
     attr(out, "file.a") <- attr(D.mat, "file.a")
     attr(out, "file.b") <- attr(D.mat, "file.b")
    }
attr(out,"D.mat.dim") <- dim(D.mat)
attr(out,"match.constr") <- match.constr
attr(out,"guessU") <- guessU
attr(out,"guessM") <- guessM
attr(out,"optim.log") <- optim.out
attr(out,"Call") <- sys.call(-1)
match.c <- switch(match.constr,
                  "1:1" = "1to1",
                  "1:N" = "1toN",
                  "N:1" = "Nto1",
                  "N:M" = "NtoM")
class(out) <- c(paste("reclink.mix.par", match.c, sep="."), "reclink.mix.par")
if ( !is.null(out[["solution"]]) ) {
     class(out) <- c(cls, class(out))
    }
if (isTRUE(print)) {
     try(print(out))
    }
if (isTRUE(plot)) {
     try(plot(out, D.mat))
    }
invisible(out)
}


plot.reclink.mix.par <- function(mix.par, D.mat){
###########################################
# Plots useful graphics illustrating the  #
# main features of the maximum likelihood #
# estimates resulting from the mixture    #
# model fit.                              # 
###########################################
if (!inherits(mix.par,"reclink.mix.par"))
    stop("Bad input: wrong class")
conformD(mix.par, D.mat)

# Check if a rare case occurred
no.trim <- identical( attr(mix.par,"optim.log"), "Empty trimmed distance matrix: only 0,1 found" )

if ( isTRUE(all.equal((Dmin <- min(D.mat)),(Dmax <- max(D.mat)))) ){
    distance <- seq(0, 1, length.out=1E3)
    }
else{
    distance <- seq(Dmin, Dmax, length.out=1E3)
    }

caption <- c("Distance Histogram\nand\nFitted Mixture Density", 
             "Distance Histogram and Fitted Mixture Density\n(very-low density region zoom)", 
             "Fitted logLikelihood Ratio at given distance", 
             "MAP linked pairs histogram")

# Plot1: Distance histogram and Fitted mixture model
plot1 <- expression(
{
h <- hist(D.mat, breaks = breaks, plot = FALSE) 
f <- mix.par[["f"]](distance)
plot(h, freq = FALSE, col = "lightgrey", xlim = c(0,1), ylim = range(0, h$density, f),
     xaxp=c(0,1,20), xlab = "distance", main = caption[1], border=FALSE)
lines(distance, f, lwd = 2, col = "red")
}
)

#Plot2: Very-low density region zoom of Plot1
plot2 <- expression(
{
h <- hist(D.mat, breaks = breaks, plot = FALSE) 
f <- mix.par[["f"]](distance)
pmfm <- mix.par$best.par["pm"] * mix.par[["f.m"]](distance)
plot(h, freq = FALSE, col = "lightgrey", xlim = c(0,1), ylim = range(0, 1.5*pmfm),
     xaxp=c(0,1,20), xlab = "distance", main = caption[2], border=FALSE)
lines(distance, f, lwd = 2, col = "red")
}
)

#Plot3: M fitted logLikelihood ratio at given distance
plot3 <- expression(
{
logLikelihood <- mix.par[["loglik"]](distance)
zero.dist <- distance[which.min(abs(logLikelihood))]
plot(distance, logLikelihood, ty="l", lwd = 2, col = "blue",
     xaxp=c(0,1,20), main = caption[3], ylab = "logLikelihood ratio")
abline(h = 0, lty = 2)
abline(v = zero.dist, lty = 2)
}
)

#Plot4: MAP linked pairs histogram
plot4 <- expression(
{
EC.M <- mix.par[["EC.M"]]
MAP.l <- mix.par[["MAP.l"]]
MAP.linked <- mapply( function(i,j) D.mat[i, j], MAP.l[,1], MAP.l[,2] )
hist(MAP.linked, breaks = breaks, col = "peachpuff", xlab = "distance", xaxp=c(0,1,20), 
     main = paste(caption[4],
                  paste("(MAP expected number of Matches given observed distances: ", EC.M, ")", sep=""),
                  sep = "\n"))
}
)

if (mix.par[["best.par"]]["pm"] > 0 ){
    if (!no.trim) {
        breaks = "scott"
        old.par <- par(mfrow = c(2,2)) 
        eval(plot1)
        breaks = "Sturges"
        eval(plot2)
        eval(plot3)
        eval(plot4)
        par(old.par)
    }
    else{
        breaks = "Sturges"
        old.par <- par(mfrow = c(3,1))
        eval(plot1)
        eval(plot2)
        eval(plot4)
        par(old.par)
    }
}	
else{
    breaks = "Sturges"
    eval(plot1)
    }
}



#############################################
# Print methods for class "reclink.mix.par" #
#############################################
print.reclink.mix.par <- function(x, ...){
cat("Fitted inflated-beta mixture model for Record Linkage\n")
cat("- MAP expected number of Matches given observed distances:", x$EC.M, "\n")
cat("- Unconditional expected number of Matches:", x$E.M, "\n")
cat("- Fitted parameters:\n\n")
print(x$best.par, digits =3)
cat("\n")
cat("- Matching constraints:", attr(x, "match.constr"), "\n")
cat("\n")
}

print.reclink.mix.par.sol.1to1 <- function(x, ...){
print.reclink.mix.par(x,...)
cat("* 1:1 RL solution found *\n")
}

print.reclink.mix.par.sol.1toN <- function(x, ...){
print.reclink.mix.par(x,...)
cat("* 1:N RL solution found *\n")
}

print.reclink.mix.par.sol.Nto1 <- function(x, ...){
print.reclink.mix.par(x,...)
cat("* N:1 RL solution found *\n")
}