EA.cluster.dedup.mix.par.sol <- function(mix.par, D.mat, ...){
############################################
# Method for class dedup.mix.par.sol: no   #
# further processing needed.               #
############################################
conformD(mix.par, D.mat)
stop("Deduplication solution already found: you don't need the Evolutionary Algorithm!")
}


mk.fitness.dedup.mix.par.1to1 <- function(mix.par, D.mat, ...){
#############################################################
# Builds the fitness function for the Clustering EA with ML #
# objective, under ONE-TO-ONE matching constraints.         #
#############################################################

# Compute (once and for all) the matrix of pairwise
# contributions to the complete logLikelihood
loglik.mat <- mix.par[["loglik"]](D.mat)

# Define actual fitness function for ML optimization
# with loglik.mat bound via lexical scoping
fitness <- function(perm){
##############################################################
# Given the matrix of pairwise contributions to the complete #
# logLikelyhood, computes the fitness of an arbitrary        #
# candidate solution to the 1:1 DD problem.                  #
# NOTE: The complete logLikelihood is calculated by dropping #
#       uninfluential constants (i.e. the contribution that  #
#       would arise by allocating all pairs to the U class). #
##############################################################
on <- perm!=0
if (sum(on)==0) return(0)
loglik1 <- function(i,j) loglik.mat[i,j]
sum( mapply(loglik1, which(on), perm[on]) )
}

return(fitness)
}


EA.cluster.dedup.mix.par <- function(mix.par, D.mat,
                                     nind = 300, ngen = 200, p.muta = 0.3, max.stallo = 50,
                                     greedy.frac = 0.05, gen.plot = 10, print = TRUE){
##############################################
# Method for class dedup.mix.par.            #
# Runs the Evolutionary Algorithm to cluster #
# Matches and Unmatches under specified      #
# matching constraints.                      #
##############################################
conformD(mix.par, D.mat)

# Make population size nind a multiple of 4
if (!identical(nind%%4,0) || nind < 4)
    nind <- as.integer(4 * (nind%/%4) + 4)

# Compute proper fitness function for ML optimization
# under given constraints
fitness <- mk.fitness(mix.par, D.mat)

monitor <- function(max.fit, avg.fit){
###############################################################
# Displays a graph to monitor the behaviour of the Clustering #
# EA while it is running.                                     #
###############################################################
plot(max.fit, ty="s", lwd=2, col="red", xlab = "Generations", ylab="Fitness maximum")
points(avg.fit, ty="l", lwd=1, col="blue")
}

# Initialize vectors to store generation statistics
max.fit <- rep(NA,ngen+1)
avg.fit <- rep(NA,ngen+1)
sd.fit  <- rep(NA,ngen+1)

# Generate initial population...
pop <- ini.pop(mix.par, D.mat, nind, greedy.frac = greedy.frac)
# ... and compute individual fitness
fit <- apply(pop, 1, fitness)

# Update statistics
max.fit[1] <- max(fit)
avg.fit[1] <- mean(fit)
sd.fit[1]  <- sd(fit)

# Initialize generations and stall counters
g <- 1
g.stallo <- 0

# Perform generation loop
while ((g < ngen+1) && (g.stallo < max.stallo)){
    # Increase generation counter
    g <- g+1
    # Selection - Reproduction - Mutation
    pop <- reprod(mix.par,pop,tourn(fit),p.muta)
    # Compute fitness of current individuals
    fit <- apply(pop, 1, fitness)
    # Update statistics
    max.fit[g] <- max(fit)
    avg.fit[g] <- mean(fit)
    sd.fit[g] <- sd(fit)
    # Check for stall and increase counter if needex
    if ( max.fit[g] <= max.fit[g-1] ) g.stallo <- (g.stallo+1) else g.stallo <- 0
    # Run monitor
	if (g %in% seq(2,ngen+1,by=gen.plot)) monitor(max.fit[1:g], avg.fit[1:g])
    }

# Shrink statistics vector
max.fit <- max.fit[!is.na(max.fit)]
avg.fit <- avg.fit[!is.na(avg.fit)]
sd.fit  <- sd.fit[!is.na(sd.fit)]
# Last monitor
if (gen.plot < Inf) monitor(max.fit, avg.fit)
# Build output list
out <- list(last=list(pop=pop, fit=fit),
            evo=list(max.fit=max.fit, avg.fit=avg.fit, sd.fit=sd.fit))

if (!is.null(attr(D.mat, "file.a"))) {
     attr(out, "file.a") <- attr(D.mat, "file.a")
     attr(out, "file.b") <- attr(D.mat, "file.b")
    }
attr(out,"D.mat.dim")    <- attr(mix.par,"D.mat.dim")
attr(out,"match.constr") <- attr(mix.par,"match.constr")
attr(out,"Call") <- sys.call(-1)
class(out) <- c("dedup.ea.out", "ea.out")
if (isTRUE(print)) {
     print(out)
    }
invisible(out)
}



ini.pop.dedup.mix.par.1to1 <- function(mix.par, D.mat, nind, greedy.frac = 0.05){
##################################################################
# Given the fitted mixture parameters and the matrix of pairwise #
# distances, generates the initial population for the Clustering #
# EA under ONE-TO-ONE matching constraints.                      #
# NOTE: A fraction greedy.frac (5% by default) of the nind       #
#       individuals is built by repairing (in diffrent ways)     #
#       'THE' greedy solution of the DD problem.                 #
#       Remaining  individuals get generated by repairing random #
#       permutations selected via a MC sampling technique which  #
#       relies on M/U posterior estimates.                       #
# VALUE: A matrix with dimension c(nind, nrow(D.mat)).           #
##################################################################
m <- nrow(D.mat)
n <- ncol(D.mat)
pop <- matrix(NA,nind,m)
n.greedy <- floor(nind*greedy.frac)
if (n.greedy > 0){
    greedy.perm <- greedy(mix.par, D.mat)
    for (i in 1:n.greedy) pop[i,] <- repair(greedy.perm)
    }
if (n.greedy < nind){
    prob.post <- prob.post(mix.par, D.mat)
    for (i in (n.greedy+1):nind) pop[i,] <- repair(sample.post(prob.post))
    }
pop
}



greedy.dedup.mix.par.1to1 <- function(mix.par, D.mat){
###############################################
# Given the observed distance matrix and the  #
# fitted parameters for the mixture, computes #
# the "most a-posteriori probable" solution   #
# for the DD problem, under ONE-TO-ONE        #
# matching constraints.                       #
# NOTE: The output permutation may (and in    #
#       general will) contain duplicated,     #
#       thus it has to be repaired before     #
#       being used.                           #
# VALUE: a permutation of length nrow(D.mat). #
###############################################
m <- nrow(D.mat)
perm <- rep(NA, m)
# first genotype allele must be always 0, due to constraints + simmetry
perm[1] <- 0
for (i in 2:m){
     # keep only meaningful columns (col < row) 
     j.legal <- 1:(i-1)
     lP0 <- sum( log( mix.par[["fu.d"]](D.mat[i, j.legal]) ) )
     lP  <- lP0 + mix.par[["loglik"]](D.mat[i, j.legal])
     # Check: if exists j such that fm.d(D.mat[i,j]) is "exactly" 1
     # this would imply lP0 is -Inf and loglik(D.mat[i,j]) is Inf;
     # thus lP[j] would be NaN (due to -Inf + Inf) and it would be
     # ignored by which.max(), uncorrectly but without errors.
     # For such (few) cases go back to the real (computationally
     # expansive) math expression (see first *if* below).
     nan.j <- which(is.nan(lP))
     if ( length(nan.j) >= 1 ) {
          lP[nan.j] <- sapply(nan.j, function(j) sum( log( mix.par[["fu.d"]](D.mat[i, 
                                                                                   j.legal[j.legal != j] ]) ) ) )
        }
     perm[i] <- as.integer(which.max(c(lP0,lP)) - 1) 
    }
perm
}


prob.post.dedup.mix.par.1to1 <- function(mix.par, D.mat){
###############################################
# Given the observed distance matrix and the  #
# fitted parameters for the mixture, computes #
# the "a-posteriori probabilities" for record #
# i of file.a being either an Unmatch or a    #
# Match with record j of file.b.              #
# NOTE: ONE-TO-ONE matching constraints are   #
#       assumed.                              #
# NOTE: The output is an (m,n+1) matrix with  #
#       first column giving the probability   #
#       of being an Unmatch.                  #
###############################################
m <- nrow(D.mat)
n <- ncol(D.mat)
prob <- matrix(NA, m, n+1)
# put to 0 probs tied to redundant pairs (diagonal + upper triangle - last column)
prob[(row(prob) <= col(prob)) & (col(prob) > 1)] <- 0
# first row is special (no M allowed)
prob[1, 1] <- 1
for (i in 2:m){
     # keep only meaningful columns (col < row) 
     j.legal <- 1:(i-1)
     lP0 <- sum( log( mix.par[["fu.d"]](D.mat[i, j.legal]) ) )
     lP  <- lP0 + mix.par[["loglik"]](D.mat[i, j.legal])
     prob[i, 1:i] <- c(exp(lP0), exp(lP))
     # Check: if exists j such that fm.d(D.mat[i,j]) is "exactly" 1
     # this would imply lP0 is -Inf and loglik(D.mat[i,j]) is Inf;
     # thus lP[j] would be NaN (due to -Inf + Inf), like prob[i, 1 + j],
     # causing an error when sampling (via *sample()*).
     # For such (few) cases go back to the real (computationally
     # expansive) math expression (see first *if* below).
	 nan.j <- which(is.nan(lP))
	 if ( length(nan.j) >= 1 ) {
          prob[i, 1 + nan.j] <- exp( sapply(nan.j, function(j) sum( log( mix.par[["fu.d"]](D.mat[i,
                                                                                j.legal[j.legal != j] ]) ) ) ) )
        }
     # Check: if exists a set of MORE THAN ONE j such that fm.d(D.mat[i,j]) is "exactly" 1
     # this means that the probability that i matches ONLY a given j
     # is ZERO, as well as the probability that i does not match any file.b record.
     # Such circumstance would imply that prob[i, ] contains ONLY ZEROS,
     # from which it would be impossible to sample (via *sample()*).
     # Solution: discard the "i matches ONLY j" constraint, see below.
     if ( sum(prob[i, ]) <= 0 ) {
          prob[i, 1 + j.legal] <- mix.par[["fm.d"]](D.mat[i, j.legal])
        }
    }
prob/rowSums(prob)
}


cross.dedup.mix.par.1to1 <- function(mix.par,perm.m,perm.f){
###############################################
# Given a pair of parent genotypes, generates #
# trough one-point crossover two new child    #
# genotypes.                                  #
# NOTE: Under ONE-TO-ONE matching constraints #
#       the childs have to be repaired to be  #
#       sure they will eventually be legal.   #
###############################################
n <- length(perm.m)
cut <- sample(1:(n-1),1)
child.1 <- repair(c(perm.m[1:cut],perm.f[(cut+1):n]))
child.2 <- repair(c(perm.f[1:cut],perm.m[(cut+1):n]))
list(child.1=child.1,child.2=child.2)
}


mutation.dedup.mix.par.1to1 <- function(mix.par,perm,p.muta){
#################################
# Mutate an individual.         #
# NOTE: For ONE-TO-ONE matching #
#       when adding a new match #
#       avoids duplicates.      #
#################################
if (runif(1) > p.muta) return(perm)
m <- length(perm)
# Get upper limit of alleles domain 
n <- attr(mix.par, "D.mat.dim")[2]
off <- perm==0
on <- perm!=0
n.on <- sum(on)
n.off <- m-n.on
# Instead of using the unconditional expected number of matches
# use the MAP estimate of the number of matches (this should
# determine some mutation which deletes a match)
plink <- min(1,mix.par[["EC.M"]]/m)
rl <- rbinom(1,m,plink)
if (n.on <= rl) {
    if (n.off==0) return(perm)
    to.on <- sample(which(off),1,rep=FALSE)
    if (to.on > 1) {
        match.on <- (1:(to.on-1))[!(1:(to.on-1) %in% unique(perm))]
        if (length(match.on) > 0) {
            to.match <- sample(match.on,1,rep=FALSE)
            }
        else to.match <- 0
        }
    else to.match <- 0
    perm[to.on] <- to.match
    }
else {
      if (n.on==0) return(perm)
      to.off <- sample(which(on),1,rep=FALSE)
      perm[to.off] <- 0
     }
perm
}



##########################################
# Print methods for class "dedup.ea.out" #
##########################################
print.dedup.ea.out <- function(x, ...){
cat("Clustering solution for Deduplication\n")
cat("- Declared Matches:", sum( getSolution(x) !=0 ), "\n")
cat("- Matching constraints:", attr(x, "match.constr"), "\n")
cat("\n")
}
