split.blocks <- function(data, block.vars){
#####################################################
# Given a set of blocking variables and a dataframe #
# returns a list containing the ids (row numbers)   #
# of the records belonging to the blocks.           #
#####################################################
interact <- interaction(data[, block.vars, drop = FALSE], drop = TRUE, lex.order = TRUE)
groups   <- .Internal(split(1:nrow(data), interact))
groups
}

blocking <- function(file.a, file.b, block.vars){
########################################################
# Given two dataframes and a set of blocking variables #
# returns the list containing the ids (row numbers) of #
# the records belonging to the smaller (S) and to the  #
# bigger (B) blocks coming from the dataframes.        #
# NOTE: The output is the dedicated data structure to  #
#       be fed into B.maerlin in order to perform a    #
#       block-wise RL job.                             # 
########################################################
g.a <- split.blocks(file.a, block.vars)
l.a <- names(g.a)
g.b <- split.blocks(file.b, block.vars)
l.b <- names(g.b)
# Only common blocks must be processed
l.ab <- intersect(l.a, l.b)
g.a <- g.a[l.a %in% l.ab]
g.b <- g.b[l.b %in% l.ab]
# Assess the dimension of a and b blocks 
n.a <- sapply(g.a, length)
n.b <- sapply(g.b, length)
a.small <- (n.a <= n.b)
# Artificial small (S) and big (B) groups (may alternate a and b blocks)
g.S <- g.a
g.S[!a.small] <- g.b[!a.small]
g.B <- g.b
g.B[!a.small] <- g.a[!a.small]
blocks <- list(S=g.S, B=g.B)
# Useful info to be stored
attr(blocks,"block.vars") <- block.vars
attr(blocks,"block.size") <- n.a*n.b
attr(blocks,"a.small") <- a.small
blocks 
}

B.maerlin <- function(file.a, file.b, block.vars,
                     match.vars = "All", dist.funs = "Levenshtein", weights = NULL,
                     match.constr = c("1:1", "1:N", "N:1", "N:M"), guessM = c(c.am = 1/2, c.bm = 2), force = TRUE,
                     verbose = FALSE, fit.plot = FALSE, print = TRUE,
                     nind = 300, ngen = 200, p.muta = 0.3, max.stallo = 50, greedy.frac = 0.05, gen.plot = Inf,
                     ...){
#########################################################
# Given two dataframes and a set of blocking variables  #
# performs a block-wise RL job.                         #
# NOTE: ... stands for further arguments to the maerlin #
#       inner functions.                                #
# NOTE: the output has the standard structure of a      #
#       global maerlin output (i.e. global info gets    #
#       reconstructed from blocks info).                #  
#########################################################
if (nrow(file.a) > nrow(file.b))
    stop("File.a (",substitute(file.a),") has more rows then file.b (",substitute(file.b),"): please exchange them")
common <- intersect(names(file.a), names(file.b))
if ( !any(block.vars %in% common) )
    stop("Blocking variables must be common to both files")

# N:M block-wise RL yet to be implemented...
# if (ntom) stop("Only 1:1 block-wise RL currently available!")

# START message
timestamp(stamp = "Block-wise Record Linkage - START", prefix = "\n###### ", suffix = "", quiet = FALSE)
timestamp(prefix = "###### ", suffix = "", quiet = FALSE)

# Timing: start
ptm <- proc.time()

# Build the blocks
cat("\n:::::: Search Space Reduction: preparing blocks\n")
blocks <- blocking(file.a, file.b, block.vars)
# Profile blocks
cat("\n:::::: Analyzing block size distribution:\n")
size <- attr(blocks,"block.size")
cat("\n> Number of blocks: ", length(size), "\n")
cat("\n> Size deciles (pairs):\n")
print(round(quantile(size, probs = seq(0, 1, by = 0.1))))
cat("\n> Size deciles (records):\n")
print(round(quantile(sqrt(size), probs = seq(0, 1, by = 0.1))))
cat("\n> Reduction Ratio: ", sum(size)/prod(nrow(file.a), nrow(file.b)), "\n")

vsep <- paste(rep("-", 0.12 * getOption("width")), collapse = "")
# Perform block-wise RL
S <- blocks[["S"]]
B <- blocks[["B"]]
a.small <-attr(blocks, "a.small")
block.vars <- attr(blocks, "block.vars")
block.tag <- paste(block.vars, collapse=".")
res <- vector(mode = "list", length = length(S))
names(res) <- paste(block.tag, " = ", names(S), sep="")
for (i.bl in seq_along(res)){
# DEBUG A BLOCK
#  if (i.bl == 51) browser()
# DEBUG A BLOCK
    s <- size[i.bl]
	n <- round(sqrt(s))
	tag <- names(res[i.bl])
    cat("\n", vsep, " Processing block ", i.bl, " (size = ", s, ", n ~ ", n, ") ", " [", tag, "] ", vsep, "\n\n", sep="")
#    cat(paste("\n**** Processing block ", i.bl, " (size = ", s, ", n ~ ", n, ") ", " [", tag, "]\n", sep=""))
    if (a.small[i.bl]) {
        r <- try( maerlin(file.a[S[[i.bl]], , drop = FALSE], file.b[B[[i.bl]], , drop = FALSE], 
                         match.vars, dist.funs, weights,
                         match.constr, guessM, force, verbose, fit.plot, print,
                         nind, ngen, p.muta, max.stallo, greedy.frac, gen.plot,
                         a.small = a.small[i.bl], ...)
                )
        }
    else{
        r <- try( maerlin(file.b[S[[i.bl]], , drop = FALSE], file.a[B[[i.bl]], , drop = FALSE], 
                         match.vars, dist.funs, weights,
                         match.constr, guessM, force, verbose, fit.plot, print,
                         nind, ngen, p.muta, max.stallo, greedy.frac, gen.plot,
                         a.small = a.small[i.bl], ...)
                )
        }

	if (!inherits(r, "try-error") ) {
        links <- attr(r, "Links")
	    r <- as.integer(getSolution(r))
	    cat("> Decoding block solution\n")
        # Find keys of positive pairs inside the block
	    np <- length( B.id <- which(r != 0) )
	    cat("\n> Matches in the block:", np, "\n")
        res[[i.bl]] <- matrix(nrow = np, ncol = 2)
	    if (a.small[i.bl]) {
            res[[i.bl]][, 1] <- S[[i.bl]][B.id]
            res[[i.bl]][, 2] <- B[[i.bl]][r[B.id]]
            }
        else{
            res[[i.bl]][, 2] <- S[[i.bl]][B.id]
            res[[i.bl]][, 1] <- B[[i.bl]][r[B.id]]
            }
        attr(res[[i.bl]], "Links") <- links 
    }
    else {
        res[[i.bl]] <- r
        cat("**** Error: block RL job failed\n")
        }
}
# Compute global RL solution
cat("\n:::::: Reconstructing global solution from blocks\n")
solution <- rep(0, nrow(file.a))
lapply(res, function(r) {
                        if (!inherits(r, "try-error")) {
                            mapply( function(i, j) solution[i] <<- j, r[, 1], r[, 2] )
                            }
                        }
      )

cat("\n:::::: Reconstructing Matches data frame\n")
Links <- lapply(res, function(el) attr(el, "Links"))
Links <- Reduce(rbind, Links[sapply(Links, is.data.frame)])
# row.names(Links) <- NULL

# Finalize output object
# attr(solution, "file.a") <- substitute(file.a)
# attr(solution, "file.b") <- substitute(file.b)
attr(solution, "block.res") <- res

# Compute for each block how many matches have been declared
  # NOTE: -Inf means that the RL job failed for that block
block.nlinks <- lapply(res, nrow)
block.nlinks[sapply(block.nlinks, is.null)] <- -Inf
block.nlinks <- unlist(block.nlinks)
attr(solution, "block.nmatch") <- block.nlinks
attr(solution, "block.size") <- size
attr(solution, "ERROR.BLOCKS.size") <- size[block.nlinks < 0]
attr(solution, "Matches") <- Links
# class(solution) <- "block.sol"

# Timing: stop
et <- proc.time() - ptm
cat("\n:::::: Elapsed time:", et[3],"\n\n")

# END message
timestamp(prefix = "###### ", suffix = "", quiet = FALSE)
timestamp(stamp = "Block-wise Record Linkage - END", prefix = "###### ", suffix = "\n\n", quiet = FALSE)

solution
}


maerlin <- function(file.a, file.b, match.vars = "All", dist.funs = "Levenshtein" , weights = NULL,
                   match.constr = c("1:1", "1:N", "N:1", "N:M"), guessM = c(c.am = 1/2, c.bm = 2),
				   force = TRUE, verbose = FALSE, plot = TRUE, print = TRUE,
                   nind = 300, ngen = 200, p.muta = 0.3, max.stallo = 50, greedy.frac = 0.05, gen.plot = 10,
                   a.small = TRUE, ...){
#############################################################################
# NOTE: Argument a.small is meaningful only when maerlin gets called in a   #
#       block-wise RL job. In these cases the semantics is that !a.small    #
#       flags blocks where file.a, though being smaller than file.b,        #
#       actually comes from the the original file B.                        #
#############################################################################
#magus <<- list(NA, NA, NA)
#names(magus) <<- c("D", "m", "g")
cat("> Computing the Distance Matrix \n")
D.ab <- Distance(file.a, file.b, match.vars, dist.funs, weights)
#magus[["D"]] <<- D.ab
cat("> Fitting the Mixture \n")
# Timing: start
#ptm <- proc.time()

m.ab <- fit.mix.par(D.ab, match.constr, guessM, force, verbose, plot, print, ...)

#magus[["m"]] <<- m.ab
if (inherits(m.ab, "reclink.mix.par.sol")) {
    cat("**** No Duplicates in MAP solution \n\n")
    out <- m.ab
    }
else {
    cat("**** Duplicates in MAP solution: One-to-One Reduction \n\n")
    cat("> Running the Clustering EA \n")
    out <- EA.cluster(m.ab, D.ab, nind, ngen, p.muta, max.stallo, greedy.frac, gen.plot, print)
									 
#    magus[["g"]] <<- out
    }
#magus <<- magus[!is.na(magus)]
# Timing: stop
#et <- proc.time() - ptm
#cat("Time:\n")
#print(et)
#attr(res,"file.a") <- substitute(file.a)
#attr(res,"file.b") <- substitute(file.b)
# class(out) <- "maerlin"

links <- try( getLinks(out, file.a, file.b, D.ab, m.ab, a.small = a.small) )
if (!inherits(links, "try-error") ) {
     attr(out, "Links") <- links}
else {
     attr(out, "Links") <- "COULD NOT FIND LINKS"}

out
}


KillLinks <- function(links, Q.th = 0.25, ...){
  nl <- nrow(links)/2
  if (nl <= 0) stop("No Matches to be processed!")
  l.a <- links[links$FILE == "A", , drop = FALSE]
  l.b <- links[links$FILE == "B", , drop = FALSE]

  # User Defined Kill Rules
  # NOTE: Must refer to existing variables
  is.toKill <- function(r.a, r.b){
      toKill <- FALSE
      # Build a string to store the rules that caused the decision to kill the link
      # NOTE: "0" is for links that are not killed
      Rule <- 0
      # --- Application Task: Relaunches on Scanner Data
      ### User Defined Kill Rules - START
      ### Rule 1: Relaunches cannot have different UNITA
      if (!is.na(r.a$UNITA) && !is.na(r.b$UNITA) && r.a$UNITA != r.b$UNITA) {
          toKill <- TRUE
          Rule <- paste(Rule, 1, sep = "|")
        }
      ### Rule 2: Relaunches cannot have a relative difference in QUANTITA
      ###         bigger than Q.th
      ### NOTE: As requested by the CPI team, relative differences are computed w.r.t. the phasing out product,  # DEPRECATED
      ###       i.e. Q.reldiff = abs(Q.in - Q.out)/Q.out                                                         # DEPRECATED
      # Qreldiff <- abs( as.numeric(r.b$QUANTITA) - as.numeric(r.a$QUANTITA) ) / as.numeric(r.a$QUANTITA)
      ### NOTE: As requested by the CPI team, relative differences are computed in a time agnostic way,          # 18/11/2019
      ###       i.e. Q.max = max(Q.in, Q.out)
      ###       i.e. Q.min = min(Q.in, Q.out)
      ###       i.e. Q.reldiff = abs(Q.max - Q.min)/Q.max
      Q.max <- pmax(as.numeric(r.b$QUANTITA), as.numeric(r.a$QUANTITA))
      Q.min <- pmin(as.numeric(r.b$QUANTITA), as.numeric(r.a$QUANTITA))
      Qreldiff <- abs( Q.max - Q.min ) / Q.max
      if ( !is.na(r.a$QUANTITA) && !is.na(r.b$QUANTITA) && (Qreldiff > Q.th) ) {
          toKill <- TRUE
          Rule <- paste(Rule, 2, sep = "|")
        }
      ### Rule 3: Relaunches cannot have missing QUANTITA
      if ( is.na(r.a$QUANTITA) || is.na(r.b$QUANTITA) ) {
          toKill <- TRUE
          Rule <- paste(Rule, 3, sep = "|")
        }
      ### User Defined Kill Rules - END
      # --- Application Task: Relaunches on Scanner Data
      if (toKill) {
         # If the link has been killed, remove the ugly initial "0|" substring from the 'Rule' string
         Rule <- substr(Rule, 3, nchar(Rule))
        }
      attr(toKill, "Rule") <- Rule
      return(toKill)
    }

  links.toKill <- sapply(1:nl, function(i) is.toKill(l.a[i, , drop = FALSE], l.b[i, , drop = FALSE]))
  Rules.toKill <- sapply(1:nl, function(i) attr( is.toKill(l.a[i, , drop = FALSE], l.b[i, , drop = FALSE]), "Rule"))

  links.toKill <- rep(links.toKill, each = 2)
  Rules.toKill <- rep(Rules.toKill, each = 2)

  links.alive <- links[!links.toKill, , drop = FALSE]
  rownames(links.alive) <- NULL

  links.dead  <- links[links.toKill, , drop = FALSE]
  links.dead$REGOLA <- Rules.toKill[links.toKill]
  rownames(links.dead) <- NULL

  attr(links.alive, "killed") <- links.dead
  return(links.alive)
}


getPRDKEY <- function(links){
  nl <- NROW(links)/2
  if (nl <= 0) stop("No Matches to be processed!")
  PRDK.a <- links[links$FILE == "A", "PRDKEY"]
  PRDK.b <- links[links$FILE == "B", "PRDKEY"]
  PRDK.ab <- data.frame(PRDKEY.A = PRDK.a, PRDKEY.B = PRDK.b)
  return(PRDK.ab)
}




PPP.dist <- function(file.a, file.b, block.vars,
                     match.vars = "All", dist.funs = "Levenshtein", weights = NULL){
#########################################################
#                                                       #
#########################################################
# if (nrow(file.a) > nrow(file.b))
#     stop("File.a (",substitute(file.a),") has more rows then file.b (",substitute(file.b),"): please exchange them")
common <- intersect(names(file.a), names(file.b))
if ( !any(block.vars %in% common) )
    stop("Blocking variables must be common to both files")

# Build the blocks
cat("\n:::::: Search Space Reduction: preparing blocks\n")
blocks <- blocking(file.a, file.b, block.vars)
# Profile blocks
cat(":::::: Analyzing block size distribution:\n")
size <- attr(blocks,"block.size")
cat("> Number of blocks: ", length(size), "\n")
cat("> Size deciles (pairs):\n")
print(round(quantile(size, probs = seq(0, 1, by = 0.1))))
cat("> Size deciles (records):\n")
print(round(quantile(sqrt(size), probs = seq(0, 1, by = 0.1))))
cat("> Reduction Ratio: ", sum(size)/prod(nrow(file.a), nrow(file.b)), "\n")

# Compute block-wise distance matrices
S <- blocks[["S"]]
B <- blocks[["B"]]
a.small <-attr(blocks, "a.small")
block.vars <- attr(blocks, "block.vars")
block.tag <- paste(block.vars, collapse=".")
res <- vector(mode = "list", length = length(S))
names(res) <- paste(block.tag, " = ", names(S), sep="")
for (i.bl in seq_along(res)){
    s <- size[i.bl]
	n <- round(sqrt(s))
	tag <- names(res[i.bl])
    cat(paste("\n**** Processing block ", i.bl, " (size = ", s, ", n ~ ", n, ") ", " [", tag, "]\n", sep=""))
    if (a.small[i.bl]) {
        file.a.i <- file.a[S[[i.bl]], , drop = FALSE]
        file.b.i <- file.b[B[[i.bl]], , drop = FALSE]
        Di <- try( Distance(file.a.i, file.b.i, match.vars, dist.funs, weights) )
        }
    else{
        file.a.i <- file.a[B[[i.bl]], , drop = FALSE]
        file.b.i <- file.b[S[[i.bl]], , drop = FALSE]
        Di <- try( Distance(file.b.i, file.a.i, match.vars, dist.funs, weights) )
        if ( !inherits(Di, "try-error") ) {
             Di <- t(Di)
            }
        }

	if (!inherits(Di, "try-error") ) {
	    # Di <- as.numeric(Di)        
        Oi <- order(Di)
        ODi <- Di[Oi]
        rowi <- row(Di)[Oi]
        coli <- col(Di)[Oi]
        res[[i.bl]] <- mapply(FUN = function(rr, cc, dd) {
                                             dfi <- rbind(file.a.i[rr, common, drop = FALSE], file.b.i[cc, common, drop = FALSE])
                                             dfi[["dist"]] <- dd
                                             dfi
                                            }, rowi, coli, ODi, SIMPLIFY = FALSE )

        res[[i.bl]] <- Reduce(f = rbind, x = res[[i.bl]]) 
	    # cat("**** Decoding block solution\n")
        # Find keys of positive pairs inside the block
	    # np <- length( B.id <- which(r != 0) )
        # res[[i.bl]] <- matrix(nrow = np, ncol = 2)
	    # if (a.small[i.bl]) {
            # res[[i.bl]][, 1] <- S[[i.bl]][B.id]
            # res[[i.bl]][, 2] <- B[[i.bl]][r[B.id]] 
            # }
        # else{
            # res[[i.bl]][, 2] <- S[[i.bl]][B.id]
            # res[[i.bl]][, 1] <- B[[i.bl]][r[B.id]] 
            # }
    }
    else {
        res[[i.bl]] <- Di
        cat("**** Error: distance computation failed for this block\n")
        }
    # attr(res[[i.bl]], "file.a") <- file.a[S[[i.bl]], , drop = FALSE]
    # attr(res[[i.bl]], "file.b") <- file.b[B[[i.bl]], , drop = FALSE]
}
# Compute global RL solution
# cat("\n:::::: Reconstructing global solution from blocks\n")
# solution <- rep(0, nrow(file.a))
# lapply(res, function(r) {
                        # if (!inherits(r, "try-error")) {
                            # mapply( function(i, j) solution[i] <<- j, r[, 1], r[, 2] )
                            # }
                        # }
      # )
# attr(solution, "block.res") <- res
# attr(solution, "block.size") <- size
# attr(solution, "file.a") <- substitute(file.a)
# attr(solution, "file.b") <- substitute(file.b)
# class(solution) <- "block.sol"
# solution
res
}