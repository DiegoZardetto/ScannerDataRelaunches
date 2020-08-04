emailCorp <- function(path, encoding = "UTF-8"){
############################################################
#                                                          #
############################################################
require(tm)
require(tm.plugin.mail)
require(proxy)
# build a vector source for tm pkg
mail.source <- DirSource(path, encoding = encoding)
# build a volatile Corpus on this source
mail.corp <- Corpus(mail.source, readerControl = list(reader = readMail))
mail.corp
}

emailBody_old <- function(email.corp, weighting = weightTfIdf){
############################################################
#                                                          #
############################################################
# compute the term-document matrix
mail.tdm <- TermDocumentMatrix(email.corp, control=list(weighting = weighting))
# extract the matrix only
mail.mat <- as.matrix(mail.tdm)
# count how many mails are being processed
nmails <- ncol(mail.mat)

# compute cosine similarity
  # first put matrix columns into a list
  l.mail  <- lapply(1:nmails, function(col) as.numeric(mail.mat[, col]))
  # build distance matrix
  d <- dist(l.mail, l.mail, "cosine")
  
# d.vect <- as.numeric(d)
# Deduplication: restrict to lower triangle of matrix d
d[row(d) <= col(d)] <- NA
d.vect <- d[row(d) > col(d)]

# Compute the standard deviation of distances for the current variable.
# If it happens to be zero (i.e. all records from both files are identical
# w.r.t. the variable, hence all distances are zero), then put sd <- 1.
s <- sd(d.vect)
# output: dist matrix
out <- d
attr(out, "stdev") <- ifelse(s > 0, s, 1)
class(out) <- "dedup.Dist"
out
}

getHFields <- function(email.corp){
nmail <- length(email.corp)
out <- data.frame(Subject = rep(NA, nmail), From = rep(NA, nmail), To = rep(NA, nmail))
i.mail <- 0
for (email in email.corp){
     i.mail <- i.mail + 1
     # get Header slot of email
     head <- attr(email, "Header")
	 # MUST WRITE DOWN SOME CHECKS ON FIELDS APPEARANCE!!!!!!!!!!!!!!!!!!!!!!!!
     # Isolate "Subject" field:
     where <- pmatch("Subject:", head)
     subj <- unlist(strsplit(head[where], "Subject:"))[2]
     # Isolate "From" field:
     where <- pmatch("From:", head)
     from <- unlist(strsplit(head[where], "From:"))[2]
     # Isolate "To" field:
     where <- pmatch("To:", head)
     to <- unlist(strsplit(head[where], "To:"))[2]
     # Insert current values inside out
     out[i.mail, ] <- c(subj, from, to)
	}
out
}


Emaildist <- function(email.corp, dist.funs = "Levenshtein", body.weighting = weightTfIdf, weights = NULL){
############################################################
############################################################
dist.funs.list <- c("Levenshtein", "JaroWinkler", "Equality", "Three.grams", "Cos.3.TfIdf")
nmail <- length(email.corp)

# Build header dataframe
Head.df <- getHFields(email.corp)
Hvars <- names(Head.df)
vars <- c(Hvars, "Body")

if (!is.null(weights)){
    if (!is.numeric(weights))
        stop("If given, weights must be numeric")
    # Avoid negative weights
    weights <- abs(weights)
    w <- rep(weights, length.out = length(vars))
    }
if (is.null(weights)){
    w <- rep(1, length.out = length(vars))
    }
names(w) <- vars

dist.funs <- rep(match.arg(dist.funs, dist.funs.list, several.ok = TRUE), length.out = length(Hvars))
names(dist.funs) <- Hvars

# Check if the function has been called from the command line
# or by another function (if so I cannot retrieve file names
# by using substitute())
directly <- !(length(sys.calls()) > 1)

dedup <- TRUE

# Compute a list with the distance matrices as components 
D.list  <- lapply(Hvars, function(var)
                          do.call(what = dist.funs[var],
                                  args = list(x = Head.df[, var], y = Head.df[, var], dedup = dedup)))

# Retrieve standard deviations vector
stdev <- sapply(D.list, function(el) attr(el, "stdev"))

# Append Body distances and stdev
Body.dist <- emailBody_old(email.corp, weighting = body.weighting)
Body.stdev <- attr(Body.dist, "stdev")
D.list <- list(D.list, Body.dist)
stdev <- c(stdev, Body.stdev)

# Convert distance matrices list into an array
D.array <- array(unlist(D.list), dim=c(nmail, nmail, length(vars)))

# If more than 1 matching variable selected, compute distance as the weighted
# average of individual distances (current version uses array algebra, this is
# ~ 2.7 times faster than apply + weighted.mean)
ww <- w/stdev
nww <- ww/sum(ww)
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
     attr(D.mat,"email.corp") <- substitute(email.corp)
    }
attr(D.mat,"match.vars") <- vars
attr(D.mat,"dist.funs") <- c(dist.funs, "Body.dist")
names(attr(D.mat,"dist.funs")) <- c(attr(D.mat,"dist.funs"), "Body")
attr(D.mat,"in.weights") <- w
attr(D.mat,"out.weights") <- ww
attr(D.mat,"cutoff") <- if (has.cutoff) cutoff
attr(D.mat,"Call") <- sys.call()
class(D.mat) <- c(if (!dedup) "reclink.Dist" else "dedup.Dist", class(D.mat))
D.mat
}