Levenshtein <- function(x, y, dedup, ...) {
#############################################################
# Normalized LEVENSHTEIN distance between character vectors #
# lower-cased with spaces removed.                          #
# Returns a Distance Matrix (either full dedup = FALSE or   #
# lower triangular if dedup = TRUE).                        #
# NOTE: Levenshtein distance d_i for each pair (x_i, y_i)   #
#       is normalized in [0,1].                             #
#############################################################
require(proxy)
# Strip spaces
x <- tolower(gsub(" ", "", x))
y <- tolower(gsub(" ", "", y))
# Compute max length for each string pair
# (first avoid NA coercion to "NA" due to function nchar)
max.len <- outer(x, y, function(x,y) {
                                      x[is.na(x)]<-""
                                      y[is.na(y)]<-""
                                      pmax(nchar(x),nchar(y))})
# Compute normalized distances:
# note that dist's weight parameter differs from the default (c(1,1,0,2))
  #29/02/2012 WARNNG: cba package sdists fun changed: 4 weights instead of previous 3
d <- dist(x, y, method="Levenshtein", weight=c(1, 1, 0, 1)) / max.len

# If x or/and y for a given pair is/are NA, put the distance
# to the (blind) average value 1/2. As this may concentrate too much
# the distance distribution, the caller function 'Distance' must be able
# to soften the importance of the variable according to its NA rate:
NAs <- is.na(d)
NAf <- sum(NAs)/prod(dim(NAs))

d[NAs] <- 0.5
d.vect <- as.numeric(d)
# If Deduplication is in order, then restrict to lower triangle of matrix d
if (dedup) {
    d[row(d) <= col(d)] <- NA
    d.vect <- d[row(d) > col(d)]
    NAs <- NAs[row(d) > col(d)]
    NAf <- sum(NAs)/length(NAs)
    }

# Compute the standard deviation of distances for the current variable.
# If it happens to be zero (i.e. all records from both files are identical
# w.r.t. the variable, hence all distances are zero), then the caller
# function 'Distance' will put the corresponding weight to zero.
s <- sd(d.vect)
# output: dist matrix
out <- list(d = d)
attr(out, "stdev") <- s
attr(out, "NAf") <- NAf
out
}


#-------------------------------------------------------------------------------

Equality <- function(x, y, dedup, ...) {
#############################################################
# Normalized EQUALITY distance between character vectors    #
# lower-cased with spaces removed.                          #
# Returns a Distance Matrix (either full dedup = FALSE or   #
# lower triangular if dedup = TRUE).                        #
# NOTE: Equality distance d_i for each pair (x_i, y_i) is   #
#       normalized in [0,1].                                #
#############################################################
# Strip spaces
x <- tolower(gsub(" ", "", x))
y <- tolower(gsub(" ", "", y))
# Compute normalized distances:
d <- 1 - outer(x, y, FUN = "==")

# If x or/and y for a given pair is/are NA, put the distance
# to the (blind) average value 1/2. As this may concentrate too much
# the distance distribution, the caller function 'Distance' must be able
# to soften the importance of the variable according to its NA rate:
NAs <- is.na(d)
NAf <- sum(NAs)/prod(dim(NAs))

d[NAs] <- 0.5
d.vect <- as.numeric(d)
# If Deduplication is in order, then restrict to lower triangle of matrix d
if (dedup) {
    d[row(d) <= col(d)] <- NA
    d.vect <- d[row(d) > col(d)]
    NAs <- NAs[row(d) > col(d)]
    NAf <- sum(NAs)/length(NAs)
    }

# Compute the standard deviation of distances for the current variable.
# If it happens to be zero (i.e. all records from both files are identical
# w.r.t. the variable, hence all distances are zero), then the caller
# function 'Distance' will put the corresponding weight to zero.
s <- sd(d.vect)
# output: dist matrix
out <- list(d = d)
attr(out, "stdev") <- s
attr(out, "NAf") <- NAf
out
}


#-------------------------------------------------------------------------------

Diffabs <- function(x, y, dedup, ...) {
#############################################################
# Normalized distance between numeric vectors (or character #
# vectors that can be interpreted as numeric).              #
#                                                           #
# The function, d_i = diffabs(x_i, y_i), is follows:        #
# - if x_i = y_i = 0       ->  d_i = 0                      #
# - if x_i != 0 & y_i = 0  ->  d_i = 1                      #
# - if x_i = 0 & y_i != 0  ->  d_i = 1                      #
# - if sign x_i is concoordant with                         #
#      sign y_i  ->  d_i = |x_i - y_i| / max(|x_i|, |y_i|)  #
# - if sign x_i is discordant with                          #
#      sign y_i  ->  d_i = 1                                #
#                                                           #
# Returns a Distance Matrix (either full dedup = FALSE or   #
# lower triangular if dedup = TRUE).                        #
# NOTE: Diffabs distance d_i for each pair (x_i, y_i) is    #
#       normalized in [0,1].                                #
#############################################################
# Strip spaces and convert to numeric
x <- as.numeric(tolower(gsub(" ", "", x)))
y <- as.numeric(tolower(gsub(" ", "", y)))

# Compute normalized distances between two numbers, a and b:
diffabs <- function(a, b) {
  sa <- sign(a)
  sb <- sign(b)
  sp <- sa * sb
  ss <- abs(sa + sb)
  r <- abs(a - b) / pmax(abs(a), abs(b))
  out <- ifelse(is.na(sp), NA, ifelse(sp > 0, r, ifelse(sp < 0, 1, ifelse(ss <= 0, 0, 1) ) ) )
  return(out)
}

d <- outer(x, y, FUN = diffabs)

# If x or/and y for a given pair is/are NA, put the distance
# to the (blind) average value 1/2. As this may concentrate too much
# the distance distribution, the caller function 'Distance' must be able
# to soften the importance of the variable according to its NA rate:
NAs <- is.na(d)
NAf <- sum(NAs)/prod(dim(NAs))

d[NAs] <- 0.5
d.vect <- as.numeric(d)
# If Deduplication is in order, then restrict to lower triangle of matrix d
if (dedup) {
    d[row(d) <= col(d)] <- NA
    d.vect <- d[row(d) > col(d)]
    NAs <- NAs[row(d) > col(d)]
    NAf <- sum(NAs)/length(NAs)
    }

# Compute the standard deviation of distances for the current variable.
# If it happens to be zero (i.e. all records from both files are identical
# w.r.t. the variable, hence all distances are zero), then the caller
# function 'Distance' will put the corresponding weight to zero.
s <- sd(d.vect)
# output: dist matrix
out <- list(d = d)
attr(out, "stdev") <- s
attr(out, "NAf") <- NAf
out
}


#-------------------------------------------------------------------------------

InRange <- function(x, y, dedup, ...) {
#######################################################################
# x: character, with the following structure 'min | max', where 'min' #
#    and 'max' can be safely converted to numeric. If either or both  #
#    'min' and 'max' are NA, then they have to appear in x as empty   #
#    strings ''.                                                      #
# y: character that can be safely converted to numeric.               #
#######################################################################
# Strip spaces and convert to numeric
x <- tolower(gsub(" ", "", x))
x <- c("-Inf|Inf", x)          #!!#

y <- as.numeric(tolower(gsub(" ", "", y)))

# Extract the range bounds from x
xl <- lapply( strsplit(x, "|", fixed = TRUE), FUN = as.numeric)
# Drop possible NAs coming from cases where either or both 'min' and 'max' are NA
xl[!is.na(xl)] <- lapply( xl[!is.na(xl)], FUN = function(el) el[!is.na(el)] )
# Convert the range bounds list to a 2 columns matrix (single values or NAs will
# be automatically replicated)
xm <- Reduce(f = rbind, xl)
xm <- xm[-1, , drop = FALSE]                #!!#
rownames(xm) <- NULL
colnames(xm) <- c("MIN", "MAX")

# Compare y values with the range bounds in xm
inrangeL <- outer( xm[, "MAX"], y, FUN = ">=" )
inrangeR <- outer( xm[, "MIN"], y, FUN = "<=" )
inrange <- inrangeL * inrangeR
# rownames(inrange) <- x
# colnames(inrange) <- y

# Map in range values to distance 0 and out of range values to distance 1
d <- 1 - inrange

# If x or/and y for a given pair is/are NA, put the distance
# to the (blind) average value 1/2. As this may concentrate too much
# the distance distribution, the caller function 'Distance' must be able
# to soften the importance of the variable according to its NA rate:
NAs <- is.na(d)
NAf <- sum(NAs)/prod(dim(NAs))

d[NAs] <- 0.5
d.vect <- as.numeric(d)
# If Deduplication is in order, then restrict to lower triangle of matrix d
if (dedup) {
    d[row(d) <= col(d)] <- NA
    d.vect <- d[row(d) > col(d)]
    NAs <- NAs[row(d) > col(d)]
    NAf <- sum(NAs)/length(NAs)
    }

# Compute the standard deviation of distances for the current variable.
# If it happens to be zero (i.e. all records from both files are identical
# w.r.t. the variable, hence all distances are zero), then the caller
# function 'Distance' will put the corresponding weight to zero.
s <- sd(d.vect)
# output: dist matrix
out <- list(d = d)
attr(out, "stdev") <- s
attr(out, "NAf") <- NAf
out
}


#-------------------------------------------------------------------------------

FormatSD <- function(x, y, dedup, ...) {
#######################################################################
# This function is domain specific, in that it has been designed to   #
# compare the 'FORMATO' field of Scanner Data (e.g. to address the    #
# relaunches identification task).                                    #
#                                                                     #
# x: character, with the following structure 'UNITA|QUANTITA', where  #
#    'UNITA' is a unit of measure and 'QUANTITA' is a quantity and    #
#    can safely be converted to numeric.                              #
# y: same as x.                                                       #
#                                                                     #
# NOTE: Both x and y can contain NAs. A NA value in x or y means that #
#       both the underlying 'UNITA' and 'QUANTITA' are missing.       #
#                                                                     #
# NOTE: pairs x_i, y_j such that either of the strings is NA are      #
#       conventionally mapped to distance 0.9: this is at odds with   #
#       other distances (e.g. Levenshtein) whose output is instead    #
#       set to 1/2.                                                   #
#######################################################################
# Strip spaces (if any) and convert to lower case
x <- tolower(gsub(" ", "", x))
x <- c("U|Q", x)          #!!#

y <- tolower(gsub(" ", "", y))
y <- c("U|Q", y)          #!!#

# Extract 'UNITA' and 'QUANTITA' from x and y
xl <- strsplit(x, "|", fixed = TRUE)

yl <- strsplit(y, "|", fixed = TRUE)

# Drop possible NAs coming from cases where either or both 'UNITA' and 'QUANTITA' are NA
xl[!is.na(xl)] <- lapply( xl[!is.na(xl)], FUN = function(el) el[!is.na(el)] )

yl[!is.na(yl)] <- lapply( yl[!is.na(yl)], FUN = function(el) el[!is.na(el)] )

# Convert the ('U', 'Q') list to a 2 columns matrix (single values or NAs will
# be automatically replicated)
xm <- Reduce(f = rbind, xl)
xm <- xm[-1, , drop = FALSE]                #!!#

colnames(xm) <- c("UNITA", "QUANTITA")
xU <- xm[, "UNITA"]
xQ <- as.numeric(xm[, "QUANTITA"])

ym <- Reduce(f = rbind, yl)
ym <- ym[-1, , drop = FALSE]                #!!#

colnames(ym) <- c("UNITA", "QUANTITA")
yU <- ym[, "UNITA"]
yQ <- as.numeric(ym[, "QUANTITA"])

# Compare U values
eqU <- 1 * outer(xU, yU, FUN = "==")

# Compare Q values
     diffabs <- function(a, b) {
     #############################################################
     # Compute normalized distances between two numeric vectors, #
     # a and b (handling a THRESHOLD):                           #
     #                                                           #
     # The function, d_i = diffabs(x_i, y_i), is as follows:     #
     # - if x_i = y_i = 0       ->  d_i = 0                      #
     # - if x_i != 0 & y_i = 0  ->  d_i = 1                      #
     # - if x_i = 0 & y_i != 0  ->  d_i = 1                      #
     # - if sign x_i is concoordant with                         #
     #      sign y_i  ->  d_i' = |x_i - y_i| / max(|x_i|, |y_i|) #
     #                    d_i'' = d_i' / d_THRESHOLD             #
     #                    d_i = sqrt(d_i'')                      #
     # - if sign x_i is discordant with                          #
     #      sign y_i  ->  d_i = 1                                #
     #                                                           #
     #############################################################
     sa <- sign(a)
     sb <- sign(b)
     sp <- sa * sb
     ss <- abs(sa + sb)
     r <- abs(a - b) / pmax(abs(a), abs(b))
     # HANDLE THE CUSTOMIZABLE THRESHOLD REQUIRED BY THE SCANNER DATA TEAM
     rTH <- 1/3
     r <- ifelse(r <= rTH, r/rTH, 1)
     r <- sqrt(r)
     out <- ifelse(is.na(sp), NA, ifelse(sp > 0, r, ifelse(sp < 0, 1, ifelse(ss <= 0, 0, 1) ) ) )
     return(out)
    }

dQ <- outer(xQ, yQ, FUN = diffabs)

# Compose U and Q distances
# If xU == yU -> d = dQ
# if xU != yU -> d = 1
d <- eqU * dQ
d[eqU <= 0] <- 1 

# If x or/and y for a given pair is/are NA, put the distance
# to the penalized value 9/10. As this may concentrate too much
# the distance distribution, the caller function 'Distance' must be able
# to soften the importance of the variable according to its NA rate:
NAs <- is.na(d)
NAf <- sum(NAs)/prod(dim(NAs))

d[NAs] <- 0.9
d.vect <- as.numeric(d)
# If Deduplication is in order, then restrict to lower triangle of matrix d
if (dedup) {
    d[row(d) <= col(d)] <- NA
    d.vect <- d[row(d) > col(d)]
    NAs <- NAs[row(d) > col(d)]
    NAf <- sum(NAs)/length(NAs)
    }

# Compute the standard deviation of distances for the current variable.
# If it happens to be zero (i.e. all records from both files are identical
# w.r.t. the variable, hence all distances are zero), then the caller
# function 'Distance' will put the corresponding weight to zero.
s <- sd(d.vect)
# output: dist matrix
out <- list(d = d)
attr(out, "stdev") <- s
attr(out, "NAf") <- NAf
out
}


#-------------------------------------------------------------------------------

JaroWinkler <- function(x, y, dedup, ...) {
##############################################################
# Normalized Jaro-Winkler distance between character vectors #
# lower-cased with spaces removed.                           #
# Returns a Distance Matrix (either full dedup = FALSE or    #
# lower triangular if dedup = TRUE).                         #
# NOTE: Levenshtein distance d_i for each pair (x_i, y_i)    #
#       is normalized in [0,1].                              #
##############################################################
require(RecordLinkage)
# Strip spaces
x <- tolower(gsub(" ", "", x))
y <- tolower(gsub(" ", "", y))
# Compute normalized distances:
d <- 1 - outer(x, y, jarowinkler)

# If x or/and y for a given pair is/are NA, put the distance
# to the (blind) average value 1/2. As this may concentrate too much
# the distance distribution, the caller function 'Distance' must be able
# to soften the importance of the variable according to its NA rate:
NAs <- is.na(d)
NAf <- sum(NAs)/prod(dim(NAs))

d[NAs] <- 0.5
d.vect <- as.numeric(d)
# If Deduplication is in order, then restrict to lower triangle of matrix d
if (dedup) {
    d[row(d) <= col(d)] <- NA
    d.vect <- d[row(d) > col(d)]
    NAs <- NAs[row(d) > col(d)]
    NAf <- sum(NAs)/length(NAs)
    }

# Compute the standard deviation of distances for the current variable.
# If it happens to be zero (i.e. all records from both files are identical
# w.r.t. the variable, hence all distances are zero), then the caller
# function 'Distance' will put the corresponding weight to zero.
s <- sd(d.vect)
# output: dist matrix
out <- list(d = d)
attr(out, "stdev") <- s
attr(out, "NAf") <- NAf
out
}


#-------------------------------------------------------------------------------

Cos.3.TfIdf <- function(x, y, dedup = FALSE, ...){
############################################################
# Normalized COS3TFIDF distance between character vectors  #
# lower-cased with spaces removed.                         #
# The distance computation requires 3 steps:               #
#   1. Extract and collect all the 3 grams from all the    #
#      strings belonging to vectors x and y.               #
#   2. Compute a TfIdf weighted term-document matrix by    #
#      using 3grams as terms and strings as documents.     #
#   3. Compute the cosine distance between columns of the  #
#      term-document matrix representing strings x_i, y_j. #
#                                                          #
# Returns a Distance Matrix (either full dedup = FALSE or  #
# lower triangular if dedup = TRUE).                       #
# NOTE: cos.3.tfidf distance d_i for each pair (x_i, y_i)  #
#       is normalized in [0,1].                            #
# NOTE: pairs x_i, y_j such that either of the strings is  #
#       NA are mapped to distance 1: this is at odds with  #
#       other distances (e.g. Levenshtein) whose output is #
#       instead 1/2.                                       #
############################################################
require(tm)
require(proxy)
# Strip spaces from input vectors
x <- tolower(gsub(" ", "", x))
y <- tolower(gsub(" ", "", y))
# Get vector lengths
m <- length(x)
n <- length(y)
# Put together strings to be compared
xy <- c(x, y)
# Build a vector source for tm pkg functions
xy.source <- VectorSource(xy)
# Build a volatile Corpus on this source
xy.corp <- VCorpus(xy.source)
# Compute the term-document matrix, where terms are 3grams
# NOTE: if some x or/and y are NA, they generate columns of
#       zeros in the term-document matrix
xy.tdm <- TermDocumentMatrix(xy.corp, control = list(tokenize =qgrams, weighting = weightTfIdf))
# Extract the matrix only
xy.mat <- as.matrix(xy.tdm)
# Following code commented: token appearing in either vectors,
# but not in both are STILL IMPORTANT since they matter in the
# final COSINE distance
# ------------------------------------------------------------
# drop tokens that appear in either vectors, but not in both  
# x.sum <- xy.mat[, 1:m, drop=FALSE] %*% rep(1, m)
# y.sum <- xy.mat[, (m+1):(m+n), drop=FALSE] %*% rep(1, n)
# xy.mat.r <- xy.mat[ (x.sum > 0) & (y.sum > 0), ,drop=FALSE]
# ------------------------------------------------------------

# Handle strings which have no common tokens with any other,
# e.g. NA strings (but not only them): their distance will be
# 1 with whatever other string but their cosine is ill defined
# (being they null vectors)
x.lonely <- (colSums(xy.mat[,         1:m, drop=FALSE]) <= 0)
y.lonely <- (colSums(xy.mat[, (m+1):(m+n), drop=FALSE]) <= 0)
# compute cosine similarity
  # first put matrix columns into lists
  l.x <- lapply(1:m, function(col) as.numeric(xy.mat[,     col]))
  l.y <- lapply(1:n, function(col) as.numeric(xy.mat[, m + col]))
  # build distance matrix
  d <- dist(l.x, l.y, "cosine")
  # asses NA frequency
  NAs <- is.na(d)
  NAf <- sum(NAs)/prod(dim(NAs))
  # recall lonely strings whose distance whith any others MUST be 1;
  # as this may concentrate too much the distance distribution, the
  # caller function 'Distance' must be able to soften the importance
  # of the variable according to its NA rate:
  d[x.lonely, ] <- 1
  d[, y.lonely] <- 1

d.vect <- as.numeric(d)
# If Deduplication is in order, then restrict to lower triangle of matrix d
if (dedup) {
    d[row(d) <= col(d)] <- NA
    d.vect <- d[row(d) > col(d)]
    NAs <- NAs[row(d) > col(d)]
    NAf <- sum(NAs)/length(NAs)
    }

# Compute the standard deviation of distances for the current variable.
# If it happens to be zero (i.e. all records from both files are identical
# w.r.t. the variable, hence all distances are zero), then the caller
# function 'Distance' will put the corresponding weight to zero.
s <- sd(d.vect)
# output: dist matrix
out <- list(d = d)
attr(out, "stdev") <- s
attr(out, "NAf") <- NAf
out
}

qgrams <- function(string, q = 3){
##################################
# Extracts qgrams from a string. #
# NOTE: by default 3grams are    #
#       computed.                #
##################################
  if (is.na(as.character(string))) return(NA)
  jolly <- paste(rep("°", q - 1), collapse = "")
  s <- paste(jolly, string, jolly, sep = "")
  sapply(1:(nchar(s) - q + 1), function(i) substr(s, i, i + q -1))
}


simCos.3.TfIdf <- function(x, y, dedup = FALSE, ...){
     dd <- Cos.3.TfIdf(x, y, dedup, ...)
     dmat <- dd$d
     one <- dmat
     one[] <- 1
     smat <- one - dmat
     d <- list(d = smat)
     attr(d, "stdev") <- attr(dd, "stdev")
     attr(d, "NAf") <- attr(dd, "NAf")
	 d
    }


#-------------------------------------------------------------------------------

Three.grams <- function(x, y, dedup = FALSE, ...){
############################################################
# Normalized THREEGRAMS distance between character vectors #
# lower-cased with spaces removed.                         #
# Returns a Distance Matrix (either full dedup = FALSE or  #
# lower triangular if dedup = TRUE).                       #
# NOTE: THREEGRAMS distance d_i for each pair (x_i, y_i)   #
#       is normalized in [0,1].                            #
############################################################

Qdist1 <- function(s1, s2, q){
######################################
# Qgrams distance between 2 strings. #
######################################
if ( ( is.na(s1) || is.na(s2) ) || ( identical(s1,"") || identical(s2,"") ) ) return(1)
q1   <- qgrams(s1, q = q)
q2   <- qgrams(s2, q = q)
int  <- intersect(q1, q2)
uni  <- union(q1, q2)
dist <- 1 - length(int)/length(uni)
dist
}

Qdistn <- function(svect1, svect2, q){
########################################
# Qgrams distance between 2 vectors of #
# strings: returns a matrix.           #
########################################
dist <- mapply(Qdist1, svect1, svect2, MoreArgs = list(q = q))
names(dist) <- NULL
dist
}

# Strip spaces
x <- tolower(gsub(" ", "", x))
y <- tolower(gsub(" ", "", y))
# Compute normalized distances:
d <- outer(x, y, Qdistn, q = 3)

# If x or/and y for a given pair is/are NA, put the distance
# to the (blind) average value 1/2. As this may concentrate too much
# the distance distribution, the caller function 'Distance' must be able
# to soften the importance of the variable according to its NA rate:
NAs <- is.na(d)
NAf <- sum(NAs)/prod(dim(NAs))

d[NAs] <- 0.5
d.vect <- as.numeric(d)
# If Deduplication is in order, then restrict to lower triangle of matrix d
if (dedup) {
    d[row(d) <= col(d)] <- NA
    d.vect <- d[row(d) > col(d)]
    NAs <- NAs[row(d) > col(d)]
    NAf <- sum(NAs)/length(NAs)
    }

# Compute the standard deviation of distances for the current variable.
# If it happens to be zero (i.e. all records from both files are identical
# w.r.t. the variable, hence all distances are zero), then the caller
# function 'Distance' will put the corresponding weight to zero.
s <- sd(d.vect)
# output: dist matrix
out <- list(d = d)
attr(out, "stdev") <- s
attr(out, "NAf") <- NAf
out
}


#-------------------------------------------------------------------------------

VSM.TfIdf <- function(x, y, dedup = FALSE, ...){
############################################################
# Normalized VSMTFIDF distance (i.e. Vector Space Model    #
# with TfIdf weights) between lower-cased character        #
# vectors.                                                 #
# The distance computation requires 3 steps:               #
#   1. Extract and collect all the 'words' from all the    #
#      strings belonging to vectors x and y.               #
#   2. Compute a TfIdf weighted term-document matrix by    #
#      using 'words' as terms and strings as documents.    #
#   3. Compute the cosine distance between columns of the  #
#      term-document matrix representing strings x_i, y_j. #
#                                                          #
# Returns a Distance Matrix (either full if dedup = FALSE  #
# or lower triangular if dedup = TRUE).                    #
# NOTE: VSM.TfIdf distance d_i for each pair (x_i, y_i)    #
#       is normalized in [0,1].                            #
# NOTE: pairs x_i, y_j such that either of the strings is  #
#       NA are mapped to distance 1: this is at odds with  #
#       other distances (e.g. Levenshtein) whose output is #
#       instead 1/2.                                       #
############################################################
require(tm)
require(proxy)
# Strip spaces from input vectors
x <- tolower(x)
y <- tolower(y)
# Get vector lengths
m <- length(x)
n <- length(y)
# Put together strings to be compared
xy <- c(x, y)
# Build a vector source for tm pkg functions
xy.source <- VectorSource(xy)
# Build a volatile Corpus on this source
xy.corp <- VCorpus(xy.source)
# Compute the term-document matrix, where terms are words
# NOTE: if some x or/and y are NA, they generate columns of
#       zeros in the term-document matrix
xy.tdm <- TermDocumentMatrix(xy.corp,
                             control = list(
                                            tokenize = NULL,
                                            weighting = weightTfIdf,
                                            removePunctuation = TRUE,
                                            stemming = FALSE,
                                            stopwords = FALSE,
                                            dictionary = NULL,
                                            minDocFreq=1,
                                            minWordLength=3
                                            )
                            )
# Extract the matrix only
xy.mat <- as.matrix(xy.tdm)
# Following code commented: token appearing in either vectors,
# but not in both are STILL IMPORTANT since they matter in the
# final COSINE distance
# ------------------------------------------------------------
# drop tokens that appear in either vectors, but not in both  
# x.sum <- xy.mat[, 1:m, drop=FALSE] %*% rep(1, m)
# y.sum <- xy.mat[, (m+1):(m+n), drop=FALSE] %*% rep(1, n)
# xy.mat.r <- xy.mat[ (x.sum > 0) & (y.sum > 0), ,drop=FALSE]
# ------------------------------------------------------------

# Handle strings which have no common tokens with any other,
# e.g. NA strings (but not only them): their distance will be
# 1 with whatever other string but their cosine is ill defined
# (being they null vectors)
x.lonely <- (colSums(xy.mat[,         1:m, drop=FALSE]) <= 0)
y.lonely <- (colSums(xy.mat[, (m+1):(m+n), drop=FALSE]) <= 0)
# compute cosine similarity
  # first put matrix columns into lists
  l.x <- lapply(1:m, function(col) as.numeric(xy.mat[,     col]))
  l.y <- lapply(1:n, function(col) as.numeric(xy.mat[, m + col]))
  # build distance matrix
  d <- dist(l.x, l.y, "cosine")
  # asses NA frequency
  NAs <- is.na(d)
  NAf <- sum(NAs)/prod(dim(NAs))
  # recall lonely strings whose distance whith any others MUST be 1;
  # as this may concentrate too much the distance distribution, the
  # caller function 'Distance' must be able to soften the importance
  # of the variable according to its NA rate:
  d[x.lonely, ] <- 1
  d[, y.lonely] <- 1

d.vect <- as.numeric(d)
# If Deduplication is in order, then restrict to lower triangle of matrix d
if (dedup) {
    d[row(d) <= col(d)] <- NA
    d.vect <- d[row(d) > col(d)]
    NAs <- NAs[row(d) > col(d)]
    NAf <- sum(NAs)/length(NAs)
    }

# Compute the standard deviation of distances for the current variable.
# If it happens to be zero (i.e. all records from both files are identical
# w.r.t. the variable, hence all distances are zero), then the caller
# function 'Distance' will put the corresponding weight to zero.
s <- sd(d.vect)
# output: dist matrix
out <- list(d = d)
attr(out, "stdev") <- s
attr(out, "NAf") <- NAf
out
}


simVSM.TfIdf <- function(x, y, dedup = FALSE, ...){
     dd <- VSM.TfIdf(x, y, dedup, ...)
     dmat <- dd$d
     one <- dmat
     one[] <- 1
     smat <- one - dmat
     d <- list(d = smat)
     attr(d, "stdev") <- attr(dd, "stdev")
     attr(d, "NAf") <- attr(dd, "NAf")
	 d
    }


#-------------------------------------------------------------------------------

emailBody <- function(x, y, dedup = FALSE, qcut = 0.99, ...){
############################################################
# This is the same as VSM.TfIdf with a modification needed #
# to scale on huge email corpora.                          #
# The trick is to restrict the term-document matrix to an  #
# 'effective dictionary' which is obtained by subsetting   #
# the original one.                                        #
# Currently the subsetting retains only the terms tied to  #
# the top-(1-qcut)% TfIdf scores.                          #
############################################################
require(tm)
require(proxy)
# Strip spaces from input vectors
x <- tolower(x)
y <- tolower(y)
# Get vector lengths
m <- length(x)
n <- length(y)
# Put together strings to be compared
xy <- c(x, y)
# Build a vector source for tm pkg functions
xy.source <- VectorSource(xy)
# Build a volatile Corpus on this source
xy.corp <- VCorpus(xy.source)
# Compute the term-document matrix, where terms are words
# NOTE: if some x or/and y are NA, they generate columns of
#       zeros in the term-document matrix
xy.tdm <- TermDocumentMatrix(xy.corp,
                             control = list(
                                            tokenize = NULL,
                                            weighting = weightTfIdf,
                                            removePunctuation = TRUE,
                                            stemming = FALSE,
                                            stopwords = FALSE,
                                            dictionary = NULL,
                                            minDocFreq=1,
                                            minWordLength=3
                                            )
                            )

# Subset the term-document matrix for scalability issues
if (qcut > 0) {
    # START subsetting 
    # Compute the qcut percentile (defaults to the 99th) of the TfIdf scores distribution
    # cut.score <- quantile(xy.tdm$v, probs = qcut)
    # Retrieve the terms tied to such top-1% scores
    # eff.dictionary <- unique(xy.tdm$dimnames$Terms[xy.tdm$i[xy.tdm$v >= cut.score]])
      # NOTE: a slower version could sort the terms by decreasing score, this
      # may be useful if one want to further distill the effective dictionary
      scores <- xy.tdm$v
      term.i <- xy.tdm$i
      o <- order(scores, decreasing = TRUE)
      o.scores <- scores[o]
      o.term.i <- term.i[o]
      cut.score <- quantile(o.scores, probs = qcut)
      eff.dictionary <- unique(xy.tdm$dimnames$Terms[o.term.i[o.scores >= cut.score]])
    # Use such terms as an effective dictionary
    xy.tdm <- TermDocumentMatrix(xy.corp,
                                 control = list(
                                                tokenize = NULL,
                                                weighting = weightTfIdf,
                                                removePunctuation = TRUE,
                                                stemming = FALSE,
                                                stopwords = FALSE,
                                                dictionary = eff.dictionary,
                                                minDocFreq=1,
                                                minWordLength=3
                                                )
                                )
    # END subsetting
    }

# Extract the matrix only
xy.mat <- as.matrix(xy.tdm)
# Following code commented: token appearing in either vectors,
# but not in both are STILL IMPORTANT since they matter in the
# final COSINE distance
# ------------------------------------------------------------
# drop tokens that appear in either vectors, but not in both  
# x.sum <- xy.mat[, 1:m, drop=FALSE] %*% rep(1, m)
# y.sum <- xy.mat[, (m+1):(m+n), drop=FALSE] %*% rep(1, n)
# xy.mat.r <- xy.mat[ (x.sum > 0) & (y.sum > 0), ,drop=FALSE]
# ------------------------------------------------------------

# Handle strings which have no common tokens with any other,
# e.g. NA strings (but not only them): their distance will be
# 1 with whatever other string but their cosine is ill defined
# (being they null vectors)
x.lonely <- (colSums(xy.mat[,         1:m, drop=FALSE]) <= 0)
y.lonely <- (colSums(xy.mat[, (m+1):(m+n), drop=FALSE]) <= 0)
# compute cosine similarity
  # first put matrix columns into lists
  l.x <- lapply(1:m, function(col) as.numeric(xy.mat[,     col]))
  l.y <- lapply(1:n, function(col) as.numeric(xy.mat[, m + col]))
  # build distance matrix
  d <- dist(l.x, l.y, "cosine")
  # asses NA frequency
  NAs <- is.na(d)
  NAf <- sum(NAs)/prod(dim(NAs))
  # recall lonely strings whose distance whith any others MUST be 1;
  # as this may concentrate too much the distance distribution, the
  # caller function 'Distance' must be able to soften the importance
  # of the variable according to its NA rate:
  d[x.lonely, ] <- 1
  d[, y.lonely] <- 1

d.vect <- as.numeric(d)
# If Deduplication is in order, then restrict to lower triangle of matrix d
if (dedup) {
    d[row(d) <= col(d)] <- NA
    d.vect <- d[row(d) > col(d)]
    NAs <- NAs[row(d) > col(d)]
    NAf <- sum(NAs)/length(NAs)
    }

# Compute the standard deviation of distances for the current variable.
# If it happens to be zero (i.e. all records from both files are identical
# w.r.t. the variable, hence all distances are zero), then the caller
# function 'Distance' will put the corresponding weight to zero.
s <- sd(d.vect)
# output: dist matrix
out <- list(d = d)
attr(out, "stdev") <- s
attr(out, "NAf") <- NAf
out
}


#-------------------------------------------------------------------------------