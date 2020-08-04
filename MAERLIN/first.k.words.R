first.k.words <- function(x, k){
require(tm)
require(proxy)
# Strip spaces from input vectors
x <- tolower(x)
y <- x
# Get vector lengths
m <- length(x)
n <- length(y)
# Put together strings to be compared
xy <- c(x, y)
# Build a vector source for tm pkg functions
xy.source <- VectorSource(xy, encoding = "UTF-8")
# Build a volatile Corpus on this source
xy.corp <- Corpus(xy.source)
# Compute the term-document matrix, where terms are 3grams
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

# Subset the term-document matrix for finding the top k score words
scores <- xy.tdm$v
term.i <- xy.tdm$i
o <- order(scores, decreasing = TRUE)
o.term.i <- term.i[o]
words <- unique(xy.tdm$dimnames$Terms[o.term.i])
k.words <- words[1:k]
k.words
}

# tapply(clust.df.mails.en$body,clust.df.mails.en$cluster.id, first.k.words, k=30)