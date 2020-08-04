cloud <- function(docs,...){
require(tm)
require(wordcloud)
# Build a vector source for tm pkg functions
source <- VectorSource(docs, encoding = "UTF-8")
# Build a volatile Corpus on this source
corp <- Corpus(source)
# Compute the term-document matrix
tdm <- TermDocumentMatrix(corp, control = list(weighting = weightTf, removePunctuation = TRUE))
# Extract the matrix only
mat <- as.matrix(tdm)
word <- rownames(mat)
freq <- rowSums(mat)
wordcloud(word,freq,...)
}
