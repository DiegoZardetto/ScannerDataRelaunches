gsub.df <- function (df, pattern, replacement, except = NULL, ...) 
{
    if (!inherits(df, "data.frame")) stop("Works only on dataframes!")
    idf <- vapply(df, is.character, NA)
    if (!is.null(except)) 
        idf[except] <- FALSE
    df[idf] <- lapply(df[idf], function(x) gsub(pattern, replacement, x, ...))
    df
}

FactorsAsStrings.df <- function(df, except = NULL)
{
    if (!inherits(df, "data.frame")) stop("Works only on dataframes!")
    idf <- vapply(df, is.factor, NA)
    if (!is.null(except)) 
        idf[except] <- FALSE
    df[idf] <- lapply(df[idf], as.character)
    df
}
