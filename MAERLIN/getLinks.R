getLinks <- function(object, file.a, file.b,
                     D.mat = NULL, mix.par = NULL, value = c("data.frame", "list"), a.small = TRUE){
#############################################################################
# Extracts record pairs that have been declared Matches by the achieved     #
# solution of a RL job.                                                     #
#                                                                           #
# NOTE: Argument a.small is meaningful only when getLinks is called in a    #
#       block-wise RL job. In these cases the semantics is that !a.small    #
#       flags blocks where file.a, though being smaller than file.b,        #
#       actually comes from the the original file B.                        #
#############################################################################

has.mix.par <- FALSE
has.D.mat <- FALSE
conformF(object, file.a, file.b)
if (inherits(object, "ea.out")){
     if ( xor(is.null(D.mat), is.null(mix.par)) )
         stop("For input objects of class 'ea.out' D.mat and mix.par must be passed both or be both NULL!")
    }
if (!is.null(D.mat)) {
     conformD(object, D.mat)
     has.D.mat <- TRUE
    }
if (!is.null(mix.par)) {
     conformF(mix.par, file.a, file.b)
     actual.mix.par <- mix.par
     has.mix.par <- TRUE
    }
if (inherits(object, "reclink.mix.par") || inherits(object, "dedup.mix.par")) {
     actual.mix.par <- object
     has.mix.par <- TRUE
     if (!is.null(mix.par)) {
         if (!identical(actual.mix.par, mix.par))
             warning("Passed mix.par does not agree with input object: ignored", immediate. = TRUE)
        }
    }
value <- match.arg(value)
common <- intersect(names(file.a), names(file.b))
mat.sol <- getMatSol(object)
n.links <- nrow(mat.sol)
if (n.links==0) {
     # cat("No Matches have been declared!\n")
	 # return(invisible(NULL))
     return("No Matches have been declared!")
    }
if (a.small) {
     a.FileName <- "A"
     b.FileName <- "B"
    }
else {
     a.FileName <- "B"
     b.FileName <- "A"
    }

a.l <- mat.sol[, 1]
fa.l <- file.a[a.l, common]
fa.l <- data.frame("FILE" = a.FileName, "ROW" = rownames(fa.l), fa.l)
rownames(fa.l) <- paste(a.FileName, rownames(fa.l), sep=".")
b.l <- mat.sol[, 2]
fb.l <- file.b[b.l, common]
fb.l <- data.frame("FILE" = b.FileName, "ROW" = rownames(fb.l), fb.l)
rownames(fb.l) <- paste(b.FileName, rownames(fb.l), sep=".")
links.df <- rbind(fa.l, fb.l)
if (has.mix.par && has.D.mat) {
    ## Add the pairwise distance? - START
     # d <- sapply(seq(along.with = a.l),
                    # function (i) {
                       # as.numeric(D.mat[a.l[i], b.l[i]])
                      # }
                    # )
     # links.df[["d"]] <- rep(d, 2)
    ## Add the pairwise distance? - END

     fm.d <- sapply(seq(along.with = a.l),
                    function (i) {
                       actual.mix.par[["fm.d"]](D.mat[a.l[i], b.l[i]])
                      }
                    )
     links.df[["fm.d"]] <- rep(fm.d, 2)
    }

if (a.small) {
     links.df <- links.df[rep(1:n.links, each = 2) + c(0, n.links),]
    }
else {
     links.df <- links.df[rep(1:n.links, each = 2) + c(n.links, 0),]
    }

if (has.mix.par && has.D.mat) {
    ## Order by decreasing posterior match probability
     links.df <- links.df[do.call(order, c(links.df["fm.d"], decreasing = TRUE)), ]
    }

switch(value, "data.frame" = links.df,
              "list"       = lapply(seq(1, nrow(links.df), by = 2), function(i) links.df[c(i, i+1), ]))
}



MAP.Links <- function(mix.par, file.a, file.b, D.mat = NULL, value = c("data.frame", "list"), a.small = TRUE){
#############################################################################
# Extracts record pairs that would be declared Matches by the MAP rule.     #
#                                                                           #
#                                                                           #
# NOTE: Argument a.small is meaningful only when getLinks is called in a    #
#       block-wise RL job. In these cases the semantics is that !a.small    #
#       flags blocks where file.a, though being smaller than file.b,        #
#       actually comes from the the original file B.                        #
#############################################################################
if (!inherits(mix.par, "reclink.mix.par") && !inherits(mix.par, "dedup.mix.par"))
     stop("Wrong input object")
has.D.mat <- FALSE
conformF(mix.par, file.a, file.b)
if (!is.null(D.mat)) {
     conformD(mix.par, D.mat)
     has.D.mat <- TRUE
    }
value <- match.arg(value)
common <- intersect(names(file.a), names(file.b))
MAP <- mix.par$MAP.l
if (nrow(MAP)==0) {
     cat("MAP rule declared no Matches!\n")
     return(invisible(NULL))
    }
n.links <- nrow(MAP)

if (a.small) {
     a.FileName <- "A"
     b.FileName <- "B"
    }
else {
     a.FileName <- "B"
     b.FileName <- "A"
    }

a.l <- MAP[, 1]
fa.l <- file.a[a.l, common]
fa.l <- data.frame("FILE" = a.FileName, "ROW" = rownames(fa.l), fa.l)
rownames(fa.l) <- paste(a.FileName, rownames(fa.l), sep=".")
b.l <- MAP[, 2]
fb.l <- file.b[b.l, common]
fb.l <- data.frame("FILE" = b.FileName, "ROW" = rownames(fb.l), fb.l)
rownames(fb.l) <- paste(b.FileName, rownames(fb.l), sep=".")
links.df <- rbind(fa.l, fb.l)
if (has.D.mat) {
    ## Add the pairwise distance? - START
     # d <- sapply(seq(along.with = a.l),
                    # function (i) {
                       # as.numeric(D.mat[a.l[i], b.l[i]])
                      # }
                    # )
     # links.df[["d"]] <- rep(d, 2)
    ## Add the pairwise distance? - END

     fm.d <- sapply(seq(along.with = a.l),
                    function (i) {
                       mix.par[["fm.d"]](D.mat[a.l[i], b.l[i]])
                      }
                    )
	 links.df[["fm.d"]] <- rep(fm.d, 2)
    }

if (a.small) {
     links.df <- links.df[rep(1:n.links, each = 2) + c(0, n.links),]
    }
else {
     links.df <- links.df[rep(1:n.links, each = 2) + c(n.links, 0),]
    }

if (has.D.mat) {
     ## Order by decreasing posterior match probability
     links.df <- links.df[do.call(order, c(links.df["fm.d"], decreasing = TRUE)), ]
    }

switch(value, "data.frame" = links.df,
              "list"       = lapply(seq(1, nrow(links.df), by = 2), function(i) links.df[c(i, i+1), ]))
}