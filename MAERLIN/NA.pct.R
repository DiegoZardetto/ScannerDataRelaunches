NA.pct <-
function(file){
vars <- names(file)
chars <- sapply(vars, function(var) is.character(file[, var]))
varchar <- vars[chars]
sapply(varchar, function(var) round(100* sum(is.na(file[,var])/nrow(file)), 2))
}