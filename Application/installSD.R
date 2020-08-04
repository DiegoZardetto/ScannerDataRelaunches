installSD <- function(path, trace = TRUE, ...){
# Source the Scanner Data code
     for (code in list.files(path, pattern = "\\.[RrSsQq]$")) {
         sys.source(file.path(path, code), envir=.GlobalEnv, ...)
         if(trace) cat(code,"\t sourced\n")
        }
     for (dat in list.files(path,pattern= "\\.RData$")) {
         load(file.path(path, dat), envir=.GlobalEnv)
         if(trace) cat(dat,"\t loaded\n")           
        }
# Verify that all the needed packages are available
     require(proxy)
     require(tm)
     # require(RecordLinkage)
}