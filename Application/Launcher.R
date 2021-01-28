# START message
release="3";
source("file.param")

.libPaths( Dir_Packages )

if (!file.exists(OutputDir_ANNO_MESE)) {
    dir.create(file.path(OutputDir_ANNO_MESE))
}

LogFile <- paste(OutputDir_ANNO_MESE, "LogRilanci.txt", sep = '\\')
MessageLogFile <- paste(OutputDir_ANNO_MESE, "LogRilanciMsg.txt", sep = '\\')

LogConn <- file(LogFile, open = "w")
sink(LogConn, type = "output")
MsgConn <- file(MessageLogFile, open = "w")
sink(MsgConn, type = "output", split = TRUE)
sink(MsgConn, type = "message")

timestamp(stamp = "Esecuzione Procedura Rilanci da Scanner Data - START", prefix = "\n###### ", suffix = "", quiet = FALSE)
timestamp(prefix = "###### ", suffix = "", quiet = FALSE)

fileinpA <- paste(InputDir_ANNO_MESE, FileA, sep = '\\')
fileinpB <- paste(InputDir_ANNO_MESE, FileB, sep = '\\')

# Build the blocks
cat("\n:::::: Parametri letti:\n")
cat("Release: ", release, "\n")
cat("InputDir_ANNO_MESE:", InputDir_ANNO_MESE, "\n")
cat("FileA:", FileA, "\n")
cat("FileB:", FileB, "\n")
cat("OutputDir_ANNO_MESE:", OutputDir_ANNO_MESE, "\n")
cat("Dir_Packages:", Dir_Packages, "\n")

cat("\n:::::: Carico i sorgenti di MAERLIN\n")
source(".\\installSD.R")
installSD("..\\Sources")

cat("\n:::::: Leggo e preparo i dataset di input\n")

# Split 'FORMATO' information into Unit of Measure and Value 
mkUnit  <- function(x) gsub(" ", "", gsub("[0-9]*[.]*", "", toupper(x)))
mkValue <- function(x) gsub(" ", "", gsub("[A-Z]", "", toupper(x)))

nast <- c("", "NA", "na", "Na", "nA", "NAN", "nan", "NaN", "nAn", "NR", "nr", "NOT IN RULE", "NOTINRULE", "not in rule","notinrule", "BADVALUE", "BAD VALUE", "Bad Value", "bad value", "BadValue", "badvalue","ND","\\N")

# Read the input files as they are
# A <- read.table(fileinpA, header = TRUE, sep = ";", na.strings = nast, stringsAsFactors = FALSE, encoding = "UTF-8", quote = "", comment.char = "")
A <- read.table(fileinpA, header = TRUE, sep = ";", na.strings = nast, stringsAsFactors = FALSE, quote = "", comment.char = "")
names(A) <- toupper(names(A))
invisible(sapply(names(A), function(col) A[[col]] <<- as.character(A[[col]])))

# B <- read.table(fileinpB, header = TRUE, sep = ";", na.strings = nast, stringsAsFactors = FALSE, encoding = "UTF-8", quote = "", comment.char = "")
B <- read.table(fileinpB, header = TRUE, sep = ";", na.strings = nast, stringsAsFactors = FALSE, quote = "", comment.char = "")
names(B) <- toupper(names(B))
invisible(sapply(names(B), function(col) B[[col]] <<- as.character(B[[col]])))

# The following columns must exist in both input files
mandatory.vars <- c("PRDKEY", "MERCATO", "MARCA", "DESCRIZIONE", "FORMATO","NUM_PEZZI")
  # Check for missing mandatory variables in file A
A.miss.vars <- mandatory.vars[!(mandatory.vars %in% names(A))]
if (length(A.miss.vars) > 0) stop("Missing mandatory variables in file A: ", paste(A.miss.vars, collapse = ", "))
  # Check for missing mandatory variables in file B
B.miss.vars <- mandatory.vars[!(mandatory.vars %in% names(B))]
if (length(B.miss.vars) > 0) stop("Missing mandatory variables in file B: ", paste(B.miss.vars, collapse = ", "))

# New var DESC: clean the "DESCRIZIONE" column, deleting the "MARCA" substring when present
# New var DESCRIZIONE1: duplicate DESCRIZIONE
A[["DESC"]] <- sapply(1:nrow(A), function(i) ifelse(!is.na(A[["MARCA"]][i]), gsub(pattern = A[["MARCA"]][i], replacement = "", x = A[["DESCRIZIONE"]][i], fixed = TRUE), A[["DESCRIZIONE"]][i]))
A[["DESCRIZIONE1"]] <- A[["DESCRIZIONE"]]
# A[["DESCRIZIONE"]] <- NULL

B[["DESC"]] <- sapply(1:nrow(B), function(i) ifelse(!is.na(B[["MARCA"]][i]), gsub(pattern = B[["MARCA"]][i], replacement = "", x = B[["DESCRIZIONE"]][i], fixed = TRUE), B[["DESCRIZIONE"]][i]))
B[["DESCRIZIONE1"]] <- B[["DESCRIZIONE"]]
# B[["DESCRIZIONE"]] <- NULL

# Split the "FORMATO" column
A[["UNITA"]] <- mkUnit(A[["FORMATO"]])
A[["QUANTITA"]] <- mkValue(A[["FORMATO"]])
# Reformat "FORMATO" column ("UNITA|QUANTITA")
A[["UNITA|QUANTITA"]] <- A[["FORMATO"]]
A[!is.na(A[["UNITA|QUANTITA"]]), "UNITA|QUANTITA"] <- paste(A[!is.na(A[["UNITA|QUANTITA"]]), "UNITA"], A[!is.na(A[["UNITA|QUANTITA"]]), "QUANTITA"], sep = "|")

B[["UNITA"]] <- mkUnit(B[["FORMATO"]])
B[["QUANTITA"]] <- mkValue(B[["FORMATO"]])
# Reformat "FORMATO" column ("UNITA|QUANTITA")
B[["UNITA|QUANTITA"]] <- B[["FORMATO"]]
B[!is.na(B[["UNITA|QUANTITA"]]), "UNITA|QUANTITA"] <- paste(B[!is.na(B[["UNITA|QUANTITA"]]), "UNITA"], B[!is.na(B[["UNITA|QUANTITA"]]), "QUANTITA"], sep = "|")

str(A)
str(B)

setwd(OutputDir_ANNO_MESE)

cat("\n:::::: Setto il seme del generatore di numeri casuali di R per assicurare la piena riproducibilita' delle elaborazioni\n")
set.seed(51713773)

#------------------------------------------------------------------------------#
cat("\n:::::: Eseguo MAERLIN\n")

res.ab <- B.maerlin(A, B, block.vars = "MERCATO", match.vars = c("MARCA", "DESCRIZIONE1", "DESC", "UNITA|QUANTITA", "NUM_PEZZI"), dist.funs = c("Equality", "VSM.TfIdf", "Three.grams", "FormatSD", "Equality"), weights = c(1, 1, 1, 1, 1))

cat("\n:::::: Scrittura file di output\n")

links.ab <- attr(res.ab, "Matches")
write.table(links.ab, file = "Links.csv", sep = ";", row.names = FALSE)

PRDK.links <- getPRDKEY(links.ab)
write.table(PRDK.links, file = "ID_Links.csv", sep = ";", row.names = FALSE, quote = FALSE, col.names = FALSE)

cat("\n:::::: Applico le regole di esclusione deterministiche:")
cat("\n:::::: REGOLA 1: I rilanci non possono avere valori diversi in UNITA")
cat("\n:::::: REGOLA 2: I rilanci non possono avere una differenza relativa in QUANTITA superiore al 25%")
cat("\n:::::: REGOLA 3: I rilanci non possono avere valori mancanti in QUANTITA")
cat("\n:::::: REGOLA 4: I rilanci non possono avere valori diversi in MARCA")
cat("\n:::::: REGOLA 5: I rilanci non possono avere valori diversi in NUM_PEZZI")
cat("\n")

alives.ab <- KillLinks(links.ab, Q.th = 0.25)
write.table(alives.ab, file = "Links_ACCETTATI.csv", sep = ";", row.names = FALSE)

PRDK.alives <- getPRDKEY(alives.ab)
write.table(PRDK.alives, file = "ID_Links_ACCETTATI.csv", sep = ";", row.names = FALSE, quote = FALSE, col.names = FALSE)


killed.ab <- attr(alives.ab, "killed")
write.table(killed.ab, file = "Links_SCARTATI.csv", sep = ";", row.names = FALSE)

PRDK.killed <- getPRDKEY(killed.ab)
write.table(PRDK.killed, file = "ID_Links_SCARTATI.csv", sep = ";", row.names = FALSE, quote = FALSE, col.names = FALSE)

# The following columns must exist in the output file of relaunches 'alives.ab' to allow downstream ELABORATIONS:
ELAB.vars <- c("FILE", "PRDKEY", "YEAR", "MONTH", "INDICODECR", "COICOP", "UNITA", "QUANTITA")
# Check for missing ELAB variables in the output file 'alives.ab'
miss.ELAB.vars <- ELAB.vars[!(ELAB.vars %in% names(alives.ab))]
# In case of missing ELAB variables, print an error-alike message to inform the users of the ELAB file 
if (length(miss.ELAB.vars) > 0) cat("\nERROR: Missing required variables in output file Links_ACCETTATI_ELAB.csv: ", paste(miss.ELAB.vars, collapse = ", "), "\n")

ELAB.file.name <- "Links_ACCETTATI_ELAB.csv"
# In case the ELAB file already exists, delete it
if (file.exists(ELAB.file.name)) file.remove(ELAB.file.name)

# Only in case *all* the ELAB vars are present, write the ELAB csv file (no header, no quotes)
if (length(miss.ELAB.vars) == 0) {
     ELAB <- alives.ab[, ELAB.vars]
     write.table(ELAB, file = ELAB.file.name, sep = ";", row.names = FALSE, quote = FALSE, col.names = FALSE)
}

#------------------------------------------------------------------------------#

cat("\n:::::: Salvo il workspace R\n")
save.image(file="MAERLIN.RData")

cat("\n\n")
cat("******************************************\n")
cat(" Prodotti uscenti:     ", nrow(A), "\n")
cat(" Nuovi in anagrafica:  ", nrow(B), "\n")
cat(" Mercati interessati:  ", length(attr(res.ab, "block.size")), "\n")
cat(" Rilanci dal modello:  ", nrow(PRDK.links), "\n")
cat(" Scarti per regole:    ", nrow(PRDK.killed), "\n")
cat("\n")
cat(" Rilanci accettati:    ", nrow(PRDK.alives), "\n")
cat("******************************************\n")
cat("\n")

timestamp(stamp = "Esecuzione Procedura Rilanci da Scanner Data - END", prefix = "\n###### ", suffix = "", quiet = FALSE)
timestamp(prefix = "###### ", suffix = "", quiet = FALSE)

sink(type = "message")
sink(type = "output")
sink(type = "output")
close(MsgConn)
close(LogConn)
closeAllConnections()
