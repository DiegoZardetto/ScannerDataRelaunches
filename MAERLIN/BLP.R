blp <- function(loglikmat){
library(lpSolve)
m <- nrow(loglikmat)
n <- ncol(loglikmat)
# initialize constraints matrix
z <- matrix(0, m+n, prod(m,n))
# handle row constraints
for (i in 1:m){
    r <- rep(0, m)
    r[i] <- 1
    z[i, ] <- rep(r,n)
    }
# handle column constraints
for (j in 1:n){
    z[m + j, (prod(m, j - 1) + 1):prod(m, j)] <- 1
    }
# handle constraint equations
const.dir <- rep("<=", m+n)
const.rhs <- rep(1, m+n)
# handle objective function coefficients
ll <- as.numeric(loglikmat)
sol <- lp (direction = "max", objective.in=ll,
           const.mat=z, const.dir=const.dir, const.rhs=const.rhs,
	       all.bin=TRUE)
loglikmat[,] <- sol$solution
loglikmat
}

mat.to.sol <- function(mat) {
# Transform result from matrix to vector (MAERLIN style)
ind.row <- function(r) if(sum(r)<=0) 0 else which(r==1)
apply(mat,1, ind.row)
}


mat.to.list <- function(loglikmat, ppostmat){
# Transform loglikelyhood + ppost matrices to "list" dataframe (RELAIS style)
m <- nrow(loglikmat)
n <- ncol(loglikmat)
i <- rep(1:m,n)
j <- rep(1:n, each=m)
patt <- rep("patt", prod(m,n))
out<-data.frame("I"=i,
                "J"=j,
				"PATTERN"=patt,
				"R"=as.numeric(loglikmat),
				"P_POST"=as.numeric(ppostmat))
out
}

list.to.sol <- function(inlist, D.mat){
m <- nrow(D.mat)
sol <- rep(0, m)
sapply(1:nrow(inlist), function(i) sol[inlist[i,1]] <<- inlist[i,2])
sol
}



o2o <- function(mix.par, D.mat, prune=TRUE, out.style=c("vec", "list")){
#carica pacchetto lpSolve
library("lpSolve")
# library("RODBC")
#lettura parametri input 
# source("paramOtoO.R")

#lettura dei dati
#connssione al db relais

#con <- odbcConnect("relais",case="mysql")
#strInput = paste(paste('select * from ', inputOtoOTableName),';')
# w <- sqlQuery(con, strInput, as.is = TRUE)
loglik.mat <- mix.par[["loglik"]](D.mat)
ppost.mat <- mix.par[["fm.d"]](D.mat)
w <- mat.to.list(loglik.mat, ppost.mat)

# seleziona solo gli elementi con peso maggiore di 1
#dat=w[w[,4]>1,] 

# seleziona solo gli elementi con p_post maggiore di 0.5
if (prune==TRUE){
     dat=w[w[,5]>0.5,]
    }
else {
     dat=w
    }


#separa le prime tre colonne dall'ultima

#indici file A
A=cbind(unique(dat[,1]),1:length(unique(dat[,1]))) 
#indici file B
B=cbind(unique(dat[,2]),1:length(unique(dat[,2]))) 
colnames(A)=c('I','A')
colnames(B)=c('J','B')

#creazione del data set di analisi con i record indicizzati
dat=merge(B,dat)
dat=merge(A,dat)
dat=t(dat)
n=nrow(A); m=nrow(B) 

#crea la matrice dei vincoli di riga e colonna e la inizializza a zero
vincoli=array(rep(0,(n+m)*ncol(dat)), dim=c(n+m,ncol(dat)))

# aggiunge gli uno nella matrice dei vincoli di riga e colonna
p=rbind(matrix(rep(dat[2,],n),c(n,ncol(vincoli)),byrow=TRUE),
        matrix(rep(as.numeric(dat[4,])+n,m),c(m,ncol(vincoli)),byrow=TRUE))
		
#matrice dei vincoli
vincoli[as.numeric(p)==row(vincoli)]=1 
#matrice delle disuguaglianze
vindir=rep('<=',m+n)  
#upper bound dei vincoli       
uno=rep(1,m+n)  
#coefficienti della funzione obiettivo             
coeff=dat[6,]    
#solutore       
soluzione=lp ("max", coeff, vincoli, vindir, uno)   
 
#unisce la colonna soluzione ai dati
dat=rbind(dat,soluzione$solution)   
#seleziona gli abbinati
coppie=t(dat[c(1,3,5,6,7),dat[8,]>0.999])   
coppieNew <- as.data.frame(coppie)
coppieNew$I <- as.numeric(as.character(coppieNew$I))
coppieNew$J <- as.numeric(as.character(coppieNew$J))
#coppieNew$R <- as.numeric(as.character(coppieNew$R))
#coppieNew$P_POST <- as.numeric(as.character(coppieNew$P_POST))
#coppieNew
list <- coppieNew[,c("I", "J")]
colnames(list) <- c("row", "col")
out.style <- match.arg(out.style)
switch(out.style,"list"=list,"vec"=list.to.sol(list,D.mat))
}
