source("scoring/operators.R")
set.seed(20261113)
x <- matrix(rnorm(400),40,10,dimnames=list(paste0("g",1:40),paste0("s",1:10)))
g <- rownames(x)[1:5]
stopifnot(max(abs(O1(x,g)$scores-as.numeric(scale(colMeans(x[g,])))))<1e-12)
z <- t(scale(t(x[g,])))
stopifnot(max(abs(O2(x,g)$scores-as.numeric(scale(colMeans(z)))))<1e-12)
for(op in c("O1","O2","O3","O7")){a<-get(op)(x,g,1:7);b<-get(op)(x[,8:10,drop=FALSE],g,reference=a$reference);stopifnot(max(abs(a$scores[8:10]-b$scores))<1e-12)}
stopifnot(!any(O1(x,c(g,"missing"))$available))
cat("PASS: synthetic mean, gene-standardized, frozen reference and missing gene checks\n")
