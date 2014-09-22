library(kernlab)

geigen <- function (A,B,C,top) {
  p       <- nrow(B) 
  q       <- nrow(C) 
  s       <- min(c(p,q))
  B       <- (B+t(B))/2
  C       <- (C+t(C))/2
  Bfac    <- chol(B)
  Cfac    <- chol(C)
  Bfacinv <- solve(Bfac)
  Cfacinv <- solve(Cfac)
  D       <- t(Bfacinv)%*%A%*%Cfacinv
  if (p >= q) {
      result <- svd(D,nu=top,nv=top)
      values <- result$d
      L      <- Bfacinv %*% result$u
      M      <- Cfacinv %*% result$v
  } else {
      result <- svd(t(D),nu=top,nv=top)
      values <- result$d
      L      <- Bfacinv %*% result$v
      M      <- Cfacinv %*% result$u
  }
  list(cor=values, xcoef=L, ycoef=M)
}

rcc <- function (X,Y,l1,l2,top) {
  geigen(cov(X,Y),var(X)+diag(l1,ncol(X)),
                  var(Y)+diag(l2,ncol(Y)),top)
}

aug <- function(x,k,type="fourier") {
  s <- sigest(x,scaled=NULL)[2]

  if(type == "linear") {
    return(function(x0) x0)
  }
  
  if(type == "nystrom") {
    w <- x[sample(1:nrow(x),k),]
    return(function(x0) kernelMatrix(rbfdot(s),x0,w))
  }
  
  if(type == "fourier") {
    w <- matrix(rnorm(ncol(x)*k,sd=sqrt(2*s)),ncol(x))
    b <- runif(k,0,2*pi)
    f <- function(x0) x0%*%w+t(matrix(b,k,nrow(x0)))
    return(function(x0) cos(f(x0)))
  }
}

rcca_fit <- function(x,y,kx,ky,type,top) {
  augx <- aug(x,kx,type)
  augy <- aug(y,ky,type)
  C    <- rcc(augx(x),augy(y),1e-10,1e-10,top)
  list(cor=sum(abs(C$cor[1:top])),a=C$xcoef,b=C$ycoef,augx=augx,augy=augy)
}

rcca_eval <- function(rcca,x,y){
  list(x=rcca$augx(x)%*%rcca$a,y=rcca$augy(y)%*%rcca$b)
}

rpca_fit <- function(x,k,type) {
  augx <- aug(x,k,type)
  list(augx=augx,pca=prcomp(augx(x)))
}

rpca_eval <- function(rpca,x) {
  predict(rpca$pca,rpca$augx(x))
}

