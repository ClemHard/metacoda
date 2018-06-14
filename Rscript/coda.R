library(rARPACK)


check_data <- function(data) {
  ## transform to column matrix if vector
  ## ensures that functions work also with one data point only
  if (is.null(dim(data))) {
    data <- t(matrix(data))
  }
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  if (any(is.na(data))) {
    warning("data with NA has been deleted")
    data <- as.matrix(data[complete.cases(data), ])
    if (nrow(data) == 0)
      return(NULL)
  }
  
  data
}

norm_data <- function(data) {
  data <- check_data(data)
  ## check presence of 0s or negative values
  if (any(data <= 0)) {
    warning("Non positive values are present in the data.")
  }
  data
}

Inner_product<-function(x,y){
  x<-clr(x)
  y<-clr(y)
  if (ncol(x) != ncol(y)){
    stop("x and y should have the same dimension")
  }
  tcrossprod(x, y) ## shortcut for x %*% t(y)
}

norm_simplex<-function(data){
  data %>% clr() %>% .^2 %>% rowSums() %>% sqrt()
}

dist_simplex <- function(data) {
  data %>% clr %>% dist
}

closure <- function(data, k=1){
  data <- norm_data(data)
  c <- rowSums(data)
  k * data / c
}

perturbation <- function(data, multiple){
  data <- norm_data(data)
  ## multiply each row in data by values in multiple
  new.data <- sweep(data, 2, multiple, "*")
  closure(new.data)
}

power <- function(data, alpha){
  data <- norm_data(data)
  new.data <- data^alpha
  closure(new.data)
}

clr <- function(data){
  data <- norm_data(data)
  ## Log transform ...
  log.data <- log(data)
  ## ... and center
  clr.data <- log.data - rowMeans(log.data)
  
  clr.data
}

mult<-function(x,t){
  x * t
}

ternary_diagram<-function(data){
  
  data <- norm_data(data)
  u0=0.2;
  v0=0.2;
  A=c(u0+0.5,v0+sqrt(3)/2)
  B=c(u0,v0)
  C=c(u0+1,v0)
  
  k=sum(data[1,])
  coord=(1./k)*(sapply(data[,1],mult,t=A)+sapply(data[,2],mult,t=B)+sapply(data[,3],mult,t=C))
  plot(coord[1,],coord[2,],xlim=c(u0-0.2,u0+1.2),ylim=c(v0-0.2,v0+1.1))
  segments(B[1],B[2],A[1],A[2])
  segments(C[1],C[2],A[1],A[2])
  segments(B[1],B[2],C[1],C[2])
  
}

ligne<-function(x0, x) {
  alpha=seq(-100, 100, length = 1000)
  perturbation(t(sapply(alpha, power, data=x)), x0)
}


simu_simplexe<-function(D, k = 1, N) {
  simplexe <- matrix(runif(D*N, max=k), N, D)
  closure(simplexe, k)
}

clr_inverse<-function(data, k = 1){
  exp(data) %>% closure(k)
}

Base_SIGMA_matrix <- function(D) {
  build_row <- function(i) {
    denom <- sqrt(i * (i+1))
    c(rep(-1 / denom, i),
      i / denom,
      rep(0, D - (i+1)))
    
  }
  sapply(1:(D-1), build_row) %>% t()
}


SIGMA_matrix <- function(sequential_binary){
  
  
  D <- ncol(sequential_binary)
  n <- nrow(sequential_binary)
  rs <- cbind(apply(sequential_binary,1,function(x){return(sum(x==1))}),apply(sequential_binary,1,function(x){return(sum(x==-1))}))
  
  mat <- matrix(0,nrow=n,ncol=D)
  
  for(i in 1:n){
    for(j in 1:D){
      if(sign(sequential_binary[i,j])==1){
        mat[i,j] <- sqrt(rs[i,1]*rs[i,2]/(rs[i,1]+rs[i,2]))/rs[i,1]
      }else if(sign(sequential_binary[i,j])==-1){
        mat[i,j] <- -sqrt(rs[i,1]*rs[i,2]/(rs[i,1]+rs[i,2]))/rs[i,2]
      }
    }
  }
  mat
}

Base_binary_matrix<-function(D){
  
  build_row <- function(i) {
    rep(c(1, -1, 0), times = c(i, 1, D - (i+1) ))
    
  }
  sapply(1:(D-1), build_row) %>% t()
  
}

ilr<-function(data, SIGMA=Base_SIGMA_matrix(ncol(data))){
  data <- norm_data(data)
  clr(data) %*% t(SIGMA)
  
}

ilr_inverse<-function(data, k=1, SIGMA=Base_SIGMA_matrix(ncol(data)+1)){
  exp(data %*% SIGMA) %>% closure(k)
  
}



variation_matrix<-function(data, norm = FALSE) {
  log.data <- log(norm_data(data))
  
  if(nrow(log.data)==1){
    stop("correlation computation is impossible with one sample only.")
  }
  
  D <- ncol(log.data)
  
  mat <- matrix(0, nrow=D, ncol=D)
  for(i in 1:D) {
    for(j in 1:i) {
      mat[i, j] <- mat[j, i] <- var(log.data[, i] - log.data[, j])
      
    }
  }
  if (norm) {
    mat <- .5 * mat
  }
  mat
}

normalised_variation_matrix <- function(data){
  variation_matrix(data, norm = TRUE)
}


center_data<-function(data) {
  
  ## geometric mean = exponential of mean in log scale
  data %>% norm_data() %>% log() %>% colMeans() %>% exp() %>% closure()
}

totvar<-function(data){
  data %>% variation_matrix(norm = TRUE) %>% rowSums() %>% mean()
}

center_scale<-function(data, center = TRUE, scale = TRUE) {
  data <- norm_data(data)
  
  texte <- dimnames(data)
  if (center) {
    data <- perturbation(data, 1. / center_data(data))
  }
  if (scale) {
    data <- power(data, 1. / sqrt(totvar(data)))
  }
  dimnames(data) <- texte
  data
}


biplot<-function(data){
  
  data <- norm_data(data)
  n <- nrow(data)
  p <- ncol(data)
  
  data <- center_scale(data, scale=FALSE)
  
  
  
  Z <- ilr(data)
  
  C <- cov(Z)
  
  if(p>n){
    C <- C %>% eigs_sym(k=min(n,p)-1)
  }else{
    C <- C %>% eigen()
  }
  
  vect <- C$vectors
  val <- C$values
  
  x <- Z%*%vect
  
  rval <- list(variance_explain=cumsum(val)/sum(val), vector=vect, values=val, coord=x)
  
  class(rval) <- "biplot"
  
  rval
}


count_to_proportion<-function(data){
  
  data <- check_data(data)
  new.data <- data/rowSums(data)
  new.data
}


MAP <- function(data, alpha=1){
  
  data <- check_data(data+alpha)
  new.data <- data/rowSums(data)
  new.data
  
}