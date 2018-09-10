library(rARPACK)

#' check if the data has NA, transform data to matrix
#'
#' @param data data to check
#'
#' @return the data transform
#'
#' @author Clement Hardy
#' @export



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

#' check if the data has NA, zero, transform data to matrix
#'
#' @param data data to check
#'
#' @return the data transform
#'
#' @author Clement Hardy
#' @export



norm_data <- function(data) {
  data <- check_data(data)
  ## check presence of 0s or negative values
  if (any(data <= 0)) {
    warning("Non positive values are present in the data.")
  }
  data
}

#' compute the Inner product of the compositional data
#'
#' @param x,y two compositional data
#'
#' @return the Inner product of the data
#' @examples 
#' \dontrun{
#' x <- matrix(c(0.1, 0.3, 0.6), nrow=1)
#' y <- matrix(c(0.4, 0.2, 0.4), nrow=1)
#' Inner_product(x,y)
#' }
#' @author Clement Hardy
#' @export
#' 
#' 
Inner_product<-function(x,y){
  x<-clr(x)
  y<-clr(y)
  if (ncol(x) != ncol(y)){
    stop("x and y should have the same dimension")
  }
  tcrossprod(x, y) ## shortcut for x %*% t(y)
}


#' compute the norm of a compositionnal data
#'
#' @param data the compositionnal data
#' 
#' @return a single number giving the norm of the compositionnal data, 
#' or a numeric vector giving the norm of each line if data is a matrix
#'
#' @examples 
#' \dontrun{
#' data <- closure(matrix(runif(30,0,1), nrow=10))
#' norm_simplex(data)
#' }
#' 
#' @author Clement Hardy
#' @export


norm_simplex<-function(data){
  data %>% clr() %>% .^2 %>% rowSums() %>% sqrt()
}


#' compute the distance in the simplex
#'
#' @param a numeric matrix containing compositionnal data
#' 
#' @return an object of class "dist"
#' @examples 
#' \dontrun{
#' data <- closure(matrix(runif(30,0,1), nrow=10))
#' norm_simplex(data)
#' }
#' @author Clement Hardy
#' @export

dist_simplex <- function(data) {
  data %>% clr %>% dist
}


#' close compositional data for a sum equal to a constant k
#'
#' @param data a matrix containing compositional data
#' @param k the sum of the closure (the default is 1)
#' @return compositionnal data (closure done)
#' @examples 
#' \dontrun{
#' x <- matrix(1:6, nrow=2)
#' closure(x)
#' closure(x, k=3)
#' }
#' @author Clement Hardy
#' @export

closure <- function(data, k=1){
  data <- norm_data(data)
  c <- rowSums(data)
  k * data / c
}


#' compute the perturbation between compositionnal data 
#'
#' @param data a matrix containing compositional data
#' @param multiple a compositional data use to perturbe data
#' 
#' @return the results of the perturbation 
#' @examples 
#' \dontrun{
#' data <- closure(matrix(1:9, nrow=3))
#' x <- closure(matrix(c(1,1.2,1.2), nrow=1))
#' perturbation(data, x)
#' }
#' @author Clement Hardy
#' @export
#' 

perturbation <- function(data, multiple){
  
  data <- norm_data(data)
  multiple <- norm_data(multiple)
  ## multiply each row in data by values in multiple
  new.data <- sweep(data, 2, multiple, "*")
  closure(new.data)
}

#' compute the power between compositionnal data and a constant
#'
#' @param data a matrix containing compositional data
#' @param alpha a numeric to power with
#' 
#' @return the results of the power  
#' @examples 
#' \dontrun{
#' data <- closure(matrix(1:9, nrow=3))
#' x <- 2
#' power(data, x)
#' }
#' @author Clement Hardy
#' @export
#' 

power <- function(data, alpha){
  data <- norm_data(data)
  new.data <- data^alpha
  closure(new.data)
}


#' compute the clr transformation
#'
#' @param data a matrix containing compositional data
#' 
#' @return the centerd log ration coordinates of the data
#' @examples 
#' \dontrun{
#' data <- closure(matrix(1:9, nrow=3))
#' clr(data)
#' }
#' @author Clement Hardy
#' @export
#' 

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

ternary_diagram<-function(data, colour="black", type_point=1, ...){
  
  data <- norm_data(data)
  u0=0.2;
  v0=0.2;
  A=c(u0+0.5,v0+sqrt(3)/2)
  B=c(u0,v0)
  C=c(u0+1,v0)
  
  k=sum(data[1,])
  coord=(1./k)*(sapply(data[,1],mult,t=A)+sapply(data[,2],mult,t=B)+sapply(data[,3],mult,t=C))
  plot(coord[1,],coord[2,], ylim=c(v0-0.03,v0+sqrt(3)/2+0.11), col=colour, pch=type_point, axes=FALSE, asp=1, xlab="", ylab="", ...)
  segments(B[1],B[2],A[1],A[2])
  segments(C[1],C[2],A[1],A[2])
  segments(B[1],B[2],C[1],C[2])
  text(A[1], A[2], "A", pos=3)
  text(B[1], B[2], "B", pos=2)
  text(C[1], C[2], "C", pos=4)
  
}

ligne<-function(x0, x) {
  alpha=seq(-100, 100, length = 1000)
  perturbation(t(sapply(alpha, power, data=x)), x0)
}


simu_simplexe<-function(D, k = 1, N) {
  simplexe <- matrix(runif(D*N, max=k), N, D)
  closure(simplexe, k)
}

#' compute the inverse transformation of the centred log ratio transformation
#'
#' @param data a matrix containing compositional data
#' 
#' @return a mtrix containing the coordinate of data in the simplex
#' @examples 
#' \dontrun{
#' data <- closure(matrix(1:9, nrow=3))
#' data1 <- clr(data)
#' clr_inverse(data1)
#' }
#' @author Clement Hardy
#' @export

clr_inverse<-function(data, k = 1){
  exp(data) %>% closure(k)
}


#' compute the SIGMA matrix of the basis (for the isometric log ration transformation)
#'
#' @param data D the size of compositional data
#' 
#' @return the  Sigma matrix of the basis (dimension of the matrix D-1,D)
#' @examples 
#' \dontrun{
#' Base_SIGMA_matrix(2)
#' Base_SIGMA_matrix(5)
#' }
#' @author Clement Hardy
#' @export


Base_SIGMA_matrix <- function(D) {
  build_row <- function(i) {
    denom <- sqrt(i * (i+1))
    c(rep(-1 / denom, i),
      i / denom,
      rep(0, D - (i+1)))
    
  }
  sapply(1:(D-1), build_row) %>% t()
}


#' compute the SIGMA matrix of a particular basis (for the isometric log ration transformation)
#'
#' @param sequential_binary a matrix of the sequence of the particular basis
#' 
#' @param data D the size of compositional data
#' 
#' @return the  binary matrix of the basis
#' @examples 
#' \dontrun{
#' mat <- matrix(c(1,1,-1,1,-1,0),byrow = TRUE, nrow=2)
#' SIGMA_matrix(mat)
#' }
#' @author Clement Hardy
#' @export

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


#' compute the binary matrix of the basis (for the isometric log ratio transformation)
#'
#' @param D a matrix of the sequence of the a particular basis
#' 
#' @return the matrix of the particular basis
#' @examples 
#' \dontrun{
#' Base_binary_matrix(2)
#' Base_binary_matrix(5)
#' }
#' @author Clement Hardy
#' @export

Base_binary_matrix<-function(D){
  
  build_row <- function(i) {
    rep(c(1, -1, 0), times = c(i, 1, D - (i+1) ))
    
  }
  -sapply(1:(D-1), build_row) %>% t()
  
}

#' compute the isometrice log ratio (ilr transformation)
#'
#' @param data a matrix containing compositional data
#' 
#' @return the isometric log ratio  coordinates of the data
#' @examples 
#' \dontrun{
#' mat <- matrix(c(1,1,-1,1,-1,0),byrow = TRUE, nrow=2)
#' data <- closure(matrix(1:9, nrow=3))
#' ilr(data)
#' ilr(data, mat)
#' }
#' @author Clement Hardy
#' @export
#' 


ilr<-function(data, base_binaire=Base_binary_matrix(ncol(data))){
  SIGMA=SIGMA_matrix(base_binaire)
  data <- norm_data(data)
  clr(data) %*% t(SIGMA)
  
}


#' compute the inverse transformation of the isometric log ratio transformation
#'
#' @param data a matrix containing compositional data
#' 
#' @return a mtrix containing the coordinate of data in the simplex
#' @examples 
#' \dontrun{
#' mat <- matrix(c(1,1,-1,1,-1,0),byrow = TRUE, nrow=2)
#' data <- closure(matrix(1:9, nrow=3))
#' ilr(data) %>% ilr_inverse()
#' ilr(data, mat) %>% ilr_inverse(base_binaire=mat)
#' }
#' @author Clement Hardy
#' @export

ilr_inverse<-function(data, k=1, base_binaire=Base_binary_matrix(ncol(data)+1)){
  
  data <- check_data(data)
  SIGMA=SIGMA_matrix(base_binaire)
  exp(data %*% SIGMA) %>% closure(k)
  
}

#' compute the additive log ratio (alr transformation)
#'
#' @param data a matrix containing compositional data
#' 
#' @return the additive log ratio  coordinates of the data
#' @examples 
#' \dontrun{
#' data <- closure(matrix(1:9, nrow=3))
#' alr(data)
#' }
#' @author Clement Hardy
#' @export

alr <- function(data){
  data <- norm_data(data)
  data <- (data[, -ncol(data)] / data[, ncol(data)]) %>% log()
  data
}


#' compute the  inverse transformation of the additive log ratio transformation
#'
#' @param data a matrix containing compositional data
#' 
#' @return a mtrix containing the coordinate of data in the simplex
#' @examples 
#' \dontrun{
#' data <- closure(matrix(1:9, nrow=3))
#' alr(data) %>% alr_inverse()
#' }
#' @author Clement Hardy
#' @export

alr_inverse <- function(data){
  ((1/sum(exp(data)) * cbind(exp(data),1))) %>% closure()
}


#' compute the covariance matrix of compositional data
#'
#' @param data a matrix containing compositional data
#' @param norm a boolean, TRUE for normalized covariance matrix, FALSE otherwise
#' 
#' @return the covariance matrix
#' @examples 
#' \dontrun{
#' data <- matrix(runif(99, 0, 1), ncol=3)
#' variation_matrix(data)
#' variation_matrix(data, norm=TRUE)
#' }
#' @author Clement Hardy
#' @export
#' 

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

#' compute the normalised covariance matrix of compositional data
#'
#' @param data a matrix containing compositional data
#' @param norm a boolean, TRUE for normalized covariance matrix, FALSE otherwise
#' 
#' @return the normalised covariance matrix
#' @examples 
#' \dontrun{
#' data <- matrix(runif(99, 0, 1), ncol=3)
#' normalised_variation_matrix(data)
#' }
#' @author Clement Hardy
#' @export
#' 

normalised_variation_matrix <- function(data){
  variation_matrix(data, norm = TRUE)
}


#' compute the center of the compositional data
#'
#' @param data a matrix containing compositional data
#' 
#' @return a composition (center of the data)
#' @examples 
#' \dontrun{
#' data <- closure(matrix(1:9, nrow=3))
#' center_data(data)
#' }
#' @author Clement Hardy
#' @export
#' 

center_data<-function(data) {
  
  ## geometric mean = exponential of mean in log scale
  data %>% norm_data() %>% log() %>% colMeans() %>% exp() %>% closure()
}


#' compute the total variance of a dataset of compositions
#'
#' @param data a matrix containing compositions data
#' 
#' @return a numeric (total variance of the dataset)
#' @examples 
#' \dontrun{
#' data <- matrix(runif(99, 0, 1), ncol=3)
#' totvar(data)
#' }
#' @author Clement Hardy
#' @export
#' 

totvar<-function(data){
  data <- norm_data(data)
  data %>% variation_matrix(norm = TRUE) %>% rowSums() %>% mean()
}

#' center and scale a dataset of compositional data
#'
#' @param data a matrix containing compositions data
#' @param center a boolean, center the dataset if TRUE
#' @param scale a boolean, scale the dataset if TRUE
#' 
#' @return the data center (if TRUE), scale (if TRUE)
#' @examples 
#' \dontrun{
#' data <- cbind(runif(100,0,1), runif(100,0,100), runif(100,20,80)) %>% closure()
#' center_scale(data, scale=FALSE)
#' center_scale(data)
#' }
#' @author Clement Hardy
#' @export
#' 

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


#' principal component analysis on compositionnal data
#'
#' @param data the dataset 
#' @param base_binaire binary sequential matrix use for the ilr transformation
#' 
#' @return variance_explain the cumulative sum of variance 
#' @return vector the eigen vector of the covariance matrix
#' @return values the eigen values of the covariance matrix
#' @return coord the coordinate of the data in the latent space
#' @examples 
#' \dontrun{
#' data <- cbind(runif(100, 0, 1), runif(100, 0, 100), runif(100, 20, 80), runif(100, 10, 11), runif(100, 23, 83)) %>% closure()
#' biplot(data)
#' }
#' @author Clement Hardy
#' @export
#' @import rARPACK
#'

biplot<-function(data, base_binaire=Base_binary_matrix(ncol(data))){
  
  data <- norm_data(data)
  n <- nrow(data)
  p <- ncol(data)
  
  data <- center_scale(data, scale=FALSE)
  
  
  
  Z <- ilr(data, base_binaire = base_binaire)

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


#' divide each line of a positive or null dataset set (no composante less than 0) by the sum of the line
#'
#' @param data a matrix containing the dataset
#' 
#' @return a matrix containing the proportion of each componant (for each line)
#' @examples 
#' \dontrun{
#' data <- matrix(1:9, nrow=3)
#' count_to_proportion(data)
#' }
#' @author Clement Hardy
#' @export
#' 


count_to_proportion<-function(data){
  
  data <- check_data(data)
  new.data <- data/rowSums(data)
  new.data
}

#' compute the maximun a posteriori estimator of a dataset of compositional data
#'
#' @param data a matrix containing compositional data
#' @param alpha the regularization parameter (default is 1)
#' 
#' @return the maximun a posteriori
#' @examples 
#' \dontrun{
#' data <- matrix(1:9, nrow=3)
#' MAP(data)
#' MAP(data, alpha=2)
#' }
#' @author Clement Hardy
#' @export

MAP <- function(data, alpha=1){
  
  data <- check_data(data+alpha)
  new.data <- data/rowSums(data)
  new.data
  
}