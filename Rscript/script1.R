library(MASS)
library(ggplot2)
library(dplyr)


norm_data <- function(data) {
  ## transform to column matrix if vector
  ## ensures that functions work also with one data point only
  if (is.null(dim(data))) {
    data <- t(matrix(data))
  }
  if(class(data)=="data.frame"){
    data <- as.matrix(data)
  }
  ## check presence of 0s or negative values
  if (any(data <= 0)) {
    stop("All values should be positive")
  }
  data
}


Inner_product<-function(x,y){
  
  x<-norm_data(x)  
  y<-norm_data(y)  

  if(ncol(x)!=ncol(y)){
    stop("x and y should have the same dimension")
  }
  if(nrow(x)!=1 && nrow(y)!=1){
    stop("wrong dimension")
  }
  D<-ncol(x)
  somme<-0
  for(i in 1:D){
    somme=somme+sum(log(x[i]/x)*log(y[i]/y))  
  }
  somme/(2*D)
}


norm_simplex<-function(x){
  
  sqrt(Inner_product(x,x))
}



closure <- function(data, k=1){
  data <- norm_data(data)
  c=rowSums(data)
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
}

mult<-function(x,t){
  
  return(x*t)
}

ternary_diagram<-function(data){
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



ligne<-function(x0,x){
  alpha=seq(-100,100,length=1000)
  perturbation(t(sapply(alpha,power,data=x)),x0)
}

simu_simplexe<-function(D,k,N){
  simplexe=matrix(runif(D*N,max=k),N,D)
  simplexe=k*apply(simplexe,2,function(x,j){x/j},j=rowSums(simplexe))
  simplexe
}



clr_inverse<-function(data,k){
  e=exp(data)
  closure(e,k)
}

Base_SIGMA_matrix<-function(D){
  mat=matrix(0,nrow=(D-1),ncol=D)
  for(i in 1:(D-1)){
    for(j in 1:(D-i)){
      mat[D-i,j]=sqrt(1/((D-i)*(D-i+1)))
    }
    mat[D-i,(D-i+1)]=-sqrt((D-i)/(D-i+1))
  }
  mat
}


Base_binary_matrix<-function(D){
  mat=matrix(0,nrow=(D-1),ncol=D)
  for(i in 1:(D-1)){
    for(j in 1:(D-i)){
      mat[D-i,j]=1
    }
    mat[D-i,(D-i+1)]=-1
  }
  mat
}

balance_coordinate=function(data,sequential_binary){
  D=dim(sequential_binary)[1]
  n=dim(data)[1]
  rs=cbind(apply(sequential_binary,1,function(x){return(sum(x==1))}),apply(sequential_binary,1,function(x){return(sum(x==-1))}))
  
  mat=matrix(0,nrow=n,ncol=D)
  
  for(i in 1:n){
    for(j in 1:D){
      mat[i,j]=sqrt(rs[j,1]*rs[j,2]/(rs[j,1]+rs[j,2]))*log((prod(data[i,which(sequential_binary[j,]==1)])^(1/rs[j,1]))/(prod(data[i,which(sequential_binary[j,]==-1)])^(1/rs[j,2])))
    }
  }
  mat
}

ilr<-function(data){
  clr(data)%*%t(Base_SIGMA_matrix(dim(data)[2])) 
}

ilr_inverse<-function(data,k){
  x=exp((data%*%Base_SIGMA_matrix(dim(data)[2]+1)))  
  
  closure(x,k)
}

variation_matrix<-function(data){
  D=dim(data)[2]
  mat=matrix(0,nrow=D,ncol=D)
  for(i in 1:D){
    for(j in 1:D){
      mat[i,j]=var(log(data[,i]/data[,j]))
    }
  }
  mat
}


cor_matrix<-function(data){
  
  if(nrow(data)==1){
    stop("correlation computation is impossible with one sample only.")
  }
  
  D <- ncol(data)
  mat <- matrix(0,nrow=D,ncol=D)
  for(i in 1:D){
    for(j in 1:D){
      mat[i,j] <- var(log(data[,i]/data[,j]))
    }
  }
  mat
}


normalised_variation_matrix <- function(data){
  log.data <- log(norm_data(data))
  ## test number of samples
  if (nrow(log.data) == 1) {
    stop("Variance computation is impossible with one sample only.")
  }
  D=ncol(log.data)
  mat=matrix(0, nrow = D, ncol = D)
  for(i in 1:D) {
    for(j in 1:i) {
      mat[i, j] <- mat[j, i] <- 0.5*var(log.data[, i] - log.data[, j])
    }
  }
  mat
}


center_data<-function(data){
  data <- norm_data(data)
  ## geometric mean = exponential of mean in log scale
  closure(exp(colMeans(log(data))))
}

totvar<-function(data){
  mean(rowSums(normalised_variation_matrix(data)))
}



center_scale<-function(data,center=TRUE,scale=TRUE){
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
  
  
  
  Z=clr(data)

  C=cov(Z)
  
  vect=eigen(C)$vector
  val=eigen(C)$values
  
  x=Z%*%vect
  
  n=dim(data)[1]
  
  stats::biplot(x[,1:2],vect[,1:2])
}


marginal_univariate_distributions<-function(data){
  
  data_ilr <- ilr(data)
  D_1 <- dim(data_ilr)[2]
  n <- dim(data_ilr)[1]
  
  u_i <- apply(data_ilr,2,mean)  
  var_i <- apply(data_ilr,2,var)
  
  q <- t(data_ilr)-u_i
  q <- t(q/var_i)
  p <- pnorm(q)
  
  z <- apply(p,2,sort)
  somme_a <- 0
  somme_c <- 0
  somme_w <- 0
  
  for(r in 1:n){
    somme_a <- somme_a+(2*r-1)*(log(z[r,])+log(1-z[n+1-r,]))
    somme_c <- somme_c+(z[r,]-(2*r-1)/(2*n))^2
    somme_w=somme_w+z[r,]-1./2
  }
  
  
  mat <- matrix(c("significance level:", " 10", " 5", " 2.5", " 1",
                  "Anderson Darling:",  0.656, 0.787, 0.918, 1.092,
                  "Cramer von Mises", 0.104, 0.126, 0.148, 0.178,
                  "Watson", 0.096, 0.116, 0.136, 0.163
                                                        ),nrow=4,byrow=TRUE)
  
  Qa <- (25./(n^2)-4./n-1)*(somme_a/n+n)
  Qc <- (somme_c+1/(12*n))*((2*n+1)/(2*n))
  Qw <- Qc-((2*n+1)/2.)*(somme_w/n)^2
  
  
  rval <- list(significance_level=mat, Anderson_Darling=Qa, Cramer_von_Mises=Qc, Watson=Qw)
  class(rval) <- "marginale univariate"
  
  rval
}





Bivariate_angle_distribution<-function(data){
  
  data <- norm_data(data)
  
  data_ilr <- ilr(data)
  D_1 <- dim(data_ilr)[2]
  n <- dim(data_ilr)[1]
  
  u_i <- apply(data_ilr,2,mean)
  var_i <- apply(data_ilr,2,var)
  
  Qa <- matrix(0, ncol = D_1, nrow = D_1)
  Qc <- matrix(0, ncol = D_1, nrow = D_1)
  Qw <- matrix(0, ncol = D_1, nrow = D_1)
  for(i in 1:D_1){
    j <- i+1
    while(j<=D_1){
     cov_ij <- cov(data_ilr[,i],data_ilr[,j]) 
     
     ur <- (1./sqrt(var_i[i]*var_i[j]-cov_ij^2))*((data_ilr[,i]-u_i[i])*sqrt(var_i[j])-(data_ilr[,j]-u_i[j])*cov_ij/sqrt(var_i[j]))
     vr <- (data_ilr[,j]-u_i[j])/sqrt(var_i[j])
     
     Or <- atan(vr/ur)+(1-sign(ur))*pi/2+(1+sign(ur))*(1-sign(vr))*pi/2
     
     z <- sort(Or/(2*pi))
     
     somme_a <- 0
     somme_c <- 0
     somme_w <- 0
     for(r in 1:n){
       somme_a <- somme_a+(2*r-1)*(log(z[r])+log(1-z[n+1-r]))
       somme_c <- somme_c+(z[r]-(2*r-1)/(2*n))
       somme_w <- somme_w+z[r]-1/2
         
     }
     
     Qa[i,j] <- -somme_a/n-n
     
     Qc[i,j] <- (somme_c-3.8/(12*n)+0.6/(n^2))*((n+1)/n)
     
     Qw[i,j] <- (somme_c-0.2/(12*n)+0.1/(n^2)-somme_w/n)*((n+0.8)/n)
     
     
     j=j+1
     }
  }
  
  mat <- matrix(c("significance level:", " 10", " 5", " 2.5", " 1",
                  "Anderson Darling:",  1.933, 2.492, 3.070, 3.857,
                  "Cramer von Mises", 0.347, 0.461, 0.581, 0.743,
                  "Watson", 0.152, 0.187, 0.221, 0.267
  ),nrow=4,byrow=TRUE)
  
  
  rval <- list(significance_level=mat, Anderson_Darling=Qa, Cramer_von_Mises=Qc, Watson=Qw)
  class(rval) <- "Bivariate angle test"
  
  rval
}


Raduis_test<-function(data){
  
  data <- norm_data(data)
  
  data_ilr <- ilr(data)
  D_1 <- ncol(data_ilr)
  n <- nrow(data_ilr)
  
  u_i <- apply(data_ilr,2,mean)
  var_ij <- var(data_ilr)
  var_ij_inv <- solve(var_ij)
  
  
  data_centrer <- t(t(data_ilr)-u_i)
  ur <- apply(data_centrer,1,function(x){return(x%*%var_ij_inv%*%x)})
  
  p <- pchisq(ur,df=D_1)
  z <- sort(p)
  
  
  somme_a <- 0
  somme_c <- 0
  somme_w <- 0
  for(r in 1:n){
    somme_a <- somme_a+(2*r-1)*(log(z[r])+log(1-z[n+1-r]))
    somme_c <- somme_c+(z[r]-(2*r-1)/(2*n))^2
    somme_w <- somme_w+z[r]-1/2
  }
  
  Qa <- -somme_a/n-n
  Qc <- (somme_c-3.8/(12*n)+0.6/(n^2))*((n+1)/n)
  Qw <- (somme_c-0.2/(12*n)+0.1/(n^2)-somme_w/n)*((n+0.8)/n)
  
  
  mat <- matrix(c("significance level:", " 10", " 5", " 2.5", " 1",
                  "Anderson Darling:",  1.933, 2.492, 3.070, 3.857,
                  "Cramer von Mises", 0.347, 0.461, 0.581, 0.743,
                  "Watson", 0.152, 0.187, 0.221, 0.267
  ),nrow=4,byrow=TRUE)
  
  
  rval <- list(significance_level=mat, Anderson_Darling=Qa, Cramer_von_Mises=Qc, Watson=Qw)
  class(rval) <- "raduis test"
  
  rval
 
}

intervalle_confiance<-function(data,alpha,case=1,moy=0,var_matrix=1){
  
  coord <- function(x,k,sigma,theta){
    u1 <- x[1]-sqrt((sigma[1,1]*sigma[2,2]-sigma[1,2]^2)/sigma[2,2])*sqrt(k)*cos(theta)-sigma[1,2]/sqrt(sigma[2,2])*sqrt(k)*sin(theta)
    u2 <- x[2]-sqrt(sigma[2,2])*sqrt(k)*sin(theta)
    
    matrix(c(u1,u2),ncol=2)
  }
  
  if(case==1){
    theta <- seq(0,2*pi,0.01)
    D <- ncol(data)
    u <- coord(moy,qchisq(1-alpha,D-1),var_matrix,theta)
    
    u
  }

  else if(case==2){
    moy <- apply(data,2,mean)
    var_matrix <- var(data)
    
    theta <- seq(0,2*pi,0.01)
    D <- dim(data)[2]
    n <- dim(data)[1]
    k <- (D-1)/(n-D+1)*qf(1-alpha,D-1,n-D+1)
    u <- coord(moy,k,var_matrix,theta)
    
    u
  }
  
}

testing<-function(data1,data2,alpha,case=1){
  
    
    u1 <- apply(data1,2,mean)
    u2 <- apply(data2,2,mean)
    
    sigma1 <-var(data1)
    sigma2 <- var(data2)

    n1 <- nrow(data1)
    n2 <- nrow(data2)
    D <- ncol(data1)+1
    
    if(case==1){
      uc <- (n1*u1+n2*u2)/(n1+n2)
      sigmac <- (n1*sigma1+n2*sigma2)/(n1+n2)+(n1*n2*(u1-u2)%*%t((u1-u2)))/((n1+n2)^2)
      
      Q <- n1*log(det(sigmac)/det(sigma1))+n2*log(det(sigmac)/det(sigma2))
      
      quant <- qchisq(1-alpha,0.5*D*(D-1))

    }
    else if(case==2){
      sigmap <- (n1*sigma1+n2*sigma2)/(n1+n2)
      
      Q <- n1*log(det(sigmap)/det(sigma1))+n2*log(det(sigmap)/det(sigma2))
      
      quant <- qchisq(1-alpha,0.5*(D-1)*(D-2))
    }
    
    else if(case==3){
      sigma1h <- sigma1
      sigma2h <- sigma2
      
      sigma_12h <- sigma1h+100
      sigma_22h <- sigma2h+100
      
      tol <- 1e-6
      nb_iter_max <- 5000000
      nb_iter <- 0
      
      while(abs(det(sigma1h)-det(sigma_12h))>tol | abs(det(sigma2h)-det(sigma_22h))>tol & nb_iter<nb_iter_max){
        
        sol_sigma1h <- solve(sigma1h)
        sol_sigma2h <- solve(sigma2h)
        
        uh <- solve(n1*sol_sigma1h+n2*sol_sigma2h)%*%(n1*sol_sigma1h%*%u1+n2*sol_sigma2h%*%u2)

        sigma_12h <- sigma1h
        sigma_22h <- sigma2h
        
        sigma1h <- sigma1+(u1-uh)%*%t(u1-uh)
        sigma2h <- sigma2+(u2-uh)%*%t(u2-uh)
        
        nb_iter <- nb_iter+1
      }
      
      Q <- n1*log(det(sigma1h)/det(sigma1))+n2*log(det(sigma2h)/det(sigma2))
      
      quant<- qchisq(1-alpha,(D-1))
    }
    
  rval<-list(statistic=Q, quantile=quant, case=case)
  class(rval)="test"
  
  rval
}



dendogram<-function(data,mat){
  
  data <- norm_data(data)
  n <- nrow(data)
  p <- ncol(data)
  
  
  data<-balance_coordinate(data,mat)
  
  
  
}



Graph_proportion_evolution<-function(data, abscisse=1:nrow(data)){
  
  data <- norm_data(data)
  time=abscisse
  D <- ncol(data)
  
  somme <- sum(data[1,])
  data <- data/somme
    
  data <- as.data.frame.table(data)
  data$Var1 <- rep(time, D)

  
  new.data <- group_by(data, Var2, Freq)
  
  g<- new.data %>% ggplot(aes(x=Var1, y=Freq))+
    geom_line(aes(group=Var2, color=Var2))+
    labs(y="proportions", x="alpha")+
    scale_color_discrete(name=NULL)+
    theme_bw()
  g
}





Graph_cumulative_evolution<-function(data, abscisse=1:nrow(data)){
  
  data <- norm_data(data)
  time=abscisse
  D <- ncol(data)
  
  data <- closure(data)
  data <- t(apply(data,1,cumsum))
  
  for(i in D:2){
    data[,i] <- data[,i]-data[,i-1]
  }
  
  data <- as.data.frame.table(data)
  data$Var1 <- rep(time, D)
  
  new.data <- group_by(data, Var2, Freq)
  
  g<- new.data %>% ggplot(aes(x=Var1, y=Freq, fill=Var2))+
    geom_area()+
    labs(y="proportions", x="alpha")+
    scale_fill_discrete(name=NULL)+
    theme_bw()
  g
  
}



simplexe=rbind(c(5,0,0),c(0,0,5),c(0,5,0))
ternary_diagram(simplexe)

test=simu_simplexe(3,10,40)

ternary_diagram(test)
ternary_diagram(perturbation(test,c(20,100,2)))
ternary_diagram(power(test,2))

x=matrix(c(1,3.8,1.2),nrow=1)
x0=matrix(c(1,4,1),nrow=1)
ternary_diagram(ligne(x0,x))

########### test ilr
x=matrix(c(1,2,-1.2,1,-6,9),nrow=3)
ternary_diagram(ilr_inverse(x,4))
x
ilr(ilr_inverse(x,4))

x=matrix(c(seq(-1,1,length=100),seq(-4,4,length=100)),ncol=2)
x1=matrix(c(seq(-4,4,length=100),seq(4,-4,length=100)),ncol=2)
ternary_diagram(ilr_inverse(x,4))
ternary_diagram(ilr_inverse(x1,4))



table_2.1=matrix(c(79.07,12.83,8.1,31.74,56.69,11.57,18.61,72.05,9.34,49.51,15.11,35.38,29.22,52.36,18.42,21.99,59.01,18.1,11.74,65.04,23.22,24.47,52.53,23.0,5.14,38.39,56.47,15.54,57.34,27.11),nrow=3)
clr(t(table_2.1))



x=matrix(c(48.29,48.83,45.61,45.5,49.27,46.53,48.12,47.93,46.96,49.16,48.41,47.9,48.45,48.98,48.74,49.61,49.2,
           2.33,2.47,1.7,1.54,3.3,1.99,2.34,2.32,2.01,2.73,2.47,2.24,2.35,2.48,2.44,3.03,2.5, 11.48,12.38,8.33,8.17,12.10,9.49,11.43,11.18,9.9,12.54,11.8,11.17,11.64,12.05,11.6,12.91,12.32,
           1.59,2.15,2.12,1.6,1.77,2.16,2.26,2.46,2.13,1.83,2.81,2.41,1.04,1.39,1.38,1.6,1.26,
           10.03,9.41,10.02,10.44,9.89,9.79,9.46,9.36,9.72,10.02,8.91,9.36,10.37,10.17,10.18,9.68,10.13, 0.18,0.17,0.17,0.17,0.17,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.17,0.18,
           13.58,11.08,23.06,23.87,10.46,19.28,13.65,14.33,18.31,10.05,12.52,14.64,13.23,11.18,12.35,8.84,10.51,
           9.85,10.64,6.98,6.79,9.65,8.18,9.87,9.64,8.58,10.55,10.18,9.58,10.13,10.83,10.45,10.96,11.05,
           1.9,2.02,1.33,1.28,2.25,1.54,1.89,1.86,1.58,2.09,1.93,1.82,1.89,1.73,1.67,2.24,2.02,
           0.44,0.47,0.32,0.31,0.65,0.38,0.46,0.45,0.37,0.56,0.48,0.41,0.45,0.8,0.79,0.55,0.48,
           0.23,0.24,0.16,0.15,0.3,0.18,0.22,0.21,0.19,0.26,0.23,0.21,0.23,0.24,0.23,0.27,0.23),
         nrow=17,
         dimnames = list(1:17,
                         paste0("v", 1:11)))
round(center_data(x),2)
totvar(x)
round(normalised_variation_matrix(x),3)

x_center_scale=center_scale(x,scale=FALSE)

x1=simu_simplexe(3,4,10)+matrix(c(rep(2,10),rep(0,20)),nrow=10)
x2=center_scale(x1)
ternary_diagram(x1)
ternary_diagram(x2)



biplot(x_center_scale)

base=matrix(c(0,0,0,1,-1,0,0,0,0,0,0,  1,0,-1,0,0,0,0,0,0,0,0,  0,1,0,0,0,0,0,0,0,0,-1,  1,-1,1,0,0,0,0,0,0,0,-1,   0,0,0,0,0,0,0,1,-1,0,0,  0,0,0,0,0,0,0,1,1,-1,0,   0,0,0,0,0,1,-1,0,0,0,0,   0,0,0,1,1,-1,-1,0,0,0,0,   0,0,0,1,1,1,1,-1,-1,-1,0,   1,1,1,-1,-1,-1,-1,-1,-1,-1,1),nrow=10,byrow = TRUE)

z=balance_coordinate(x,base)
m=matrix(0,nrow =10,ncol=10)
m[upper.tri(m)]=cor(z)[upper.tri(cor(z))]
m[lower.tri(m)]=var(z)[lower.tri(var(z))]
diag(m)=diag(var(z))
m


########################## simulation multivarie
simplexe=simu_simplexe(3,4,50000)
Sigma <- matrix(c(1,0.4,0.4,1),2,2)
Sigma1 <- matrix(c(1.2,0.8,0.8,1),2,2)
normale=mvrnorm(n = 5000, c(0,0.8), Sigma1)
normale1=mvrnorm(n = 4000, c(0,2), Sigma1)
ternary_diagram(ilr_inverse(rbind(normale,normale1),4))

marginal_univariate_distributions(simplexe)
marginal_univariate_distributions(ilr_inverse(normale,4))

Bivariate_angle_distribution(simplexe)
Bivariate_angle_distribution(ilr_inverse(normale,4))

Raduis_test(simplexe)
Raduis_test(ilr_inverse(normale,4))


u=intervalle_confiance(normale,0.05,1,c(0,0.8),Sigma1)
plot(normale)
lines(u,col='red')

u=intervalle_confiance(normale,0.05,2)
plot(normale)
lines(u,col='red')


normale=mvrnorm(n = 50000, c(0,0.8), Sigma1)
normale1=mvrnorm(n = 40000, c(0,0.81), Sigma1)
testing(normale,normale1,0.05,3)

