library(MASS)


regularisation <- function(var1){

  if(class(try(solve(var1), silent=TRUE))=="matrix"){
    var1
  }else{
    D <- ncol(var1)
    ID <- diag(1, nrow=D)

    l <- 1e-24
    while(!(class(try(solve(var1+l*ID), silent=TRUE))=="matrix")){
      l <- l*10
    }

    var1 <- var1+l*ID
    var1
  }
}

solve_regularisation <- function(var1) {
  var1 %>% regularisation() %>% solve()
}


#' realize three univariate normality test (Anderson-Darling statistic, Cramer-von Mises statistic, Watson statistic) 
#' on each variable of compositionnal dataset
#'
#' @param data a compositionnal dataset (composition on line)
#'
#' @return significance_level: critical values for the test
#' @return Anderson_Darling statistic
#' @return Cramer-von Mises statistic
#' @return Watson statistic
#' @return marginale_normale the results for the watson test (confidence 5%)
#' @author Clement Hardy
#' @export
#' @import stats

marginal_univariate_distributions<-function(data){

  data_ilr <- ilr(data)
  D_1 <- dim(data_ilr)[2]
  n <- dim(data_ilr)[1]

  u_i <- apply(data_ilr,2,mean)
  var_i <- apply(data_ilr,2,var)

  q <- t(data_ilr)-u_i
  q <- t(q/sqrt(var_i))
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


  rval <- list(significance_level=mat, Anderson_Darling=Qa, Cramer_von_Mises=Qc, Watson=Qw, marginale_normale=which(Qw<0.116))
  class(rval) <- "marginale univariate"

  rval
}



#' realize three bivariate normality test (Anderson-Darling statistic, Cramer-von Mises statistic, Watson statistic) 
#' on each couple of variable of a compositionnal dataset (via a raduis angle test)
#'
#' @param data a compositionnal dataset (composition on line)
#'
#' @return significance_level: critical values for the test
#' @return Anderson_Darling statistic
#' @return Cramer-von Mises statistic
#' @return Watson statistic
#' @return couple_normal the results for the watson test, the couple of variables following a normal (confidence 5%)
#' @author Clement Hardy
#' @export
#' @import stats

Bivariate_angle_distribution<-function(data){

  data <- norm_data(data)

  data_ilr <- ilr(data)
  D_1 <- dim(data_ilr)[2]
  n <- dim(data_ilr)[1]

  u_i <- apply(data_ilr,2,mean)
  var_i <- apply(data_ilr,2,var)

  Qa <- matrix(NA, ncol = D_1, nrow = D_1)
  Qc <- matrix(NA, ncol = D_1, nrow = D_1)
  Qw <- matrix(NA, ncol = D_1, nrow = D_1)
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


  rval <- list(significance_level=mat, Anderson_Darling=Qa, Cramer_von_Mises=Qc, Watson=Qw, couple_normal=which(Qw<0.187, arr.ind=TRUE))
  class(rval) <- "Bivariate angle test"

  rval
}


#' realize three bivariate normality test (Anderson-Darling statistic, Cramer-von Mises statistic, Watson statistic) 
#'
#' @param data a compositionnal dataset (composition on line)
#'
#' @return significance_level: critical values for the test
#' @return Anderson_Darling statistic
#' @return Cramer-von Mises statistic
#' @return Watson statistic
#' 
#' @author Clement Hardy
#' @export
#' @import stats


Raduis_test<-function(data){

  data <- norm_data(data)

  data_ilr <- ilr(data)
  D_1 <- ncol(data_ilr)
  n <- nrow(data_ilr)

  u_i <- apply(data_ilr,2,mean)
  var_ij <- var(data_ilr)
  var_ij_inv <- solve_regularisation(var_ij)


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


#' give an confidence ellipse  for 2D compositionnal data
#' 
#' @param data les donnees compositionnelles
#' @param alpha risk level
#' @param case the type of ellipse givent (1 mean, variance unknown, probabilty ellipse)
#' (2 mean unknown, variance unknown, confidence ellipse for the mean)
#' (2 mean unknown, variance unknown, probability ellipse)
#' @param moy the mean (if known)
#' @param var_matrix the variance (if known)
#' @return significance level
#' @author Clement Hardy
#' @export
#' @import stats


intervalle_confiance<-function(data, alpha=0.05, case=3, moy=NULL, var_matrix=NULL){

  coord <- function(x,k,sigma,theta){
    u1 <- x[1]-sqrt((sigma[1,1]*sigma[2,2]-sigma[1,2]^2)/sigma[2,2])*sqrt(k)*cos(theta)-sigma[1,2]/sqrt(sigma[2,2])*sqrt(k)*sin(theta)
    u2 <- x[2]-sqrt(sigma[2,2])*sqrt(k)*sin(theta)

    return(matrix(c(u1,u2),ncol=2))
  }

  if(case==1){
    theta <- seq(0,2*pi,0.05)
    D <- ncol(data)
    u <- coord(moy,qchisq(1-alpha,D-1),var_matrix,theta)

    return(u)
  }

  else if(case==2){
    moy <- apply(data,2,mean)
    var_matrix <- var(data)

    theta <- seq(0,2*pi,0.01)
    D <- ncol(data)
    n <- nrow(data)
    k <- (D-1)/(n-D+1)*qf(1-alpha,D-1,n-D+1)
    u <- coord(moy,k,var_matrix,theta)
    return(u)

  }else if(case==3){

    moy <- apply(data,2,mean)
    var_matrix <- var(data)

    theta <- seq(0,2*pi,0.01)
    D <- ncol(data)
    n <- nrow(data)
    k <- (n-1)*(D-1)/(n-D+1)*((n+1)/n)*qf(1-alpha,D-1,n-D+1)
    u <- coord(moy,k,var_matrix,theta)

    return(u)

  }

}


#' test the equality of mean or covariance matrix of tow group of compositionnal dataset
#' 
#' @param data1 the first dataset
#' @param data2 the second dataset
#' @param alpha the risk level
#' @param case type of test (1 equal mean and equal covariance matrix)
#' (2 equal covariance matrix)
#' (3 equal mean)
#' 
#' @return statistic=Q, quantile=quant, result=(Q<quant), case=case
#' @author Clement Hardy
#' @export
#' @import stats



testing<-function(data1, data2, alpha=0.05, case=1){
  
    data1 <- data1 %>% norm_data() %>% ilr()
    data2 <- data2 %>% norm_data() %>% ilr()

    u1 <- apply(data1,2,mean)
    u2 <- apply(data2,2,mean)

    sigma1 <- data1 %>% var() %>% regularisation()
    sigma2 <- data2 %>% var() %>% regularisation()

    n1 <- nrow(data1)
    n2 <- nrow(data2)
    D <- ncol(data1)+1

    if(case==1){
      uc <- (n1*u1+n2*u2)/(n1+n2)
      sigmac <- (n1*sigma1+n2*sigma2)/(n1+n2)+(n1*n2*(u1-u2)%*%t((u1-u2)))/((n1+n2)^2) %>% regularisation()

      Q <- n1*log(det(sigmac)/det(sigma1))+n2*log(det(sigmac)/det(sigma2))
      quant <- qchisq(1-alpha,0.5*D*(D-1))

    }
    else if(case==2){
      sigmap <- (n1*sigma1+n2*sigma2)/(n1+n2) %>% regularisation()

      Q <- n1*log(det(sigmap)/det(sigma1))+n2*log(det(sigmap)/det(sigma2))

      quant <- qchisq(1-alpha,0.5*(D-1)*(D-2))
    }

    else if(case==3){
      sigma1h <- sigma1
      sigma2h <- sigma2

      sigma_12h <- sigma1h+1
      sigma_22h <- sigma2h+1


      tol <- 1e-12
      nb_iter_max <- 500
      nb_iter <- 0

      while(abs(det(sigma1h)-det(sigma_12h))>tol | abs(det(sigma2h)-det(sigma_22h))>tol & nb_iter<nb_iter_max){

        sol_sigma1h <- solve_regularisation(sigma1h)
        sol_sigma2h <- solve_regularisation(sigma2h)

        uh <- solve_regularisation(n1*sol_sigma1h+n2*sol_sigma2h)%*%(n1*sol_sigma1h%*%u1+n2*sol_sigma2h%*%u2)

        sigma_12h <- sigma1h
        sigma_22h <- sigma2h

        sigma1h <- (sigma1+(u1-uh)%*%t(u1-uh)) %>% regularisation()
        sigma2h <- (sigma2+(u2-uh)%*%t(u2-uh)) %>% regularisation()

        nb_iter <- nb_iter+1
      }

      Q <- n1*log(det(sigma1h)/det(sigma1))+n2*log(det(sigma2h)/det(sigma2))
      quant<- qchisq(1-alpha,(D-1))
    }

  rval<-list(statistic=Q, quantile=quant, result=(Q<quant), case=case)
  class(rval)="test"

  rval
}



