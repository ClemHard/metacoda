library(MASS)
library(ggplot2)
library(dplyr)


check_data <- function(data) {
  ## transform to column matrix if vector
  ## ensures that functions work also with one data point only
  if (is.null(dim(data))) {
    data <- t(matrix(data))
  }
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  if(any(is.na(data))){
    warning("NA in data has been deleted")
    data <- as.matrix(data[complete.cases(data),])
    
    if(nrow(data)==0) data <- NULL
  }
  
  data
}


norm_data <- function(data) {
  ## transform to column matrix if vector
  ## ensures that functions work also with one data point only
  if (is.null(dim(data))) {
    data <- t(matrix(data))
  }
  if(is.data.frame(data)){
    data <- as.matrix(data)
  }
  ## check presence of 0s or negative values
  if (any(data <= 0)) {
    stop("All values should be positive")
  }
  if(any(is.na(data))){
    warning("NA in data has been deleted")
    data <- as.matrix(data[complete.cases(data),])
    
    if(nrow(data)==0) data <- NULL
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
  
  clr.data
}

mult<-function(x,t){
  
  return(x*t)
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
      mat[D-i,j]=-sqrt(1/((D-i)*(D-i+1)))
    }
    mat[D-i,(D-i+1)]=sqrt((D-i)/(D-i+1))
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
  
  data <- norm_data(data)
  clr(data)%*%t(Base_SIGMA_matrix(ncol(data))) 
}

ilr_inverse<-function(data,k=1){
  
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
  
  data <- center_scale(data, scale=FALSE)
  
  Z=ilr(data)

  C=cov(Z)
  
  vect=eigen(C)$vector
  val=eigen(C)$values
  
  x=Z%*%vect
  
  n=dim(data)[1]
  rval <- list(variance_explain=cumsum(val)/sum(val), vector=vect, values=val, coord=x)
  
  class(rval) <- "biplot"
  
  rval
}


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

solve_regularisation<-function(var1){
  var1 %>% regularisation() %>% solve()
}

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

intervalle_confiance<-function(data,alpha,case=3,moy=NULL,var_matrix=NULL){
  
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

testing<-function(data1,data2,alpha,case=1){
    data1 <- ilr(norm_data(data1))
    data2 <- ilr(norm_data(data2))
    
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


count_to_proportion<-function(data){
  K <- ncol(data)
  alpha <- rep(2,K)
    
  #new.data <- apply(data, 2, sum) +alpha
  #somme <- sum(new.data)
  #new.data <- new.data/somme

  new.data <- sweep(data, 2, alpha, "+")
  somme <- apply(new.data, 1, sum)
  new.data <- sweep(new.data, 1, somme, "/")
  
  new.data
}


MAP <- function(data){
  
  K <- ncol(data)
  alpha <- rep(2,K)
  
  #new.data <- apply(data, 2, sum) +alpha
  #somme <- sum(new.data-K)
  #new.data <- (new.data-1)/somme
  
  new.data <- sweep(data, 2, alpha-1, "+")
  somme <- apply(new.data, 1, sum)
  new.data <- sweep(new.data, 1, somme, "/")
  
  new.data
  
}


graph_biplot_normale <- function(data, metadata_group, nb_graph=1, title=NULL, legend_title="group"){
  
  data_MAP <- MAP(data)
  
  b_data <- biplot(data_MAP);
  
  if(ncol(b_data$coord)<nb_graph-1){
    stop("number of graph incorrect")
  }
  
  metadata_group <- metadata_group %>% as.factor()
  
  
  m <- data.frame(b_data$coord, group=metadata_group)
  
  name_group <- nth(summarise(group_by(m,group)),1)
  
  ellipse_confiance <- list()
  for(i in 1:nb_graph) ellipse_confiance[[i]] <- data.frame()
  
  for(i in name_group){
    for(j in 1:nb_graph){
    
      if(((filter(m,group==i)[,j:(j+1)]) %>% as.matrix() %>%ilr_inverse() %>% Raduis_test())$Watson<0.187){
        temp <- as.matrix(filter(m,group==i)[,j:(j+1)])
        el <- intervalle_confiance(temp, alpha=0.05,1,moy=apply(temp,2,mean),var=var(temp))
      
        ellipse_confiance[[j]] <- rbind(ellipse_confiance[[j]], cbind(el, i))
    
      }
    }
    
  }
  
  check_type_data <- function(data_list){
    data_list[,1] <- data_list[,1] %>% as.character() %>% as.numeric()
    data_list[,2] <- data_list[,2] %>% as.character() %>% as.numeric()
    data_list[,3] <- data_list[,3] %>% as.character() %>% as.factor()
    
    data_list
  }
  
  for(i in 1:length(ellipse_confiance)){
    if(length(ellipse_confiance[[i]])!=0){
      ellipse_confiance[[i]] %<>% check_type_data()
    }
  }
  
  create_graph <- function(data, ellipse){
 
    graph_bc <- list()
    for(i in 1:nb_graph){
      if(length(ellipse_confiance[[i]])!=0){
        g <- eval(substitute( ggplot()+geom_point(aes(x=m[,i], y=m[,(i+1)],col=group),m)
                              +labs(title=title, y=paste("comp",i+1), x=paste("comp",i))
                              +scale_color_discrete(name=legend_title)
                              +geom_path(aes(x=ellipse_confiance[[i]][,1], y=ellipse_confiance[[i]][,2], group=ellipse_confiance[[i]][,3],col=ellipse_confiance[[i]][,3]),ellipse_confiance[[i]])
                              +theme(plot.title = element_text(hjust=0.5))
        ))
      }else{
        g <- eval(substitute( ggplot()+geom_point(aes(x=m[,i], y=m[,(i+1)],col=group),m)
                              +scale_color_discrete(name=legend_title)
                              +labs(title=title, y=paste("comp",i+1), x=paste("comp",i))
                              +theme(plot.title = element_text(hjust=0.5))
        ))
      }
      graph_bc[[i]] <- g
    }
    graph_bc
  }
  
  create_graph(m,ellipse_confiance)
  
}


comparaison_k_means <- function(data, metadata, nb_cluster=2, nb_graph=1, nb_start=50){
  
  data_ilr <- data %>% MAP() %>% ilr()
  
  k_data_ilr <- kmeans(data_ilr, nb_cluster, nstart=nb_start)
  k_data <- kmeans(data, nb_cluster, nstart=nb_start)
  
  table1 <- table(k_data$cluster, metadata)
  table2 <- table(k_data_ilr$cluster, metadata)
  
  grob <- c(graph_biplot_normale(data, k_data$cluster, title="comptage", nb_graph = nb_graph), graph_biplot_normale(data, metadata, title = "correct", nb_graph = nb_graph),  graph_biplot_normale(data, k_data_ilr$cluster, title = "ilr", nb_graph = nb_graph))
  
  list( comptage_table=table1, ilr_table=table2, graphics=grob)
}


comparaison_hclust <- function(data, metadata, nb_cluster,nb_graph){
  
  data_ilr <- data %>% MAP() %>% ilr()
  
  hclust_data <- data %>% dist() %>% as.dist() %>% hclust(method="ward.D") %>% cutree(k=nb_cluster)
  hclust_data_ilr <- data_ilr %>% dist() %>% as.dist() %>% hclust(method="ward.D") %>% cutree(k=nb_cluster)
  
  table1 <- table(hclust_data, metadata)
  table2 <- table(hclust_data_ilr, metadata)
  
  grob <- c(graph_biplot_normale(data, hclust_data, title="comptage", nb_graph = nb_graph), graph_biplot_normale(data, metadata, title = "correct", nb_graph = nb_graph),  graph_biplot_normale(data, hclust_data_ilr, title = "ilr", nb_graph = nb_graph))
  
  list( comptage_table=table1, ilr_table=table2, graphics=grob)
  
}