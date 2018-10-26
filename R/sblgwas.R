#' @title Sparse Bayesian Learning for QTL Mapping and Genome-Wide Association Studies
#' @description The sparse Bayesian learning (SBL) method for quantitative trait
#' locus (QTL) mapping and genome-wide association studies (GWAS)
#' deals with a linear mixed model. This is also
#' a multiple locus model that includes all markers (random effects) in a single
#' model and detect significant markers simultaneously. SBL method adopts
#' coordinate descent algorithm to update parameters by estimating one parameter
#' at a time conditional on the posterior modes of all other parameters. The parameter
#' estimation process requires multiple iterations and the final estimated
#' parameters take the values when the entire program converges.
#' @details The multiple locus hierarchical linear mixed model of SBL is
#' \deqn{y=X\beta+Z\gamma+\epsilon}
#' where y is an \eqn{n*1} vector of response variables (observations of the trait);
#' X is an \eqn{n*p} design matrix for fixed effects; \eqn{\beta} is a \eqn{p*1}
#' vector of fixed effect; Z is an \eqn{n*m} genotype indicator matrix;
#' \eqn{\gamma} is an \eqn{m*1} vector of marker effects and \eqn{\epsilon} is
#' an \eqn{n*1} vector of residual errors with an aassumed
#' \eqn{\epsilon~N(0,\Sigma)} distribution. Each marker effect, \eqn{\gamma[k]}
#' for marker \emph{k}, is treated as a random variable following \eqn{N(0,\Phi[k])}
#' distribution, where \eqn{\Phi[k]} is the prior variance. The estimate of
#' \eqn{\gamma[k]} is best linear unbiased prediciton (BLUP). The estimate of
#' \eqn{\Phi[k]} is miximum likelihood estimate (MLE).
#' @param x a design matrix for fixed effects
#' @param y a vector of response variables
#' @param z a design matrix for random effects
#' @param t a number between [-2,0] to control sparseness of the model,
#' default is -1.
#' @param max.iter maximum number of iterations set by user to stop the program,
#' default is 200.
#' @param min.err minimum threshold of mean squared error of random effects estimated
#' from the current and the previous iteration to stop the program, default is 1e-6.
#' @return
#' \item{iteration}{a matrix storing intermediate results of each iteration before
#' the entire program converges, including \cr
#'
#' "\code{iter}" iteration indicator \cr
#'
#' "\code{error}" mean squared error of random effects estimated from the current and
#' the previous iteration \cr
#'
#' "\code{s2}" estimated variance of residual error \cr
#'
#' "\eqn{beta[1]\dots beta[p]}" estimates of fixed effects \cr
#'
#' "\eqn{gamma[1]\dots gamma[m]}" estimates of random effects}
#'
#' \item{parm}{a vector containing 5 elements: "\code{iter}", "\code{error}",
#' "\code{s2}", "\code{beta}" and "\code{df}" \cr
#'
#' "\code{iter}" the number of iterations required by program to stop \cr
#'
#' "\code{error}" mean square error of random effect estimated from the last iteration
#' before the program stops \cr
#'
#' "\code{s2}" estimated variance of residual error \cr
#'
#' "\code{beta}" estimate of fixed effect \cr
#'
#' "\code{df}" the effective degree of freedom from total random effects}
#'
#' \item{blup}{a matrix containing 4 columns: "\code{gamma}", "\code{vg}",
#' "\code{wald}" and "\code{p_wald}" \cr
#'
#' "\code{gamma}" estimate of random effect \cr
#'
#' "\code{vg}" estimated variance of random effect \cr
#'
#' "\code{wald}" Wald statistic calculated as \eqn{\gamma^2}/\eqn{\Phi}
#'
#' "\code{p_wald}" the \emph{p}-value of Wald statistic following Chi-squared
#' distribution with 1 degree of freedom}
#'
#' @author Meiyue Wang and Shizhong Xu \cr
#' Maintainer: Meiyue Wang \email{mwang024@@ucr.edu}
#' @examples
#' # Load example data from sbl package
#' data(gen)
#' data(phe)
#' data(intercept)
#'
#' # Run sblgwas() to perform association study of example data
#' # setting t = 0 leads to the most sparse model
#' fit<-sblgwas(x=intercept, y=phe, z=gen, t=0)
#' my.blup<-fit$blup
#'
#' # setting t = -2 leads to the least sparse model
#' fit<-sblgwas(x=intercept, y=phe, z=gen, t=-2)
#' my.blup<-fit$blup
#' @importFrom stats pchisq
#' @export

sblgwas<-function(x,y,z,t=-1,max.iter=200,min.err=1e-6){
  x<-as.matrix(x)
  y<-as.matrix(y)
  z<-as.matrix(z)
  n<-length(y)
  q<-ncol(x)
  m<-ncol(z)
  b0<-solve(t(x)%*%x,tol=1e-50)%*%(t(x)%*%y)
  s2<-sum((y-x%*%b0)^2)/(n-q)
  b0<-matrix(0,q,1)
  b<-b0
  g0<-matrix(0,m,1)
  g<-g0
  lambda<-matrix(0,m,1)
  tau<-g0
  v<-g0
  xx<-NULL
  xy<-NULL
  for(i in 1:q){
    xx<-c(xx,sum(x[,i]^2))
    xy<-c(xy,sum(x[,i]*y))
  }
  zz<-NULL
  zy<-NULL
  for(k in 1:m){
    zz<-c(zz,sum(z[,k]^2))
    zy<-c(zy,sum(z[,k]*y))
  }
  d<-numeric(m)
  a<-matrix(0,n,1)
  iter<-0
  err<-1e8
  my.iter<-NULL
  while(iter < max.iter & err > min.err){
    for(i in 1:q){
      a<-a-x[,i]*b0[i]
      ai<-sum(x[,i]*a)
      b[i]<-(xy[i]-ai)/xx[i]
      a<-a+x[,i]*b[i]
    }
    df<-0
    for(k in 1:m){
      a<-a-z[,k]*g0[k]
      ak<-sum(z[,k]*a)
      c1<- -(t+3)*zz[k]^2
      c2<- -(2*t+5)*zz[k]+(zy[k]-ak)^2
      c3<- -(t+2)
      if( ((c2^2-4*c1*c3) < 0) | (c2 < 0) ){
        tau[k]<-0
      } else {
        tau[k]<-(-c2-sqrt(c2^2-4*c1*c3))/(2*c1)
      }
      lambda[k]<-tau[k]/s2
      g[k]<-lambda[k]*(zy[k]-ak)-lambda[k]^2*zz[k]*(zy[k]-ak)/(lambda[k]*zz[k]+1)
      d[k]<-lambda[k]*(zz[k]-lambda[k]*zz[k]^2/(lambda[k]*zz[k]+1))
      v[k]<-tau[k]-tau[k]*d[k]
      df<-df+d[k]
      a<-a+z[,k]*g[k]
    }

    if((n-q-df) > 0){s2<-sum((y-a)^2)/(n-q-df)
    }else{
      s2<-sum((y-a)^2)/(n-q)
    }

    iter<-iter+1
    err<-sum((g-g0)^2)/m
    g0<-g
    b0<-b
    my.iter<-rbind(my.iter,cbind(iter,err,s2,t(b),t(g)))
  }
  my.parm<-data.frame(iter,err,s2,b,df)
  names(my.parm)<-c("iter","error","s2","beta","df")

  posv<-which(v!=0)
  m<-length(g)
  wald<-c(rep(0,m))
  gg<-g[posv]
  vv<-v[posv]
  wald[posv]<-gg^2/vv
  p<-pchisq(wald,1,lower.tail=FALSE)

  my.blup<-data.frame(g,v,wald,p)
  names(my.blup)<-c("gamma","vg","wald","p_wald")

  var.beta<-NULL
  for(i in 1:q){
    var.beta<-c(var.beta,paste("beta",i,sep=""))
  }
  var.gamma<-NULL
  for(k in 1:m){
    var.gamma<-c(var.gamma,paste("gamma",k,sep=""))
  }
  var.names<-c(c("iter","error","s2"),var.beta,var.gamma)
  my.iter<-data.frame(my.iter)
  names(my.iter)<-var.names

  out<-list(my.iter,my.parm,my.blup)
  names(out)<-c("iteration","parm","blup")
  return(out)
}


