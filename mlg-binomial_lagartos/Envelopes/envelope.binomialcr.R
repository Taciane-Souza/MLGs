envelope.binomialcr <- function(form=form, Fam=Fam, k=k, alfa=alfa, total = total){ 
  
  d<-length(total)
  
  alfa1<-ceiling(k*alfa)
  
  alfa2<-ceiling(k*(1-alfa)) 
  
  glm1<-glm(formula=form,family=Fam,maxit=50) 
  
  X<-model.matrix(glm1)
  
  w<-glm1$weights
  
  W<-diag(w)
  
  H<-sqrt(W)%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%sqrt(W)
  
  h<-diag(H)
  
  rd<-resid(glm1,type="deviance")   
  
  phi<-1
  
  td<-rd*sqrt(phi/(1-h)) 
  
  res<-matrix(0,d,k) 
  
  y<-numeric(d) 
  
  for(i in 1:k){ 
    
    for(j in 1:d){ 
      
      dif<-runif(total[j])-fitted(glm1)[j] 
      
      dif[dif>=0]<-0 
      
      dif[dif<0]<-1 
      
      y[j]<-sum(dif)
      
    } 
    
    xmat<-cbind(y,total-y) 
    
    fit<-glm(xmat~-1+X,family=Fam,maxit=50)
    
    w<-fit$weights
    
    W<-diag(w)
    
    H<-sqrt(W)%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%sqrt(W)
    
    h<-diag(H)
    
    res[,i]<-sort(resid(fit,type="pearson")*sqrt(phi/(1-h)))
    
  } 
  
  e1<-numeric(d) 
  e2<-numeric(d) 
  e11 <- numeric(d)
  e21 <- numeric(d)
  e12 <- numeric(d)
  e22 <- numeric(d)
  
  for(i in 1:d){ 
    
    eo<-sort(res[i,]) 
    
    e1[i]<-eo[ceiling(k*0.01)]
    e2[i]<-eo[ceiling(k*(1 - 0.01))]
    e11[i] = eo[ceiling(k*0.05)]
    e21[i] = eo[ceiling(k*(1 - 0.05))]
    e12[i] = eo[ceiling(k*0.1)]
    e22[i] = eo[ceiling(k*(1 - 0.1))]
    
  } 
  
  xb<-apply(res,1,mean)
  
  faixa<-range(td,e1,e2)
  
  par(pty="s")
  #--------------------------
  a <- qqnorm(e1, plot.it = FALSE)$x
  
  #par(new=TRUE)
  
  #qqnorm(e2,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=2)
  a2 <- qqnorm(e2, plot.it = FALSE)$x
  #par(new=TRUE)
  
  #qqnorm(xb,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=1)
  rxb <- qqnorm(xb, plot.it = FALSE)$x
  
  #par(new=TRUE)

  r <- qqnorm(td, plot.it = FALSE)$x
  #'@importFrom graphics par points polygon
  #'@importFrom stats qqnorm
  #'@import ggplot2
  #' @export
  
  df <- data.frame(r=r,xab=a,emin=cbind(e1,e11,e12),emax=cbind(e2,e21,e22),xb=xb,td=td,rxb=rxb)
  ggplot(df,aes(r,td)) + 
    geom_ribbon(aes(x=xab, ymin=e1, ymax=e2), fill= "#cf3266",alpha=0.5) +
    geom_ribbon(aes(x=xab, ymin=e11, ymax=e21),fill="#cf3266",alpha=0.5) +
    geom_ribbon(aes(x=xab, ymin=e12, ymax=e22),fill="#cf3266",alpha=0.5) +
    scale_fill_gradient(low = "#a30b3e", high = "#fa87ad") + 
    geom_point() + 
    geom_line(aes(rxb,xb),lty=2) +
    xlab("Quantis Teóricos") + 
    ylab("Quantis Empíricos") + 
    theme(text=element_text(size=15,family="serif"))
}




  
  
  