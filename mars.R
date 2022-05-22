library(rpart)
library(earth)

# main function
mars = function(formula, data, control){
  fwd_out = fwd_stepwise()
  bwd_out = bwd_stepwise(bwd_in = fwd_out)
  return(bwd_out)
}


#sub functions
fwd_stepwise = function(y, x, Mmax){
  if(Mmax<2) {
    warning("Mmax must be >= 2. Mmax will be set to 2")
    Mmax = 2
  } 
  
  N = length(y)
  n = ncol(x)
  B = matrix(1, nrow=N, ncol=1)
  splits = list(data.frame(m=0, v=0, s=NA, t=NA))
  
  M = 2
  while(!(M>Mmax)){
    lof_best = Inf
    
    for(m in 1:(M-1)){
      svars = setdiff(1:n, splits[[m]]$v)
      
      for(v in svars){
        tt = split_points(x[,v], B[,m])
        for(t in tt){
          Bnew = data.frame(B, 
                            Btem1=B[,m]*h(x[,v],1,t), 
                            Btem2=B[,m]*h(x[,v],-1,t))
          gdat = data.frame(y=y, Bnew)
          lof = LOF(y~., gdat)
          if(lof < lof_best){
            lof_best = lof 
            m_best = m
            v_best = v
            t_best = t
          }
        }
      }
    }
    
    left_split = rbind(splits[[m_best]], 
                       c(m_best, v_best, -1, t_best))
    right_split = rbind(splits[[m_best]], 
                        c(m_best, v_best, 1, t_best))
    
    splits = c(splits, list(left_split), list(right_split))
    
    B = cbind(B, 
              B[,m_best]*h(x[,v_best],1,t_best), 
              B[,m_best]*h(x[,v_best],-1,t_best))
    M = M + 2
    
  }
  
  colnames(B) = paste0("B", (0:(ncol(B)-1)))
  return(list(y=y, B=B, splits=splits))
  
}


bwd_stepwise = function(bwd_in){
  Mmax = ncol(bwd_in$B)
  Jstar = 2:Mmax
  Kstar = Jstar
  dat = data.frame(y=bwd_in$y, Jstar)
  lofstar = LOF(y~., dat)
  
  for(M in Mmax:2){
    b = Inf
    L = Kstar
    
    for(m in L){
      K = setdiff(L, m)
      dat = data.frame(y=bwd_in$y, bwd_in$B[,K])
      lof = LOF(y~., dat)
      if(lof < b){
        b = lof
        Kstar = K
      }
      if(lof < lofstar){
        lofstar = lof
        Jstar = K
      }
    }
  }
  Jstar = c(1,Jstar)
  return(list(y=bwd_in$y, B=bwd_in$B[,Jstar], splits=bwd_in$splits[Jstar]))
}

LOF = function(form, data, d=3){
  ff = lm(form,data)
  RSS = sum(residuals(ff)^2)
  N = nrow(data)
  M = length(coef(ff))-1
  Ctilde = sum(diag(hatvalues(ff))) + d*M
  return(RSS * N/(N-Ctilde)^2)
}

init_B <- function(N,Mmax) {
  B <- data.frame(matrix(NA,nrow=N,ncol=1))
  B[,1] <- 1
  names(B) <- c("B0",paste0("B",1:Mmax))
  return(B)
}

h = function(x, s, t){
  return(pmax(0,s*(x-t)))
}

split_points = function(xv, Bm){
  out = sort(unique(xv[Bm>0]))
  return(out[-length(out)])
}

mars.control <- function(){
  control <- list()
  
  return(control)
}


