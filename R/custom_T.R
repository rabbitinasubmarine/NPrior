# S: samples from prior
custom_T = function(S,act_type=1,plot=T){
  N = length(S)
  target_sam = sort(S)
  qq = c(0.01,0.1,0.2,0.4,0.6,0.8,0.9,0.99)
  q0 = stats::quantile(target_sam,qq)
  wgt_q = c(1.2,0.7,0.5,0.4,0.4,0.5,0.7,1.2)^5
  a1_cand = seq(0,0.5, length.out = 10)
  a2_cand = seq(-2,2, length.out = 20)
  a3_cand = seq(-2, 2, length.out = 20)
  if(act_type==1){
    act_fun = function(x, a1, a2, a3){
      exp(a1*x*abs(x) + a2*x + a3) 
    }
  }
  if(act_type == 2){
    act_fun = function(x, a1, a2, a3){
      exp(a1*x*x + a2*x + a3) 
    }
  }
  alpha = stats::rnorm(N)
  w = stats::rnorm(N)
  for(k in 1:10){
    if(k>1){
      i1 = ind[1,1]
      i2 = ind[1,2]
      i3 = ind[1,3]
      if(i1 == 1 | i1 == length(a1_cand)){
        m = a1_cand[i1] - 0.1
        if(m < 0) m = 0
        a1_cand = seq(m, a1_cand[i1]+0.1,length.out = 10)
      }else{
        a1_cand = seq(a1_cand[i1-1]-0.02, a1_cand[i1+1]+0.02,length.out = 10)
      }
      if(i2 == 1 | i2 == length(a2_cand)){
        a2_cand = seq(a2_cand[i2]-0.1, a2_cand[i2]+0.1,length.out = 10)
      }else{
        a2_cand = seq(a2_cand[i2-1]-0.02, a2_cand[i2+1]+0.02,length.out = 10)
      }
      if(i3 == 1 | i3 == length(a3_cand)){
        a3_cand = seq(a3_cand[i3]-0.1, a3_cand[i3]+0.1,length.out = 10)
      }else{
        a3_cand = seq(a3_cand[i3-1]-0.02, a3_cand[i3+1]+0.02,length.out = 10)
      }
      
    }
  D = array( 0, dim=c(length(a1_cand),length(a2_cand),length(a3_cand))) 
  i1 = 1
  for(i1 in 1:length(a1_cand)){
    for(i2 in 1:length(a2_cand)){
      for(i3 in 1:length(a3_cand)){
        a1 = a1_cand[i1]
        a2 = a2_cand[i2]
        a3 = a3_cand[i3]
        propose_sam = act_fun(alpha, a1, a2, a3)*w
        propose_sam = sort(propose_sam)
        a = abs(target_sam - propose_sam)#*wgt
        q1 = stats::quantile(propose_sam,qq)
        D[i1,i2,i3] = mean(abs(q0-q1)/wgt_q) + mean(a)
      }
      #print(i2)
    }
  }  
  m = min(D)
  ind = which(D == m, arr.ind = T)
  a1 = a1_cand[ind[1,1]]
  a2 = a2_cand[ind[1,2]]
  a3 = a3_cand[ind[1,3]]
  propose_sam = act_fun(alpha, a1, a2, a3)*w
  #print(ind)
  print(paste("lam1: ",round(a1,4)))
  print(paste("lam2: ",round(a2,4)))
  print(paste("lam3: ",round(a3,4)))
  print("#####################")
  }
  qq = c(0.6,0.7,0.8,0.95)
  print(paste("Target: "))
  print(quantile(target_sam,qq))
  print(paste("Neuronized: "))
  print(quantile(propose_sam,qq))
  if(plot==T){
    stats::qqplot(propose_sam, target_sam,xlab="Neuronized",ylab="Target",cex=0.8,lwd=0.8)
    graphics::abline(0,1,col="red",lwd=1.5)
  }
return(list(lam1=a1,lam2=a2,lam3=a3, act_type=act_type))

}



