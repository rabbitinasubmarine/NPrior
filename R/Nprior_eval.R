NPrior_run = function(X, y, N=10000, BURN=2000, prior="SpSL-L", sig = NULL, eta = NULL, alpha0 = NULL, alpha = NULL, w = NULL,
                  sig_update = T, alpha0_update = T, method = "exact",
                  prior_prop = 0.01, s=2, K=10, B_size = NULL, a0 = 1, b0 = 1, prior_sig_type = 0, verbose = T){
# SpSL-g
# SpSL-c
# BL
# HS
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

n = nrow(X)
p = ncol(X)
if(prior == "SpSL-L"){type=3L}
if(prior == "SpSL-C"){type=4L}
if(prior == "SpSL-G"){type=5L}
if(prior == "HS"){type=2L}
if(prior == "BL"){type=1L}

if( is.null(alpha0)  == T ){
   if(type == 1 | type == 2){
     alpha0 = 0
     
   }else{
     alpha0 = -stats::qnorm(1/p)
   }
  }
if( is.null(alpha)  == T ){alpha = stats::rnorm(p)*0.1}
if( is.null(w)  == T ){w = stats::rnorm(p)*0.1}
if( is.null(sig)  == T ){sig = stats::rnorm(1)^2; sig1 = matrix(sig,1,1)}
if( is.null(B_size)  == T ){
  if(p > 50){
    B_size = 50
  }else{
    B_size = p
  }
}
if(is.null(eta)){
   if(prior == "HS"){
      #eta = p^-2
     tau_cand = c(-4,-3,seq(-2.5,-0.5,length.out=100))*log(p)
     kappa_target = 1-min(0.01, 0.1*n/p)
     N = 10000
     C = abs(stats::rcauchy(N))
     KAPPA = matrix(0,N,length(tau_cand))
     for(k in 1:length(tau_cand)){
       KAPPA[,k] =  1/(1+exp(tau_cand[k]+2*log(C)))
     }
     kappa = apply(KAPPA,2,mean)
     id = which.min(abs(kappa - kappa_target))
     eta = exp(tau_cand[id])
   }else{
      eta = 1
   }
  if(prior == "BL"){
    eta = p^-1
  }
}

if(B_size > p){ B_size = p; print("B_size cannot be larger than p!") }
id = which(prior == c("SpSL-L", "SpSL-C", "BL", "HS", "SpSL-G"))
if(length(id) == 0){
  prior = "SpSL-L"
  print("The prior type is not correct! The default (SpSL-L) is going to be used.")
}
if(sig_update == T){ s_u = 1}else{ s_u = 0 }
if(alpha0_update == T){ a0_u = 1}else{ a0_u = 0 }
if(type == 1 | type == 2){a0_u = 0}
if(verbose == T){verbose = 1}else{verbose = 0}
pmt = proc.time()[3]
#alpha0 = -1.0*stats::qnorm(prior_prop)
  if(type == 3 & method == "exact"){
    res = Neuro_ReLU_Linear( y, X, N, BURN = BURN, alpha,  w,
                                sig1, n1=n,  p1=p, eta = eta, alpha0 = alpha0,
                                a0=a0, b0=b0, alpha0_update=a0_u,
                                eta_update = 0L,  sig_update = s_u,
                                B_size = B_size, prior_sig_type, verbose = verbose)
  }
  if(method == "RWMH" | type == 1 | type == 2 | type == 4 | type == 5 ){
    res = Neuro_Linear( y=y, X=X, N=N, BURN=BURN, alpha=alpha,
                      w=w, sig=sig1, n1=n,  p1=p, eta = eta,
                      alpha0 = alpha0, size_a = s,
                      type = type, eta_update = 0L, alpha0_update = a0_u,
                      sig_update = s_u, K = K, B_size=B_size, a0=a0, b0=b0, prior_sig_type, verbose = verbose)
  }
time0 = proc.time()[3] - pmt
print(time0)
setting = list(type = type, eta = eta, alpha0 = alpha0, alpha0_update = alpha0_update, a0 = a0, b0 = b0, sig_update = sig_update, prior_sig_type = prior_sig_type)
return(list( THETA = res$THETA, SIG = res$sig, ALPHA0 = res$ALPHA0, ALPHA = res$ALPHA, W = res$W, POST = res$POST, setting = setting))
}










