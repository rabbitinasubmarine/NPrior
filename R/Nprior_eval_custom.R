NPrior_run_custom = function(X, y, N=10000, BURN=2000, lam1, lam2, lam3, act_type, SpSL_type, 
                  sig = NULL, eta = 1, alpha0 = NULL, alpha = NULL, w = NULL,
                  sig_update = T, alpha0_update = T, 
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
if( is.null(alpha0)  == T ){
  if(SpSL_type == 0){
    alpha0 = 0
    alpha0_update = F
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

if(B_size > p){ B_size = p; print("B_size cannot be larger than p!") }
if(sig_update == T){ s_u = 1}else{ s_u = 0 }
if(alpha0_update == T){ a0_u = 1}else{ a0_u = 0 }
if(verbose == T){verbose = 1}else{verbose = 0}
pmt = proc.time()[3]
res = Neuro_Linear_custom( y=y, X=X, N=N, BURN=BURN, alpha=alpha,
                      w=w, sig=sig1, n1=n,  p1=p, eta = eta,
                      alpha0 = alpha0, size_a = s,
                      act_type = act_type, SpSL_type = SpSL_type, 
                      lam1=lam1, lam2=lam2, lam3=lam3,
                      eta_update = 0L, alpha0_update = a0_u,
                      sig_update = s_u, K = K, B_size=B_size, a0=a0, b0=b0, prior_sig_type, verbose = verbose)
time0 = proc.time()[3] - pmt
print(time0)
setting = list(act_type = act_type, SpSL_type= SpSL_type, 
               lam1 = lam1, lam2 = lam2, lam3 = lam3,
               eta = eta, alpha0 = alpha0, alpha0_update = alpha0_update, a0 = a0, b0 = b0, sig_update = sig_update, prior_sig_type = prior_sig_type)
return(list( THETA = res$THETA, SIG = res$sig, ALPHA0 = res$ALPHA0, ALPHA = res$ALPHA, W = res$W, POST = res$POST, setting = setting))
}










