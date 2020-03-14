NPrior_run = function(X, y, N=10000, BURN=2000, prior="SpSL-g", sig = NULL, eta = NULL, alpha0 = 0, alpha = NULL, w = NULL,
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
   if(prior == "HS" | prior == "BL"){
      eta = p^-2
   }else{
      eta = 1
   }
}

if(B_size > p){ B_size = p; print("B_size cannot be larger than p!") }
if(prior == "SpSL-g"){type=3L}
if(prior == "SpSL-c"){type=4L}
if(prior == "HS"){type=2L}
if(prior == "BL"){type=1L}
id = which(prior == c("SpSL-g", "SpSL-c", "BL", "HS"))
if(length(id) == 0){
  prior = "SpSL-g"
  print("The prior type is not correct! The default (SpSL-g) is going to be used.")
}
if(sig_update == T){ s_u = 1}else{ s_u = 0 }
if(alpha0_update == T){ a0_u = 1}else{ a0_u = 0 }
if(verbose == T){verbose = 1}else{verbose = 0}
pmt = proc.time()[3]
  if(type == 3 & method == "exact"){
    alpha0 = -1.0*stats::qnorm(prior_prop)
    res = Neuro_ReLU_Linear( y, X, N, BURN = BURN, alpha,  w,
                                sig1, n1=n,  p1=p, eta = eta, alpha0 = alpha0,
                                a0=a0, b0=b0, alpha0_update=a0_u,
                                eta_update = 0L,  sig_update = s_u,
                                B_size = B_size, prior_sig_type, verbose = verbose)
  }
  if(method == "RWMH" | type == 1 | type == 2 | type == 4 ){
    res = Neuro_Linear( y=y, X=X, N=N, BURN=BURN, alpha=alpha,
                      w=w, sig=sig1, n1=n,  p1=p, eta = eta,
                      alpha0 = alpha0, size_a = s,
                      type = type, eta_update = 0L, alpha0_update = a0_u,
                      sig_update = s_u, K = K, B_size=B_size, a0=a0, b0=b0, prior_sig_type, verbose = verbose)
  }
time0 = proc.time()[3] - pmt
print(time0)
setting = list(type = type, eta = eta, alpha0 = alpha0, alpha0_update = alpha0_update, a0 = a0, b0 = b0, sig_update = sig_update, prior_sig_type = prior_sig_type)
return(list( THETA = res$THETA, SIG = res$sig, ALPHA0 = res$ALPHA0 , POST = res$POST, setting = setting))
}










