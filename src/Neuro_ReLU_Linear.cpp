//[[Rcpp::depends( RcppArmadillo )]]
#include <RcppArmadillo.h>
//#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp/Benchmark/Timer.h>
#include <cmath>

using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List Neuro_ReLU_Linear(const arma::vec & y, const arma::mat & X,
                              int  N, int  BURN, arma::colvec alpha,  arma::colvec w,
                              arma::colvec sig, double n1, double p1, double  eta,  double  alpha0,
                              double a0, double b0, int alpha0_update,
                              const int & eta_update, const int & sig_update, const int & B_size , const int & prior_sig_type, const int & verbose){
  int n =  X.n_rows, p = X.n_cols;
  int k;
  int k1;
  double lam = 0.0;
  double C = 0.0;
  double k0 = 0.0;
  double alpha_cand = 0.0;
  double weight = 1.0;
  double acc_a1 = 0.0;
  double acc_a2 = 0.0;
  
  
  int J0 = 5;
  int J = 10;
  int jj;
  double aux=0.0;
  double aux_cand=0.0;
  double m0;
  double acc_a0 = 0.0;
  double cand1, curr1;
  arma::colvec res2(1);
  NumericVector sam_prob(J);
  IntegerVector sam_J(J);
  IntegerVector s(1);
  NumericVector CAND(J);
  arma::colvec alpha0_CAND(J);
  
  arma::colvec r(B_size);
  arma::colvec X_norm(p);
  double tau0, par0, inv_par0, S, mu1;
  arma::colvec ALPHA0(N);
  //arma::mat C_S(p,p);
  arma::mat G_alpha(n,p);
  arma::colvec theta(p); theta = alpha % w;
  arma::colvec Act_alpha(p); Act_alpha.zeros();
  arma::colvec res(n);
  arma::colvec rx(1), obj(1), b(1);
  double aa, bb, cc, a, c, d;
  double sig0 = 1.0;
  IntegerVector sam(p), sam1(p);

  arma::mat SAVE_theta(p,N); SAVE_theta.zeros();
  arma::mat SAVE_alpha(p,N); SAVE_alpha.zeros();
  arma::mat SAVE_w(p,N); SAVE_w.zeros();
  arma::colvec OBJ(N), SIG(N), theta_save(p), gam_save(p), gam(p);
  theta_save.zeros();
  gam_save.zeros();
  gam.zeros();
  arma::uvec ind;
  Timer timer;
  int i;
  int ii;
  int j;
  //int start, ii_B, iii_B;
  double wt0, wt1,  a_t, sig_t,  prop ;
  double tau_a = 1.0;
  //NumericVector prob0(p);
  for(ii=0; ii<p; ii++){
    sam1(ii) = ii;
    sam(ii) = 1;
    X_norm(ii) = sum(X.col(ii) % X.col(ii));
   // prob0(ii) = 1/p1;
  }
  wt0 = R::pnorm( alpha0 , 0.0, 1.0, 1, 1);
  for(i=0; i<(N+BURN); i++){
    for(ii=0; ii<p; ii++){
      G_alpha.col(ii) = X.col(ii)*Act_alpha(ii);
    }
    
    //if(type == 3 || type == 4 || type == 5){
      if(alpha0_update == 1){
        if(a0 == 1.0 & b0 == 1.0){
          //if(i % 10 == 0){
          //a = sum(gam) + a0;
          //d = p - sum(gam) + b0;
          //c = R::rbeta(a,d);
          //alpha0 = R::qnorm(c, 0.0, 1.0, 0, 0);
          //}
          a = -1.0*(alpha0 + sum(alpha))/(1.0 + p1);
          d = 1.0/(1.0 + p1);
          aux = a + R::rnorm(0.0, 1.0)*sqrt(d); 
          alpha0 += aux;
          alpha += aux;
          //Rcpp::Rcout << "#######" << std::endl;
        }else{
          a = -1.0*(alpha0 + sum(alpha))/(1.0+p1);
          d = 1.0/(1.0+p1);
          for(jj=0; jj<J0; jj++){
            curr1 = 0.0;
            curr1 += -0.5*(alpha0 + aux)*(alpha0 + aux) + (a0-1.0)*R::pnorm(alpha0+aux, 0.0, 1.0, 0, 1);
            res2 = (alpha+aux).t()*(alpha+aux);
            curr1 += -0.5*res2(0) + (b0-1.0)*R::pnorm(alpha0+aux, 0.0, 1.0, 1, 1);
            curr1 += 0.5*d*(aux-a)*(aux-a);
            for(ii=0; ii<J; ii++){
              aux_cand = a + R::rnorm(0.0, 1.0)*sqrt(d); 
              alpha0_CAND(ii) = aux_cand;
              cand1 = 0.0;
              cand1 += -0.5*(alpha0 + aux_cand)*(alpha0 + aux_cand) + (a0-1.0)*R::pnorm(alpha0+aux_cand, 0.0, 1.0, 0, 1);
              res2 = (alpha+aux_cand).t()*(alpha+aux_cand);
              cand1 += -0.5*res2(0) + (b0-1.0)*R::pnorm(alpha0+aux_cand, 0.0, 1.0, 1, 1);
              cand1 += 0.5*d*(aux_cand-a)*(aux_cand-a);
              CAND(ii) = cand1;
            }
            m0 = max(CAND);
            sam_prob = exp(CAND-m0);//+0.00000000001;
            s = Rcpp::sample(sam_J, 1, FALSE, sam_prob);
            aux_cand = alpha0_CAND(s(0));
            bb = sum(sam_prob);
            aa = bb - sam_prob(s(0)) + exp(curr1-m0);
            //curr1 = log(aa+0.00000000001);
            //cand1 = log(bb+0.00000000001);
            cc = log(bb/aa);//cand1 - curr1;
            aa = R::runif( 0.0, 1.0);
            if(cc > log(aa) ){
              aux = aux_cand;
              acc_a0 += 1.0;
            }
          }
          alpha += aux;
          alpha0 += aux;
        }
      }
    //}
    
    res = y - X*theta;
      for(j=0; j<p; j++){
      res = res +  X.col(j)*theta(j);
      //if(gam(j) == 0.0){
      //   w(j) = R::rnorm(0.0, 1.0)*sqrt(eta);
      if(alpha(j)>alpha0){
        rx = res.t()*X.col(j);
        S = 1.0 / (X_norm(j)*Act_alpha(j)*Act_alpha(j) + sig0/eta);
        mu1 = rx(0)*Act_alpha(j)*S;
        r(0) = R::rnorm( 0.0 , 1.0 );
        w(j) = mu1 + sqrt(sig(0)*S)*r(0);
        theta(j) = w(j)*Act_alpha(j);
      }else{
        w(j) = sqrt(eta*sig0)*R::rnorm( 0.0 , 1.0 );
        theta(j) = 0.0;
      }
        res = res -  X.col(j)*theta(j);
      }



    //sam = RcppArmadillo::sample(sam1, p, FALSE, prob);
    sam = Rcpp::sample(sam1, p, FALSE, NumericVector::create());
      res = y - X*theta;
      for(ii=0; ii<p; ii++){
        j = sam(ii);
        res = res +  X.col(j)*theta(j);
        rx = res.t()*X.col(j);
        a_t = (rx(0)*w(j) + X_norm(j)*alpha0*w(j)*w(j))/(X_norm(j)*w(j)*w(j) + sig(0));
        sig_t = sig(0)/(X_norm(j)*w(j)*w(j) + sig(0));
        wt0 = R::pnorm( alpha0, 0, 1.0, 1, 1 );
        //wt1 = 0.5*log( sig(0) ) - 0.5*log( X_norm(j)*w(j)*w(j) + sig(0)/tau_a  ) + R::pnorm( alpha0, a_t, sqrt(sig_t), 0, 1 );
        wt1 = 0.5*log( sig(0) ) - 0.5*log( X_norm(j)*w(j)*w(j) + sig(0)  ) + R::pnorm( alpha0, a_t, sqrt(sig_t), 0, 1);
        wt1 = wt1 +  0.5*a_t*a_t / sig_t ;
        wt1 = wt1 - 0.5*( X_norm(j)*alpha0*alpha0*w(j)*w(j) + 2*rx(0)*alpha0*w(j) )/sig(0);
        aa = R::runif( 0.0, 1.0);
        prop = 1/(1 + exp(wt0 - wt1) );
        for(k1 =0; k1 < 1; k1++){  
          if(  prop > aa ){
            for(k=0; k < 20; k++){
              if((alpha0 - a_t)/sqrt(sig_t) > 0.0){  
                k0 = (alpha0- a_t)/sqrt(sig_t);
                lam = (k0+sqrt(k0*k0+4.0))/2.0;
                alpha_cand = R::rexp(lam);
                alpha_cand += k0;
                C = 0.5*(lam*lam - 2.0*lam*k0) - log(lam);
                weight = -0.5*alpha_cand*alpha_cand - C - log(lam) + lam*alpha_cand; 
                aa = R::runif( 0.0, 1.0);
                if(log(aa) < weight){
                  alpha(j) = a_t + alpha_cand*sqrt(sig_t);
                  acc_a1 += 1.0;
                }
              }else{
                alpha_cand = a_t + R::rnorm(0.0, 1.0)*sqrt(sig_t);
                if(alpha_cand > alpha0){
                  alpha(j) = alpha_cand;
                  acc_a1 += 1.0;
                }
              }
            }
            Act_alpha(j) = alpha(j) - alpha0;
            theta(j) = (alpha(j) - alpha0)*w(j);
            gam(j) = 1.0;
          }else{
            if(alpha0 < 0.0){
              k0 = -1.0*alpha0;
              lam = (k0+sqrt(k0*k0+4.0))/2.0;
              alpha_cand = R::rexp(lam);
              alpha_cand += k0;
              C = 0.5*(lam*lam - 2.0*lam*k0) - log(lam);
              weight = -0.5*alpha_cand*alpha_cand - C - log(lam) + lam*alpha_cand; 
              aa = R::runif( 0.0, 1.0);
              if(log(aa) < weight){
                alpha(j) = -1.0*alpha_cand;
                acc_a2 += 1.0;
              }
            }else{
              alpha_cand = R::rnorm(0.0, 1.0);
              if(alpha_cand < alpha0){
                alpha(j) = alpha_cand;
                acc_a2 += 1.0;
              }
            }
            
            Act_alpha(j) = 0.0;
            theta(j) = 0.0;
            gam(j) = 0.0;
          }
        }
        res = res - X.col(j)*theta(j);
      }
      if(eta_update == 1){
      tau0 = eta;
      par0 = sum(w % w);
      inv_par0 = 1.0/(tau0);
      aa = 1.0/(inv_par0 + 1.0);
      //cc = R::runif( 0.0, aa);
      if(aa > 0.00000001){
        cc = R::runif( 0.0, aa);
      }else{
        cc = R::runif(0.0, 0.00000001);
      }
      aa = 1.0/cc - 1.0;
      bb = R::pgamma(aa, (1 + p1)/2, (2*sig(0))/par0, 1, 0);
      cc = R::runif( 0.0, bb);
      inv_par0 = R::qgamma(cc, (1 + p1)/2, (2*sig(0))/par0, 1, 0);
      eta = 1.0/inv_par0;
    }
    res = y - X*theta;
    if(prior_sig_type == 0){
      a = n1*0.5 + 1.0;
      b = res.t()*res*0.5 + 1.0;
    }else{
      a = n1*0.5 + p1*0.5 + 1.0;
      b = res.t()*res*0.5 + w.t()*w*0.5/eta + 1.0;
    }
    if( sig_update == 1 ){
      sig(0) = R::rgamma(a,1.0/b(0));
      sig(0) = 1.0/sig(0);
    }
    if(prior_sig_type == 0){ sig0 = sig(0);}else{sig0 = 1.0;}
    obj = -0.5*res.t()*res/sig(0) - w.t()*w*0.5/eta/sig0 - alpha.t()*alpha*0.5  - a*log(sig(0));
    if((i+1) > BURN){
      timer.step("");
      OBJ(i-BURN) = obj(0);
      ALPHA0(i-BURN) = alpha0;
      theta_save = theta_save + theta/N;
      gam_save = gam_save + gam/N;
      SAVE_theta.col(i-BURN) = theta;
      SAVE_alpha.col(i-BURN) = alpha;
      SAVE_w.col(i-BURN) = w;
      SIG(i-BURN) = sig(0);
    }
    ind = arma::find(gam == 1.0);

    if(verbose == 1){
      if((i+1) % 3000 == 0){
        Rcpp::Rcout << "ReLU: #iterations = " << i+1 << std::endl;
        Rcpp::Rcout << "Selected Variables = " << (ind+1).t() << std::endl;
      }
    }
    //yy = as<double>(Act(j,alpha0));
    //j = j + 1.0;
  }
  return Rcpp::List::create(Rcpp::Named("CompTime") = timer,
                            Rcpp::Named("THETA") = SAVE_theta,
                            Rcpp::Named("ALPHA0") = ALPHA0,
                            Rcpp::Named("theta") = theta_save,
                            Rcpp::Named("ALPHA") = SAVE_alpha,
                            Rcpp::Named("gam") = gam_save,
                            //Rcpp::Named("alpha") = SAVE_alpha,
                            Rcpp::Named("W") = SAVE_w,
                            Rcpp::Named("POST") = OBJ,
                            Rcpp::Named("acc_a0") = acc_a0/(N+BURN),
                            
                            //Rcpp::Named("wt1") = wt1,
                            //Rcpp::Named("wt0") = wt0,
                            //Rcpp::Named("a_t") = a_t,
                            //Rcpp::Named("sig_t") = sig_t,
                            //Rcpp::Named("j") = j,
                            //Rcpp::Named("alpha0") = alpha0,
                            Rcpp::Named("sig") = SIG);
}
