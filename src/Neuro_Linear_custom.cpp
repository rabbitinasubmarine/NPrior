//[[Rcpp::depends( RcppArmadillo )]]
#include <RcppArmadillo.h>
//#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp/Benchmark/Timer.h>
using namespace Rcpp;


double T_fun(const double& x, const double& alpha0, 
             const double& lam1, const double& lam2, const double& lam3,
             const int& act_type, const int & SpSL_type){
  double act;
  if(SpSL_type == 0){
    if(act_type == 1){
      act = exp( lam1*x*fabs(x)  + lam2*x + lam3);
    }
    if(act_type == 2){
      //  act = exp( 0.37*x*fabs(x)  + 0.89*x +0.08 );
     act = exp( lam1*x*x  + lam2*x + lam3);
    }
  }else{
    if(act_type == 1){
      //  act = exp( 0.37*x*fabs(x)  + 0.89*x +0.08 );
      if(x > alpha0){
        act = exp( lam1*x*fabs(x)  + lam2*x + lam3); 
      }else{
        act = 0.0;
      }
    }
    if(act_type == 2){
      //  act = exp( 0.37*x*fabs(x)  + 0.89*x +0.08 );
      if(x > alpha0){
        act = exp( lam1*x*x  + lam2*x + lam3); 
      }else{
        act = 0.0;
      }
    }
  }
  return (act);
}

// type 1 : Linear activation (Approximately Bayesian LASSO); t/alpha0
// type 2 : Exponetial activation (Approximately horseshoe prior); exp( 0.5*t*|t| + alpha0*t )
// type 3 : ReLU; max(t,alpha0)

// [[Rcpp::export]]
Rcpp::List Neuro_Linear_custom(const arma::vec & y, const arma::mat & X, int  N, int BURN, arma::colvec alpha,  arma::colvec w, arma::colvec sig, double n1, double p1, double  eta, double  alpha0, double size_a,  
                               const int & act_type, const int & SpSL_type,
                        const double & lam1, const double & lam2, const double & lam3,
                        const int & eta_update, const int & alpha0_update, const int & sig_update,const int & K, const int & B_size, 
                        const double & a0, const double & b0, const int & prior_sig_type, const int & verbose){
  int n =  X.n_rows, p = X.n_cols;
  int B = p / B_size;
//double a0 = 1.0, b0 = 1.0;
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
  
  arma::mat X_sub(n,B_size);
  arma::mat XtX(B_size,B_size);
  arma::colvec Xty(B_size);
  arma::mat D(B_size,B_size);  D.eye();
  arma::mat G_alpha2(n,B_size);
  arma::colvec alpha_sub(B_size);
  arma::colvec mu(B_size);
  arma::uvec indx(B_size);
  arma::colvec r(B_size);

  arma::mat G_alpha(n,p);
  arma::colvec theta(p); theta = alpha % w;
  arma::colvec prob(p), Act_alpha(p); Act_alpha.ones();
  //arma::colvec sig(1); sig(0) = 1.0;
  arma::colvec res(n), cand(1), curr(1);
  arma::colvec rx(1), obj(1), b(1);
  double alpha_cand,  theta_cand; theta_cand = 0.0;
  double aa, bb, cc, a,d,c, S, mu1;
  double act; act = 0.0;
  double act_cand; act_cand = 0.0;
  double sig0 = 1.0;
  //arma::colvec prob0(p);
  //double theta_curr;
  IntegerVector sam(p), sam1(p);
  arma::colvec X_norm(p);
  double tau0, par0, inv_par0;
  arma::mat SAVE_theta(p,N); SAVE_theta.zeros();
  arma::mat SAVE_alpha(p,N); SAVE_alpha.zeros();
  arma::uvec ind;
  arma::mat SAVE_w(p,N); SAVE_w.zeros();
  arma::colvec OBJ(N), SIG(N), ETA(N), theta_save(p), gam_save(p), gam(p);
  theta_save.zeros();
  gam_save.zeros();
  gam.zeros();
  Timer timer;
  arma::colvec ALPHA0(N);
//if(type < 3){
//    alpha0=0;
//  }
  //NumericVector prop(p);
  //double res;
  int i;
  int ii;
  int j;
  int ii_B, iii_B;
  int start;
  double pq = 0.0;

  for(ii=0; ii<p; ii++){
    sam1(ii) = ii;
    sam(ii) = 1;
    X_norm(ii) = sum(X.col(ii) % X.col(ii));
    act = T_fun(alpha(ii), alpha0, lam1, lam2, lam3,
                act_type, SpSL_type);
    //prob0(ii) = 1/p1;
  }
  for(i=0; i< (N+BURN); i++){
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
    
    for(ii=0; ii<p; ii++){
      Act_alpha(ii) = T_fun(alpha(ii), alpha0, lam1, lam2, lam3,
                act_type, SpSL_type);
      G_alpha.col(ii) = X.col(ii)*Act_alpha(ii);
    }
    if(eta_update == 1){
      par0 = sum(w % w);
      tau0 = eta;
      inv_par0 = 1.0/(tau0);
      aa = 1.0/(inv_par0 + 1.0);
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
if( i % 2 == 0){
    if( SpSL_type==1 ){
      for(j=0; j<p; j++){
        res = res +  X.col(j)*theta(j);
        rx = res.t()*X.col(j);
        if(gam(j) == 0.0){
          w(j) = R::rnorm(0.0, 1.0)*sqrt(eta*sig0);
        }else{
          S = 1.0 / (X_norm(j)*Act_alpha(j)*Act_alpha(j) + sig0/eta);
          mu1 = rx(0)*Act_alpha(j)*S;
          r(0) = R::rnorm( 0.0 , 1.0 );
          w(j) = mu1 + sqrt(sig(0)*S)*r(0);
        }
        theta(j) = w(j)*Act_alpha(j);
        res = res -  X.col(j)*theta(j);
      }
    }else{
     //sam = RcppArmadillo::sample(sam1, p, FALSE, prob0);
      for(ii_B = 1; ii_B <  (B+1); ii_B++ ){
       r.randn();
       start = (ii_B - 1)*B_size;
       for(iii_B = 0;iii_B <  B_size; iii_B++ ){
         //indx(iii_B) = sam(start + iii_B);
         indx(iii_B) = start + iii_B;
       }
      X_sub = X.cols(indx);
      alpha_sub = alpha.rows(indx);
      G_alpha2 = G_alpha.cols(indx);
      res  = res +  X.cols(indx)*theta.rows(indx);
      //Rcpp::Rcout << "#####################"<< std::endl;
      XtX = G_alpha2.t()*G_alpha2 + sig0*D/eta;
      Xty = G_alpha2.t()*res;
      XtX = XtX.i();
      mu = XtX*Xty;
      XtX = arma::chol(XtX);
      w.rows(indx)  = mu + sqrt(sig(0))*(XtX.t()*r);
      theta.rows(indx) = Act_alpha.rows(indx) % w.rows(indx);
      res = res - X.cols(indx)*theta.rows(indx);
    }
  }
  }

  //  sam = RcppArmadillo::sample(sam1, p, FALSE, prob0);
    res = y - X*theta;
    for(ii=0; ii<p; ii++){
    //  j = sam(ii);
      j = ii;
      res = res +  X.col(j)*theta(j);
      for(ii_B = 0; ii_B<K; ii_B++){
      rx = res.t()*X.col(j);

      r(0) = R::rnorm( 0.0 , 1.0 );
      r(1) = R::rnorm( 0.0 , 1.0 );
      alpha_cand = alpha(j) + size_a*r(0);
      act = T_fun(alpha(j), alpha0, lam1, lam2, lam3,
                  act_type, SpSL_type);
      act_cand =   T_fun(alpha_cand, alpha0, lam1, lam2, lam3,
                         act_type, SpSL_type);
      cand = -0.5*log(X_norm(j)*act_cand*act_cand + sig0/eta) - 0.5*alpha_cand*alpha_cand + 0.5*rx*rx*act_cand*act_cand/(X_norm(j)*act_cand*act_cand + sig0/eta)/sig(0);
      curr = -0.5*log(X_norm(j)*act*act + sig0/eta) - 0.5*alpha(j)*alpha(j) + 0.5*rx*rx*act*act/(X_norm(j)*act*act + sig0/eta)/sig(0);
      aa = R::runif( 0.0, 1.0);
      cc = cand(0) - curr(0);
      if(cc > log(aa) ){
        alpha(j) = alpha_cand;
        theta(j) = theta_cand;
        act = act_cand;
        pq = pq + 1;
      }
      if(alpha(j)>alpha0){
        gam(j) = 1.0;
      }else{
        gam(j) = 0.0;
      }
      S = 1.0 / (X_norm(j)*act*act + sig0/eta);
      mu1 = rx(0)*act*S;
      r(0) = R::rnorm( 0.0 , 1.0 );
      w(j) = mu1 + sqrt(sig(0)*S)*r(0);
      theta(j) = w(j)*act;
      }
    res = res - X.col(j)*theta(j);
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
      OBJ(i-BURN) = obj(0);
      timer.step("");
      theta_save = theta_save + theta/N;
      gam_save = gam_save + gam/N;
      SAVE_theta.col(i-BURN) = theta;
      SAVE_w.col(i-BURN) = w;
      SAVE_alpha.col(i-BURN) = alpha;
      SIG(i-BURN) = sig(0);
      ETA(i-BURN) = eta;
      ALPHA0(i-BURN) = alpha0;
    }
    prob = 1.0-1.0/(1.0 + Act_alpha%Act_alpha);

    ind = arma::find(prob > 0.5);

    if(verbose == 1){
     if((i+1) % 2000 == 0){
       Rcpp::Rcout << "#iterations = " << i+1 << std::endl;
       //Rcpp::Rcout << "sig0 = " << sig0 << std::endl;
       Rcpp::Rcout << "Selected Variables = " << (ind+1).t() << std::endl;
     }
    }
  }
  return Rcpp::List::create(Rcpp::Named("CompTime") = timer,
                            //Rcpp::Named("X_norm") = X_norm,
                            Rcpp::Named("THETA") = SAVE_theta,
                            Rcpp::Named("ALPHA") = SAVE_alpha,
                            Rcpp::Named("ALPHA0") = ALPHA0,
                            Rcpp::Named("theta") = theta_save,
                            Rcpp::Named("gam") = gam_save,
                            Rcpp::Named("W") = SAVE_w,
                            Rcpp::Named("acc_a0") = acc_a0/(N+BURN),
                            Rcpp::Named("POST") = OBJ,
                            Rcpp::Named("sig") = SIG,
                            Rcpp::Named("eta") = ETA
                            );
}
