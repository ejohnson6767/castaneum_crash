
     functions {
     
  real f(real ntp1, real nt, real ft, real ovi_period, real mu_alpha, real beta0, real beta1, real eta_alpha, real eta_beta) {
  
  real beta = beta0 + nt * beta1;
  real focal_exponent = ft * ovi_period * beta;
  real mu_zt = mu_alpha * (1-exp(-focal_exponent)) / beta;
    
  vector[12] ls;
      ls[1] = log(2 * pow(mu_alpha,2) * pow(beta,2)) - 2 * focal_exponent;
      ls[2] = log(2 * pow(mu_alpha,2) * pow(beta,2)) + ft * ovi_period * (pow(eta_beta,2) - 2 * beta);
      ls[3] = log(pow(beta,3) * pow(eta_alpha,2));
      ls[4] = log(pow(beta,3) * pow(eta_alpha,2)) + ft * ovi_period * (pow(eta_beta,2) - 2 * beta);
      ls[5] = log( pow(mu_alpha,2) * beta * pow(eta_beta,2));
      ls[6] =  log(3 * pow(mu_alpha,2) * beta * pow(eta_beta,2))  - 2 * focal_exponent;
      ls[7] =  log(4 * pow(mu_alpha,2) * beta * pow(eta_beta,2))  - focal_exponent;
      ls[8] =  log( pow(beta,2) * pow(eta_alpha,2) * pow(eta_beta,2));
      ls[9] =  log( pow(beta,2) * pow(eta_alpha,2) * pow(eta_beta,2)) + ft * ovi_period * (pow(eta_beta,2) - 2 * beta);
      ls[10] = log(pow(mu_alpha,2) * pow(eta_beta,4));
      ls[11] = log(pow(mu_alpha,2) * pow(eta_beta,4)) - 2 * focal_exponent;
      ls[12] = log( 2* pow(mu_alpha,2) * pow(eta_beta,4)) - focal_exponent;

 real sigma_zt = sqrt( pow(pow(beta,2)*(2*pow(beta,2) - 3*beta*pow(eta_beta,2) + pow(eta_beta,4)),-1) *(-exp(ls[1]) +exp(ls[2]) +exp(ls[3]) -exp(ls[4]) +exp(ls[5]) +exp(ls[6]) -exp(ls[7]) -exp(ls[8]) +exp(ls[9]) -exp(ls[10]) -exp(ls[11]) +exp(ls[12])));

  real log_prob_nt_to_ft = normal_lpdf(ft | 0.5*nt, 0.5*sqrt(nt)) - normal_lccdf(0 | 0.5*nt, 0.5*sqrt(nt));
  real log_prob_ft_to_ntp1 = normal_lpdf(ntp1 | 0.91*mu_zt, 0.91*sigma_zt) - normal_lccdf(0 | 0.91*mu_zt, 0.91*sigma_zt);
  return(log_prob_nt_to_ft + log_prob_ft_to_ntp1);
  }
  
  real get_simp_add_on(int index, int M)
  {
  real res;
  if(index == 0 || index == M)
  {
  res = 0;
  } else if(index % 2 == 0)
  {
  res = 0.6931472;
  } else{
  res = 1.386294;
  }
  return(res);
  }
  
  real log_lik_Simpson_f(real ntp1, real nt, real ovi_period, real mu_alpha, real beta0, real beta1, real eta_alpha, real eta_beta, int M) {
    vector[M+1] lp; // M intervals implies M+1 partition points
    real h; // delta x, i.e. interval length
      real li = 0.5*nt - 4 * 0.5 * sqrt(nt);
      real ui = 0.5*nt + 4 * 0.5 * sqrt(nt);
      if(li < 0.5)
      {
        li = 0.5;
      }
    h = (ui - li)/M;
    
    for(m in 0:M)
  {
    lp[m + 1] = get_simp_add_on(m, M) + f(ntp1, nt, li + m * h, ovi_period, mu_alpha, beta0, beta1, eta_alpha, eta_beta);
  }

    return(log(h/3) + log_sum_exp(lp));
  }
  
     }
  
      data {
      int <lower = 1> N; // number of total data
      vector<lower = 0>[N] nt;
      vector<lower = 0>[N] ntp1;
      vector<lower = 0>[N] ovi_period;
      

      real <lower = 0> m_scale;
      real <lower = 0> beta0_loc;
      real <lower = 0> beta0_scale;
      real <lower = 0> beta1_loc;
      real <lower = 0> beta1_scale;
      real <lower = 0> eta_alpha_loc;
      real <lower = 0> eta_alpha_scale;
      real <lower = 0> eta_beta_loc;
      real <lower = 0> eta_beta_scale;
      }
      
      parameters {
      real <lower = 0> m_raw; 
      real <lower = 0> beta0_raw; // cannibalism rate intercept (per adult, per egg)
      real <lower = 0> beta1_raw; 
      real <lower = 0> eta_alpha_raw; // noise in ntp1
      real <lower = 0> eta_beta_raw; // noise in ntp1
      vector<lower = 0>[N] ft;
      }
      
      transformed parameters {
      // declarations
      real <lower = 0> m; // parameter m = mu_alpha / beta0. The parameters mu_alpha and beta0 are highly correlated in posterior sampling, leading to divergences. The reparameterization to m prevents this
      real <lower = 0> mu_alpha; // the population-level mean oviposition rate (per adult)
      vector<lower = 0>[N] beta; // cannibalism rate (per adult, per egg)
      real <lower = 0> beta0; 
      real <lower = 0> beta1; 
      real <lower = 0> eta_alpha; // noise in ntp1
      real <lower = 0> eta_beta;
      vector<lower = 0>[N] mu_zt;
      vector<lower = 0>[N] sigma_zt;


      // rescale non-centered parameters   
      m = m_scale * m_raw;
      beta0 = beta0_loc + beta0_scale * beta0_raw; 
      mu_alpha = m * beta0; 
      beta1 = beta1_loc + beta1_scale * beta1_raw; 
      eta_alpha = eta_alpha_loc + eta_alpha_scale * eta_alpha_raw;
      eta_beta = eta_beta_loc + eta_beta_scale * eta_beta_raw;

      // true transformed params
      beta = beta0 + nt * beta1;

      // larvae model
      vector[N] focal_exponent = ft.*ovi_period.*beta;
      
      mu_zt = mu_alpha * (1-exp(-focal_exponent)) ./ beta;
      
    for(i in 1:N)
    {
      vector[12] ls;
      ls[1] = log(2 * pow(mu_alpha,2) * pow(beta[i],2)) - 2 * focal_exponent[i];
      ls[2] = log(2 * pow(mu_alpha,2) * pow(beta[i],2)) + ft[i] * ovi_period[i] * (pow(eta_beta,2) - 2 * beta[i]);
      ls[3] = log(pow(beta[i],3) * pow(eta_alpha,2));
      ls[4] = log(pow(beta[i],3) * pow(eta_alpha,2)) + ft[i] * ovi_period[i] * (pow(eta_beta,2) - 2 * beta[i]);
      ls[5] = log( pow(mu_alpha,2) * beta[i] * pow(eta_beta,2));
      ls[6] =  log(3 * pow(mu_alpha,2) * beta[i] * pow(eta_beta,2))  - 2 * focal_exponent[i];
      ls[7] =  log(4 * pow(mu_alpha,2) * beta[i] * pow(eta_beta,2))  - focal_exponent[i];
      ls[8] =  log( pow(beta[i],2) * pow(eta_alpha,2) * pow(eta_beta,2));
      ls[9] =  log( pow(beta[i],2) * pow(eta_alpha,2) * pow(eta_beta,2)) + ft[i] * ovi_period[i] * (pow(eta_beta,2) - 2 * beta[i]);
      ls[10] = log(pow(mu_alpha,2) * pow(eta_beta,4));
      ls[11] = log(pow(mu_alpha,2) * pow(eta_beta,4)) - 2 * focal_exponent[i];
      ls[12] = log( 2* pow(mu_alpha,2) * pow(eta_beta,4)) - focal_exponent[i];

 sigma_zt[i] = sqrt( pow(pow(beta[i],2)*(2*pow(beta[i],2) - 3*beta[i]*pow(eta_beta,2) + pow(eta_beta,4)),-1) *(-exp(ls[1]) +exp(ls[2]) +exp(ls[3]) -exp(ls[4]) +exp(ls[5]) +exp(ls[6]) -exp(ls[7]) -exp(ls[8]) +exp(ls[9]) -exp(ls[10]) -exp(ls[11]) +exp(ls[12])  ));
    }

      }
      
      model {
      // priors  
      m_raw ~ std_normal();
      beta0_raw ~ std_normal();
      beta1_raw ~ std_normal();
      eta_alpha_raw ~ exponential(1);
      eta_beta_raw ~ exponential(1);
      ft ~ normal(0.5 * nt, 0.5 * sqrt(nt));

      // likelihood
      for(i in 1:N)
      {
        ntp1[i] ~ normal(0.91*mu_zt[i], 0.91*sigma_zt[i]) T[0,];
      }
      
      //ntp1 ~ normal(mu_zt, sigma_zt);
      
      }
      
       generated quantities{
       int M  = 100;
      vector[N] log_lik;
      for(i in 1:N)
      {
        log_lik[i] = log_lik_Simpson_f(ntp1[i], nt[i], ovi_period[i], mu_alpha, beta0, beta1, eta_alpha, eta_beta, M); 
      }
       }
      
      
