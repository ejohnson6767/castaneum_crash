
     functions {
     
  real f(real ntp1, real nt, real ft, real ovi_period, real mu_alpha, real beta0,  real sigma_alpha) {
  

  real focal_exponent = ft * ovi_period * beta0;
  real mu_ntp1 = mu_alpha * (1-exp(-focal_exponent)) / beta0;
      real t1 = log(2*pow(sigma_alpha, 2)) - focal_exponent;
      real t2 = log(pow(sigma_alpha, 2)) -2 * focal_exponent;
      real t3 = log(pow(sigma_alpha, 2));
    
      real sigma_ntp1 = sqrt((-exp(t1) + exp(t2) + exp(t3)) / (ft * pow(beta0,2)));
  real log_prob_nt_to_ft = normal_lpdf(ft | 0.5*nt, 0.5*sqrt(nt)) - normal_lccdf(0 | 0.5*nt, 0.5*sqrt(nt));
  real log_prob_ft_to_ntp1 = normal_lpdf(ntp1 | mu_ntp1, sigma_ntp1) - normal_lccdf(0 | mu_ntp1, sigma_ntp1);
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
  
  real log_lik_Simpson_f(real ntp1, real nt, real ovi_period, real mu_alpha, real beta0,  real sigma_alpha, int M) {
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
    lp[m + 1] = get_simp_add_on(m, M) + f(ntp1, nt, li + m * h, ovi_period, mu_alpha, beta0,  sigma_alpha);
  }

    return(log(h/3) + log_sum_exp(lp));
  }
  
     }
  
      data {
      int <lower = 1> N; // number of total data
      vector<lower = 0>[N] nt;
      vector<lower = 0>[N] ntp1;
      vector<lower = 0>[N] ovi_period;
      
      real <lower = 0> mu_alpha_loc;
      real <lower = 0> mu_alpha_scale;
      real <lower = 0> beta0_loc;
      real <lower = 0> beta0_scale;


      real <lower = 0> sigma_alpha_loc;
      real <lower = 0> sigma_alpha_scale;
      }
      
      parameters {
      real <lower = 0> mu_alpha_raw; // the population-level mean oviposition rate (per adult)
      real <lower = 0> beta0_raw; // cannibalism rate intercept (per adult, per egg)

      real <lower = 0> sigma_alpha_raw; // noise in ntp1
      vector<lower = 0>[N] ft;
      }
      
      transformed parameters {
      // declarations
      real <lower = 0> mu_alpha; // the population-level mean oviposition rate (per adult)

      real <lower = 0> beta0; 

      real <lower = 0> sigma_alpha;
      vector<lower = 0>[N] mu_ntp1;
      vector<lower = 0>[N] sigma_ntp1;


      // rescale non-centered parameters   
      mu_alpha = mu_alpha_loc + mu_alpha_scale * mu_alpha_raw;
      beta0 = beta0_loc + beta0_scale * beta0_raw; 

      sigma_alpha = sigma_alpha_loc + sigma_alpha_scale * sigma_alpha_raw;
      //ft = nt ./ 2; 

      // true transformed params


      // larvae model
      vector[N] focal_exponent = ft.*ovi_period*beta0;
      
      mu_ntp1 = mu_alpha * (1-exp(-focal_exponent)) / beta0;
      
      vector[N] t1 = log(2*pow(sigma_alpha, 2)) - focal_exponent;
      vector[N] t2 = log(pow(sigma_alpha, 2)) -2 * focal_exponent;
      real t3 = log(pow(sigma_alpha, 2));
    
      sigma_ntp1 = sqrt((-exp(t1) + exp(t2) + exp(t3)) ./ (ft * pow(beta0,2)));

      }
      
      model {
      // priors  
      mu_alpha_raw ~ std_normal();
      beta0_raw ~ std_normal();

      sigma_alpha_raw ~ exponential(1);
      ft ~ normal(0.5 * nt, 0.5 * sqrt(nt));
      
      // likelihood
      for(i in 1:N)
      {
        ntp1[i] ~ normal(mu_ntp1[i], sigma_ntp1[i]) T[0,];
      }
      
      //ntp1 ~ normal(mu_ntp1, sigma_ntp1);
      
      }
      
       generated quantities{
       int M  = 100;
      vector[N] log_lik;
      for(i in 1:N)
      {
        log_lik[i] = log_lik_Simpson_f(ntp1[i], nt[i], ovi_period[i], mu_alpha, beta0,  sigma_alpha, M); 
      }
       }
      
      
