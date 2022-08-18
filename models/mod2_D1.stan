
     functions {
     
  real f(real ntp1, real nt, real ft, real ovi_period, real mu_alpha, real beta0,  real eta_alpha, real sigma_alpha) {
  

  real focal_exponent = ft * ovi_period * beta0;
  real mu_zt = mu_alpha * (1-exp(-focal_exponent)) / beta0;
      real t1 = log(ft * beta0 *  pow(eta_alpha,2));
      real t2 = log(ft * beta0 *  pow(eta_alpha,2)) -2 * focal_exponent;
      real t3 = log(2 * pow(sigma_alpha, 2));
      real t4 = log(2 * pow(sigma_alpha, 2)) - 2* focal_exponent;
      real t5 = log(4 * pow(sigma_alpha,2)) - focal_exponent;
    
      real sigma_zt = sqrt((exp(t1) - exp(t2) + exp(t3) + exp(t4) - exp(t5)) ./ (2 * ft * pow(beta0,2)));
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
  
  real log_lik_Simpson_f(real ntp1, real nt, real ovi_period, real mu_alpha, real beta0,  real eta_alpha, real sigma_alpha, int M) {
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
    lp[m + 1] = get_simp_add_on(m, M) + f(ntp1, nt, li + m * h, ovi_period, mu_alpha, beta0,  eta_alpha, sigma_alpha);
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


      real <lower = 0> eta_alpha_loc;
      real <lower = 0> eta_alpha_scale;
      real <lower = 0> sigma_alpha_loc;
      real <lower = 0> sigma_alpha_scale;
      }
      
      parameters {
      real <lower = 0> m_raw; 
      real <lower = 0> beta0_raw; // cannibalism rate intercept (per adult, per egg)

      real <lower = 0> eta_alpha_raw; // noise in ntp1
      real <lower = 0> sigma_alpha_raw; // noise in ntp1
      vector<lower = 0>[N] ft;
      vector<lower = 0>[N] alpha_raw;
      }
      
      transformed parameters {
      // declarations
      real <lower = 0> m; // parameter m = mu_alpha / beta0. The parameters mu_alpha and beta0 are highly correlated in posterior sampling, leading to divergences. The reparameterization to m prevents this
      real <lower = 0> mu_alpha; // the population-level mean oviposition rate (per adult)

      real <lower = 0> beta0; 

      real <lower = 0> eta_alpha; // noise in ntp1
      real <lower = 0> sigma_alpha;
      vector<lower = 0>[N] mu_zt;
      vector<lower = 0>[N] sigma_zt;
      vector<lower = 0>[N] alpha;


      // rescale non-centered parameters   
      m = m_scale * m_raw;
      beta0 = beta0_loc + beta0_scale * beta0_raw; 
      mu_alpha = m * beta0; // the oviposition rate
      eta_alpha = eta_alpha_loc + eta_alpha_scale * eta_alpha_raw;
      sigma_alpha = sigma_alpha_loc + sigma_alpha_scale * sigma_alpha_raw;

      // true transformed params

      alpha = mu_alpha + (sigma_alpha./sqrt(ft)) .* alpha_raw; 

      // larvae model
      vector[N] focal_exponent = ft.*ovi_period*beta0;
      
      mu_zt = alpha .* (1-exp(-focal_exponent)) / beta0;
      
      vector[N] t1 = log(ft * beta0 *  pow(eta_alpha,2));
      vector[N] t2 = log(ft * beta0 *  pow(eta_alpha,2)) -2 * focal_exponent;
      real t3 = log(2 * pow(sigma_alpha, 2));
      vector[N] t4 = log(2 * pow(sigma_alpha, 2)) - 2* focal_exponent;
      vector[N] t5 = log(4 * pow(sigma_alpha,2)) - focal_exponent;
    
      sigma_zt = sqrt((exp(t1) - exp(t2) + exp(t3) + exp(t4) - exp(t5)) ./ (2 * ft * pow(beta0,2)));

      }
      
      model {
      // priors  
      m_raw ~ std_normal();
      beta0_raw ~ std_normal();

      eta_alpha_raw ~ exponential(1);
      sigma_alpha_raw ~ exponential(1);
      ft ~ normal(0.5 * nt, 0.5 * sqrt(nt));
      alpha_raw ~ std_normal();

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
        log_lik[i] = log_lik_Simpson_f(ntp1[i], nt[i], ovi_period[i], mu_alpha, beta0,  eta_alpha, sigma_alpha, M); 
      }
       }
      
      
