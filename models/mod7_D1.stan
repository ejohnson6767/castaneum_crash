
      data {
      int <lower = 1> N; // number of total data
      vector<lower = 0>[N] nt;
      vector<lower = 0>[N] ntp1;
      vector<lower = 0>[N] ovi_period;
      

      real <lower = 0> m_scale;
      real <lower = 0> beta0_loc;
      real <lower = 0> beta0_scale;


      real <lower = 0> sigma_alpha_loc;
      real <lower = 0> sigma_alpha_scale;
      }
      
      parameters {
      real <lower = 0> m_raw; 
      real <lower = 0> beta0_raw; // cannibalism rate intercept (per adult, per egg)

      real <lower = 0> sigma_alpha_raw; // noise in ntp1
      }
      
      transformed parameters {
      // declarations
      real <lower = 0> m; // parameter m = mu_alpha / beta0. The parameters mu_alpha and beta0 are highly correlated in posterior sampling, leading to divergences. The reparameterization to m prevents this
      real <lower = 0> mu_alpha; // the population-level mean oviposition rate (per adult)

      real <lower = 0> beta0; 

      real <lower = 0> sigma_alpha;
      vector<lower = 0>[N] mu_zt;
      vector<lower = 0>[N] sigma_zt;
      vector<lower = 0>[N] ft;


      // rescale non-centered parameters   
      m = m_scale * m_raw;
      beta0 = beta0_loc + beta0_scale * beta0_raw; 
      mu_alpha = m * beta0; // the oviposition rate
      sigma_alpha = sigma_alpha_loc + sigma_alpha_scale * sigma_alpha_raw;
      ft = nt ./ 2; 

      // true transformed params


      // larvae model
      vector[N] focal_exponent = ft.*ovi_period*beta0;
      
      mu_zt = mu_alpha * (1-exp(-focal_exponent)) / beta0;
      
      vector[N] t1 = log(2*pow(sigma_alpha, 2)) - focal_exponent;
      vector[N] t2 = log(pow(sigma_alpha, 2)) -2 * focal_exponent;
      real t3 = log(pow(sigma_alpha, 2));
    
      sigma_zt = sqrt((-exp(t1) + exp(t2) + exp(t3)) ./ (ft * pow(beta0,2)));

      }
      
      model {
      // priors  
      m_raw ~ std_normal();
      beta0_raw ~ std_normal();

      sigma_alpha_raw ~ exponential(1);

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
        log_lik[i] = normal_lpdf(ntp1[i] | 0.91*mu_zt[i], 0.91*sigma_zt[i]) - normal_lccdf(0 | 0.91*mu_zt[i], 0.91*sigma_zt[i]); 
      }
       }
      
      
