

  
      data {
      int <lower = 1> N; // number of total data
      vector<lower = 0>[N] nt;
      vector<lower = 0>[N] ntp1;
      vector<lower = 0>[N] ovi_period;
      
      real <lower = 0> mu_alpha_loc;
      real <lower = 0> mu_alpha_scale;
      real <lower = 0> beta0_loc;
      real <lower = 0> beta0_scale;


      real <lower = 0> eta_alpha_loc;
      real <lower = 0> eta_alpha_scale;
      real <lower = 0> sigma_alpha_loc;
      real <lower = 0> sigma_alpha_scale;
      }
      
      parameters {
      real <lower = 0> mu_alpha_raw; // the population-level mean oviposition rate (per adult)
      real <lower = 0> beta0_raw; // cannibalism rate intercept (per adult, per egg)

      real <lower = 0> eta_alpha_raw; // noise in ntp1
      real <lower = 0> sigma_alpha_raw; // noise in ntp1
      vector<lower = 0>[N] alpha_raw;
      }
      
      transformed parameters {
      // declarations
      real <lower = 0> mu_alpha; // the population-level mean oviposition rate (per adult)

      real <lower = 0> beta0; 

      real <lower = 0> eta_alpha; // noise in ntp1
      real <lower = 0> sigma_alpha;
      vector<lower = 0>[N] mu_ntp1;
      vector<lower = 0>[N] sigma_ntp1;
      vector<lower = 0>[N] alpha;
      vector<lower = 0>[N] ft;

      // rescale non-centered parameters   
      mu_alpha = mu_alpha_loc + mu_alpha_scale * mu_alpha_raw;
      beta0 = beta0_loc + beta0_scale * beta0_raw; 

      eta_alpha = eta_alpha_loc + eta_alpha_scale * eta_alpha_raw;
      sigma_alpha = sigma_alpha_loc + sigma_alpha_scale * sigma_alpha_raw;

      // true transformed params

      ft = nt ./ 2;
      alpha = mu_alpha + (sigma_alpha./sqrt(ft)) .* alpha_raw; 

      // larvae model
      vector[N] focal_exponent = ft.*ovi_period*beta0;
      
      mu_ntp1 = alpha .* (1-exp(-focal_exponent)) / beta0;
      
      vector[N] t1 = log(ft * beta0 *  pow(eta_alpha,2));
      vector[N] t2 = log(ft * beta0 *  pow(eta_alpha,2)) -2 * focal_exponent;
      real t3 = log(2 * pow(sigma_alpha, 2));
      vector[N] t4 = log(2 * pow(sigma_alpha, 2)) - 2* focal_exponent;
      vector[N] t5 = log(4 * pow(sigma_alpha,2)) - focal_exponent;
    
      sigma_ntp1 = sqrt((exp(t1) - exp(t2) + exp(t3) + exp(t4) - exp(t5)) ./ (2 * ft * pow(beta0,2)));

      }
      
      model {
      // priors  
      mu_alpha_raw ~ std_normal();
      beta0_raw ~ std_normal();

      eta_alpha_raw ~ exponential(1);
      sigma_alpha_raw ~ exponential(1);
      alpha_raw ~ std_normal();

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
        log_lik[i] = normal_lpdf(ntp1[i] | mu_ntp1[i], sigma_ntp1[i]) - normal_lccdf(0 | mu_ntp1[i], sigma_ntp1[i]); 
      }
       }
      
      