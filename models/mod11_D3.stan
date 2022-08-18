

      data {
      int <lower = 1> N; // number of total data
      vector<lower = 0>[N] nt;
      vector<lower = 0>[N] ntp1;
      vector<lower = 0>[N] ovi_period;

      real <lower = 0> mu_alpha_loc;
      real <lower = 0> mu_alpha_scale;
      real <lower = 0> beta0_loc;
      real <lower = 0> beta0_scale;


      real <lower = 0> gamma_alpha_loc;
      real <lower = 0> gamma_alpha_scale;
      }

      parameters {
      real <lower = 0> mu_alpha_raw; // the population-level mean oviposition rate (per adult)
      real <lower = 0> beta0_raw; // cannibalism rate intercept (per adult, per egg)

      real <lower = 0> gamma_alpha_raw; // noise in ntp1
      }

      transformed parameters {
      // declarations
      real <lower = 0> mu_alpha; // the population-level mean oviposition rate (per adult)

      real <lower = 0> beta0;

      real <lower = 0> gamma_alpha; // noise in ntp1
      vector<lower = 0>[N] mu_ntp1;
      vector<lower = 0>[N] sigma_ntp1;
      vector<lower = 0>[N] ft;


      // rescale non-centered parameters
      mu_alpha = mu_alpha_loc + mu_alpha_scale * mu_alpha_raw;
      beta0 = beta0_loc + beta0_scale * beta0_raw;

      gamma_alpha = gamma_alpha_loc + gamma_alpha_scale * gamma_alpha_raw;

      // true transformed params

      ft = nt ./ 2;
      // larvae model
      vector[N] focal_exponent = ft.*ovi_period*beta0;
      
      vector[N] t1 = -2 * focal_exponent;
      vector[N] t2 = log(2) - focal_exponent;

      mu_ntp1 = mu_alpha * (1-exp(-focal_exponent)) / beta0;
      sigma_ntp1 = sqrt( pow(gamma_alpha / beta0, 2) * (1 + exp(t1) - exp(t2)));

      }

      model {
      // priors
      mu_alpha_raw ~ std_normal();
      beta0_raw ~ std_normal();

      gamma_alpha_raw ~ exponential(1);

      // likelihood
      for(i in 1:N)
      {
        ntp1[i] ~ normal(mu_ntp1[i], sigma_ntp1[i]) T[0,];
      }


      }

       generated quantities{
      vector[N] log_lik;
      for(i in 1:N)
      {
        log_lik[i] = normal_lpdf(ntp1[i] | mu_ntp1[i], sigma_ntp1[i]) - normal_lccdf(0 | mu_ntp1[i], sigma_ntp1[i]); 
      }
       }

      
