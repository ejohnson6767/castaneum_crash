

      data {
      int <lower = 1> N; // number of total data
      vector<lower = 0>[N] nt;
      vector<lower = 0>[N] ntp1;
      vector<lower = 0>[N] ovi_period;
      real<lower = 0, upper = 1> theta_L;
      

      real <lower = 0> m_scale;
      real <lower = 0> beta0_loc;
      real <lower = 0> beta0_scale;


      real <lower = 0> eta_alpha_loc;
      real <lower = 0> eta_alpha_scale;
      }
      
      transformed data{
      real <lower = 0> epsilon_L = sqrt(theta_L*(1-theta_L));
      }
      
      parameters {
      real <lower = 0> m_raw; 
      real <lower = 0> beta0_raw; // cannibalism rate intercept (per adult, per egg)

      real <lower = 0> eta_alpha_raw; // noise in ntp1
      vector<lower = 0>[N] ft;
      vector<lower = 0>[N] zt;

      }
      
      transformed parameters {
      // declarations
      real <lower = 0> m; // parameter m = mu_alpha / beta0. The parameters mu_alpha and beta0 are highly correlated in posterior sampling, leading to divergences. The reparameterization to m prevents this
      real <lower = 0> mu_alpha; // the population-level mean oviposition rate (per adult)

      real <lower = 0> beta0; 

      real <lower = 0> eta_alpha; // noise in ntp1
      vector<lower = 0>[N] mu_zt;
      vector<lower = 0>[N] sigma_zt;
      //vector<lower = 0>[N] ft;


      // rescale non-centered parameters   
      m = m_scale * m_raw;
      beta0 = beta0_loc + beta0_scale * beta0_raw; 
      mu_alpha = m * beta0; // the oviposition rate
      eta_alpha = eta_alpha_loc + eta_alpha_scale * eta_alpha_raw;

      // true transformed params

      //ft = nt ./ 2;
      // larvae model
      vector[N] focal_exponent = ft.*ovi_period*beta0;
      
      mu_zt = mu_alpha * (1-exp(-focal_exponent)) / beta0;
      sigma_zt = eta_alpha * sqrt( (1- exp(-2 * focal_exponent)) / (2 * beta0));
      
      }
      
      model {
      // priors  
      m_raw ~ std_normal();
      beta0_raw ~ std_normal();
      eta_alpha_raw ~ std_normal();
      //ft ~ normal(0.5 * nt, 0.5 * sqrt(nt));
      
      
      // likelihood
      for(i in 1:N)
      {
        ft[i] ~ normal(0.5 * nt[i], 0.5 * sqrt(nt[i]))T[0,];
        zt[i] ~ normal(mu_zt[i], sigma_zt[i]) T[0,];
        ntp1[i] ~ normal(theta_L * zt[i], epsilon_L * sqrt(zt[i])) T[0,];
      }
      

      }
      
      
