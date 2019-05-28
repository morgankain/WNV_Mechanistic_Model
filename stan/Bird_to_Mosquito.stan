data {
    int<lower=0> N;                     	// Number of Data Points
    int<lower=0> N_VS;                  	// Number of Vector Species
    int<lower=0> N_CIT;                 	// Number of Citations
    int<lower=1, upper=N_CIT> CIT[N];   	// Citation
    int<lower=1, upper=N_VS> VS[N];     	// Vector Species
    real LD[N];                         	// Log Dose
    real Temp[N];				// Temperature
    int<lower=0> Samp_Max;              	// Max number of samples
    int<lower=0> Inf_Max;               	// Max number of infected
    int<lower=0, upper=Samp_Max> N_Samp[N];     // Number of samples
    int<lower=0, upper=Inf_Max> N_Inf[N];       // Number Infected
}

parameters {  
  // Fixed Effects: Slopes
     real beta_LD;				// slope across LD
     real beta_Temp;				// slope across temp            
  
  // Random Effects
    real alpha_CIT_r[N_CIT];      		// Citation       
    real beta_LD_CIT_r[N_CIT];                 
    real alpha_VS_r[N_VS];                 	// Vector Species
    real beta_LD_VS_r[N_VS];         
  
  // Random Effect Varaiances
    real<lower=0> sigmasq_alpha_CIT;       
    real<lower=0> sigmasq_beta_LD_CIT;   
    real<lower=0> sigmasq_alpha_VS;      
    real<lower=0> sigmasq_beta_LD_VS;   
    
}

transformed parameters {
  // Square root for convenience
    real<lower=0> sigma_alpha_CIT;       
    real<lower=0> sigma_beta_LD_CIT;
    real<lower=0> sigma_alpha_VS;       
    real<lower=0> sigma_beta_LD_VS;
    
  // linear predictor
    vector[N] lin_pred;
    
  // 1 / Sqaure root
    sigma_alpha_CIT = 1.0 / sqrt(sigmasq_alpha_CIT);
    sigma_beta_LD_CIT = 1.0 / sqrt(sigmasq_beta_LD_CIT);
    sigma_alpha_VS = 1.0 / sqrt(sigmasq_alpha_VS);
    sigma_beta_LD_VS = 1.0 / sqrt(sigmasq_beta_LD_VS);
    
    for (i in 1:N){
    lin_pred[i] = 	
     alpha_CIT_r[CIT[i]] +   		        // random intercept for each Citation
     alpha_VS_r[VS[i]] + 		        // random intercept for each Vector Species  
     (beta_LD +		   		        // base slope
     beta_LD_CIT_r[CIT[i]] + 		        // linear component of Log Dose for each Citation
     beta_LD_VS_r[VS[i]]) * 		        // linear component of Log Dose for each Vector Species
     LD[i] + 				        // across Log Dose
     beta_Temp * Temp[i];       	        // temp
    }
}

model {
   beta_LD ~ normal(0.0, 1.0E3);                // diffuse priors fixed effects
   beta_Temp ~ normal(0.0, 1.0E3);
         
   sigmasq_alpha_CIT ~ gamma(1.0E-3, 1.0E-3);	// Citation
   sigmasq_beta_LD_CIT ~ gamma(1.0E-3, 1.0E-3);
   sigmasq_alpha_VS ~ gamma(1.0E-3, 1.0E-3);	// Vector Species
   sigmasq_beta_LD_VS ~ gamma(1.0E-3, 1.0E-3);
   
  for (l in 1:N_CIT) {
    alpha_CIT_r[l] ~ normal(0, sigma_alpha_CIT);
    beta_LD_CIT_r[l] ~ normal(0, sigma_beta_LD_CIT);
  }
  
  for (m in 1:N_VS) {
    alpha_VS_r[m] ~ normal(0, sigma_alpha_VS);
    beta_LD_VS_r[m] ~ normal(0, sigma_beta_LD_VS);
  }

  for (i in 1:N) {
     N_Inf[i] ~ binomial_logit(N_Samp[i], lin_pred[i]);
  }
  
}
