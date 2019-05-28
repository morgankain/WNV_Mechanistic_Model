data {
    int<lower=0> N;                     		// Number of Data Points
    int<lower=0> N_IE;                  		// Number of Unique Lines
    int<lower=0> N_CIT;                 		// Number of Citations
    int<lower=0> N_VS;                  		// Number of Vector Species
    int<lower=1, upper=N_IE> IE[N];     		// Unique Line
    int<lower=1, upper=N_CIT> CIT[N];   		// Citation
    int<lower=1, upper=N_VS> VS[N];     		// Vector Species    
    real LD[N];                         		// Log Dose
    real<lower=0> Temp_X_Day[N];                        // Temp_Day Interaction
    real<lower=0> Temp[N];               	        // Temperature
    real<lower=0> Day[N];               	        // Day
    int<lower=0> Samp_Max;              		// Max number of samples
    int<lower=0> Trans_Max;                             // Max number transmitting
    int<lower=0, upper=Samp_Max> N_Samp[N];       	// Number of samples
    int<lower=0, upper=Trans_Max> N_Trans[N];     	// Number transmitting   
}

parameters {
 // Fixed Effects: Intercepts
    real alpha;       		  	    // Virus Lineage 1            
    
 // Fixed Effects: Slopes
    real beta_Day;          	   	    // Day      		     
    real beta_LD;                           // LD       
    real beta_Temp;                         // Temp      
    real beta_Temp_X_Day;	            // Temp x Day
  
 // Random Effects
    real alpha_IE_r[N_IE];		    // Infection Experiment
    real beta_Day_IE_r[N_IE];  
    real alpha_CIT_r[N_CIT];         	    // Citation     
    real beta_Day_CIT_r[N_CIT];             
    real alpha_VS_r[N_VS];                  // Vector Species
    real beta_Day_VS_r[N_VS];    
    real beta_LD_VS_r[N_VS]; 
    real beta_Temp_VS_r[N_VS]; 

 // Random Effect Varaiances
    real<lower=0> sigmasq_alpha_IE;
    real<lower=0> sigmasq_beta_Day_IE;
    real<lower=0> sigmasq_alpha_CIT;      
    real<lower=0> sigmasq_beta_Day_CIT;    
    real<lower=0> sigmasq_alpha_VS;      
    real<lower=0> sigmasq_beta_Day_VS;  
    real<lower=0> sigmasq_beta_LD_VS;
    real<lower=0> sigmasq_beta_Temp_VS;
    
}

transformed parameters {

 // Square root for convenience
    real<lower=0> sigma_alpha_IE;
    real<lower=0> sigma_beta_Day_IE;
    real<lower=0> sigma_alpha_CIT;       
    real<lower=0> sigma_beta_Day_CIT;    
    real<lower=0> sigma_alpha_VS;      
    real<lower=0> sigma_beta_Day_VS;  
    real<lower=0> sigma_beta_LD_VS;
    real<lower=0> sigma_beta_Temp_VS;
    
 // Linear Predictor
    vector[N] lin_pred;
    
 // 1 / Sqaure root
    sigma_alpha_IE = 1.0 / sqrt(sigmasq_alpha_IE);
    sigma_beta_Day_IE = 1.0 / sqrt(sigmasq_beta_Day_IE);
    sigma_alpha_CIT = 1.0 / sqrt(sigmasq_alpha_CIT);
    sigma_beta_Day_CIT = 1.0 / sqrt(sigmasq_beta_Day_CIT);
    sigma_alpha_VS = 1.0 / sqrt(sigmasq_alpha_VS);
    sigma_beta_Day_VS = 1.0 / sqrt(sigmasq_beta_Day_VS);
    sigma_beta_LD_VS = 1.0 / sqrt(sigmasq_beta_LD_VS);
    sigma_beta_Temp_VS = 1.0 / sqrt(sigmasq_beta_Temp_VS);
    
 // Linear Predictor
     for (i in 1:N) {
     lin_pred[i] = 
     alpha +       				  // intercept for each virus lineage
     alpha_IE_r[IE[i]] +   			  // random intercept for each Infection Experiment
     alpha_CIT_r[CIT[i]] +   			  // for each Citation
     alpha_VS_r[VS[i]] + 			  // for each Vector Species
     (beta_Day + 				  // slope over days
     beta_Day_IE_r[IE[i]] +		          // random component of day for each infection experiment
     beta_Day_CIT_r[CIT[i]] +			  // for each Citation
     beta_Day_VS_r[VS[i]]) * 		          // for each vector species
     Day[i] +  					  // across days
     (beta_LD +		 		          // Log Dose slope for each virus lineage
     beta_LD_VS_r[VS[i]]) *			  // random component of Log Dose for each Vector Species
     LD[i] +					  // across Log Dose
     (beta_Temp_VS_r[VS[i]] +			  // Random component of temperature for each Vector Species
     beta_Temp) *				  // Temperature slope for each Virus Lineage
     Temp[i] +					  // across Temperature
     beta_Temp_X_Day * Temp_X_Day[i];	          // Temp_Day Continuous Interaction

  }
}

model {
  
 // Intercept and Slope diffuse priors for Virus Lineage Parameters (Over Day and Titer)
    alpha ~ normal(0.0, 1.0E3);
    beta_Day ~ normal(0.0, 1.0E3);
    beta_LD ~ normal(0.0, 1.0E3);
    beta_Temp ~ normal(0.0, 1.0E3);
    beta_Temp_X_Day ~ normal(0.0, 1.0E3);
  
 // Diffuse Priors on Random Effect variance
    sigmasq_alpha_IE ~ gamma(1.0E-3, 1.0E-3);
    sigmasq_beta_Day_IE ~ gamma(1.0E-3, 1.0E-3);

    sigmasq_alpha_CIT ~ gamma(1.0E-3, 1.0E-3);
    sigmasq_beta_Day_CIT ~ gamma(1.0E-3, 1.0E-3);    

    sigmasq_alpha_VS ~ gamma(1.0E-3, 1.0E-3);      
    sigmasq_beta_Day_VS ~ gamma(1.0E-3, 1.0E-3);  
    sigmasq_beta_LD_VS ~ gamma(1.0E-3, 1.0E-3);
    sigmasq_beta_Temp_VS ~ gamma(1.0E-3, 1.0E-3);

  for (n in 1:N_IE) {
    alpha_IE_r[n] ~ normal(0, sigma_alpha_IE);
    beta_Day_IE_r[n] ~ normal(0, sigma_beta_Day_IE);
  }

  for (l in 1:N_CIT) {
    alpha_CIT_r[l] ~ normal(0, sigma_alpha_CIT);
    beta_Day_CIT_r[l] ~ normal(0, sigma_beta_Day_CIT);
  }

  for (m in 1:N_VS) {
    alpha_VS_r[m] ~ normal(0, sigma_alpha_VS);
    beta_Day_VS_r[m] ~ normal(0, sigma_beta_Day_VS);
    beta_LD_VS_r[m] ~ normal(0, sigma_beta_LD_VS);
    beta_Temp_VS_r[m] ~ normal(0, sigma_beta_Temp_VS);
  }

  for (i in 1:N) {
     N_Trans[i] ~ binomial_logit(N_Samp[i], lin_pred[i]);
  }
  
}
