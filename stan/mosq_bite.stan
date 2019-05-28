data {

	int<lower=0> N;                       			// number of observations in proportion data
	int<lower=0> prop[N];         				// proportion data of length N
	int<lower=0> bites[N];         				// bite data of length N
	real prev_obs_b;					// prior on biting preference (conceivably could be previous observations)
	vector<lower=0>[N] flat_alpha_p;			// default flat prior
	vector<lower=0>[N] prev_obs_p;				// scaling from previous data

}

transformed data {
	vector<lower=0>[N] alpha_p;
	alpha_p = flat_alpha_p .* prev_obs_p;			// sets up strength of prior
}

parameters { 
	simplex[N] theta_p;	 				// vector of length = K dimensions that sums to 1
	vector<lower=0>[N] bite_pref;
} 

model { 
        bite_pref ~ gamma(prev_obs_b, prev_obs_b);		// prior on bite preference
	theta_p ~ dirichlet(alpha_p);				// prior on bird proportions
	prop ~ multinomial(theta_p);				// bird proportions
	bites ~ multinomial((to_vector(theta_p) .* bite_pref) /
		 sum(to_vector(theta_p) .* bite_pref));		// number of bites on each bird
} 
