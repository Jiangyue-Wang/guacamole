
		data{
		
		int<lower=1> n;
		int<lower=1> max_locs;
		int<lower=1> n_locs[n];
		matrix[n, max_locs] time_step;
		vector<lower=0, upper=1>[n] delta;
		int<lower=1> n_knots;
		int cell_mat[n, max_locs];
		int ind_cell_effect;
		vector[2] beta_prior;
		vector[2] llambda_prior;
		vector[2] alpha_prior;
		int<lower=0> num_hab_covs; 
		matrix[n_knots, num_hab_covs] hab_cov; 
		int<lower=0> num_indv_covs; 		
		matrix[n, num_indv_covs] z;
		int ind_month;
		int n_months;
		vector[2] month_prior;
		int month_mat[n, max_locs];
	
		
		}
		
	
		parameters{
		
		real llambda;
		vector[num_hab_covs] alpha; 
		vector[num_indv_covs] beta; 
    vector[n_knots*ind_cell_effect] eta;
    vector[(n_months - 1)*ind_month] month_raw;

}
		
		transformed parameters{
		
		  vector[n_knots] cell_effect;
		  vector[n] nonspat;
		  real hab_spat;
		  vector[max_locs] haz;
		  vector[n] log_h;
		  vector[n] log_S;
		  vector[n] log_lik;
		  vector[1 + (n_months - 1)*ind_month] month = append_row(0, month_raw);
		  vector[n_months] month_effect;

		
		  cell_effect = rep_vector(0, n_knots);

	  if(num_indv_covs > 0){
		  nonspat = llambda + z * beta;
	  }else{
	    nonspat = rep_vector(llambda, n);
	  }
	  
	   for(mon in 1:n_months){
		    if(ind_month > 0){
		      		month_effect[mon] = month[mon];
		    }else{
		          month_effect[mon] = 0;
		    }
		 }
	  
		
		for(i in 1:n){
		  
		  haz = rep_vector(0, max_locs);
		  
		  for(j in 1:n_locs[i]){
		    
		    if(num_hab_covs > 0){
		      hab_spat = dot_product(alpha, hab_cov[cell_mat[i,j], 1:num_hab_covs]);
		    }else{
		      hab_spat = 0;
		    }
		    
		    haz[j] = exp(nonspat[i] + hab_spat + cell_effect[cell_mat[i,j]] + month_effect[month_mat[i, j]]); 
          
		  }
      
      log_S[i] = -dot_product(haz[1:n_locs[i]], time_step[i, 1:n_locs[i]]);
			log_h[i] = log(haz[n_locs[i]]);
		}
		
		log_lik = delta.*(log_h+log_S)	+ (1-delta).*log_S; 
}
		
		model{
		
		llambda ~ normal(llambda_prior[1], llambda_prior[2]);
		if(num_hab_covs > 0){
				alpha ~ normal(alpha_prior[1], alpha_prior[2]);
		}
		if(num_indv_covs > 0){
		    beta ~ normal(beta_prior[1], beta_prior[2]); 
		}
    if(ind_cell_effect == 1){
        eta ~ std_normal();
    }
    if(ind_month > 0){
        month_raw ~ normal(month_prior[1], month_prior[2]);
    }
    target += sum(log_lik); 
    
    
}
		

