function [ Sigma, intercept] = unwrap_NL_IRF( params,epsilon_vec,setup)
%function returns intercept and array of Sigma matrices - setup.size_obs by setup.size_obs by setup.lag_length 
%with Sigma_0 ordered first
%parametrs are ordered as follows: [intercepts;alpha_diags;beta_diags;alpha_gen;beta_gen;b_gen;c_gen]
%epsilons are ordered from most recent to epsilon with largest lag
[ params ] = params_mod( params,setup );

%unwrapping parameters


Sigma=zeros(setup.size_obs,setup.size_obs,setup.lags+1);


Sigma(:,:,1)=setup.store_responses(:,:,1) ;



if epsilon_vec(setup.index_unrestricted,1)<0
Sigma(:,:,1)=params.beta_diag_neg;
elseif epsilon_vec(setup.index_unrestricted,1)>0
Sigma(:,:,1)=params.beta_diag_neg;
Sigma(1:end,setup.index_unrestricted,1)=params.beta_diag_pos;

end





for oo=1:setup.size_obs
for jj=1:setup.lags
    
   Sigma(:,oo,jj+1)=sum(SL_NL_3(params.beta_gen_neg(:,:,oo),params.b_gen_neg(:,:,oo),params.c_gen_neg(:,:,oo),params.beta_gen_pos(:,:,oo),params.b_gen_pos(:,:,oo),params.c_gen_pos(:,:,oo),epsilon_vec(setup.index_unrestricted,jj+1),jj,setup.size_obs),2);
   
end
end

intercept=params.intercepts;
end

