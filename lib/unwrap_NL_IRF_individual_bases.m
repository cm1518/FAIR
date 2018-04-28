function [ Sigma] = unwrap_NL_IRF_individual_bases( params,epsilon_vec,setup)
%function returns individual Gaussian basis functions for inspection

%unwrapping parameters
[ params ] = params_mod( params,setup );







for oo=1:setup.size_obs
for jj=1:setup.lags
    
   Sigma(:,:,oo,jj+1)=(SL_NL_3(params.beta_gen_neg(:,:,oo),params.b_gen_neg(:,:,oo),params.c_gen_neg(:,:,oo),params.beta_gen_pos(:,:,oo),params.b_gen_pos(:,:,oo),params.c_gen_pos(:,:,oo),epsilon_vec(setup.index_unrestricted,jj+1),jj,setup.size_obs));
   
end
end


end

