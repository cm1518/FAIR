function [ Sigma] = unwrap_NL_lagged2(  beta_gen_neg,b_gen_neg,c_gen_neg,beta_gen_pos,b_gen_pos,c_gen_pos,epsilon_vec,setup )
%function returns intercept and array of  lagged Sigma matrices - setup.size_obs by setup.size_obs by setup.lag_length 
%with Sigma_0 ordered first
%parametrs are ordered as follows: [intercepts;alpha_diags;beta_diags;alpha_gen;beta_gen;b_gen;c_gen]
%epsilons are ordered from most recent to epsilon with largest lag


%unwrapping parameter vector


Sigma=zeros(setup.size_obs,setup.size_obs,setup.lags+1);


% Sigma(:,:,1)=tril(SL_NL( alpha_gen,beta_gen,b_gen,c_gen,epsilon_vec(:,1),0,setup.size_obs ),-1) ;
% 
% for kk=1:setup.size_obs
%    Sigma(kk,kk,1)=exp(alpha_diag(kk)*abs(epsilon_vec(kk,1))+beta_diag(kk)); 
%     
% end



% for jj=1:setup.lags
%     for nn=1:max(setup.num_gaussian)
%    Sigma(:,setup.index_unrestricted,jj+1)=SL_NL_2( alpha_gen_neg(:,nn),beta_gen_neg(:,nn),b_gen_neg(:,nn),c_gen_neg(:,nn),gamma_gen_neg(:,nn),alpha_gen_pos(:,nn),beta_gen_pos(:,nn),b_gen_pos(:,nn),c_gen_pos(:,nn),gamma_gen_pos(:,nn),epsilon_vec(setup.index_unrestricted,jj+1),jj,setup.size_obs,setup.threshold_vec(setup.index_unrestricted),indicator_vec(end-jj+1) ) +Sigma(:,setup.index_unrestricted,jj+1);
%     end
% end
for oo=1:setup.size_obs
for jj=1:setup.lags
    
   Sigma(:,oo,jj+1)=sum(SL_NL_3(beta_gen_neg(:,:,oo),b_gen_neg(:,:,oo),c_gen_neg(:,:,oo),beta_gen_pos(:,:,oo),b_gen_pos(:,:,oo),c_gen_pos(:,:,oo),epsilon_vec(setup.index_unrestricted,jj+1),jj,setup.size_obs),2);
   
end
end

Sigma=Sigma(:,:,2:end);

end

