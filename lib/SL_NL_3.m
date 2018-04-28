function [ Sigma_lagged ] = SL_NL_3( beta_gen_neg,b_gen_neg,c_gen_neg,beta_gen_pos,b_gen_pos,c_gen_pos,epsilon_vec,lag,size_obs )
%computes one lagged Sigma matrix for the non-linear case
ind=epsilon_vec<0*ones(size(beta_gen_neg));



beta_gen=beta_gen_neg.*(ind)+beta_gen_pos.*((1-ind));

b_gen=b_gen_neg.*(ind)+b_gen_pos.*((1-ind));
c_gen=c_gen_neg.*(ind)+c_gen_pos.*((1-ind));


Sigma_lagged=exp(-((lag-b_gen).^2)./c_gen)...
    .*(beta_gen);



end
