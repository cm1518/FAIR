function [ Sigma_lagged ] = SL_NL_symmetric( beta_gen_neg,b_gen_neg,c_gen_neg,lag )
%computes one lagged Sigma matrix 




beta_gen=beta_gen_neg;

b_gen=b_gen_neg;
c_gen=c_gen_neg;


Sigma_lagged=exp(-((lag-b_gen).^2)./c_gen)...
    .*(beta_gen);



end
