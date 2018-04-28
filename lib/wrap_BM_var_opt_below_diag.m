function [ MA_matrices ] = wrap_BM_var_opt_below_diag( params, setup )
%function that maps a parameter vector into the objects needed to evaluate
%the likelihood function

%this version of the code is for the case of an IRF situated below the main
%diagonal



j_matrix=kron(0:1:setup.lags,ones(1,1));



%parameters for negative shocks
a_neg=params(1);
b_neg=params(2);

c_neg=params(3);





MA_matrices=zeros(1,1,setup.lags+1);

%parameters for positive shocks
% a_pos=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);
% ind_for_loop=ind_for_loop+setup.size_obs^2;
% b_pos=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);
% ind_for_loop=ind_for_loop+setup.size_obs^2;
% c_pos=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);

a_pos=a_neg;
b_pos=b_neg;
c_pos=c_neg;


MA_matrices(:,:,1:end)=a_neg.*exp((-(j_matrix-b_neg).^2)./c_neg);


