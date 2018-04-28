function [ MA_matrices ] = wrap_BM_var_opt_11_3_gaussian( params, setup )
%function that maps a parameter vector into the objects needed to evaluate
%the likelihood function

%transforming parameters back to the original parameter space (NOT
%NECESSARY FOR IRFS, FINAL DRAWS ALREADY IN UN-TRANSFORMED VARIABLES)
%[ params ] = inv_transform( params,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );



current_matrices=params(1);



MA_matrices=zeros(1,1,setup.lags+1);
MA_matrices(:,:,1)=current_matrices;
j_matrix=kron(1:1:setup.lags,ones(1,1));



%parameters for negative shocks
a_neg=params(2);
b_neg=params(3);

c_neg=params(4);

a_neg2=params(5);
b_neg2=params(6);

c_neg2=params(7);

a_neg3=params(8);
b_neg3=params(9);

c_neg3=params(10);

%parameters for positive shocks
% a_pos=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);
% ind_for_loop=ind_for_loop+setup.size_obs^2;
% b_pos=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);
% ind_for_loop=ind_for_loop+setup.size_obs^2;
% c_pos=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);



MA_matrices(:,:,2:end)=a_neg.*exp((-(j_matrix-b_neg).^2)./c_neg)+a_neg2.*exp((-(j_matrix-b_neg2).^2)./c_neg2)...
    +a_neg3.*exp((-(j_matrix-b_neg3).^2)./c_neg3);


