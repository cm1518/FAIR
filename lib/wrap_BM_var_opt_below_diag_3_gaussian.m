function [ MA_matrices ] = wrap_BM_var_opt_below_diag_3_gaussian( params, setup )
%function that maps a parameter vector into the objects needed to evaluate
%the likelihood function

%this version of the code is for the case of an IRF situated below the main
%diagonal



j_matrix=kron(0:1:setup.lags,ones(1,1));



%parameters for negative shocks
a_neg=params(1);
b_neg=params(2);

c_neg=params(3);

a_neg2=params(4);
b_neg2=params(5);

c_neg2=params(6);

a_neg3=params(7);
b_neg3=params(8);

c_neg3=params(9);

MA_matrices=zeros(1,1,setup.lags+1);


MA_matrices(:,:,1:end)=a_neg.*exp((-(j_matrix-b_neg).^2)./c_neg)+a_neg2.*exp((-(j_matrix-b_neg2).^2)./c_neg2)...
    +a_neg3.*exp((-(j_matrix-b_neg3).^2)./c_neg3);



