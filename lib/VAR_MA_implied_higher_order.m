function [ A ] = VAR_MA_implied_higher_order( setup,j )
%calculates the MA coefficient of order j implied by a higher order VAR
%(coefficients for which are part of setup)

temp=[setup.VARsymA;kron(eye(setup.VARsym_order-1),eye(setup.size_obs)) zeros(setup.size_obs*(setup.VARsym_order-1),setup.size_obs)];
%companion form of VAR
%eig(temp)
A=temp^(j);
A=A(1:setup.size_obs,1:setup.size_obs)*setup.VARsymchol;
end

