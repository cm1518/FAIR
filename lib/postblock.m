function [ lnp ] = postblock( params,setup,data,block, params_estimated )
%computes the posterior as a function of the parameters in one block,
%keeping the parameters in all other blocks fixed
params_total=params_estimated;
params_total(setup.index_block{block})=params;

lnp=prior(params_total,setup)+likelihood_wrap( params_total,setup,data );

end

