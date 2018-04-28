function [ params_new ] = param_zero_restriction( params,setup )
%function takes restricted parameter vector (missing parameter values that are fixed to be zero) and pads it with zeros in the
%correct places to return a parameter vector that is of full length

params_new=zeros(setup.length_param_vector,1);
params_new(setup.index_non_zero)=params;

end

