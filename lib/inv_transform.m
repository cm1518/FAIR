function [ param_vector ] = inv_transform( param_vector,index_log, index_logit,index_logit_general, length_log,length_logit,length_logit_general,logit_general_lb,logit_general_ub )
%transforms unbounded parameter values into bounded transformed ones.
%if parameter value has to be positive, exp are taken
%if parameter is bounded between 0 and 1, an inverse logit transformation is used

if length_log>0
param_vector(index_log)=exp(param_vector(index_log));
end

if length_logit>0
param_vector(index_logit)=1./(1+exp(-param_vector(index_logit)));
end

if length_logit_general>0
 param_vector(index_logit_general)=1./(1+exp(-param_vector(index_logit_general)));
%  param_vector
%  param_vector(index_logit_general)
%  logit_general_ub - logit_general_lb
 param_vector(index_logit_general)=param_vector(index_logit_general).*(logit_general_ub - logit_general_lb)+logit_general_lb;
end


end

