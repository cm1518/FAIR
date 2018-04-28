function [ param_vector ] = transform( param_vector,index_log, index_logit,index_logit_general, length_log,length_logit,length_logit_general,logit_general_lb,logit_general_ub )
%transforms bounded parameter values into unbounded transformed ones.
%if parameter value has to be positive, logs are taken
%if parameter is bounded between 0 and 1, a logit transformation is used

if length_log>0
param_vector(index_log)=log(param_vector(index_log));
end

if length_logit>0
param_vector(index_logit)=log(param_vector(index_logit)./(1-param_vector(index_logit)));
end

if length_logit_general>0
 param_vector(index_logit_general)= (param_vector(index_logit_general) -logit_general_lb)./(logit_general_ub -logit_general_lb);  
 param_vector(index_logit_general)=log(param_vector(index_logit_general)./(1-param_vector(index_logit_general)));   
end

end

