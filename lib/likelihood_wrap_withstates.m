function [lnpost, x, add_matrices]=likelihood_wrap_withstates( params,setup,data )


[ params ] = inv_transform( params,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );

try
[ lnpost, ~, x, ~] = likelihood(data, params,setup );
add_matrices=x; %epsilon shocks
x=0;
catch ME
  lnpost=-1e100;
  x=0;
  add_matrices=zeros(setup.size_obs,setup.lags+size(data,2));
end
if isnan(lnpost)
    lnpost=-1e100;
end
end

