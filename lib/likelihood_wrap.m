function [lnpost]=likelihood_wrap( params,setup,data )

[ params ] = inv_transform( params,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );
try
[ lnpost, ~, ~, ~] = likelihood(data, params,setup );
catch ME
   lnpost=-1e100;
end
    
if isnan(lnpost)
    lnpost=-1e100;
end


end

