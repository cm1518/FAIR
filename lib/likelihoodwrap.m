function [lnpost x add_matrices]=likelihood_wrap_withstates( params,setup,data )
[ param ] = inv_transform( params,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );   


current_matrices=zeros(1,1,2);
current_matrices(1,1,1)=param(2);
current_matrices(1,1,1)=param(3);
[ lnpost indicator x uvec] = likelihood(data, param(1), current_matrices,[],[], setup.lags );
add_matrices=indicator;

end

