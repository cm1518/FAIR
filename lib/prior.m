function [ lnp ] = prior( param,setup )
%computes log prior for normal and gamma prior distributions
[ param ] = inv_transform( param,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );

if ~isempty(setup.index_normal)
log_normal=sum(log(normpdf(param(setup.index_normal),setup.normal_prior_means,setup.normal_prior_std)));
else
log_normal=0;
end


if ~isempty(setup.index_gamma)
log_gamma=sum(log(pdf('Gamma',param(setup.index_gamma),setup.gamma_prior_shape,setup.gamma_prior_scale)));
else
log_gamma=0;
end


lnp=log_normal+log_gamma;


end

