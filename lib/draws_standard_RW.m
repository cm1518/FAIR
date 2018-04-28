function [ acc_rate ,param_draws,log_posteriors,x_draws,add_matrices_draws] = draws_standard_RW( posterior_draw,cov_matrix,scaling,old_posterior,setup,data  )
%returns posterior draws for the standard MH algorithm with a RW propsal
%density

acceptances=zeros(setup.number_blocks,1);
param_draws=zeros(setup.length_param_vector,setup.number_of_draws/setup.keep_draw);
temp_draw=inv_transform( posterior_draw,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );   


if setup.add_matrices
add_matrices_draws=zeros(setup.dim_add_matrices(1),setup.dim_add_matrices(2),setup.number_of_draws/setup.keep_draw);
else
 add_matrices_draws=[];   
end




log_posteriors=zeros(1,setup.number_of_draws/setup.keep_draw);

x_draws=zeros(setup.state_size,setup.sample_size+1,setup.number_of_draws/setup.keep_draw);

    
%    if setup.likelihood==1
%       % ll_function=SSKF_wrap( params,setup,data );
%        posterior=@(params,setup,data) (prior(params,setup)+SSKF_wrap( params,setup,data ));
%        
%    elseif setup.likelihood==2
%        %ll_function=KF_wrap;
%         posterior=@(params,setup,data) (prior(params,setup)+KF_wrap( params,setup,data ));
%      
%    else
%   
%    end

[~, xdraw, add_matrix]=posteriorwithstates(posterior_draw,setup,data);

for nn=1:setup.number_of_draws
   if (nn/setup.disp_iter)==floor(nn/setup.disp_iter)
    nn
acceptances/nn
   end
   
   for bb=1:setup.number_blocks
   
param_prop=proposal_draw(posterior_draw(setup.index_block{bb}),cov_matrix{bb}, scaling(bb),setup);
param_prop=param_prop';
posterior_draw_prop=posterior_draw;
posterior_draw_prop(setup.index_block{bb})=param_prop;
[post_prop, x_prop, add_matrix_prop]=posteriorwithstates(posterior_draw_prop,setup,data);
%post_prop=posterior(param_prop,setup,data);
alpha=min(1,exp(post_prop-old_posterior+adjustment( posterior_draw,posterior_draw_prop,setup )));

if rand<alpha
   posterior_draw=posterior_draw_prop;
   old_posterior=post_prop;
   acceptances(bb)=acceptances(bb)+1;
   temp_draw=inv_transform( posterior_draw,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );
   xdraw=x_prop;
   add_matrix=add_matrix_prop;
end



   end
if ((nn)/setup.keep_draw)==floor(((nn)/setup.keep_draw))
param_draws(:,(nn)/setup.keep_draw)=temp_draw;
log_posteriors((nn)/setup.keep_draw)=old_posterior;
x_draws(:,:,(nn)/setup.keep_draw)=xdraw;



if setup.add_matrices
add_matrices_draws(:,:,(nn)/setup.keep_draw)=add_matrix;
end


end

end


acc_rate=acceptances/setup.number_of_draws;

end
