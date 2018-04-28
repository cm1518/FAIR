function [ acc_rate ,last_param_draw,new_posterior] = acceptance_rate( posterior_draw,cov_matrix,scaling,number_of_draws,old_posterior,setup, data,block )
%calculates the acceptance rate of running the standard MH algorithm for
%number_of_draws
acceptances=0;


   
   if setup.likelihood==1
      % ll_function=SSKF_wrap( params,setup,data );
       posterior=@(params,setup,data) (prior(params,setup)+SSKF_wrap( params,setup,data ));
       
   elseif setup.likelihood==2
       %ll_function=KF_wrap;
        posterior=@(params,setup,data) (prior(params,setup)+KF_wrap( params,setup,data ));
     
   elseif setup.likelihood==3
       %ll_function=KF_wrap;
        posterior=@(params,setup,data) (prior(params,setup)+likelihood_wrap( params,setup,data ));
  
   end

for nn=1:number_of_draws
  
param_prop=proposal_draw(posterior_draw(setup.index_block{block}),cov_matrix, scaling,setup);
param_prop=param_prop';
posterior_draw_prop=posterior_draw;
posterior_draw_prop(setup.index_block{block})=param_prop;

post_prop=posterior(posterior_draw_prop,setup,data);
alpha=min(1,exp(post_prop-old_posterior+adjustment( posterior_draw,posterior_draw_prop,setup )));
if rand<alpha
   posterior_draw=posterior_draw_prop;
   old_posterior=post_prop;
   acceptances=acceptances+1;
end

end


acc_rate=acceptances/number_of_draws;
last_param_draw=posterior_draw;
new_posterior=old_posterior;
end

