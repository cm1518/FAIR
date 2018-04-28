function [ draws, acc_rate, log_posteriors, statedraws, add_matrices] = sampling_MH_MC_GB( setup )
%main function for the MH algorithm
data=load(setup.data);
data=data.data;
data=data(:,end-setup.sample_size+1:end);



warning off;

options = optimset('Display','iter','TolFun',1e-4,'TolX',1e-4,'MaxIter',1000, 'MaxFunEvals',50000);

if setup.proposal==1
   temp=setup.scaling_draws/setup.check_scaling;
   %pr_function=str2func(setup.prior_function);
   
   if setup.likelihood==1
      % ll_function=SSKF_wrap( params,setup,data );
       post=@(params,setup,data) (prior(params,setup)+SSKF_wrap( params,setup,data ));
       postopt=@(params,setup,data) -(prior(params,setup)+SSKF_wrap( params,setup,data ));
       postopt2=@(params) -(prior(params,setup)+SSKF_wrap( params,setup,data ));
       %since we use a minimizer
   elseif setup.likelihood==2
       %ll_function=KF_wrap;
        post=@(params,setup,data) (prior(params,setup)+KF_wrap( params,setup,data ));
       postopt=@(params,setup,data) -(prior(params,setup)+KF_wrap( params,setup,data )); %since we use a minimizer
       postopt2=@(params) -(prior(params,setup)+KF_wrap( params,setup,data ));
   elseif setup.likelihood==3
    post=@(params,setup,data) (prior(params,setup)+likelihood_wrap( params,setup,data ));
       postopt=@(params,setup,data) -(prior(params,setup)+likelihood_wrap( params,setup,data )); %since we use a minimizer
       postopt2=@(params) -(prior(params,setup)+likelihood_wrap( params,setup,data ));
        postopt3=@(params) -(prior(param_zero_restriction(params,setup),setup)+likelihood_wrap( param_zero_restriction(params,setup),setup,data ));
   end
if setup.skip_opt==0
   %initial minimization
   [ x0 ] = transform( setup.initial_parameter,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub  );
 
   disp('finding starting value using simplex optimizer')
   [xestimate1,functionvalue1]=fminsearch(postopt3,x0(setup.index_non_zero),options);
  disp('finding starting value using derivative-based optimizer')
 [xestimate1,functionvalue1,exitflag,output,grad,hessian]=fminunc(postopt3,xestimate1,options);
   
 
       xestimate=param_zero_restriction(xestimate1,setup);
       functionvalue=functionvalue1;
 
 save temp_max xestimate
    elseif setup.skip_opt==1
load temp_max

functionvalue=postopt2(xestimate);
end
  
  
    

 for bb=1:setup.number_blocks  
     try
        Hess_temp = hessian3(func2str(@postblock),xestimate(setup.index_block{bb}),setup,data,bb,xestimate);

   Hess_temp=reshape(-Hess_temp,length(setup.index_block{bb}),length(setup.index_block{bb}));
  Hess_temp=eye(size(Hess_temp,1))/Hess_temp; 
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   %making sure the Hessian is positive definite
  [V,D]=eig(Hess_temp);

   d=diag(D);
   d(d<=0)=eps;

   Hess_temp= V*diag(d)*V';
   Hess{bb}=Hess_temp;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     catch ME
         disp('Hessian computation in the following block did not work:')
         bb
Hess{bb}=.5*eye(length(setup.index_block{bb}));
     end
end
  %end
   scaling=ones(setup.number_blocks,1).*setup.initial_scaling;
   old_posterior=-functionvalue;
   posterior_draw=xestimate;
   acc_rate=zeros(setup.number_blocks,1);
   for jj=1:temp %loop to adjust scaling matrix
       disp('finding scaling for propsal densities.... progress:')
       jj/temp
    for bb=1:setup.number_blocks
       [ acc_rate(bb) ,posterior_draw,old_posterior] = acceptance_rate( posterior_draw,Hess{bb},scaling(bb),setup.check_scaling,old_posterior,setup,data,bb  );  
     
     if acc_rate(bb)<.2
       scaling(bb)=.75*scaling(bb);
       elseif acc_rate(bb)>.5
           scaling(bb)=1.5*scaling(bb);
     end
    end 
   end
    if jj==temp
    disp('scaling:')
    scaling
    disp('acceptance rate:')
    acc_rate
    end
    
   %actual draws
    
    [ acc_rate ,draws,log_posteriors, statedraws, add_matrices] = draws_standard_RW( posterior_draw,Hess,scaling,old_posterior,setup,data  );
    
 else %adaptive MH
     
     disp('this code only uses the standard MH algorithm for now')

end



end

