function [ lnp x add_matrices] = posteriorwithstates( params,setup,data )
%function that returns log posterior and unobserved states
if setup.likelihood==1
      % ll_function=SSKF_wrap( params,setup,data );
      [lnpost x add_matrices]=SSKF_wrap_withstates( params,setup,data );
       lnp=(prior(params,setup)+lnpost);
       
   elseif setup.likelihood==2
       %ll_function=KF_wrap;
       [lnpost x add_matrices]=KF_wrap_withstates( params,setup,data );
        lnp=(prior(params,setup)+lnpost);
     
   elseif setup.likelihood==3
       %ll_function=KF_wrap;
       [lnpost x add_matrices]=likelihood_wrap_withstates( params,setup,data );
        lnp=(prior(params,setup)+lnpost);
  
   end

end

