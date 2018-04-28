function [ mod_params ] = params_mod( params,setup )
%transforms parameter vector into matrices used to evaluate likelihood function




mod_params.intercepts=params(1:setup.size_obs);
counter=setup.size_obs+1;
poly_coefficients=params(counter:counter+setup.polynomials*setup.size_obs-1);
mod_params.poly_coefficients=reshape(poly_coefficients,setup.size_obs,setup.polynomials);
counter=counter+setup.polynomials*setup.size_obs;
mod_params.beta_diag_neg=reshape(params(counter:counter+setup.size_obs*(setup.size_obs)-1),setup.size_obs,setup.size_obs);
counter=counter+setup.size_obs*(setup.size_obs);



mod_params.beta_gen_neg=zeros(setup.size_obs,max(setup.num_gaussian),setup.size_obs);
for oo=1:setup.size_obs
for jj=1:max(setup.num_gaussian)
    mod_params.beta_gen_neg(setup.num_gaussian>=jj,jj,oo)=params(counter:(counter+sum(setup.num_gaussian>=jj)-1));
    counter=counter+(sum(setup.num_gaussian>=jj));
    
end
end

mod_params.b_gen_neg=zeros(setup.size_obs,max(setup.num_gaussian),setup.size_obs);
for oo=1:setup.size_obs
for jj=1:max(setup.num_gaussian)
    mod_params.b_gen_neg(setup.num_gaussian>=jj,jj,oo)=params(counter:(counter+sum(setup.num_gaussian>=jj)-1));
    counter=counter+(sum(setup.num_gaussian>=jj));
    
end
end

mod_params.c_gen_neg=ones(setup.size_obs,max(setup.num_gaussian),setup.size_obs); %note that here the default value is one to avoid division by zero
for oo=1:setup.size_obs
for jj=1:max(setup.num_gaussian)
    mod_params.c_gen_neg(setup.num_gaussian>=jj,jj,oo)=params(counter:(counter+sum(setup.num_gaussian>=jj)-1));
    counter=counter+(sum(setup.num_gaussian>=jj));
    
end
end


mod_params.beta_diag_pos=params(counter:counter+setup.size_obs-1);
counter=counter+setup.size_obs;





mod_params.beta_gen_pos=mod_params.beta_gen_neg;
for jj=1:max(setup.num_gaussian)
    mod_params.beta_gen_pos(setup.num_gaussian>=jj,jj,setup.index_unrestricted)=params(counter:(counter+sum(setup.num_gaussian>=jj)-1));
    counter=counter+(sum(setup.num_gaussian>=jj));
    
end


mod_params.b_gen_pos=mod_params.b_gen_neg;
for jj=1:max(setup.num_gaussian)
    mod_params.b_gen_pos(setup.num_gaussian>=jj,jj,setup.index_unrestricted)=params(counter:(counter+sum(setup.num_gaussian>=jj)-1));
    counter=counter+(sum(setup.num_gaussian>=jj));
    
end

mod_params.c_gen_pos=mod_params.c_gen_neg; %note that here the default value is one to avoid division by zero
for jj=1:max(setup.num_gaussian)
    mod_params.c_gen_pos(setup.num_gaussian>=jj,jj,setup.index_unrestricted)=params(counter:(counter+sum(setup.num_gaussian>=jj)-1));
    counter=counter+(sum(setup.num_gaussian>=jj));
    
end







end

