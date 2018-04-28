%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code plots the IRFs from a VAR and the fitted GMA impulse responses
% (along with the corresponding Gaussian basis functions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Plot fitted and VAR IRFs
epsilon_vec=zeros(setup.size_obs,setup.lags+1);
[ Sigma, intercept] = unwrap_NL_IRF( setup.initial_parameter,epsilon_vec,setup);
%for now only plot the response to the unrestricted shock

for kk=1:setup.size_obs
    figure;
    for jj=1:setup.size_obs
        subplot(setup.size_obs,1,jj)
        plot(1:size(Sigma,3),squeeze(Sigma(jj,kk,:)),1:size(Sigma,3),squeeze(setup.store_responses(jj,kk,:)),'LineWidth',2)
        legend('fitted GMA response','VAR response')
        grid on
        str = sprintf('Variable %d , shock %d',jj,kk);
        title(str);
        
        
    end
end

% Plot the individual Gaussian basis functions
[ Sigma_ind] = unwrap_NL_IRF_individual_bases( setup.initial_parameter,epsilon_vec,setup);

for kk=1:setup.size_obs
    figure;
    for jj=1:setup.size_obs
        subplot(setup.size_obs,1,jj)
        for bb=1:max(setup.num_gaussian) %note that if one irf has less than the maximum number of gaussians, some zero lines will be plotted
            hold on, plot(1:size(Sigma,3),squeeze(Sigma_ind(jj,bb,kk,:)),'LineWidth',2)
        end
        
        grid on
        str = sprintf('Variable %d , shock %d',jj,kk);
        ylabel(str);
        if jj==1
            title('Individual Gaussian basis functions')
        end
        
    end
end

