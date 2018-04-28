
% load (([fileroot,'/Tools/Data/data_file',num2str(setup.size_obs),'.mat']))

%build X matrix of RHS variables
X=[];
for jj=1:setup.polynomials+1
    % X=[X [0:setup.sample_size-1]'.^(jj-1)];
    X=[X (2*[0:size(data,2)-1]'/(size(data,2)-1)-1).^(jj-1)]; %rescale polynomials to make algorithm more stable
    % see https://www.rssd.esa.int/SP/LISAPATHFINDER/docs/Data_Analysis/DA_Nine/Annex7b.pdf
end

%run regressions to get trend
intercepts=zeros(setup.size_obs,1);
poly_coefficients=zeros(setup.size_obs,setup.polynomials);
for kk=1:setup.size_obs
    params_temp=X\data(kk,:)';
    intercepts(kk)=params_temp(1);    
    poly_coefficients(kk,:)=params_temp(2:end);   
end

trends=[intercepts poly_coefficients]*X';
data_original=data;
data=data-[intercepts poly_coefficients]*X';

save(['data_file_detrended',num2str(setup.size_obs),'.mat'])

figure;
for jj=1:setup.size_obs
    subplot(setup.size_obs,1,jj)
    plot(1:size(data,2),data_original(jj,:)',1:size(data,2),trends(jj,:));
    title('data and trends')
end

