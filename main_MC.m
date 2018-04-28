%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Master file for GMA estimation of bivariate system
% data simulated from DGP based on one Gaussian basis function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;
tic

addpath('lib');

%we will use 2 obs
size_obs=2;
%length of simulation
length_simul=200;
%number of lags in moving average
lags=40;

% Generate data from DGP with 1 Gaussian basis function
%storing across simulations
Atrue(:,:)=[2 1;-.5 2];
B(:,:)=[5 10; 7 -10];
C(:,:)=[300/2 200/2; 200/2 900/2];
INIT(:,:)=[1 0;-.1 1];


shocks=zeros(size_obs,length_simul+lags);

shocks(:,lags+1:end)=randn([size_obs length_simul]); %initial shocks set to 0

SHOCKS(:,:)=shocks;

data=zeros(size_obs,length_simul);


for ll=1:lags
Sigma(:,:,ll)=SL_NL_symmetric( Atrue(:,:),B(:,:),C(:,:),ll );
end



A_NEG(:,:)=Atrue(:,:);
B_NEG(:,:)=B(:,:);
C_NEG(:,:)=C(:,:);
INIT_NEG(:,:)=INIT(:,:);

%creating asymmetry in DGP
A_POS(:,:)=Atrue(:,:);
A_POS(:,2)=3*A_POS(:,2);
B_POS(:,:)=B(:,:);
B_POS(:,:)=B_POS(:,:);
C_POS(:,:)=C(:,:);
INIT_POS(:,:)=INIT(:,:);
INIT_POS(:,2)=3*INIT_POS(:,2);


for ll=1:lags
Sigma_pos(:,:,ll)=SL_NL_symmetric( A_POS(:,:),B_POS(:,:),C_POS(:,:),ll);
end
SIGMA_POS(:,:,:)=Sigma_pos; 

for ll=1:lags
Sigma_neg(:,:,ll)=SL_NL_symmetric(A_NEG(:,:),B_NEG(:,:),C_NEG(:,:),ll );
end
SIGMA_NEG(:,:,:)=Sigma_neg;
data=zeros(size_obs,length_simul);



figure;
for pp=1:4
    ind_vec=[1 3 2 4];
subplot(2,2,pp);
[I, J]=ind2sub([2,2],ind_vec(pp));
plot([INIT_POS(ind_vec(pp)) squeeze(SIGMA_POS(I,J,:))'],'LineWidth',2)
grid on
title('true response to positive shocks')
end

figure;
for pp=1:4
subplot(2,2,pp);
[I, J]=ind2sub([2,2],ind_vec(pp));
plot([INIT_NEG(ind_vec(pp)) squeeze(SIGMA_NEG(I,J,:))'],'LineWidth',2)
grid on
title('true response to negative shocks')
end


for kk=1:length_simul
    if shocks(2,kk+lags)>0
data(:,kk)=INIT_POS(:,:)*shocks(:,kk+lags);
    else
     data(:,kk)=INIT_NEG(:,:)*shocks(:,kk+lags);   
    end
for jj=1:lags
    if shocks(2,kk+lags-jj)>0
data(:,kk)=Sigma_pos(:,:,jj)*shocks(:,kk+lags-jj)+data(:,kk);
    else
     data(:,kk)=Sigma_neg(:,:,jj)*shocks(:,kk+lags-jj)+data(:,kk);   
    end
end
end


true_pos=zeros(2,2,lags+1);
true_pos(:,:,2:end)=SIGMA_POS(:,:,:);
true_pos(:,:,1)=INIT_POS(:,:);

true_neg=zeros(2,2,lags+1);
true_neg(:,:,2:end)=SIGMA_NEG(:,:,:);
true_neg(:,:,1)=INIT_NEG(:,:);

save ('data/true_asymmetry.mat', 'true_pos', 'true_neg')

save(['data/data_file2.mat'],'data');

%estimation
Principal_MC;

plots_irfs_MC


toc

save results/MC_results_GB

