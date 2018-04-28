%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main prg to run the BM asymmetric GMA estimation of a bivariate system using simulated data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
clc

%Initialization/Parametrization of routine
setup_NL_MC;

%To find starting values (or, in our case, priors), match IRFs with those of VAR:
VAR_resp_match_NL;

%Set priors
%In practical applications, the resulting parameter values can be used as starting values
%setup.initial_parameter=[zeros(setup.size_obs,1);poly_coefficients(:);beta_diag;beta2(:);b2(:);c2(:);beta_diag;beta2(:);b2(:);c2(:)];
%Here we just use them as a prior
setup.normal_prior_means=[0;0;poly_coefficients(:);beta_diag2;beta2(:);b2(:);c2(:);beta_diag;beta22(:);b22(:);c22(:)];

setup.normal_prior_std=abs([2*setup.normal_prior_means(1:end)']');
% we will fix some parameter values (the intercepts and the element above the main diagonal in the contemporaneous response) to be basically zero, this is the standard recursive identifying restriction
setup.normal_prior_std([1 2 5 19])=1e-10;

setup.index_normal=1:length(setup.normal_prior_means);
setup.index_gamma=[];
close all;

% Alternative starting value for optimizers (which in turn return starting value for MCMC): use symmetric GBs used for the simulation
% setup.initial_parameter=[0 0 INIT_NEG(:)' A_NEG(:)' B_NEG(:)' C_NEG(:)' INIT_POS(:,setup.index_unrestricted)' A_NEG(:,setup.index_unrestricted)' B_NEG(:,setup.index_unrestricted)' C_NEG(:,setup.index_unrestricted)']';

setup.length_param_vector=length(setup.initial_parameter);

%the code uses, among others, a derivative-based optimizer to find the starting value for
%the MCMC algorithm. With a tight prior (the zero restriction used for the recursive ordering), this optimizer will not work, so
%we will restrict the parameter vector for that optimizer. To do so, we use
%the following index  of unrestricted parameters.

setup.index_non_zero=ones(setup.length_param_vector,1);
setup.index_non_zero([1 2 5 19])=0; %index of upper-diagonal coefficients of impact matrix
setup.index_non_zero=logical(setup.index_non_zero); %we turn this into a logical array so that it can be used for indexing

%estimation
[ draws, acc_rate, log_posteriors, statedraws, add_matrices] = sampling_MH_MC( setup ); %add_matrices are the estimated shocks



