%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main prg to run the FAIR asymmetric estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

tic

addpath(['lib']);
addpath(['results']);

%Initialization/Parametrization of routine
setup_NL;

%To find starting values (or, in our case, priors), match IRFs with those of VAR:
VAR_resp_match_NL;
setup.length_param_vector=length(setup.initial_parameter);

close all
pause(.1)


%% Set parameter restrictions

% Impose lower-bnd / upper-bnd on parameters
% logit - parameter has to be between 0 and 1
% log - parameter has to be non-negative
% logit_general - upper and lower bounds are set by hand

% Impose -5<b<25 and 4<c<5000 (for illustration)
setup.index_logit_general=[22:30 46:48 31:39 49:51];
setup.length_logit_general=length(setup.index_logit_general);
setup.logit_general_lb=[-5*ones(1,12) 4*ones(1,12)]';
setup.logit_general_ub=[25*ones(1,12) 5000*ones(1,12)]';

% % Alternative identification scheme: sign restriction to (set) identify shock
% % beta0_pi<0 and a_pi<0 for shock to third variable (monetary shock)
% % a_r>0 (and beta0_r>0, which is always imposed in setup_NL) 
% % also impose -25<b<25 and 8<c<4000 (for illustration)
% setup.index_logit_general=[11 41 20 21 44 45 22:30 46:48 31:39 49:51];
% setup.length_logit_general=length(setup.index_logit_general);
% setup.logit_general_lb=[-10  -10  -10 -1e-9 -10 -1e-9 -25*ones(1,12) 8*ones(1,12)]';
% setup.logit_general_ub=[1e-9 1e-9 1e-9 10 1e-9 10 25*ones(1,12) 5000*ones(1,12)]';

%Verify that initial guess is inside the bounds:
if max(setup.initial_parameter(setup.index_logit_general)<=setup.logit_general_lb)==1 | max(setup.initial_parameter(setup.index_logit_general)>=setup.logit_general_ub)==1
    'Initial Guess inside the bounds!!'
    stop
end

%% Set prior parameters (and identifying restrictions imposed through priors)
% %scaling for adaptive MCMC (see handbook of MCMC, page 104) ADAPTED
% setup.scaling_adaptive=.03^2/setup.length_param_vector;

setup.index_gamma=[];
setup.index_normal=[];

%In our examples, we only use Normal priors
setup.index_normal=[1:setup.length_param_vector]';
setup.normal_prior_means=[setup.initial_parameter;];
setup.normal_prior_std=2*abs([setup.initial_parameter;]);

%contemporaneous matrix:
setup.normal_prior_std(setup.index_block{1})=10;
% a:
setup.normal_prior_std(setup.index_block{2})=10;
% b:
setup.normal_prior_std(setup.index_block{3})=10;
% c:
setup.normal_prior_std(setup.index_block{4})=20^2;

% Example of short-run pointwise restriction: tight prior at 0 for upper diagonal elements of last column of contemp impact matrix)
setup.normal_prior_std([10 11 40 41])=1e-3;


%scaling for adaptive MCMC (see handbook of MCMC, page 104) ADAPTED
setup.scaling_adaptive=.03^2/setup.length_param_vector;


%% Estimation
[ draws, acc_rate, log_posteriors, statedraws, add_matrices] = sampling_MH( setup );
%add_matrices store the estimated shocks

mkdir results

save results/test_NL

toc

% Get shock series
inds=setup.number_of_draws/setup.keep_draw;
indices_for_draws=unidrnd(inds,1000,1);
shocks=add_matrices(:,:,indices_for_draws);
median_shock=squeeze(prctile(shocks,50,3));
lower_shock=squeeze(prctile(shocks,5,3));
upper_shock=squeeze(prctile(shocks,95,3));

% Plot IRFs
plots_irfs

