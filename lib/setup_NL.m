%creates cell with objects needed for estimation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%General Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of draws after burn-in 
setup.number_of_draws=50000;%500000
%number of draws for choosing the scaling matrix in case of standard RW proposal
setup.scaling_draws=4000;%10000
%scaling is recomputed after check_scaling draws
setup.check_scaling=100;%200
%display every disp iter draw
setup.disp_iter=100;
% keep every keep_draw th draw
setup.keep_draw=10;
%proposal=1 ->standard RW
%proposal=2 ->adaptive RW
setup.proposal=1;
%log likelihood computation
%likelihood=3 -> user-supplied logL (only option available in this code)
setup.likelihood=3;
setup.skip_opt=0; %skip optimization and go directly to MCMC

% Assign blocks
setup.number_blocks=4;
setup.index_block{1}=[1:12 40:42]'; %psi_0
setup.index_block{2}=[13:21 43:45]'; %a
setup.index_block{3}=[22:30 46:48]'; %b
setup.index_block{4}=[31:39 49:51]'; %c
% setup.number_blocks=1;
% setup.index_block{1}=[1:51]';
%initial scaling for standard RW proposal
setup.initial_scaling=[.01 .01 .1 .1]'; %one for each block

% Basic setup for data input
setup.lags=40; %lag length of IRF
setup.size_obs=3; %number of observables
% setup.freq=2;% frequency of data: 1 for monthly, 2 for quarterly
setup.shocks=0; %0-> initial shocks are zero, 1-> initial shocks from VAR (reduces sample size)
setup.polynomials=0; %degree of polynomial detrending
% setup.symmetry=1000; %making sure non of the options are turned on
setup.VARsym_order=4; %order of symmetric VAR used for starting values

%setting up symmetry/asymmetry restrictions 
setup.index_restricted=[1 2]; %shocks with symmetric IRFs
setup.index_unrestricted=3; %shocks with asymmetric IRFs

%number of gaussians per IR (same order as observables)
setup.num_gaussian=[1 1 1]'; %[nb for IRs to shock1, nb for IRs to shock2, .. ]

% Options for parameter restrictions:
% log - parameter has to be non-negative
% logit_general - upper and lower bounds are set by hand
% logit - parameter has to be between 0 and 1

% impose diagonal coefficients of impact matrix to be positive
% beta_1>0 and beta_2>0 and beta_3>0
setup.index_log=[4 8 12 42];
setup.length_log=length(setup.index_log);

setup.length_logit_general=0;
setup.index_logit_general=[];
setup.logit_general_lb=[]';
setup.logit_general_ub=[]';

setup.length_logit=0;
setup.index_logit=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the data file, add polynomial trend if needed and estimate VAR
load (['data/data_file',num2str(setup.size_obs),'.mat'])
run (['data/detrending']); % include polynomial detrending if needed
run (['data/symVAR_estimation']); %estimate VAR
save ((['data/VAR',num2str(setup.size_obs),'.mat']), 'A', 'sigma','errors') %save dataset
setup.data=['data/data_file',num2str(setup.size_obs),'.mat'];

%setting up initial shocks - either form VAR or set to zero
if setup.shocks==0
    setup.initial_eps=zeros(setup.size_obs,setup.lags); %intial shocks set to 0
    setup.sample_size=length(data);
else
    setup.initial_eps=errors(:,1:setup.lags);
    setup.sample_size=length(data)-setup.VARsym_order-setup.lags;
end

%size of state vector
setup.state_size=1;
%epsilon for adaptive RW
setup.eps_adaptive=.005^2;

%initial parameter value for optimization
load ((['data/VAR',num2str(setup.size_obs),'.mat']))

setup.VARsymA=A; % A matrix for option 3 above
setup.VARsymcov=sigma; %covariance matrix of residuals from estimated VAR (for option 3 above)
setup.VARsymchol=chol(sigma,'lower'); 

% Store VAR IRFs
setup.horizon=setup.lags; %horizon up to which IRFS are matched (for initialization and setting of data-driven priors)
for kk=1:setup.horizon+1
setup.store_responses(:,:,kk)=VAR_MA_implied_higher_order( setup,kk-1); %store IRFs from VAR
end

%should additional matrices be stored
setup.add_matrices=1;
%dimension of those variables
setup.dim_add_matrices=[setup.size_obs setup.sample_size+setup.lags];

%this option is needed because the code would in general allow the
%likelihood function to return estimates of the state if the Kalman filter
%was used
setup.state_size=1;
