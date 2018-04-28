%creates cell with objects needed for estimation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%General Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of draws after burn-in 
setup.number_of_draws=10000;
%number of draws for choosing the scaling matrix in case of standard RW proposal
setup.scaling_draws=2000;
%scaling is recomputed after check_scaling draws
setup.check_scaling=100;
%display every disp iter draw (After scaling for proposal is found)
setup.disp_iter=1000;
% keep every keep_draw th draw
setup.keep_draw=20;

%proposal=1 ->standard RW
%proposal=2 ->adaptive RW
setup.proposal=1;
%log likelihood computation
%likelihood=3 -> user-supplied logL (only option available in this code)
setup.likelihood=3;
setup.skip_opt=0; %skip optimization and go directly to MCMC

%number of gaussians per IRF (same order as observables)
setup.num_gaussian=[1 1]';

% Assign blocks for multiple block version of MH algorithm
setup.number_blocks=5;
setup.index_block{1}=[1:2 5 19]';
setup.index_block{2}=[3 4 6 20]';
setup.index_block{3}=[7:10 21:22]';
setup.index_block{4}=[11:14 23:24]';
setup.index_block{5}=[15:18 25:26]';

%initial scaling for standard RW proposal
setup.initial_scaling=[1 200 2 20 20]'; 

% Basic setup for data input
setup.lags=40; %lag length of IRF
setup.size_obs=2; %number of observables
% setup.freq=1;% frequency of data: 1 for monthly, 2 for quarterly
setup.shocks=0; %0-> initial shocks are zero, 1-> initial shocks from VAR (reduces sample size)
setup.polynomials=0; %degree of polynomial detrending in estimation
setup.VARsym_order=6; %order of symmetric VAR used for starting values

%setting up restrictions (which shock's resposnes are restricted to be symmetric> Here only the responses to shock 2 are asymmetric) 
setup.index_restricted=1;
setup.index_unrestricted=2;

% Options for parameter restrictions:
%log - parameter has to be non-negative
%logit_general - upper and lower bounds are set by hand
%logit - parameter has to be between 0 and 1

% impose diagonal coefficients of impact matrix to be positive, beta_ii>0
setup.length_log=3;
setup.index_log=[3 6 20];

setup.length_logit_general=0;
setup.index_logit_general=[];
 setup.logit_general_lb=[]';
setup.logit_general_ub=[]';

setup.length_logit=0;
setup.index_logit=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data input: load the data file, add polynomial trend if needed and estimate VAR
run (['data/detrending']);
run (['data/symVAR_estimation']);
save ((['data/VAR',num2str(setup.size_obs),'.mat']), 'A', 'sigma','errors')
setup.data=['data/data_file',num2str(setup.size_obs),'.mat'];

%setting up initial shocks - either form VAR or set to zero
if setup.shocks==0
  setup.initial_eps=zeros(setup.size_obs,setup.lags); %intial shocks set to 0 
  setup.sample_size=length(data);
else
    setup.initial_eps=errors(:,1:setup.lags);
    setup.sample_size=length(data)-setup.VARsym_order-setup.lags;
end


%initial parameter value for optimization
load ((['data/VAR',num2str(setup.size_obs),'.mat']))
setup.VARsymA=A; % A matrix for option 3 above
setup.VARsymcov=sigma; %covariance matrix of residuals from estimated VAR (for option 3 above)
setup.VARsymchol=chol(sigma,'lower'); 


%should additional matrices be stored
setup.add_matrices=1;
%dimension of those variables
setup.dim_add_matrices=[setup.size_obs setup.sample_size+setup.lags];


%this option is needed because the code would in general allow the
%likelihood function to return estimates of the state if the Kalman filter
%was used
setup.state_size=1;

setup.horizon=setup.lags; %horizon up to which IRFS are matched (for initialization and setting of data-driven priors)

for kk=1:setup.horizon+1
setup.store_responses(:,:,kk)=VAR_MA_implied_higher_order( setup,kk-1);
end

