%--------------------------------------------------------------------------
% Bayesian estimation, prediction and impulse response analysis in VAR
% models. Dependent on your choice of forecasting, the VAR model is:
%
% Iterated forecasts:
%     Y(t) = A0 + Y(t-1) x A1 + ... + Y(t-p) x Ap + e(t)
%
% so that in this case there are p lags of Y (from 1 to p).
%
% Direct h-step ahead foreacsts:
%     Y(t+h) = A0 + Y(t) x A1 + ... + Y(t-p+1) x Ap + e(t+h)
%
% so that in this case there are also p lags of Y (from 0 to p-1).
%
% In any of the two cases, the model is written as:
%
%                   Y(t) = X(t) x A + e(t)
%
% where e(t) ~ N(0,SIGMA), and A summarizes all parameters. Note that we
% also use the vector a which is defined as a=vec(A). 
%
%--------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HOW TO USE THIS CODE TO SET PRIOR FOR BM PROJECT:
%1. SET OPTION CONSTANT BELOW (PROBABLY EASIEST TO RUN ON DEMEANED DATA AND
%SET CONSTANT TO 0)
%2. SAVE DATA IN VARIABLE YRAW, WITH THE TIME DIMENSION BEING THE ROWS
%3. RUN CODE (WE ONLY NEED TO RUN TO LINE 140)
%4. THE A MATRIX WE NEED IS THE TRANSPOSE (!!) OF A_OLS (IF YOU ESTIMATE A
%CONSTANT THEN YOU NEED TO TAKE OUT THE FIRST ROW OF A_OLS (WHICH CONTAINS
%THE INTERCEPTS).
%5.THE SIGMA MATRIX WE NEED IS SIGMA_OLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------LOAD DATA-----------------------------------
if setup.polynomials>0
    load ((['data_file_detrended',num2str(setup.size_obs),'.mat']),'data')
else
    load ((['data_file',num2str(setup.size_obs),'.mat']), 'data')
end

Yraw=data';

% In any case, name the data you load 'Yraw', in order to avoid changing the
% rest of the code. Note that 'Yraw' is a matrix with T rows by M columns,
% where T is the number of time series observations (usually months or
% quarters), while M is the number of VAR dependent macro variables.
%----------------------------PRELIMINARIES---------------------------------
% Define specification of the VAR model
constant = 1;        % 1: if you desire intercepts, 0: otherwise 
p = setup.VARsym_order;  % Number of lags on dependent variables
forecasting = 0;     % 1: Compute h-step ahead predictions, 0: no prediction
forecast_method = 0; % 0: Direct forecasts 
                     % 1: Iterated forecasts
h = 0;               % Number of forecast periods

% Set prior for BVAR model:
prior = 1;  % prior = 1 --> Noninformative Prior
            % prior = 2 --> Minnesota Prior
            % prior = 3 --> Natural conjugate Prior
  
%--------------------------DATA HANDLING-----------------------------------
% Get initial dimensions of dependent variable
[Traw M] = size(Yraw);

% The model specification is different when implementing direct forecasts,
% compared to the specification when computing iterated forecasts.
if forecasting==1
    if h<=0
        error('You have set forecasting, but the forecast horizon choice is wrong')
    end    

    % Now create VAR specification according to forecast method
    if forecast_method==0       % Direct forecasts
        Y1 = Yraw(h+1:end,:);
        Y2 = Yraw(2:end-h,:);
        Traw = Traw - h - 1;
    elseif forecast_method==1   % Iterated forecasts
        Y1 = Yraw;
        Y2 = Yraw;
    else
        error('Wrong choice of forecast_method')
    end
else
   Y1 = Yraw;
   Y2 = Yraw;
end
        
% Generate lagged Y matrix. This will be part of the X matrix
Ylag = mlag2(Y2,p); % Y is [T x M]. ylag is [T x (Mp)]

% Now define matrix X which has all the R.H.S. variables (constant, lags of
% the dependent variable and exogenous regressors/dummies)
if constant
    X1 = [ones(Traw-p,1) Ylag(p+1:Traw,:)];
else
    X1 = Ylag(p+1:Traw,:);  %#ok<UNRCH>
end

% Get size of final matrix X
[Traw3 K] = size(X1);

% Create the block diagonal matrix Z
Z1 = kron(eye(M),X1);

% Form Y matrix accordingly
% Delete first "LAGS" rows to match the dimensions of X matrix
Y1 = Y1(p+1:Traw,:); % This is the final Y matrix used for the VAR

% Traw was the dimesnion of the initial data. T is the number of actual 
% time series observations of Y and X
T = Traw - p;

%========= FORECASTING SET-UP:
% Now keep also the last "h" or 1 observations to evaluate (pseudo-)forecasts
if forecasting==1
    if forecast_method==0  % Direct forecasts, we only need to keep the 
        Y = Y1(1:end-1,:);                             % last observation
        X = X1(1:end-1,:);
        Z = kron(eye(M),X);
        T = T - 1;
    else              % Iterated forecasts, we keep the last h observations
        Y = Y1(1:end-h,:);
        X = X1(1:end-h,:);
        Z = kron(eye(M),X);
        T = T - h;
    end
else
    Y = Y1;
    X = X1;
    Z = Z1;
end

 
%--------------------------------PRIORS------------------------------------
% First get Ordinary Least Squares (OLS) estimators
A_OLS = inv(X'*X)*(X'*Y); % This is the matrix of regression coefficients
a_OLS = A_OLS(:);         % This is the vector of coefficients, i.e. it holds
                          % that a_OLS = vec(A_OLS)
SSE = (Y - X*A_OLS)'*(Y - X*A_OLS);
SIGMA_OLS = SSE./(T-K);
errors=(Y - X*A_OLS)';

% ATTENTION: WATCH OUT IF YOU WANT A CONSTANT OR NOT IN THE VAR
if constant==1
    A=A_OLS(2:end,:)';
else
    A=A_OLS';
end

sigma=SIGMA_OLS;

save ((['VAR',num2str(setup.size_obs),'.mat']), 'A', 'sigma','errors')

