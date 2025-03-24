function y = simvar(T,Bet,Sigma,y0)
% SIMVAR - simulate observations from a VAR 
%
% Underlying VAR:
%   yt = nu0 + B1 yt-1 + ... + Bp yt-p + ut,
%       where yt, its lags, and ut are [nvar X 1] vectors, and 
%       V(ut) = Sigma
%
%   INPUT:
%       T   -   Number of observations (rows of y)
%       Bet -   [nvar X 1+nvar*p] matrix of slope coefficients (each rows refers to
%               one of the LHS variables, structure: Bet = { nu0 B1 ... Bp } )  
%       Sigma - [nvar x nvar] VC-matrix, sufficient to provide upper triangular
%               matrix
%     [ y0  -   [p x nvar] matrix of initial observations ]
%
%   OUTPUT:
%     y     - [T x N] matrix generated from VAR(L) model
%
% AUTHOR: Frieder Mokinski

%% Input parsing and basic processing

arguments
    T {mustBeScalarOrEmpty,mustBePositive} = 100
    Bet = [zeros(3,1), .45*eye(3)+0.15*ones(3,3)];
    Sigma = .5*ones(3,3)+0.5*eye(3);
    y0 = [];
end

nvar = size(Sigma,1); % # endogenous variables

p = (size(Bet,2)-1)/nvar; % # lags

if istriu(Sigma), Sigma = Sigma + Sigma' - diag(diag(Sigma)); end

if ~isPositiveInteger(p), error("Bet has inappropriate dimensions: must have nvar rows and 1+nvar*p columns."), end

if ~checkVarianceMatrix(Sigma), error("Sigma is not an admissible variance matrix"), end

if ~isempty(y0), if ~(size(y0,1)==p) || ~(size(y0,2)==nvar)
    error("y0 has inappropriate dimensions: must have ""p"" rows and ""nvar"" columns."),
end, end 
  

%% Simulate observations

% if starting values are not provided, generate from random uniform, then
% simulate a few burn-in draws to discard in final output
if isempty(y0)
    y0 = rand(p,nvar);
    tBurnIn = 100;
else 
    tBurnIn = 0;
end

% init y 
y = y0;

cholSigma = chol(Sigma);
for nn = 1:tBurnIn+T
    u = cholSigma'*randn(nvar,1);
    ylag = y( end:(-1):(end-p+1), : )';
    y(p+nn,:) = ( Bet*[1;ylag(:)] + u )';
end

% remove initial + burn-in observations
y = y(end-T+1:end,:);

