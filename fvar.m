% Estimate frequentist VAR
%--------------------------------------------------------------------------
% AUTHORS: Frieder Mokinski 
%--------------------------------------------------------------------------  

function [Y,X,AHat,SigmaHat] = fvar(YRaw,p,constant)
% function [Y,X,AHat,SigmaHat] = fvar(Yraw,p,constant)
%
% Model is:
%
% 	Y(t) = a0 + A1 x Y(t-1) + ... + Ap x Y(t-p) + e(t)
%
%
% where Y(t) is Nx1 and e(t) ~ N(0,SIGMA) 
%
%  A := [a0,A1,...,Ap]' [K x N]         where K = N * p + constant
%
% Input:  Yraw          [Traw x N]      raw endogenous series
%         p             [1 x 1]         lag length
%         constant      [1 x 1]         =0: without
%                                       =1: with constant
%        
% Output: Y             [T x N]         left hand side variables in VAR
%         X             [T x K]         right hand side variables in VAR
%         AHat          [K x N]         draws coefficients
%         SigmaHat      [N x N]         covariance matrix

% Get initial dimensions of dependent variable
[TRaw, N] = size(YRaw);

% Generate lagged Y matrix. This will be part of the X matrix
Ylag = mlag2(YRaw,p); % Y is [Traw x N]. Ylag is [Traw x (N*p)]

% X contains lags of Y and exogenous variables (constant and/or trend)
if constant == 0 % no constant
    X = Ylag(p+1:TRaw,:);
elseif constant == 1 % constant 
    X = [ones(TRaw-p,1) Ylag(p+1:TRaw,:)];
end    

% Form Y matrix accordingly
% Delete first "LAGS" rows to match the dimensions of X matrix
Y = YRaw(p+1:TRaw,:); % This is the final Y matrix used for the VAR

% Traw was the dimension of the initial data. T is the number of actual 
% time series observations of Y and X (we lose the p-lags)
T = TRaw - p;

AHat = ((X'*X)\(X'*Y))';
epsilonHat = Y-X*AHat'; 
SigmaHat = (epsilonHat'*epsilonHat) / (T-N*p-constant);
