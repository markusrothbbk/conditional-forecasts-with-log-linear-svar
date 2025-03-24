rng(1);
clear;
clc;
close all;

%% Simulate data: gen. obs from unit root var, where levels resemble typical macro
% variables in levels / log levels. Concept
% 1st variable --> to resemble interest rate / unemployment rate 
% 2nd variable --> to resemble log(price index), log(oil price)
% 3rd variable --> to resemble log(real GDP)

Bet = [ [0;0.0025;0.0025] , [ [0.5 0 0]; [0.1 0.5 0.1] ; [0.1 0.1 0.5] ], diag([.1 .1 .1])]; % coefficients of VAR with 2 lags and constant
R =     [ 1  0.5  0.5;   
          0   1   0.5;
          0   0    1  ]; % correlation matrix of disturbances
std_dev = [ 0.001 0.002 0.002 ]; % standard deviations of disturbances
Sigma = diag(std_dev) * (R+triu(R,1)') * diag(std_dev); % variance-covariance matrix
Sigma = (Sigma+Sigma')/2; % ensure symmetry of Sigma 
dy = simvar(200,Bet,Sigma,zeros(2,3)); % simulate obs, thought to be in first differences or log first differences
y0 =  [0.05, log(100), log(1e12)] ; % initial (log) levels
y = cumsum([y0;dy]);   % (log) levels

% estimate parameters of (log) levels var and (log) differences var 
p = 2; % number of lags in estimated VAR
[yy,xx,bb,ss] = fvar(y,p,1);
[dyy,dxx,dbb,dss] = fvar(dy,p,1);

%% Experiment with levels VAR (part of note)

conditions4 = table();
conditions6 = table();
resultsLevels = [];
for vv = [ -.5 -.3 -.1 0 0.1 .3 .5 ]
    conditions4{1,["type","variable","h0","h1","value"]} = [4 3 1 4 log(1+vv)];
    conditions6{1,["type","variable","h0","h1","value"]} = [6 3 1 4 vv];
    [c4_array,~,~,~] = conditionalForecast(bb,chol(ss)',3,p,yy,xx,4,conditions4); 
    [c6_array,~,~,~] = conditionalForecast(bb,chol(ss)',3,p,yy,xx,4,conditions6,showIterationOutput=false); 
    resultsLevels = ...
        [ resultsLevels;
            [   ...
                vv, ...
                mean(exp(c4_array(2:5,3))) / mean(exp(yy(end-3:end,3)))-1, ...
                mean(exp(c6_array(2:5,3))) / mean(exp(yy(end-3:end,3)))-1 ...
            ]
        ];
end
resultsLevelTbl = array2table(resultsLevels,"VariableNames",["Condition","Value without iteration", "Value with iteration"])


%% Experiment with first-differences VAR (not part of note)

conditions5 = table();
conditions7 = table();
resultsDifferences = [];
for vv = [ -.5 -.3 -.1 0 0.1 .3 .5 ]
    conditions5{1,["type","variable","h0","h1","value"]} = [5 3 1 4 log(1+vv)];
    conditions7{1,["type","variable","h0","h1","value"]} = [7 3 1 4 vv];
    [c5_array,~,~,~] = conditionalForecast(dbb,chol(dss)',3,p,dyy,dxx,4,conditions5,showIterationOutput=false); 
    [c7_array,~,~,~] = conditionalForecast(dbb,chol(dss)',3,p,dyy,dxx,4,conditions7,showIterationOutput=false); 
    c5_levels = cumsum( [yy(end,:);c5_array(2:end,1:3)] );
    c7_levels = cumsum( [yy(end,:);c7_array(2:end,1:3)] );
    resultsDifferences = ...
        [ resultsDifferences;
            [   ...
                vv, ...
                mean(exp(c5_levels(2:5,3))) / mean(exp(yy(end-3:end,3)))-1, ...
                mean(exp(c7_levels(2:5,3))) / mean(exp(yy(end-3:end,3)))-1 ...
            ]
        ];
end
resultsDifferencesTbl = array2table(resultsDifferences,"VariableNames",["Condition","Conditional forecast without iteration", "... with iteration"])




