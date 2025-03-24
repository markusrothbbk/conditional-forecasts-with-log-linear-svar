rng(1);
clear;
clc;
close all;

%% Simulate data: gen. obs from unit root var, where levels resemble typical macro
% variables in levels / log levels. Concept
% 1st variable --> to resemble interest rate / unemployment rate 
% 2nd variable --> to resemble log(price index), log(oil price)
% 3rd variable --> to resemble log(real GDP)

Bet = [ [0;0.0025;0.0025] , [ [0.5 0 0]; [0.1 0.5 0.1] ; [0.1 0.1 0.5] ], diag([.1 .1 .1])];
R =     [ 1 -0.25  -0.5;   
          0   1     0.25;
          0   0      1  ]; % correlations
std_dev = [ 0.001 0.002 0.002 ]; % standard deviations
Sigma = diag(std_dev) * (R+triu(R,1)') * diag(std_dev); % variance-covariance matrix
Sigma = (Sigma+Sigma')/2; % ensure symmetry of Sigma 
dy = simvar(200,Bet,Sigma,zeros(2,3)); % simulate obs, thought to be in first differences or log first differences
y0 = -sum(dy) + [0.05, log(100), log(1e12)] ; % initial level such that terminal level is second term in sum (to allow sensible value restrictions)
y = cumsum([y0;dy]);   % (log) levels

% estimate parameters of (log) levels var and (log) differences var 
p = 2; % number of lags in estimated VAR
[yy,xx,bb,ss] = fvar(y,p,1);
[dyy,dxx,dbb,dss] = fvar(dy,p,1);

%% (1) Selection of conditions for levels models according to conditionsTable.xlsx.
% - example combines cases 1.a, 1.b, 2, 3.a, 3.b, 4.a, 5.a, 6.a.
% - example also uses conditionalForecast to obtain an unconditional
% forecast

conditions1 = table();
conditions1{1:8,["type","variable","h0","h1","value"]} = [...
    [1 1  4  4  0.05];  % (1.a) variable 1 = 0.05 in period T+4
    [9 1  4  8  0.005]; % (2)   variable 1 to rise by 0.005 between T+4 and T+8 
    [3 1  9 12  0.06];  % (3.a) variable 1 to average 0.06 between T+9 and T+12 
    [4 1 13 16 -0.02];  % (4.a) variable 1 to obtain 0.02 lower average between T+13 and T+16 than
                        %       between T+9 and T+12  
    [5 2 4 4 log(100)]; % (1.b) exp of variable 2 = 100 in period T+4
    [8 2 5 8 110];      % (3.b) exp of variable 2 to average 110 between periods T+5 and T+8
    [6 3 1 4 0.1]       % (5.a) exp of variable 3 to obtain 10% higher average between T+1 and T+4 than
                        %       between T-3 and T
    [9 3 4 8 log(1.08)] % (6.a) exp of variable 3 to grow by 8% between T+4 and T+8
];

% obtain conditional forecast
[c1,c1_draws,c1_table,c1_checks] = conditionalForecast(bb,chol(ss)',3,p,yy,xx,16,conditions1,1000); 
c1_draws_p05 = quantile(c1_draws,0.05,3);
c1_draws_p95 = quantile(c1_draws,0.95,3);

% for exposition, check for c1_array if conditions hold, note that similar
% information is collected automatically in c1_checks

    % (1.a) variable 1 = 0.05 in period T+4
    tmp1a = c1(5,1)  
    
    % (2) variable 1 to rise by 0.005 between T+4 and T+8 
    tmp2 = c1(9,1)-c1(5,1) 
    
    % (3.a) variable 1 to average 0.06 between T+9 and T+12 
    tmp3a = mean(c1(10:13,1)) 
    
    % (4.a) variable 1 to obtain 0.02 lower average between T+13 and T+16 than
    % between T+9 and T+12 
    tmp4a = mean(c1(14:17,1))-mean(c1(10:13,1)) 
    
    % (1.b) exp of variable 2 = 100 in period T+4
    tmp1b = exp(c1(5,2))
    
    % (3.b) exp of variable 2 to average 110 between periods T+5 and T+8
    % [case with condition approximation and approximation error iteration]
    tmp3b = mean(exp(c1(6:9,2)))
    
    % (5.a) exp of variable 3 to obtain 10% higher average between T+1 and T+4
    % than between T-3 and T
    % [case with condition approximation and approximation error iteration]
    tmp5a = mean(exp(c1(2:5,3))) / mean(exp(yy(end-3:end,3))) - 1 
    
    % (6.a) exp of variable 3 to grow by 8% between T+4 and T+8
    tmp6a = exp(c1(9,3)) / exp(c1(5,3)) - 1 

% unconditional forecast
[u1_array,u1_draws,u1_table] = conditionalForecast(bb,chol(ss)',3,p,yy,xx,16,[],1000); 
u1_draws_p05 = quantile(u1_draws,0.05,3);
u1_draws_p95 = quantile(u1_draws,0.95,3);

% illustration
figure('Name','Example (1)'); 
HOR = (0:16)';
for jj=1:3
subplot(1,3,jj);
hold on 
plot(HOR,c1(:,jj),'r-')
plot(HOR,c1_draws_p05(:,jj),'r--')
plot(HOR,c1_draws_p95(:,jj),'r--')
plot(HOR,u1_array(:,jj),'k-')
plot(HOR,u1_draws_p05(:,jj),'k--')
plot(HOR,u1_draws_p95(:,jj),'k--')
title(strcat("variable ",string(jj)))
xlabel("h")
end

%% (2) Selection of conditions for 1st differences models according to conditionsTable.xlsx.
% - example combines cases 4.b, 5.b, 6.b.
% - example also uses conditionalForecast to obtain an unconditional
% forecast

conditions2 = table();
conditions2{1:4,["type","variable","h0","h1","value"]} = [...
    [5 1 13 16 -0.02];    % (4.b) level of variable 1 to obtain 0.02 lower average between T+13 and T+16 than
                          %       between T+9 and T+12  
    [7 3 1 4 0.1]         % (5.b) exp of level of variable 3 to obtain 10% higher average between T+1 and T+4 than
                          %       between T-3 and T
    [3 3 5 8 log(1.08)/4] % (6.b) exp of level of variable 3 to grow by 8% between T+4 and T+8
                          % NOTE: h0=5 NOT h0=4 <=> condition is imposed on 
                          % d.ln(T+5) + d.ln(T+6) + d.ln(T+7) + d.ln(T+8) = ln(T+8) - ln(T+4)
    [2 3 1 4 0]           % (7) 3rd shock = 0 for every period between T+1 and T+4,  
                          % NOTE: 2nd argument of conditionalForecast is
                          % impulse response function @ horizon zero - in
                          % the case below the impulse response function is
                          % identified through a recursive ordering, in
                          % which the 3rd shock only affects the 3rd
                          % variable on impact.
];


% obtain conditional forecast
[c2,c2_draws,c2_table,c2_checks] = conditionalForecast(dbb,chol(dss)',3,p,dyy,dxx,16,conditions2,1000); 
c2_draws_p05 = quantile(c2_draws,0.05,3);
c2_draws_p95 = quantile(c2_draws,0.95,3);

% for exposition, check for c2_array if conditions hold, note that similar
% information is collected automatically in c2_checks
    
    % obtain forecasts for (log) levels of variables
    c2_levels = cumsum([yy(end,:);c2(2:end,1:3)]);
    
    % (4.b) (level of) variable 1 to obtain 0.02 lower average between T+13 and T+16 than
    % between T+9 and T+12 
    tmp4b = mean(c2_levels(14:17,1))-mean(c2_levels(10:13,1)) 
    
    % (5.b) exp of (level of) variable 3 to obtain 10% higher average between T+1 and T+4
    % than between T-3 and T
    % [case with condition approximation and approximation error iteration]
    tmp5b = mean(exp(c2_levels(2:5,3))) / mean(exp(yy(end-3:end,3))) - 1 
    
    % (6.b)  exp of (level of) variable 3 to go down by 10% between T+4 and T+8
    tmp6b = exp(c2_levels(9,3)) / exp(c2_levels(5,3)) - 1 
    
% unconditional forecast
[u2_array,u2_draws,u2_table] = conditionalForecast(dbb,chol(dss)',3,p,dyy,dxx,16,[],1000); 
u2_draws_p05 = quantile(u2_draws,0.05,3);
u2_draws_p95 = quantile(u2_draws,0.95,3);

% illustration
figure('Name','Example (2): forecasts in (log-)differences'); 
HOR = (0:16)';
for jj=1:3
subplot(1,3,jj);
hold on 
plot(HOR,c2(:,jj),'r-')
plot(HOR,c2_draws_p05(:,jj),'r--')
plot(HOR,c2_draws_p95(:,jj),'r--')
plot(HOR,u2_array(:,jj),'k-')
plot(HOR,u2_draws_p05(:,jj),'k--')
plot(HOR,u2_draws_p95(:,jj),'k--')
title(strcat("variable ",string(jj)))
xlabel("h")
end


% transform forecasts from (log-)differences to (log-)levels by summing up
% - the log-levels in the last period before the forecast horizon
% [yy(end,:)], and
% - the cumulated log-differences over the periods of the forecast horizon
% [cumsum(c2_array(:,1:3),1)]
% logic: y_T+h = y_T + (y_T+1-y_T) + ... + (y_T+h-y_T+h-1) 
c2_levels = [ yy(end,:) ; yy(end,:) + cumsum(c2(2:end,1:3),1) ]; 
c2_draws_levels = cat( 1, repmat(yy(end,:),[1 1 1000]), ...
                        yy(end,:) + cumsum(c2_draws(2:end,1:3,:),1)); 
c2_draws_p05_levels = quantile(c2_draws_levels,0.05,3);
c2_draws_p95_levels = quantile(c2_draws_levels,0.95,3);

u2_levels = yy(end,:)+cumsum(u2_array(:,1:3),1) - u2_array(1,1:3);
u2_draws_levels = yy(end,:)+cumsum(u2_draws(:,1:3,:),1) - u2_array(1,1:3);
u2_draws_p05_levels = quantile(u2_draws_levels,0.05,3);
u2_draws_p95_levels = quantile(u2_draws_levels,0.95,3);

figure('Name','Example (2): forecasts in (log-)levels'); 
HOR = (0:16)';
for jj=1:3
subplot(1,3,jj);
hold on 
plot(HOR,c2_levels(:,jj),'r-')
plot(HOR,c2_draws_p05_levels(:,jj),'r--')
plot(HOR,c2_draws_p95_levels(:,jj),'r--')
plot(HOR,u2_levels(:,jj),'k-')
plot(HOR,u2_draws_p05_levels(:,jj),'k--')
plot(HOR,u2_draws_p95_levels(:,jj),'k--')
title(strcat("variable ",string(jj)))
xlabel("h")
end

