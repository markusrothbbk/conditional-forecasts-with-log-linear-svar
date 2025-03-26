function [ cFcastPoint, cFcastDraws, cFcastPointTbl, varargout ] = conditionalForecast( Bet , A0inv , k , p , yy , xx , hFcst , conditions, nDraws, options )
%CONDITIONALFORECAST Obtain conditional forecasts from a Structural VAR,
%   [ cFcastPoint, cFcastDraws, cFcastPointTbl [, value of restricted variables in conditional forecast] ] = conditionalForecast( Bet , A0inv , k , p , yy , xx , hFcst , conditions, nDraws )
%
%   Code based on
%   (*) Camba-Mendez, 2012, Conditional Forecasts on SVAR models using
%   the Kalman filter, Economics Letters, and
%   (*) Mokinski / Roth, 2025, Forecasting with Log-Linear (S)VAR Models:
%   Incorporating Annual Growth Rate Conditions
%
%   Underlying VAR:
%   yt = nu0 + B1 yt-1 + ... + Bp yt-p + A0inv et,
%   where yt, its lags, and et are kx1 vectors
%
%   INPUTS (unless implied by above equation):
%
%   > Bet = { nu0     B1 ... Bp }
%
%   > A0_inv = [S]IRF(h=0)
%       use A0inv = chol(Sigma)' if model is not structurally identified
%       where Sigma = Var[A0inv et]
%
%   > yy = { y1 ... yT-1 yT }', where
%       T is the point in time of the last observation to be used in the forecast
%
%   > xx = { 1    yt-q-1' ... yt-q-p'
%            .     .          .
%            .     .          .
%            .     .          .
%            1    yT-1'  ...  YT-p' };
%
%   > hFcst is the forecast horizon, i.e. the latest forecast refers to YT+hFcst
%
%   > conditions is either
%       (i) a table with columns ["type","variable","h0","h1","value"], or
%       (ii) a matrix with typical row {type variable h0 h1 value},
%
%       where below for the ease of exposition we abbreviate "variable"
%       with "j"
%
% vvvvvvvvvv how to find the correct "type? vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
%
%  >>>>>>>> [simple condition on variables' single-period values]
%           type-1 condition demands that yT+h0(j) = ... = yT+h1(j) = value,
%           where (j) denotes the j-th element; and
%
%  >>>>>>>> [simple condition on shocks' single-period values]
%           type-2 condition demands that eT+h0(j) = ... = eT+h1(j) = value,
%
%  >>>>>>>> [condition on multiple-period (e.g. annual) average (or sum) if variables are in levels]
%  >>>>>>>> [condition on multiple-period growth rate if variables are in log-first-differnces]
%           type-3 condition demands that the average of the observables
%               1/n  * ( yT+h0(j) + ... + yT+h1(j) ) = value
%               (*) where n = h1 - h0 + 1
%               (*) type-3 can also handle conditions zT+h1(j)/zT+h0(j)-1 = g 
%               when zt(j) enters the VAR in first differences of logs,
%               i.e. yt(j) = d.ln(zt(j)) = ln(zt(j)-ln(zt-1(j))):
%               ==> zT+h1(j)/zT+h0(j)-1 = g
%               <=> ln(zT+h1(j))-ln(zT+h0(j)) = ln(1+g)
%               <=> d.ln(zT+h1(j)) + ... + d.ln(zT+h0+1(j)) = ln(1+g)
%               <=> 1/(h1-h0) * [d.ln(zT+h1(j)) + ... + d.ln(zT+h0+1(j))] = 1/(h1-h0+1) * ln(1+g)
%
%  >>>>>>>> [condition on average (e.g. annual) change]
%           type-4 condition demands that
%               1/n * { (yT+h0(j) + ... + yT+h1(j)) - (yT+h0-n(j) + ... + yT+h1-n(j)) } = value
%               (*) where n = h1 - h0 + 1
%               (*) i.e. the condition refers to the change in the average level
%               of a variable in one period (from T+h0 to T+h1) compared to the
%               equally long period directly preceding this period (from
%               T+h0-n to T+h1-n)
%               (*) example: suppose yt(j) is the level of the unemployment
%               rate in period t, then the condition {4 j 1 4 0.05}
%               demands that the average unemployment rate in periods T+1
%               to T+4 shall be by 0.05 higher than the average
%               unemployment rate in the periods T-3 to T
%
%  >>>>>>>> [condition on average (e.g. annual) change if variables are in first-differnces]
%           type = 5 condition demands same as type-4 condition but with
%               the  difference that the VAR is specified in
%               first differences, i.e. condition is
%               1/n * { (yT+h0(j) + ... + yT+h1(j)) - (YT+h0-n(j) + ... + YT+h1-n(j)) } = value
%               and VAR is specified as
%               d.yt = nu0 + B1 d.yt-1 + ... + Bp d.yt-p + A0inv et,
%               where d.yt = yt - yt-1
%               (*) where n = h1 - h0 + 1
%               (*) example: {5 variable 1 4 0.05} would impose
%               the same condition as in the previous example, but in
%               case of a VAR where the unemployment rate enters  the VAR
%               in first differences (instead of levels)
%
%  >>>>>>>> [condition on average growth rate if variables are in log-levels]
%           type-6 condition demands that
%               { (zT+h0(j) + ... + zT+h1(j)) / (zT+h0-n(j) + ... + zT+h1-n(j)) } - 1 = value
%               (*) where n = h1 - h0 + 1, and
%               (*) where zt enters the VAR in log levels, i.e. yt(j) = ln(zt)
%               (*) example: suppose yt(j) is the log level of real GDP in
%               period t, then the condition {6 j 1 4 0.05} demands that
%               the average level of real GDP in periods T+1 to T+4 shall
%               be 5 percent higher than the average level in the periods
%               T-3 to T
%
%  >>>>>>>> [condition on average growth rate if variables are in log-levels]
%           type = 7 conditions demands same as type-6 condition but 
%               (*) where zt enters the VAR in log first differences, i.e.
%               yt(j) = ln(zt) -  ln(zt-1)
%               (*) example: {7 j 1 4 0.05} would impose the same
%               condition as in the previous example, but in case of a VAR
%               where real GDP enters in first differences of logs (instead of
%               log levels)
%
%  >>>>>>>> [condition on multiple-period (e.g. annual) average (or sum) if variables are in log-levels]
%           type-8 condition demands the same as type-3 conditions, i.e.
%               1/n  * ( zT+h0(j) + ... + zT+h1(j) ) = value
%               (*) BUT where zt enters the VAR in log levels, i.e. yt(j) =
%               ln(zt) [instead of levels, yt(j) = zt, for type-3 conditions ]  
%               (*) example: suppose yt(j) = ln(zt) is the log oil price in
%               period t, then the condition {8 j 1 4 50} demands that
%               the average oil price level in periods T+1 to T+4 shall
%               be 50 
%               (*) where n = h1 - h0 + 1     
%
%  >>>>>>>> [condition on multiple-period change if variables are in levels]
%           type-9 condition demands that yT+h1(j) - yT+h0(j) = value,
%               where (j) denotes the j-th element.
%               (*) Example: suppose yt(j) is the unemployment rate
%               measured in %, then the condition {9 j 0 4 1.2} demands
%               that the unemployment rate increases by 1.2 percentage
%               points between T and T+4 
%
% ^^^^^^^^^^ how to find the correct "type? ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%
%   > nDraws = number of draws from conditional predictive density        
%
%   OUTPUT:
%
%   > cFcastPoint = {   yT'               cf_eT+1'
%                     cf_yT+1'            cf_eT+2'
%                         .                   .
%                         .                   .
%                         .                   .
%                    cf_yT+hFcst'   cf_eT+hFcst+1' }
%
%       --> Array of conditional point forecasts, note that the last sample
%           period's observation YT is included in the first row
%
%   > cFcastDraws = {   yT'               cf_eT+1(dd)'
%                     cf_yT+1(dd)'        cf_eT+2(dd)'
%                         .                   .
%                         .                   .
%                         .                   .
%                   cf_yT+hFcst(dd)'     cf_eT+hFcst+1(dd)' }
%
%       --> 3D-Array of draws from the conditional predictive density
%
%   > cFcastPointTbl --> cFcastPoint converted into a table with row and
%       column names
%
%   > varargout --> diagnostic table with imposed conditions and the mean
%       and standard deviations of the variables on which these conditions
%       are impose; Purpose: check if imposed conditions are satisfied
%
% AUTHOR: Frieder Mokinski

%% 1 INPUT PARSING + BASIC PROCESSING OF INPUTS
% ----------------------------------------------
% ----------------------------------------------
arguments
    Bet 
    A0inv
    k
    p 
    yy
    xx
    hFcst
    conditions
    nDraws = 0
    options.showIterationOutput = true
end

% conditions can be provided as table or array;
% > if table: check that the right columns are included
% > if array: convert to table
noConditions = false;
if ~isempty(conditions)
    if istable(conditions)
        tmpVariablesInTable = conditions.Properties.VariableNames;
        if ~all(ismember(["type","variable","h0","h1","value"],tmpVariablesInTable))
            error("If 'conditions' is provided as a table, then it must contain the columns [type,variable,h0,h1,value].")
        end
        conditionsTbl = conditions(:,["type","variable","h0","h1","value"]);
    else
        conditionsTbl = array2table(conditions,"VariableNames",["type","variable","h0","h1","value"]);
    end
else 
    noConditions = true;
end

if not((size(Bet,1)==k) && size(Bet,2)==1+k*p), error("Bet does not have dimensions (k X 1+kp). Either Bet is incorrect or dimensions k or p are incorrect."); end
if not((size(A0inv,1)==k) && (size(A0inv,2)==k)), error("A0inv does not have dimensions (k X k). Either A0inv is incorrect or dimension k is incorrect."); end
if not(size(yy,2)==k), error("yy does not have dimensions (? X k). Either yy is incorrect or dimension k is incorrect."); end
if not(size(xx,2)==1+k*p), error("xx does not have dimensions (? X 1+kp). Either xx is incorrect or dimensions k or p are incorrect."); end
if not(size(yy,1)==size(xx,1)), error("xx and yy must have the same number of rows."); end
if ~noConditions
    if max(conditionsTbl.variable) > k, error("conditionsTbl.variable has elements that exceed k."); end
    if min(conditionsTbl.type)<1 || max(conditionsTbl.type)>9, error("conditionsTbl.type has elements outside the valid range from 1 to 9."); end
end
warning('off','MATLAB:table:RowsAddedExistingVars')
alpha = 0.5; 

% case where draws are not allocated to an output--> no drawing from
% conditional predictive density
if nargout == 1
    nDraws = 0;
end

% init pStar for case when no conditions are specified (pStar refers to the
% number of lags of the endogenous variables of the VAR to be included in
% the state equation
pStar = p;

% convert information about conditions into "cases" that are processed
% one-after-another in the code below
% > case variable: type 1 / type 2 / type 3 with summation over horizon
% range 1 / type 3 with summation over horizon range 2 / ... / type 4 with
% horizon range 1 / ...
% > why case variable? each case needs own treatment in terms of porting of
% conditions and - for type-3... conditions - in terms of the augmentation of
% matrices C and Y
% > pStar variable: in case of type-3... conditions, the state
% equation may need additional lags of yt (beyond the lag order p / reason
% in case of type-3: summation extends over more periods than p)
if ~isempty(conditions)
    conditionsTbl.tmpHRange = (conditionsTbl.type >= 3) .* (conditionsTbl.h1 - conditionsTbl.h0 + 1);
    conditionsTbl = sortrows(conditionsTbl,["type","tmpHRange"]);
    conditionsTbl.tmpLagsRequired = ...
        ismember(conditionsTbl.type,[1,2]) * 1 + ...
        ismember(conditionsTbl.type,[3,8,9]) .* conditionsTbl.tmpHRange + ...
        ismember(conditionsTbl.type,[4,6]) .* conditionsTbl.tmpHRange * 2 + ...
        ismember(conditionsTbl.type,[5,7]) .* (conditionsTbl.tmpHRange * 2 - 1); % lags required in state equation for conditions
    pStar = max([conditionsTbl.tmpLagsRequired;p]);
    [~,~,conditionsTbl.case] = unique(conditionsTbl{:,["type","tmpHRange"]},"rows","stable");
%     conditionsTbl = conditionsTbl(:,["type","variable","h0","h1","value","case"]);
end

if exist("conditionsTbl","var")

    % this is a copy of the original conditionsTbl
    if ~exist("conditionsTblOriginal","var")
        conditionsTblOriginal = conditionsTbl;
    end

    % the type-6 or type-7 conditions are assumed to hold for the level of
    % variable z
    % { (zT+h0 + ... + zT+h1) / (zT+h0-n + ... + zT+h1-n) - 1 = value
    % BUT the variable  enters the model either in log-levels, yt(j) = ln(zt),
    % or in log differences, yt(j) = ln(zt)-ln(zt-1)
    % --> because the condition is applied to a log approximation of the
    % condition, i.e.
    % [ (ln(zT+h0) + ... + ln(zT+h1)) - (ln(zT+h0-n) + ... + ln(zT+h1-n)) ]
    % we also replace the value "value" of the condition with log(1+"value")
    % REMARK: this replacement of condition values is the reason why we
    % keep a copy of the original conditionsTbl (called conditionsTblOriginal)
    % - we will use this copy of the original value of the exact
    % condition, when the code will check how large the approximation
    % error to the exact condition is, and if this error is to big, it
    % will adjust the value of the approximate condition from log(1+g) to
    % log(1+g)+d until the approximation error is sufficiently small.
    for ii=1:size(conditionsTbl,1)
        if ismember(conditionsTbl.type(ii),[6 7])
            conditionsTbl.value(ii) = log(1+conditionsTbl.value(ii));
        end
    end

    % similarly, type-8 conditions are assumed to hold for the level of
    % variable z
    %  1/n  * ( zT+h0(j) + ... + zT+h1(j) ) = value
    %  BUT the variable  enters the model either in log-levels, yt(j) = ln(zt),
    % --> because the condition is applied to a log approximation of the
    % condition, i.e.
    % 1/n  * ( ln(zT+h0(j)) + ... + ln(zT+h1(j)) ) = ln(value)
    % and we also replace the value "value" of the condition with ln("value")
    % REMARK: this replacement of condition values is the reason why we
    % keep a copy of the original conditionsTbl (called conditionsTblOriginal)
    % - we will use this copy of the original value of the exact
    % condition, when the code will check how large the approximation
    % error to the exact condition is
    for ii=1:size(conditionsTbl,1)
        if ismember(conditionsTbl.type(ii),8)
            conditionsTbl.value(ii) = log(conditionsTbl.value(ii));
        end
    end

end

% tmpGoOn is a dummy indicator that will be set to 'false' as soon as the
% approximation error to the exact condition(s) [type-6/type-7] is
% sufficiently small
tmpGoOn=true;
tmpRun=1;
optimizerOutput=table();

while tmpGoOn

%% 2. SET UP STATE-SPACE SYSTEM
% ------------------------------
% ------------------------------

    % A from state equation st+1 = A st + ut, where st' = [1 et+1' yt' ...
    % yt-pStar+1']
    A =[
        [1 zeros(1,k*(pStar+1))];
        zeros(k,1+k*(pStar+1));
        [Bet(:,1) A0inv Bet(:,2:end) zeros(k,k*(pStar-p))];
        [zeros(k*(pStar-1),1) zeros(k*(pStar-1),k) eye(k*(pStar-1)) zeros(k*(pStar-1),k)]
        ];

    % V(u)
    Sigma_u = zeros(1+k+k*pStar,1+k+k*pStar);
    Sigma_u(2:k+1,2:k+1) = eye(k);

    % C from measurement equation [yt' et+1']' = C * st (must be extended in
    % case of type-3 ... conditions)
    C =[                                            % C =[
        [zeros(k,1+k) eye(k) zeros(k,k*(pStar-1))]; %     [zeros(k,1)  zeros(k,k)  eye(k)      zeros(k,k*(pStar-1)) ] 
        [zeros(k,1)   eye(k) zeros(k,k*pStar)]      %     [zeros(k,1)  eye(k)      zeros(k,k)  zeros(k,k*(pStar-1)) ]
        ];                                          %     ];

    % extract structural shock until T
    res_draw = yy-xx*Bet';
    shocks_draw = res_draw*inv(A0inv)';

    % Y is like Y_smoo without the forecasted elements (see OUTPUT above), i.e.
    %   Y     = [   [YT-1'          eT']
    %               [YT'            NaN]  ]
    Y = NaN(hFcst+2,2*k);
    Y(1:2,1:k) = yy(end-1:end,:);
    Y(1,k+1:2*k) = shocks_draw(end,:);

    % case-by-case: extract respective conditions, augment matrices C and Y if
    % needed, and port all conditions into Y matrix
    if ~isempty(conditions)

        for iCase = unique(conditionsTbl.case,"rows","sorted")'

            tmpConditionsCase = conditionsTbl(conditionsTbl.case == iCase,:);
            tmpTypeCase = tmpConditionsCase.type(1);

            % type-3 to type-9 conditions: generate block that will augment C
            % (measurment equation) / generate empty block that will augment Y
            % (measurement equation) / port conditions into block
            if ismember(tmpTypeCase ,[3,8])
                tmpHRangei = tmpConditionsCase.h1(1)-tmpConditionsCase.h0(1)+1;
                tmpCAug = [ zeros(k,k+1), repmat(eye(k),1,tmpHRangei)/tmpHRangei , zeros(k,k*(pStar-tmpHRangei)) ];
                tmpYAug = NaN(hFcst+2,k);
            elseif ismember(tmpTypeCase ,[4,6])
                tmpHRangei = tmpConditionsCase.h1(1)-tmpConditionsCase.h0(1)+1;
                tmpCAug = [ zeros(k,k+1), repmat(eye(k),1,tmpHRangei)/tmpHRangei, -repmat(eye(k),1,tmpHRangei)/tmpHRangei , zeros(k,k*(pStar-2*tmpHRangei)) ];
                tmpYAug = NaN(hFcst+2,k);
            elseif ismember(tmpTypeCase ,[5,7])
                tmpHRangei = tmpConditionsCase.h1(1)-tmpConditionsCase.h0(1)+1;
                tmpCAugSubBlock = [];
                for ii=1:2*tmpHRangei-1
                    if ii<=tmpHRangei
                        tmpCAugSubBlock = [tmpCAugSubBlock,ii*eye(k)/tmpHRangei];
                    else
                        tmpCAugSubBlock = [tmpCAugSubBlock,(2*tmpHRangei-ii)*eye(k)/tmpHRangei];
                    end
                end
                tmpCAug = [ zeros(k,k+1), tmpCAugSubBlock , zeros(k,k*(pStar-2*tmpHRangei+1)) ];
                tmpYAug = NaN(hFcst+2,k);
            elseif ismember(tmpTypeCase ,9)
                tmpHRangei = tmpConditionsCase.h1(1)-tmpConditionsCase.h0(1)+1;
                tmpCAug = [ zeros(k,k+1), eye(k), zeros(k,k*(tmpHRangei-2)), -eye(k), zeros(k,k*(pStar-tmpHRangei)) ];
                tmpYAug = NaN(hFcst+2,k);                
            end

            for iCondition=1:size(tmpConditionsCase,1)

                tmpCndi = tmpConditionsCase(iCondition,:);

                % type-1 or type-2 conditions: port conditions into
                % Y:
                if tmpTypeCase == 1
                    Y(2+tmpCndi.h0:2+tmpCndi.h1,tmpCndi.variable) = tmpCndi.value;
                elseif tmpTypeCase == 2
                    Y(1+tmpCndi.h0:1+tmpCndi.h1,tmpCndi.variable+k) = tmpCndi.value;
                % type=3...9 conditions: port conditions into
                % augmentation block:
                elseif ismember(tmpTypeCase,3:9)
                    tmpYAug(2+tmpCndi.h1,tmpCndi.variable) = tmpCndi.value;
                end

            end

            % type=3...9 conditions: after processing all conditions of
            % the case append the augmentation blocks to C and Y
            if ismember(tmpTypeCase,3:9)
                C = [C;tmpCAug];
                Y = [Y,tmpYAug];
            end

        end
    end

%% 3. OBTAIN CONDITIONAL POINT FORECAST THROUGH FIXED INTERVAL SMOOTHER
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

    % initialize
    s_tmin1_filt = [1; shocks_draw(end,:)'; vec(yy(end-1:(-1):end-pStar,:)')]; % i.e. [ 1 eT' yT-1'  ...  YT-p' ]
    P_tmin1_filt = zeros(size(A));
    s_tmin1_plus = s_tmin1_filt; % stored to be used in simulation smoother

    % to hold the estimates
    Y_filt = Y;
    Y_smoo = Y;
    s_filt = NaN(size(Y_filt,1),size(s_tmin1_filt,1));
    s_filt(1,:) = s_tmin1_filt';
    s_smoo = NaN(size(Y_filt,1),size(s_tmin1_filt,1));
    P_filt = NaN(size(P_tmin1_filt,1),size(P_tmin1_filt,1),hFcst+2);
    P_filt(:,:,1) = P_tmin1_filt;
    P_fcast = NaN(size(P_tmin1_filt,1),size(P_tmin1_filt,1),hFcst+2);
    P_smoo = NaN(size(P_tmin1_filt,1),size(P_tmin1_filt,1),hFcst+2);

    % KALMAN FILTER
    for t=2:size(Y,1)
        J_t = diag(~isnan(Y(t,:)'));
        P_t_fcast = A*P_tmin1_filt*A' + Sigma_u;
        K_t = P_t_fcast*C'*pinv(J_t*C*P_t_fcast*C'*J_t);
        Y_t = Y(t,:)'; Y_t(isnan(Y_t)) = 0;
        s_t_filt = A*s_tmin1_filt + K_t*(Y_t-C*A*s_tmin1_filt);
        P_t_filt = P_t_fcast - K_t*C*P_t_fcast;

        % save
        Y_filt(t,:) = (C*s_t_filt)';
        s_filt(t,:) = s_t_filt';
        P_filt(:,:,t) = P_t_filt;
        P_fcast(:,:,t) = P_t_fcast;

        % update
        s_tmin1_filt = s_t_filt;
        P_tmin1_filt = P_t_filt;
    end

    Y_smoo(end,:) = Y_filt(end,:);
    s_smoo(end,:) = s_filt(end,:);
    P_smoo(:,:,end) = P_filt(:,:,end);

    s_tplus1_smoo = s_t_filt;
    P_tplus1_smoo = P_t_filt;

    % FIXED INTERVAL SMOOTHER
    for t=size(Y,1)-1:-1:1
        s_t_filt = s_filt(t,:)';
        P_t_filt = P_filt(:,:,t);
        P_tplus1_fcast = P_fcast(:,:,t+1);
        H_t = P_filt(:,:,t)*A'*pinv(P_fcast(:,:,t+1)); %!!!: pinv
        s_t_smoo = s_t_filt + H_t*(s_tplus1_smoo-A*s_t_filt);
        P_t_smoo = P_t_filt - H_t*(P_tplus1_fcast-P_tplus1_smoo)*H_t';

        % save
        Y_smoo(t,:) = (C*s_t_smoo)';
        s_smoo(t,:) = s_t_smoo';
        P_smoo(:,:,t) = P_t_smoo;

        s_tplus1_smoo = s_t_smoo;
        P_tplus1_smoo = P_t_smoo;
    end

    cFcastPoint = Y_smoo(2:end,:);

    %% 3.1 COMPUTE APPROXIMATION ERROR TO THE EXACT RESTRICTION(S) [TYPE-6...8]
    % ADJUST CONSTRAINT ON APPROXIMATION UNTIL APPROXIMATION ERROR IS
    % SUFFICIENTLY SMALL
    % ----------------------------------------------------------------------
    % ----------------------------------------------------------------------
    tmpGoOn = false;
    tmpExact=[];
    tmpExactApproximated=[];

    if ~isempty(conditions)

        if any(ismember(conditionsTbl.type,6:8))

            tmpFcstsWithHistory = [ yy(1:end,:); cFcastPoint(2:end,1:k) ];
            tmpH =  (-size(yy,1)+1:hFcst)'; % forecast horizons corresponding to tmpFcstsWithHistory

            for iCase = unique(conditionsTbl.case,"rows","sorted")'

                tmpConditionsCase = conditionsTbl(conditionsTbl.case == iCase,:);
                tmpConditionsCaseOriginal = conditionsTblOriginal(conditionsTblOriginal.case == iCase,:);
                tmpTypeCase = tmpConditionsCase.type(1);

                if ismember(tmpTypeCase,6:7)

                    for iCondition=1:size(tmpConditionsCase,1)

                        tmpCndi = tmpConditionsCase(iCondition,:);
                        tmpCndOriginali = tmpConditionsCaseOriginal(iCondition,:);
                        tmpHRangei = tmpConditionsCase.h1(1)-tmpConditionsCase.h0(1)+1;

                        if tmpTypeCase==6
                            tmpLevels = exp(tmpFcstsWithHistory);
                        elseif tmpTypeCase==7
                            tmpLevels = exp(cumsum(tmpFcstsWithHistory));
                        end

                        tmpForNom = tmpLevels( tmpH >= tmpCndi.h0 & tmpH <= tmpCndi.h1, tmpCndi.variable );
                        tmpForDenom = tmpLevels( tmpH >= tmpCndi.h0-tmpHRangei & tmpH <= tmpCndi.h1-tmpHRangei, tmpCndi.variable );
                        tmpValueExactrRstrctn = mean(tmpForNom)/mean(tmpForDenom)-1;
                        tmpValueApprxRstrctnNew = tmpCndi.value - alpha*(tmpValueExactrRstrctn-tmpCndOriginali.value); % new value = old value - alpha * ( current value exact condition - target value exact restrction)
                        if tmpGoOn == false
                            tmpValueExactrRstrctnAbsDeviation = abs(tmpValueExactrRstrctn-tmpCndOriginali.value); 
                            tmpGoOn = abs(tmpValueExactrRstrctnAbsDeviation) > 1e-6; % stop condition expessed as absolute deviation from target
                        end
                        tmpExact=[tmpExact,tmpCndOriginali.value];
                        tmpExactApproximated=[tmpExactApproximated,tmpValueExactrRstrctn];
                        conditionsTbl{prod(conditionsTbl{:,:} == tmpCndi{:,:},2)==1,"value"} = tmpValueApprxRstrctnNew;

                    end

                elseif tmpTypeCase==8

                    for iCondition=1:size(tmpConditionsCase,1)

                        tmpCndi = tmpConditionsCase(iCondition,:);
                        tmpCndOriginali = tmpConditionsCaseOriginal(iCondition,:);
                        
                        tmpLevels = exp(tmpFcstsWithHistory);
                        tmpLevelsh0h1 = tmpLevels( tmpH >= tmpCndi.h0 & tmpH <= tmpCndi.h1, tmpCndi.variable );   
                        tmpValueExactrRstrctn = mean(tmpLevelsh0h1);
                        tmpValueApprxRstrctnNew = log( exp(tmpCndi.value) - alpha*(tmpValueExactrRstrctn-tmpCndOriginali.value) ); % new value = old value - alpha * ( current value exact condition - target value exact restrction) )
                                                                                                                                   % note: inside the log() function, we compute the new candidate value in levels 
                        if tmpGoOn == false
                            tmpValueExactrRstrctnPctDeviation = 100*(tmpValueExactrRstrctn-tmpCndOriginali.value)/tmpCndOriginali.value; % percent deviation of value of exact condition from target
                            tmpGoOn = abs(tmpValueExactrRstrctnPctDeviation) > 1e-3; % stop condition expessed as percent (!) deviation from target
                        end
                        tmpExact=[tmpExact,tmpCndOriginali.value];
                        tmpExactApproximated=[tmpExactApproximated,tmpValueExactrRstrctn];
                        conditionsTbl{prod(conditionsTbl{:,:} == tmpCndi{:,:},2)==1,"value"} = tmpValueApprxRstrctnNew;

                    end

                end

            end

            optimizerOutput{tmpRun,"Iteration"} = tmpRun;
            optimizerOutput{tmpRun,"Exact conditions"} = tmpExact;
            optimizerOutput{tmpRun,"Values obtained through approximation"} = tmpExactApproximated;
            tmpRun=tmpRun+1;

            if tmpRun == 1000
                disp(optimizerOutput)
                error("Iteration did not converge after 1000 iteration steps: operation aborted.")
            end

        end

    end

end


%% 4: DRAW FROM CONDITIONAL PREDICTVE DENSITY 
% --------------------------------------------
% --------------------------------------------
% uses the simulation smoother with mean corrections as proposed by DK ed.
% 2, Sec. 4.9.2

cFcastDraws = double.empty;

if nDraws>0
    cFcastDraws_ = NaN(size(Y_smoo,1),size(Y_smoo,2),nDraws);
    drawsS = NaN(size(s_smoo,1),size(s_smoo,2),nDraws);
    for ii=1:nDraws
        % simulated alternative path of Y and S from unconditional distribution
        R = Sigma_u(:,2:k+1); % selection matrix for state disturbances such that s_t=A*s_t-1+R*e_t with u_t=R*e_t
        Y_plus = Y; % init
        s_plus = NaN(size(Y,1),size(s_tmin1_plus,1));
        s_plus(1,:) = s_tmin1_plus';
        s_t_plus = s_tmin1_plus; % init
        for t=2:size(Y,1)
            s_t_plus = A*s_t_plus  + R*randn(k,1);
            s_plus(t,:) = s_t_plus';
            Y_plus(t,:) = (C*s_t_plus)';
        end

        % replace the cells in Y_plus with missings, on which there are no conditions
        Y_plus_o = Y_plus;
        Y_plus_o(isnan(Y)) = NaN(1,1);

        % compute conditional forecast based on the simulated conditions
        % KALMAN FILTER
        P_tmin1_filt = zeros(size(A));
        s_tmin1_filt = s_tmin1_plus;

        Y_plus_filt = Y_plus_o;
        Y_plus_smoo = Y_plus_o;
        s_plus_filt = NaN(size(Y_plus_filt,1),size(s_tmin1_filt,1));
        s_plus_filt(1,:) = s_tmin1_filt';
        s_plus_smoo = NaN(size(Y_plus_filt,1),size(s_tmin1_filt,1));
        P_plus_filt = NaN(size(P_tmin1_filt,1),size(P_tmin1_filt,1),hFcst+2);
        P_plus_filt(:,:,1) = P_tmin1_filt;
        P_plus_fcast = NaN(size(P_tmin1_filt,1),size(P_tmin1_filt,1),hFcst+2);
        for t=2:size(Y,1)
            J_t = diag(~isnan(Y_plus_o(t,:)'));
            P_t_fcast = A*P_tmin1_filt*A' + Sigma_u;
            K_t = P_t_fcast*C'*pinv(J_t*C*P_t_fcast*C'*J_t);
            Y_t = Y_plus_o(t,:)'; Y_t(isnan(Y_t)) = 0;
            s_t_filt = A*s_tmin1_filt + K_t*(Y_t-C*A*s_tmin1_filt);
            P_t_filt = P_t_fcast - K_t*C*P_t_fcast;
            % update
            s_tmin1_filt = s_t_filt;
            P_tmin1_filt = P_t_filt;
            % save
            Y_plus_filt(t,:) = (C*s_t_filt)';
            s_plus_filt(t,:) = s_t_filt';
            P_plus_filt(:,:,t) = P_t_filt;
            P_plus_fcast(:,:,t) = P_t_fcast;
        end

        Y_plus_smoo(end,:) = Y_plus_filt(end,:);
        s_plus_smoo(end,:) = s_plus_filt(end,:);
        s_tplus1_smoo = s_t_filt;
        P_tplus1_smoo = P_t_filt;

        % FIXED INTERVAL SMOOTHER
        for t=size(Y,1)-1:-1:1
            s_t_filt = s_plus_filt(t,:)';
            P_t_filt = P_plus_filt(:,:,t);
            P_tplus1_fcast = P_plus_fcast(:,:,t+1);
            H_t = P_plus_filt(:,:,t)*A'*pinv(P_plus_fcast(:,:,t+1)); %!!!: pinv
            s_t_smoo = s_t_filt + H_t*(s_tplus1_smoo-A*s_t_filt);
            P_t_smoo = P_t_filt - H_t*(P_tplus1_fcast-P_tplus1_smoo)*H_t';
            % update
            s_tplus1_smoo = s_t_smoo;
            P_tplus1_smoo = P_t_smoo;
            % save
            Y_plus_smoo(t,:) = (C*s_t_smoo)';
            s_plus_smoo(t,:) = s_t_smoo';
        end

        % DK show that if (s_T+1',...,s_T+h')'|Y_n ~  N(s_smoo,P_smoo)
        % ==> s_plus-s_plus_smoo ~ N(0,P_smoo) such that:
        %   (s_smoo+s_plus-s_plus_smoo) ~ N(s_smoo,P_smoo)
        cFcastDraws_(:,:,ii) = (C*(s_smoo+s_plus-s_plus_smoo)')';
        drawsS(:,:,ii) = ((s_smoo+s_plus-s_plus_smoo)')';
    end

    cFcastDraws = cFcastDraws_(2:end,:,:);
  
end

%% 5 COLLECT FUNCTION OUTPUTS
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

if exist("optimizerOutput","var") && options.showIterationOutput
    disp(optimizerOutput)
end

% check if conditions are met with precision of at least 1e-5
% note: this checks if the conditions can be set for the forecast and does
% not evaluate the error of the approximation procedure
check = cFcastPoint - Y(2:end,:);
if max(abs(round(check(:),5))) > 0  
    error("Conditions are not met with sufficient precision. A possible explanation is that conditions are not well defined, see Camba-Mendez, 2012, Section 2.2.")
end

cFcastPoint = cFcastPoint(:,1:2*k);

tmpRowNames = strcat( repmat("T",hFcst+1,1) , ["";repmat("+",hFcst,1)], ["";(1:hFcst)']);
tmpColNames = strcat( [repmat("y_t(",k,1);repmat("eps_{t+1}(",k,1)] , string([(1:k)';(1:k)']), repmat(")",2*k,1) );
cFcastPointTbl = array2table(cFcastPoint(:,1:2*k),'RowNames',tmpRowNames,'VariableNames',tmpColNames); % 20240131 MR: get an error with M;atlab 2019b version issue? 
% cFcastPointTbl = array2table(cFcastPoint(:,1:2*k),RowNames=tmpRowNames,VariableNames=tmpColNames); 

% value of restricted variables in conditional forecast --> put into extra colums of conditionsTblOriginal
% purpose: check if / how well conditions are satisfies 
% ==> optional function output
if ~isempty(conditions)

    tmpFcstsWithHistory = [ yy(1:end,:); cFcastPoint(2:end,1:k) ];
    tmpFcstsWithHistoryCumulated = cumsum([ yy(1:end,:); cFcastPoint(2:end,1:k) ]);
    tmpH =  (-size(yy,1)+1:hFcst)'; % forecast horizons corresponding to tmpFcstsWithHistory

    conditionsTblOriginal = addvars(conditionsTblOriginal,strings(size(conditionsTblOriginal,1),1),strings(size(conditionsTblOriginal,1),1),'NewVariableNames',["cndFcast:avgValue","cndFcast:stdValue"]); % <=== 20240202 MR avoids warning messages in Matlab 2019b
    for ii=1:size(conditions,1)

        tmpCndi = conditionsTblOriginal(ii,:);

        if conditionsTblOriginal.type(ii) == 1

            tmpValues = tmpFcstsWithHistory( tmpH >= tmpCndi.h0 & tmpH <= tmpCndi.h1, tmpCndi.variable );
            conditionsTblOriginal{ii,"cndFcast:avgValue"} = mean(tmpValues,1);
            conditionsTblOriginal{ii,"cndFcast:stdValue"} = std(tmpValues,1);

        elseif conditionsTblOriginal.type(ii) == 2

            tmpValues = cFcastPoint(tmpCndi.h0:tmpCndi.h1,k+tmpCndi.variable);
            conditionsTblOriginal{ii,"cndFcast:avgValue"} = mean(tmpValues,1);
            conditionsTblOriginal{ii,"cndFcast:stdValue"} = std(tmpValues,1);            

        elseif conditionsTblOriginal.type(ii) == 3

            tmpValues = tmpFcstsWithHistory( tmpH >= tmpCndi.h0 & tmpH <= tmpCndi.h1, tmpCndi.variable );
            conditionsTblOriginal{ii,"cndFcast:avgValue"} = mean(tmpValues,1);
            conditionsTblOriginal{ii,"cndFcast:stdValue"} = std(tmpValues,1);                        
           
        elseif conditionsTblOriginal.type(ii) == 4

            tmpValues1 = tmpFcstsWithHistory( tmpH >= tmpCndi.h0 & tmpH <= tmpCndi.h1, tmpCndi.variable );
            tmpValues2 = tmpFcstsWithHistory( tmpH >= tmpCndi.h0-tmpCndi.tmpHRange & tmpH <= tmpCndi.h1-tmpCndi.tmpHRange, tmpCndi.variable );
            conditionsTblOriginal{ii,"cndFcast:avgValue"} = mean(tmpValues1,1)-mean(tmpValues2,1);
            conditionsTblOriginal{ii,"cndFcast:stdValue"} = nan;                        

        elseif conditionsTblOriginal.type(ii) == 5

            tmpValues1 = tmpFcstsWithHistoryCumulated( tmpH >= tmpCndi.h0 & tmpH <= tmpCndi.h1, tmpCndi.variable );
            tmpValues2 = tmpFcstsWithHistoryCumulated( tmpH >= tmpCndi.h0-tmpCndi.tmpHRange & tmpH <= tmpCndi.h1-tmpCndi.tmpHRange, tmpCndi.variable );
            conditionsTblOriginal{ii,"cndFcast:avgValue"} = mean(tmpValues1,1)-mean(tmpValues2,1);
            conditionsTblOriginal{ii,"cndFcast:stdValue"} = nan;            
        
        elseif conditionsTblOriginal.type(ii) == 6

            tmpValues1 = exp(tmpFcstsWithHistory( tmpH >= tmpCndi.h0 & tmpH <= tmpCndi.h1, tmpCndi.variable ));
            tmpValues2 = exp(tmpFcstsWithHistory( tmpH >= tmpCndi.h0-tmpCndi.tmpHRange & tmpH <= tmpCndi.h1-tmpCndi.tmpHRange, tmpCndi.variable ));
            conditionsTblOriginal{ii,"cndFcast:avgValue"} = mean(tmpValues1,1)/mean(tmpValues2,1)-1;
            conditionsTblOriginal{ii,"cndFcast:stdValue"} = nan; 

        elseif conditionsTblOriginal.type(ii) == 7

            tmpValues1 = exp(tmpFcstsWithHistoryCumulated( tmpH >= tmpCndi.h0 & tmpH <= tmpCndi.h1, tmpCndi.variable ));
            tmpValues2 = exp(tmpFcstsWithHistoryCumulated( tmpH >= tmpCndi.h0-tmpCndi.tmpHRange & tmpH <= tmpCndi.h1-tmpCndi.tmpHRange, tmpCndi.variable ));
            conditionsTblOriginal{ii,"cndFcast:avgValue"} = mean(tmpValues1,1)/mean(tmpValues2,1)-1;
            conditionsTblOriginal{ii,"cndFcast:stdValue"} = nan; 
                      
        elseif conditionsTblOriginal.type(ii) == 8

            tmpValues = exp(tmpFcstsWithHistory( tmpH >= tmpCndi.h0 & tmpH <= tmpCndi.h1, tmpCndi.variable ));
            conditionsTblOriginal{ii,"cndFcast:avgValue"} = mean(tmpValues,1);
            conditionsTblOriginal{ii,"cndFcast:stdValue"} = std(tmpValues,1);     

        elseif conditionsTblOriginal.type(ii) == 9

            tmpValues = tmpFcstsWithHistory( tmpH == tmpCndi.h1, tmpCndi.variable ) - ...
                     tmpFcstsWithHistory( tmpH == tmpCndi.h0 , tmpCndi.variable );
            conditionsTblOriginal{ii,"cndFcast:avgValue"} = mean(tmpValues,1);
            conditionsTblOriginal{ii,"cndFcast:stdValue"} = nan;              
        
        end

    end

    varargout{1} = conditionsTblOriginal(:,["type","variable","h0","h1","value","cndFcast:avgValue","cndFcast:stdValue"]);

end

% varargout{1} = optimizerOutput;
