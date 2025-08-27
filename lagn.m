% Copyright (C) 2025 Frieder Mokinski and Markus Roth
%  
% This program is free software: you can redistribute it and/or modify  
% it under the terms of the GNU General Public License as published by  
% the Free Software Foundation, either version 3 of the License, or  
% (at your option) any later version.  
%  
% This program is distributed in the hope that it will be useful,  
% but WITHOUT ANY WARRANTY; without even the implied warranty of  
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  
% GNU General Public License for more details.  
%  
% You should have received a copy of the GNU General Public License  
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

function lagged = lagn(data, n)  
%LAGN Create lagged (and/or leaded) versions of a data matrix or table for given lags.  
%  
%   lagged = lagn(data, n)  
%  
%   INPUTS:  
%       data : [T x N] numeric matrix or table (non-empty)  
%           The input time series data, where T is the number of observations  
%           and N is the number of variables (columns).  
%       n    : vector of integers  
%           The lags to be applied. Positive values correspond to standard lags  
%           (shift data down), negative values to leads (shift data up), and zero  
%           returns the original data.  
%  
%   OUTPUT:  
%       lagged : [T x (N*numel(n))] matrix or table  
%           The lagged data, same type as input.  
%  
%   EXAMPLE:  
%       A = (1:5)';  
%       lagn(A, [0 1 2])  
%       % Returns a 5x3 matrix with columns: original, lag 1, lag 2  
%  
%       T = array2table(rand(5,2), 'VariableNames', {'GDP','CPI'});  
%       lagn(T, [0 1])  
%       % Returns a table with columns: GDP_lag0, GDP_lag1, CPI_lag0, CPI_lag1  
%  
%   NOTES:  
%       - For lags that exceed the available data, NaN values are used for padding.  
%       - Negative lags correspond to leads (future values).  
%       - Argument validation requires MATLAB R2019b or newer.  
  
arguments  
    data {mustBeNonempty}  
    n {mustBeInteger, mustBeVector}  
end  
  
% Convert table to array if needed, and get variable names  
if istable(data)  
    dataArray = table2array(data);  
    varNames = data.Properties.VariableNames;  
else  
    dataArray = data;  
    % Generate default variable names for matrix input  
    varNames = matlab.lang.makeUniqueStrings("Var" + (1:size(data,2)));  
end  
  
n = n(:)'; % Ensure n is a row vector  
ncol_ = size(dataArray,2);  
T = size(dataArray,1);  
laggedArray = NaN(T, ncol_ * numel(n));  
laggedNames = strings(1, ncol_ * numel(n));  
  
for i_ = 1:numel(n)  
    n_ = n(i_);  
    if n_ >= 0  
        lagged_ = [NaN(n_, ncol_); dataArray(1:end-n_, :)];  
    else  
        lagged_ = [dataArray(1-n_:end, :); NaN(-n_, ncol_)];  
    end  
    laggedArray(:, (i_-1)*ncol_+1 : i_*ncol_) = lagged_;  
    % Create variable names  
    for j = 1:ncol_  
        laggedNames((i_-1)*ncol_+j) = varNames{j} + "_lag" + n_;  
    end  
end  
  
% Return as table if input was table, else as array  
if istable(data)  
    lagged = array2table(laggedArray, 'VariableNames', cellstr(laggedNames));  
else  
    lagged = laggedArray;  
end  
  
end  