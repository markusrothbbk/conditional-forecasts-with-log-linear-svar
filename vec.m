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

% This program create a vector from a matrix
% V = vec(X,rc)
% X must be a matrix 
% stacking rows rc=1 or columns rc=2.
function V = vec(X,rc)

if nargin == 1
    rc = 2;
end

if isempty(X)
    V = [];
else
    [a, b] = size(X);
    if rc==1
        for i=1:a
            V(1,(i-1)*b+1:i*b) = X(i,:);  %% stack rows
        end
        V = V';
    elseif rc==2
        for i=1:b
            V((i-1)*a+1:i*a,1)=X(:,i); %% stack columns
        end
    end
end
