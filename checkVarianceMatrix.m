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

function is_admissible = checkVarianceMatrix(matrix)
    % Check if matrix is square
    if size(matrix,1) ~= size(matrix,2)
        is_admissible = false;
        return;
    end
    
    % Check if matrix is symmetric
    if ~isequal(matrix, matrix')
        is_admissible = false;
        return;
    end
    
    % Check if matrix is positive semi-definite (all eigenvalues are non-negative)
    if any(eig(matrix) < 0)
        is_admissible = false;
        return;
    end
    
    is_admissible = true;
end