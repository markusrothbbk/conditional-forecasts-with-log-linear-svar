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