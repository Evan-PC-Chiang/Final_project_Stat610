function Get_rates4inv(dict,Assign_slip)
    keysArray = keys(dict);
    
    for i = 1:length(keysArray)
        key = keysArray{i}; % Current key
        cellArray = dict(key); % Get the cell array of matrices
        
        % Initialize a sum matrix of the same size as the first matrix in the cell array
        matrixSize = size(cellArray{1});
        sumMatrix = zeros(matrixSize);
        
        % Loop through each matrix in the cell array and sum them
        for j = 1:length(cellArray)
            sumMatrix = sumMatrix + Assign_slip(j)*cellArray{j};
        end
        assignin('caller', key, sumMatrix);
    end
end