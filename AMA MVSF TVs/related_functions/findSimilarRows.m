function similarGroups = findSimilarRows(matrix)
    % Number of rows in the matrix
    numRows = size(matrix, 1);  
    visited = false(numRows, 1); % Array to keep track of visited rows
    similarGroups = {}; % Cell array to store groups of similar rows
    
    % Loop through each row in the matrix
    for i = 1:numRows
        if visited(i)
            continue; % Skip if already grouped
        end
        
        % Start a new group with the current row
        currentGroup = i;
        
        % Compare the current row with all subsequent rows
        for j = i+1:numRows
            if ~visited(j) && isequal(matrix(i, :), matrix(j, :))
                currentGroup = [currentGroup, j]; % Add similar row to the group
                visited(j) = true; % Mark the row as visited
            end
        end
        
        % Mark the current row as visited
        visited(i) = true; 
        % Add the group to the list of similar groups
        similarGroups{end+1} = currentGroup; 
    end
end
