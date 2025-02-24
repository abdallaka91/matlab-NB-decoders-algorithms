function [bestOrder, sp_frwd] = simulatedAnnealing(n, costFunction, maxIter, initialTemp, coolingRate, currentOrder, wght)

    % n: Number of elements (64 in your case)
    % costFunction: Handle to your function that evaluates an order
    % maxIter: Maximum iterations
    % initialTemp: Starting temperature
    % coolingRate: Factor to reduce temperature each step

    % Generate initial random order
%     currentOrder = randperm(n);
    [currentCost, sp_frwd] = costFunction(currentOrder, wght);

    % Store the best solution found
    bestOrder = currentOrder;
    bestCost = currentCost;

    % Set initial temperature
    T = initialTemp;

    for iter = 1:maxIter
        % Swap two random elements to generate a new order
        newOrder = currentOrder;
        swapIdx = randperm(n, 2); % Pick two indices to swap
        newOrder([swapIdx(1), swapIdx(2)]) = newOrder([swapIdx(2), swapIdx(1)]);

        % Compute cost of the new order
        [newCost, sp_frwd] = costFunction(newOrder, wght);

        % Acceptance probability (Metropolis criterion)
        if newCost >= currentCost || rand < exp((newCost - currentCost) / T)
            % Accept the new solution
            currentOrder = newOrder;
            currentCost = newCost;
            tab(iter) = newCost;

            % Update best solution
            if newCost >= bestCost
                bestOrder = newOrder;
                bestCost = newCost;
                disp(sp_frwd);
            end
        end
        tab_best(iter) = bestCost;


        % Cool down temperature
        T = T * coolingRate;

        % Display progress
        if mod(iter,100)==0
        fprintf('Iteration %d: Best Cost = %.4f\n', iter, bestCost);
        
        end
    end
    stem(tab); hold on; plot(tab_best); grid on; hold off;
end
