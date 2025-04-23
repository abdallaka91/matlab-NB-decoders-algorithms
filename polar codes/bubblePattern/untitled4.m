clc; clear; close all;

N = 50; % Size of each matrix (NxN)
num_levels = 6; % Number of levels (cells)

% Generate random matrices (1s and 0s) for demonstration
matrices = cell(1, num_levels);
for i = 1:num_levels
    num_matrices = 2^(i-1); % Number of matrices at this level
    matrices{i} = randi([0 1], N, N, num_matrices); % Random binary matrices
end

plot2D_bubble(matrices, N, num_levels);