function [p_Si_given_y, py_given_Si, p_y] = Prb_Si_given_y(y, N0, norm_const)


M = numel(norm_const);
N = length(y);
p_Si_given_y = zeros(M,N);
py_given_Si = zeros(M, N);
p_y = zeros(M, 1);
p_Si = 1/M;

for j = 1 : N
    for i = 1:M
        py_given_Si(i,j) = (1/(pi * N0)) * exp(-abs(y(j) - norm_const(i))^2 / N0);  % Likelihood
        
    end
    py_given_Si(:,j) = py_given_Si(:,j)/sum(py_given_Si(:,j));
    p_y(j) = sum(py_given_Si(:,j) * p_Si);
end
for j = 1 : N
    p_Si_given_y(:,j) = (py_given_Si(:,j) * p_Si) / p_y(j);
end

