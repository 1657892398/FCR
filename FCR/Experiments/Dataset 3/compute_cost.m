function cost = compute_cost(x, x_t, y, w, b, mu, m, gamma, lambda)
% Calculate the cost function value for the Lasso model
% Input:
% x: Independent variables
% x_t: Variables for regression
% y: Dependent variable
% w: Membership
% b: Model parameter variables
% mu: Cluster center matrix
% m: Membership exponent
% gamma: Gamma parameters
% lambda: Penalty parameters
% Output:
% cost: Cost function value

[N, ~] = size(x);
[~, C] = size(b);
cost = 0;

for j = 1:C
    for i = 1:N
        L = gamma(1) * (y(i) - x_t(i, :) * b(:, j))^2 + gamma(2) * sum((x(i, :) - mu(:, j)').^2); 
        cost = cost + w(i, j)^m * L;
    end
    ma = max(abs(b(:, j)));
    cost = cost + lambda(j) * ma(1);
end

cost = cost / N;
end
