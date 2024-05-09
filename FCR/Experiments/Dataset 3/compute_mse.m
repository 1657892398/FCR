function mse = compute_mse(x, y, w, b)
% Calculate the loss function value for the Lasso model
% Input:
% x: Independent variables
% y: Dependent variable
% w: Membership
% b: Model parameter variables
% Output:
% mse: Mean squared error

[N, ~] = size(x);
[~, C] = size(b);
mse = 0;

for i = 1:N
    yy = 0;
    for j = 1:C
        yy = yy + w(i, j) * x(i, :) * b(:, j);
    end
    mse = mse + (y(i) - yy).^2;
end

mse = mse / N;
end
