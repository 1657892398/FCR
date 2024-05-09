function mu = iterative_mu(w,x,m)
% w: Menbership matrix N*C
% x: Data set N*p
% m: Membership exponent
[~,p] = size(x);
[~,C] = size(w);
mu = zeros(p,C);
for j = 1:C
    mu(:,j) = x' * w(:,j).^m / sum( w(:,j).^m);
end
end

