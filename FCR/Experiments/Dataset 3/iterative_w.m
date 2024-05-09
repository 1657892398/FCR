function w = iterative_w(x_tilde,x,y,b,mu,m,ga)
% x_tilde: Variables for regression N*p+1
% x: Independent variables N*p
% y: Dependent variable N*1
% b: Model parameter variables N*C
% mu: Cluster center matrix p*C
% m: Membership exponent
% ga: Gamma parameter

[N,~] = size(x);
[~,C] = size(b);
w = zeros(N,C);
for i = 1:N
    ww = zeros(1,C);
    for k = 1:C
        ww(k) = ga(1)*( y(i) - x_tilde(i,:)*b(:,k) )^2 + ga(2)*sum((x(i,:)'-mu(:,k)).^2);
    end  
        ww = 1./ww;
        ww = ww.^(1/(m-1));
        ww = sum(ww);
    for j = 1:C
        w(i,j) = 1 / ( ga(1)*(y(i) - x_tilde(i,:)*b(:,j))^2 + ga(2)*sum((x(i,:)'-mu(:,j)).^2) )^(1/(m-1))/ww;
    end
end
end

