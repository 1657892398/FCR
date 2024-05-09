function b1 = iterative_b(x,y,w,b,m,lamda)
% x: Independent variables N*p
% y: Dependent variable N*1
% w: Menbership matrix N*C
% b: Model parameter variables N*C
% m: Membership exponent
% lambda: Lasso parameter

w = w.^m;
b1 = b;
b2 = 100.*ones(size(b));
k = 0;
[~,n] = size(x);% n = p+1
[~,C] = size(b);
while norm(b1-b2)>0.0001
    b2 = b1;
for p = 1:C
    for t = 1:n
        b1(t,p) = 0;
        fenmu = w(:,p)'*x(:,t).^2;
        fenzi1 = sum(w(:,p).*y.*x(:,t));
        xx = x*b1;
        fenzi2 = sum( xx(:,p).*w(:,p).*x(:,t) );
        fenzi = fenzi1-fenzi2;
        if fenzi - lamda(p)/2 >0
            b1(t,p) = ( fenzi - lamda(p)/2 )/fenmu;
        elseif fenzi + lamda(p)/2 <0
            b1(t,p) = ( fenzi + lamda(p)/2 )/fenmu;
        else
            b1(t,p) = 0;
        end
    end
end
k = k+1;
if k > 10000
    print('error')
    break;
end
end

