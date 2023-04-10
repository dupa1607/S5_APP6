function [m,b] = moindreCarre(x,y)
    X = [sum(x.^0) sum(x.^1); sum(x.^1) sum(x.^2)];
    Y = [sum(y); sum(y.*x)];
    A = inv(X)*Y;
    b = A(1);
    m = A(2);
end