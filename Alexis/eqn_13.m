function f = eqn_13(t,x)
Kp = 13;
Kd = 3.4;
xdes = 0;
u = Kp * (xdes - x(1)) - Kd * x(2) + x(1)^2;
f(1) = x(2);
f(2) = -0.6*x(2) - 3*x(1) - x(1)^2 + u;


f = f(:);