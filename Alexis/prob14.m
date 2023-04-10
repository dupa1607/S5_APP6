function f = prob14(t,y)
  f(1) = y(2);
  f(2) = -0.6*y(2)-3*y(1) - y(1)^2;
  
  f = f(:);
