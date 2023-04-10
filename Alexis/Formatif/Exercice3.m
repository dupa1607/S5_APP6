clear all
close all
clc

zeta = (2^0.5)/2;
ts = 4;
wn = 4/(zeta*ts);

pdes = -zeta*wn + i*wn*(1-zeta^2)^0.5;

value = conv([1 -(-1+1i)],[1 -(-1-1i)]);
value2 = conv(value, [1 -(-3)]);
