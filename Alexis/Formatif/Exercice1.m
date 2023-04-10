

%% y1m
clear all
close all
clc
load('APP6e_Formatif_Probleme1.mat');

Y1 = log(y1m);
X1 = xm;
[m,b] = moindreCarre(X1,Y1);
tau = -1/m;
K = exp(b);
yapprox = K*exp(-xm/tau);

figure()
hold on
plot(xm, y1m, 'o');
plot(xm, yapprox)

[R2,RMSEa, RMSEr] = ErrorCalculator(yapprox,y1m);

%% y2m
clear all
close all
clc
load('APP6e_Formatif_Probleme1.mat');

Y1 = y2m;
X1 = xm.*sin(xm./7);
[m,b] = moindreCarre(X1,Y1);
beta = b;
alpha = m;
yapprox = alpha.*xm.*sin(xm./7)+beta;

figure()
hold on
plot(xm, y2m, 'o');
plot(xm, yapprox)

[R2,RMSEa, RMSEr] = ErrorCalculator(yapprox,y2m);