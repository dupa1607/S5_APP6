clear all
close all
clc


epsilon = 10e-8;
x = 5;
F = 80*exp(-x/12)-3*x*sin(x/7)-8;
D = (-80/12)*exp(-x/12)-(3*x/7)*cos(x/7)-3*sin(x/7);
i = 0;
while epsilon < abs(F) && i<501
    x = x-F/D;
    F = 80*exp(-x/12)-3*x*sin(x/7)-8;
    D = (-80/12)*exp(-x/12)-(3*x/7)*cos(x/7)-3*sin(x/7);
    i=i+1;
end
finalX = x;
figure
hold on
x = [0:0.01:20];
plot(x, 80*exp(-x./12))
plot(x, 3*x.*sin(x./7)+8)

figure
hold on
x = [0:0.01:20];
plot(x, 80*exp(-x./12)-3*x.*sin(x./7)-8)
