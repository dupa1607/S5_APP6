clear all
close all
clc


u = 1.0/(82.45);
z0_1 = [-3 7.831];
z0_2 = [-3 7.832];
tspan = [0 20];

reltol1 = 1e-03;
options = odeset('abstol' ,1e-07, 'reltol', reltol1);
[t, z1] = ode45('prob14', tspan, z0_1, options);
reltol1 = 1e-03;
options = odeset('abstol' ,1e-07, 'reltol', reltol1);
[t, z2] = ode45('prob14', tspan, z0_2, options);

figure()
hold on
plot(t, z1(:,1));
plot(t, z2(:,1));

figure()
hold on
plot(t, z1(:,2));
plot(t, z2(:,2));