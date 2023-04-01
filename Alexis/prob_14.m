clc
clear all
close all


z0 = [-3 7.831];
tspan = [0 20];
options = odeset('abstol' ,1e-02);
[t,z] = ode45('eqn_13',tspan,z0, options);

figure(1)
plot(t,z(:,1));
hold on
figure(2)
plot(t,z(:,2));
hold on

figure(3)
plot(z(:,1), z(:,2),'.g', 'Markersize',10)
hold on

z0 = [-3, 7.831];
tspan = [0 20];
[t,z] = ode45('eqn_13',tspan,z0, options);

figure(3)
plot(z(:,1), z(:,2),'.b', 'Markersize',10)
legend('1', '2')


figure(2)
plot(t,z(:,2));
figure(1)
plot(t,z(:,1));
legend('7.831', '7.832')

% Kp = 13;
% Kd = 3.4;

% ft = tf([Kp Kd], [1 0.6+Kd 3+Kp])
% Kp = 16;
% Kd = 4;
% ft2 = tf([Kp Kd], [1 0.6+Kd 3+Kp])
% figure()
% step(ft)
% hold on
% step(ft2)
% legend('lin juste x^2', 'lin tout')
% 
