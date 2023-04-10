clear all
close all
clc



ci = [-10 -5.0 10.0];
t = [0:0.01:10];


options = odeset('RelTol',1e-3,'AbsTol',1e-6);
[t,x] = ode45('eq4a',t,ci, options);

figure(1)
hold on
plot(t, x(:,1));
plot(t, x(:,2));
plot(t, x(:,3));
hold off

figure(2)
plot3(x(1,1), x(1,2), x(1,3),'go','Markersize',10) % début
hold on
plot3(x(:,1), x(:,2), x(:,3),'Linewidth',2) % trajectoire
plot3(x(end,1), x(end,2), x(end,3),'rx','Markersize',10) % fin
grid on
xlabel('Position', 'Fontsize',15)
ylabel('Vitesse' , 'Fontsize',15)
zlabel('Accélération', 'Fontsize',15)
legend('Début', 'Trajectoire', 'Fin')