% S4-APP4 JdeLafontaine 17 juin 2008
% S5-APP6 JdeLafontaine 16 novembre 2019 (r�vision)
% 
% Int�grateur num�rique Runge Kutta

% (1) Effet de la 'reltol' de l'int�grateur sur le nombre de pas et la
%     pr�cision de la trajectoire (qui doit se refermer sur elle-m�me)
% (2) Variation du pas d'int�gration selon l'amplitude des deriv�es (pr�s de la Terre et de la Lune)
% (3) Utilisation de ODE45

% Etat: z = [x, y, dx/dt, dy/dt]

  clc
  clear
  close all

% Conditions initiales et temps final

  u = 1.0/(82.45);
  z0 = [1.2, 0.0, 0.0, -1.04935751];
  tspan = [0, 6.19];
  
% Place la Lune et la Terre sur le graphique

  plot((1-u),0,'rp', 'Markersize',5)
  hold on, plot(-u,0,'ro', 'Markersize',5)

% 1�re int�gration avec param�tres nominaux

  reltol1 = 1e-03;
  options = odeset('abstol' ,1e-06, 'reltol', reltol1);
  [t, z] = ode45('apollo', tspan, z0, options);

  plot(z(:,1), z(:,2),'.b', 'Markersize',10)
  Lz = length(z);
  plot(z(Lz,1), z(Lz, 2),'xb', 'Markersize',10)
  
% 2�me int�gration avec param�tres reltol = 1e-06

  reltol2 = 1e-06;
  options = odeset('abstol', 1e-06, 'reltol', reltol2);
  [t, z] = ode45('apollo', tspan, z0, options);
  
  plot(z(:,1), z(:,2),'.g', 'Markersize',10)
  Lz = length(z);
  plot(z(Lz,1), z(Lz, 2),'xg', 'Markersize',10)
  
  xlabel('Position X', 'Fontsize',15)
  ylabel('Position Y', 'Fontsize',15)
  title('Capsule Apollo', 'Fontsize',15)

  legend('Lune', 'Terre', ['Pr�cision: ', num2str(reltol1)], 'Fin trajectoire', ['Pr�cision: ', num2str(reltol2)], 'Fin trajectoire')
