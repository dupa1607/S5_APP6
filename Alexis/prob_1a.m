% S4-APP4 JdeLafontaine 17 juin 2008
% S4-APP6 JdeLafontaine 18 juillet 2013 (r�vision)
% S5-APP6 JdeLafontaine 16 novembre 2019 (r�vision)
% 
% Probl�me no 1: Int�grateur num�rique,  stabilit� et convergence
% 
% Objectifs: 
% (a) savoir programmer et utiliser un int�grateur num�rique sur MATLAB
% (b) conna�tre les m�thodes d'int�gration explicite et pr�diction-correction
% (c) savoir reconna�tre une instabilit� num�rique et en trouver la cause
% (d) conna�tre et savoir utiliser le lien entre l'erreur de discr�tisation 
%     et le pas d'int�gration.

% Probl�me no 1a : �quation diff�rentielle marginalement stable qui est d�stabilis�e
% 				   par les erreurs de discr�tisation.
%
% Dans la partie 1a, on trouve la solution num�rique d'un probl�me dynamique
% avec 3 diff�rents int�grateurs num�riques:
% (1) Euler explicite
% (2) Euler implicite
% (3) Runge-Kutta
% On compare ensuite leur stabilit� et pr�cision en fonction du pas d'int�gration.

% NOTE: Les int�grateurs num�riques de type Euler sont rarement utilis�s en
% pratique. Ils ne sont utilis�s ici que pour des raisons p�dagogiques: leur
% mauvaise performance permet d'illuster plus facilement des probl�mes de
% stabilit� et de convergence parfois rencontr�s avec des int�grateurs
% num�riques plus performants quand les �quations diff�rentielles �
% int�grer comportent des probl�mes particuler ('stiffness' par exemple). 
% --------------------------------------------------------------------------------
% 
  clc
  close all  
  clear 

  figure(1)		% Commence nouvelle figure
  N = 10;		% Nombre de diff�rents pas d'int�gration
% dt = 2.0;		% Pas d'int�gration 
  dt = 1.0;		% Pas d'int�gration 
  y0 = [1 0];	% Conditions initiales
% y0 = [1];	    % Conditions initiales
  
% Boucle pour calculer les r�sultats num�riques pour chaque pas d'int�gration
  
  for n=1:N,
     
%    Int�gration avec m�thode de Euler     
     
     pas(n) = dt;									% pas d'int�gration � utiliser
     [t, yex]=eulerex('eqn1a',0, dt, 4.5, y0);		% int�gration avec Euler explicite
     disp([n, max(t)])
     [t, ypc]=eulerpc('eqn1a',0, dt, 4.5, y0);		% int�gration avec Euler pr�diction-correction
     disp([n, max(t)])

%    Integration avec m�thode Runge-Kutta 
%    (On utilise une grande tol�rance pour avoir un pas fixe �gal � dt)
     
     options = odeset('abstol',1000,'reltol',1000,'maxstep',dt);
     [t, yrk] = ode45('eqn1a',t,y0,options);
     
%    Solution exacte et indice de la valeur finale f
     
     y = cos(3*t);
%    y = y0*exp(2*t);
     f = length(t);
     
%    Graphique des erreurs de position et calcul de l'erreur 

     subplot(3,1,1)
     title('Comparaison: solutions num�rique et exacte')
     plot(t, [yex(:,1), y])
     ylabel('Y EX')
     hold on
     errex(n) = abs(yex(f,1)-y(f));
%    errex(n) = sqrt((yex(:,1)-y)'*(yex(:,1)-y));
     
     subplot(3,1,2)
     plot(t, [ypc(:,1), y])
     ylabel('Y PC')
     hold on
     errpc(n) = abs(ypc(f,1)-y(f));
%    errpc(n) = sqrt((ypc(:,1)-y)'*(ypc(:,1)-y));
     
     subplot(3,1,3)
     plot(t, [yrk(:,1), y])
     ylabel('Y RK')
     hold on
     errrk(n) = abs(yrk(f,1)-y(f));
%    errrk(n) = sqrt((yrk(:,1)-y)'*(yrk(:,1)-y));
     
%    R�duction du pas d'int�gration pour la prochaine boucle
       
     dt = dt/2;

  end

  figure(1)
  subplot(3,1,1), xlabel('Temps (s)'), grid, axis([0.0 4.5 -2 2])
  hold off
  subplot(3,1,2), xlabel('Temps (s)'), grid, axis([0.0 4.5 -2 2])
  hold off
  subplot(3,1,3), xlabel('Temps (s)'), grid, axis([0.0 4.5 -2 2])
  hold off
  
% Graphique des erreurs sur 2eme figure
  
  figure(2)

  subplot(3,1,1)
  plot(pas',errex') 
  title('Erreur en fonction du pas d''int�gration')
  ylabel('EULER EX')
  grid on
  
  subplot(3,1,2)
  plot(pas',errpc')
  ylabel('EULER PC')
  grid on
  
  subplot(3,1,3)
  plot(pas',errrk')
  ylabel('RK')
  xlabel('Pas d''int�gration')
  grid on

% Calcul de l'ordre des m�thodes avec les "vraies" erreurs, utilisant la solution analytique exacte  
  
% ordre_ex = log(errex(2:N)./errex(1:N-1))./log(pas(2:N)./pas(1:N-1));
% ordre_pc = log(errpc(2:N)./errpc(1:N-1))./log(pas(2:N)./pas(1:N-1));
% ordre_rk = log(errrk(2:N)./errrk(1:N-1))./log(pas(2:N)./pas(1:N-1));

  ordre_ex = log(errex(1:N-1)./errex(N))./log(pas(1:N-1)./pas(N));
  ordre_pc = log(errpc(1:N-1)./errpc(N))./log(pas(1:N-1)./pas(N));
  ordre_rk = log(errrk(1:N-1)./errrk(N))./log(pas(1:N-1)./pas(N));

%   ordre_ex = log(errex(1:N-1)./errex(2:N))./log(pas(1:N-1)./pas(2:N));
%   ordre_pc = log(errpc(1:N-1)./errpc(2:N))./log(pas(1:N-1)./pas(2:N));
%   ordre_rk = log(errrk(1:N-1)./errrk(2:N))./log(pas(1:N-1)./pas(2:N));
  
  figure(3)
  subplot(3,1,1)
  plot(pas(1:N-1)',ordre_ex')
  ylabel('EULER EX')
  title('Calcul de l''ordre de l''int�grateur: d�tection de convergence')
  axis([0.0 0.8 0 2])
  grid on
  
  subplot(3,1,2)
  plot(pas(1:N-1)',ordre_pc')
  ylabel('EULER PC')
  axis([0.0 0.8 0 3])
  grid on
  
  subplot(3,1,3)
  plot(pas(1:N-1)',ordre_rk')
  ylabel('ODE45')
  xlabel('Pas d''int�gration')
  axis([0.0 0.8 4.5 5.5])
  grid on

  figure(4)
  subplot(3,1,1)
  plot(log(pas(1:N-1)'/pas(N)),log(errex(1:N-1)'/errex(N)))
  ylabel('EULER EX')
  title('Erreur (rapport logarithmique)')
  grid on
  
  subplot(3,1,2)
  plot(log(pas(1:N-1)'/pas(N)),log(errpc(1:N-1)'/errex(N)))
  ylabel('EULER PC')
  grid on
  
  subplot(3,1,3)
  plot(log(pas(1:N-1)'/pas(N)),log(errrk(1:N-1)'/errex(N)))
  ylabel('ODE45')
  xlabel('Pas d''int�gration (rapport logarithmique)')
  grid on


