% S4-APP4 JdeLafontaine 17 juin 2008
% S4-APP6 JdeLafontaine 18 juillet 2013 (r�vision)
% S5-APP6 JdeLafontaine 16 novembre 2019 (r�vision)

% Probl�me no 1: Int�grateur num�rique,  stabilit� et convergence
% 
% Objectifs: 
% (a) savoir programmer et utiliser un int�grateur num�rique sur MATLAB
% (b) conna�tre les m�thodes d'int�gration explicite et pr�diction-correction
% (c) savoir reconna�tre une instabilit� num�rique et en trouver la cause
% (d) conna�tre et savoir utiliser le lien entre l'erreur de discr�tisation 
%     et le pas d'integration.

% Probl�me no 1b : �quation diff�rentielle stable mais "stiff" (ses deux valeurs propres
%                  sont tr�s diff�rentes). L'int�grateur num�rique est 
%                  conditionnellement stable (d�pend du pas d'int�gration).
%
% Dans cette partie 1b, on trouve la solution num�rique d'un probl�me dynamique de 
% nature "stiff" i.e. avec des valeurs propres tr�s diff�rentes. On utilise 5 diff�rents
% pas d'int�gration et on observe ensuite la stabilit� en fonction du pas d'int�gration.

% � noter que pour le pas = 0.002, la stabilit� est marginale (p�les sur cercle unitaire)

% NOTE: Les int�grateurs num�riques de type Euler sont rarement utilis�s en
% pratique. Ils ne sont utilis�s ici que pour des raisons p�dagigiques: leur
% mauvaise performance permet d'illuster plus facilement des probl�mes de
% stabilit� et de convergence parfois rencontr�s avec des int�grateurs
% num�riques plus performants quand les �quations diff�rentielles �
% int�grer comportent des probl�mes particuler ('stiness' par exemple). 
  
  clear errex errpc
  close all
  clc

  pas = [0.010, 0.005, 0.002, 0.0015, 0.0010];
  for n=1:5,
     [t, yex]=eulerex('eqn1b',0, pas(n), 0.05,[1 1]);
     [t, ypc]=eulerpc('eqn1b',0, pas(n), 0.05,[1 1]);
     y1 = 4*exp(-t) - 3*exp(-1000*t);
     y2 =-2*exp(-t) + 3*exp(-1000*t);

     figure(1)							% Sortie 1
     subplot(5,2,2*n-1)
     plot(t, [yex(:,1), y1])
     ylabel(['Pas = ', num2str(pas(n))])
     if(n==1), title('EULEREX'), end
     subplot(5,2,2*n)
     plot(t, [ypc(:,1), y1])
     ylabel(['Pas = ', num2str(pas(n))])
     if (n==1), title('EULERPC'), end
   
     figure(2)							% Sortie 2
     subplot(5,2,2*n-1)
     plot(t, [yex(:,2), y2])
     ylabel(['Pas = ', num2str(pas(n))])
     if (n==1), title('EULEREX'), end
     subplot(5,2,2*n)
     plot(t, [ypc(:,2), y2])
     ylabel(['Pas = ', num2str(pas(n))])
     if (n==1), title('EULERPC'), end
 
  end

  figure(1), subplot(5,2,1), legend('EulerEX', 'exacte')
  figure(1), subplot(5,2,2), legend('EulerPC', 'exacte')

  figure(2), subplot(5,2,1), legend('EulerEX', 'exacte')
  figure(2), subplot(5,2,2), legend('EulerPC', 'exacte')



