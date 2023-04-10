% S4-APP4 JdeLafontaine 17 juin 2008
% S4-APP6 JdeLafontaine 18 juillet 2013 (r�vision)
% S5-APP6 JdeLafontaine 16 novembre 2019 (r�vision)
% 
% Probleme no 1: Int�grateur num�rique,  stabilit� et convergence
% 
% Objectifs: 
% (a) savoir programmer et utiliser un int�grateur num�rique sur MATLAB
% (b) conna�tre les m�thodes d'int�gration explicite et pr�diction-correction
% (c) savoir reconna�tre une instabilit� num�rique et en trouver la cause
% (d) conna�tre et savoir utiliser le lien entre l'erreur de discr�tisation 
%     et le pas d'integration.

% Probl�me no 1c : �quation diff�rentielle instable.
%
% Dans cette partie 1c, on trouve la solution num�rique d'un probl�me dynamique de 
% nature instable. On utilise 3 diff�rents pas d'int�gration. 

% NOTE: Les int�grateurs num�riques de type Euler sont rarement utilis�s en
% pratique. Ils ne sont utilis�s ici que pour des raisons p�dagigiques: leur
% mauvaise performance permet d'illuster plus facilement des probl�mes de
% stabilit� et de convergence parfois rencontr�s avec des int�grateurs
% num�riques plus performants quand les �quations diff�rentielles �
% int�grer comportent des probl�mes particuler ('stiness' par exemple). 

  clc
  close all
  clear
  
  pas = [0.002, 0.001, 0.0005];
  for n=1:3,
     [t, yex]=eulerex('eqn1c',0, pas(n), 3.0 ,0.08);
     [t, ypc]=eulerpc('eqn1c',0, pas(n), 3.0 ,0.08);
     
     y = t.^2 + 0.4*t + 0.08;										% Solution exacte

     figure(1)							
     subplot(3,2,2*n-1)
     plot(t, [yex, y])
     if(n==1), title('EULEREX'), end
     ylabel(['Pas = ', num2str(pas(n))])

     subplot(3,2,2*n)
     plot(t, [ypc, y])
     if (n==1), title('EULERPC'), end
     ylabel(['Pas = ', num2str(pas(n))])
  end
  figure(1), subplot(3,2,1), legend('EulerEX', 'Exacte')
  figure(1), subplot(3,2,2), legend('EulerPC', 'Exacte')
  figure(1), subplot(3,2,5), xlabel('Temps (s)'), subplot(3,2,6), xlabel('Temps (s)')




