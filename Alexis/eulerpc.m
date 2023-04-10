  function [t,y]=eulerpc(fct,t0,dt,tf,y0)
%
% Fonction qui int�gre les �quations diff�rentielles "fct" selon la
% m�thode d'Euler pr�diction-correction, de t0 a tf avec pas dt.
%
% Utilisation: [t,y]=eulerpc(fct,y0,t0,dt,tf)
%
% Entr�es:
% fct : variable alpha-num�rique avec le nom de la fonction f(y,t) � int�grer
% y0  : condition initiale
% t0  : temps initial
% dt  : pas d'int�gration (fixe)
% tf  : temps final
%
% Sorties:
% y   : solution num�rique avec une rang�e par pas d'int�gration et 
%       une colonne par variables de sortie (nbre de col = ordre du syst�me).
% t   : matrice-colonne temps de t0 � tf par intervals dt
%--------------------------------------------------------------------
%
  t = [t0:dt:tf]';		% Calcul de la mat-col temps
  nstep=size(t)-1 ;     % Calcul du nombre de pas
  y(1,:)=y0(:)';			% Initialisation
  for n=1:nstep,			% Boucle pour propager l'�tat du syst�me
    fn = feval(fct,t(n),y(n,:))';
    yp = y(n,:)+dt*fn;			% Pr�diction
    fp = feval(fct,t(n+1),yp)';
    y(n+1,:)=y(n,:)+dt*(fn+fp)/2;	% Correction
  end
%








