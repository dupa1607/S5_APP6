clc
clear all
close all

format short G
load('Accelero_Data_from_NASA\Accelero_Data_from_NASA.mat')
constantes

acc_mes = -acc_mes;

% Méthode des trapèzes pour faire l'intégrale de l'accélération
N = length(t);
h = t(2) - t(1);
v_mes(1,1) = v_ini;
for i = 2:N
    v_mes(i,1) = v_ini + ((acc_mes(1) + acc_mes(i) + 2*sum(acc_mes(2:i-1))) * h/2);
end
% Erreur de la méthode des trapèzes
df_a = (acc_mes(2) - acc_mes(1))/h;
df_b = (acc_mes(N-1) - acc_mes(N))/h;
Err_vit = (h^2/12)*(df_b - df_a);

%méthode de simpson pour faire l'intégrale de la vitesse
h_mes(1,1) = h_ini;
ts(1) = t(1);
for i = 3:2:N
    h_mes((i-1)/2 + 1,1) = h_ini - ((v_mes(1) + v_mes(i) + 4*sum(v_mes((2:2:i-1))) + 2*sum(v_mes(3:2:i-1))) * (h/3));
    ts((i-1)/2 + 1) = t(i);
end
ts = ts';
% Erreur de la méthode simpson
d3f_a = (v_mes(4) - 3*v_mes(3) + 3*v_mes(2) - v_mes(1))/h^3;
d3f_b = (v_mes(N) - 3*v_mes(N-1) + 3*v_mes(N-2) - v_mes(N-3))/h^3;
Err_pos = h^4/180*(d3f_b - d3f_a);

% Identification des paramètres (boite grise) en trouvant la forme linéaire
gr_sin_gamma = -MUmars./((Rmars + h_mes).^2);
Daero = -m.*(acc_mes(1:2:end) + gr_sin_gamma);
Pdyn = Daero./(S * CD0);

X = h_mes;
Y = log((2.*Pdyn)./(v_mes(1:2:end).^2));

mat1 = [length(ts)  sum(X);
        sum(X)      sum(X.^2)];

mat2 = [sum(Y);
        sum(Y.*X)];

params = inv(mat1) * mat2;
hs = -1/params(2);
rho0 = exp(params(1));

P_obt = 0.5*(rho0*exp(-h_mes./hs)).*(v_mes(1:2:end).^2);
Daero_obt = P_obt * S *CD0;
acc_approx = -Daero_obt/m - gr_sin_gamma;
%Les données semblent être shiftés d'une case ???


% code pour calcul de l'erreur copié du laboratoire
E = sum((P_obt - Pdyn).^2);
N2 = length(h_mes);
RMS = sqrt(E/N2);
moy = sum(Pdyn)/N2;
R_2 = sum((P_obt - moy).^2)/sum((Pdyn - moy).^2);




%% figures
% ts = ts + 1; % Il semble y avoir une erreur de shift de 1 dans les données
figure()
plot(t, acc_mes, 'o')
hold on
plot(ts, acc_approx)

grid minor
xlabel('temps (s)')
ylabel('Accélération (m/s^2)')
legend('acc mesuré', 'acc obtenu')
title('Accélératione du test NASA en fonction du temps')


figure()
plot(t, v_mes, 'o')

grid minor
xlabel('temps (s)')
ylabel('Vitesse (m/s)')
legend('vit mesuré')

figure()
plot(ts, h_mes, 'o')

grid minor
xlabel('temps (s)')
ylabel('altitude (m)')
legend('altitude mesuré')

%Les données semblent être shiftés d'une case ???
% P_obt = circshift( P_obt , 1 )
figure()
plot(ts, Pdyn, 'o')
hold on
plot(ts, P_obt)
grid minor
legend('Pdyn calc. avec Daero', 'Pdyn calc avec hs et rho0')
title('Pression dynamique comparaison')
