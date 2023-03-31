clc
clear all
close all

format short G
load('accelero_Data_from_NASA')
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
Err_vit = h^2/12*(df_b - df_a);

%méthode de simpson pour faire l'intégrale de la vitesse
h_mes(1,1) = h_ini;
xs(1) = t(1);
for i = 3:2:N
    h_mes((i-1)/2 + 1,1) = h_ini + ((v_mes(1) + v_mes(i) + 4*sum(v_mes((2:2:i-1))) + 2*sum(v_mes(3:2:i-1))) * h/3);
    xs((i-1)/2 + 1) = t(i);
end
% Erreur de la méthode des trapèzes
d3f_a = (v_mes(4) - 3*v_mes(3) + 3*v_mes(2) - v_mes(1))/h^3;
d3f_b = (v_mes(N) - 3*v_mes(N-1) + 3*v_mes(N-2) - v_mes(N-3))/h^3;
Err_pos = h^4/180*(d3f_b - d3f_a);

% Identification des paramètres (boite grise) en trouvant la forme linéaire
gr_sin_gamma = -MUmars./(Rmars + h_mes).^2;
Daero = -m.*(acc_mes(1:2:end) + gr_sin_gamma);
Pdyn = Daero/(S * CD0);

X = h_mes;
Y = log((2.*Pdyn)./(v_mes(1:2:end).^2));

mat1 = [length(h_mes)       sum(X);
        sum(X)  sum(X.^2)];

mat2 = [sum(Y);
        sum(Y.*X)];

params = inv(mat1) * mat2;
hs = -1/params(2);
rho0 = exp(params(1));

P_obt = 0.5*(rho0*exp(-h_mes./hs)).*v_mes(1:2:end).^2;

% code pour calcul de l'erreur copié du laboratoire
E = sum((P_obt - Pdyn).^2);
N2 = length(h_mes);
RMS = sqrt(E/N2);
moy = sum(Pdyn)/N2;
R_2 = sum((P_obt - moy).^2)/sum((Pdyn - moy).^2);




%% figures
figure()
plot(t, acc_mes, 'o')

grid minor
xlabel('temps (s)')
ylabel('Accélération (m/s^2)')
legend('acc mesuré')

figure()
plot(t, v_mes, 'o')

grid minor
xlabel('temps (s)')
ylabel('Vitesse (m/s)')
legend('vit mesuré')

figure()
plot(xs, h_mes, 'o')

grid minor
xlabel('temps (s)')
ylabel('altitude (m)')
legend('altitude mesuré')

figure()
plot(t(1:2:end), Pdyn, 'o')
hold on
plot(t(1:2:end), P_obt)
grid minor
legend('Pdyn calc. avec Daero', 'Pdyn calc avec hs et rho0')