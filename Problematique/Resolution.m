clc
clear all
close all

format short G
load('Accelero_Data_from_NASA\Accelero_Data_from_NASA.mat')
constantes

%% Partie 1, voir photo tutorat 1

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
Err_alt = h^4/180*(d3f_b - d3f_a);

% Identification des paramètres (boite grise) en trouvant la forme linéaire
gr_sin_gamma = -MUmars./((Rmars + h_mes).^2);
% Daero = -m.*(acc_mes(1:2:end) + gr_sin_gamma);
Daero = -m.*(acc_mes(1:2:end));
Pdyn = Daero./(S * CD0);

X = h_mes;
% Y = log((2.*Pdyn)./(v_mes(1:2:end).^2));
Y = log((2.*Pdyn)./(decimate(v_mes,2).^2));
% Y = log((2*m*acc_mes(1:2:end))./(S*CD0.*vcc_mes(1:2:end).^2))
Acceleration = -acc_mes(1:2:end);
Vitesse = -v_mes(1:2:end);
Y = log((2.*m.*Acceleration)./(S.*CD0.*(Vitesse.^2)));



mat1 = [length(ts)  sum(X);
        sum(X)      sum(X.^2)];

mat2 = [sum(Y);
        sum(Y.*X)];

params = inv(mat1) * mat2;
hs = -1/params(2);
rho0 = exp(params(1));

P_obt = 0.5*(rho0*exp(-h_mes./hs)).*(v_mes(1:2:end).^2);
Daero_obt = P_obt * S *CD0;
acc_approx = -Daero_obt/m;

var = -1/hs;
b = log(rho0)
g = var * X + b;
E = sum((g - Y).^2);
N2 = length(Y);
RMS = sqrt(E/N2);
RMS_abs = sqrt(mean((acc_mes(1:2:end) - acc_approx).^2))
RMS_rel = sqrt(mean(((acc_mes(1:2:end) - acc_approx)./acc_mes(1:2:end)).^2))
moy = sum(Y)/N2;
R_2 = sum((g - moy).^2)/sum((Y - moy).^2);






% code pour calcul de l'erreur copié du laboratoire
% E = sum((P_obt - Pdyn).^2);
% N2 = length(h_mes);
% RMS = sqrt(E/N2);
% moy = sum(Pdyn)/N2;
% R_2 = sum((P_obt - moy).^2)/sum((Pdyn - moy).^2);


%% Partie 2

B = CD0*S/m;
rho_ini = rho0 * exp(-h_ini/hs);
v_RAA = v_ini*exp(-0.5*B*hs*(rho0*exp(-h_mes/hs) - rho_ini));



close all
figure()
plot(ts', v_RAA, 'o')
grid minor

%% Partie 3
close all
clc
B = CD0*S/m;
h_fin = 10000;
rho_ini = rho0 * exp(-h_ini/hs);
rho_fin = rho0 * exp(-h_fin/hs);
rho =  rho0 * exp(-h_ini./hs);
rfin = 10000 + Rmars;
r = Rmars + h_ini;
dv_aero_300 = 300 - sqrt(v_ini.^2 + 2 * MUmars.*(1/rfin - 1./r));



dv_aero_250 = 250 - sqrt(v_ini.^2 + 2 * MUmars.*(1/rfin - 1./r));
gamma_ref_300 = asind(0.5 * B * hs * (rho_fin-rho)./(log(1+dv_aero_300./v_ini)));
% gamma_ref_300 = asind(0.5*B*hs*(rho_fin-rho)./(log(300/v_ini)))
gamma_ref_250 = asind(0.5 .* B .* hs .* (rho_fin-rho)./(log(1+dv_aero_250./v_ini)));
% gamma_ref_250 = asind(0.5*B*hs*(rho_fin-rho)./(log(250/v_ini)))
h_simul = [120000:-1000:10000]';
% h_simul = 46112.5;
rho_simul = rho0 .* exp(-h_simul./hs);

v_simul_300 = v_ini*exp(0.5*B*hs*(rho_simul - rho_ini)/sind(gamma_ref_300));
Pdyn_simul_300 = 0.5*rho_simul .* v_simul_300.^2;

v_simul_250 = v_ini*exp(0.5*B*hs*(rho_simul - rho_ini)/sind(gamma_ref_250));
Pdyn_simul_250 = 0.5*rho_simul .* v_simul_250.^2;
% 
% figure()
% plot(h_simul, Pdyn_simul_300)
% grid minor
% title("pression dynamique")
% 
% figure()
% plot(h_simul, Pdyn_simul_300*S*CD0)
% grid minor
% title("Daero")
% 
% figure()
% plot(h_simul, v_simul_300)
% grid minor
% title('Vitesse simuler avec RAA et gravité')
Pdyn_max_300 = max(Pdyn_simul_300);
Pdyn_max_250 = max(Pdyn_simul_250);


%Newton-Raphson 300m/s
h_raph = 16000;

rho_raph = rho0 * exp(-h_raph/hs);
F = 0.5*rho0*exp(-h_raph./hs).*(v_ini*exp(0.5*B*hs*(rho_raph - rho_ini)./sind(gamma_ref_300))).^2*S*CD0 - 2650;
tol = 1e-10;
a = v_ini^2*S*CD0*0.5*rho0;
D = a*((-B*rho0.*exp(-h_raph./hs))./sind(gamma_ref_300) - 1/hs) .* exp(-h_raph./hs + (B*hs*(rho0.*exp(-h_raph./hs) - rho_ini))./sind(gamma_ref_300));
nb_iter = 0;
while (abs(F)) > tol & nb_iter < 501
    h_raph = h_raph - F/D;
    rho_raph = rho0 * exp(-h_raph/hs);
    F = 0.5*rho0*exp(-h_raph./hs).*(v_ini*exp(0.5*B*hs*(rho_raph - rho_ini)./sind(gamma_ref_300))).^2*S*CD0 - 2650;
    D = a*((-B*rho0.*exp(-h_raph./hs))./sind(gamma_ref_300) - 1/hs) .* exp(-h_raph./hs + (B*hs*(rho0.*exp(-h_raph./hs) - rho_ini))./sind(gamma_ref_300));
    nb_iter = nb_iter + 1;
end
disp(['Valeur de h1 (300): ', num2str(h_raph), ' en ', num2str(nb_iter),  ' itérations'])
rho_1_300 = rho0 .* exp(-h_raph./hs);
v1_300 = v_ini*exp(0.5*B*hs*(rho_1_300 - rho_ini)/sind(gamma_ref_300));
h1_300 = h_raph;
disp(['Valeur de v1 (300): ', num2str(v1_300)])


h_raph = 51000;
rho_raph = rho0 * exp(-h_raph/hs);
F = 0.5*rho0*exp(-h_raph./hs).*(v_ini*exp(0.5*B*hs*(rho_raph - rho_ini)./sind(gamma_ref_300))).^2*S*CD0 - 2650;
tol = 1e-10;
a = v_ini^2*S*CD0*0.5*rho0;
D = a*((-B*rho0.*exp(-h_raph./hs))./sind(gamma_ref_300) - 1/hs) .* exp(-h_raph./hs + (B*hs*(rho0.*exp(-h_raph./hs) - rho_ini))./sind(gamma_ref_300));
nb_iter = 0;
while (abs(F)) > tol & nb_iter < 501
    h_raph = h_raph - F/D;
    rho_raph = rho0 * exp(-h_raph/hs);
    F = 0.5*rho0*exp(-h_raph./hs).*(v_ini*exp(0.5*B*hs*(rho_raph - rho_ini)./sind(gamma_ref_300))).^2*S*CD0 - 2650;
    D = a*((-B*rho0.*exp(-h_raph./hs))./sind(gamma_ref_300) - 1/hs) .* exp(-h_raph./hs + (B*hs*(rho0.*exp(-h_raph./hs) - rho_ini))./sind(gamma_ref_300));
    nb_iter = nb_iter + 1;
end
disp(['Valeur de h (300): ', num2str(h_raph), ' en ', num2str(nb_iter),  ' itérations'])

rho_2_300 = rho0 .* exp(-h_raph./hs);
v2_300 = v_ini*exp(0.5*B*hs*(rho_2_300 - rho_ini)/sind(gamma_ref_300));
h2_300 = h_raph;
v_moy_300 = (v2_300 + v1_300)/2;
dist_300 = (h2_300 - h1_300)/sind(-gamma_ref_300);
dt_300 = dist_300 / v_moy_300;
disp(['Valeur de v1 (300): ', num2str(v2_300)])
fprintf('\n')
disp(['Temps limite obtenu: ', num2str(dt_300)])
fprintf('\n')







%Newton-Raphson 250m/s
h_raph = 16000;
rho_raph = rho0 * exp(-h_raph/hs);
F = 0.5*rho0*exp(-h_raph./hs).*(v_ini*exp(0.5*B*hs*(rho_raph - rho_ini)./sind(gamma_ref_250))).^2*S*CD0 - 2650;
tol = 1e-10;
a = v_ini^2*S*CD0*0.5*rho0;
D = a*((-B*rho0.*exp(-h_raph./hs))./sind(gamma_ref_250) - 1/hs) .* exp(-h_raph./hs + (B*hs*(rho0.*exp(-h_raph./hs) - rho_ini))./sind(gamma_ref_250));
nb_iter = 0;
while (abs(F)) > tol & nb_iter < 501
    h_raph = h_raph - F/D;
    rho_raph = rho0 * exp(-h_raph/hs);
    F = 0.5*rho0*exp(-h_raph./hs).*(v_ini*exp(0.5*B*hs*(rho_raph - rho_ini)./sind(gamma_ref_250))).^2*S*CD0 - 2650;
    D = a*((-B*rho0.*exp(-h_raph./hs))./sind(gamma_ref_250) - 1/hs) .* exp(-h_raph./hs + (B*hs*(rho0.*exp(-h_raph./hs) - rho_ini))./sind(gamma_ref_250));
    nb_iter = nb_iter + 1;
end
disp(['Valeur de h (250): ', num2str(h_raph), ' en ', num2str(nb_iter),  ' itérations'])
rho_1_250 = rho0 .* exp(-h_raph./hs);
h1_250 = h_raph;
v1_250 = v_ini*exp(0.5*B*hs*(rho_1_250 - rho_ini)/sind(gamma_ref_250));
disp(['Valeur de v1 (250): ', num2str(v1_250)])


h_raph = 51000;
rho_raph = rho0 * exp(-h_raph/hs);
F = 0.5*rho0*exp(-h_raph./hs).*(v_ini*exp(0.5*B*hs*(rho_raph - rho_ini)./sind(gamma_ref_250))).^2*S*CD0 - 2650;
tol = 1e-10;
a = v_ini^2*S*CD0*0.5*rho0;
D = a*((-B*rho0.*exp(-h_raph./hs))./sind(gamma_ref_250) - 1/hs) .* exp(-h_raph./hs + (B*hs*(rho0.*exp(-h_raph./hs) - rho_ini))./sind(gamma_ref_250));
nb_iter = 0;
while (abs(F)) > tol & nb_iter < 501
    h_raph = h_raph - F/D;
    rho_raph = rho0 * exp(-h_raph/hs);
    F = 0.5*rho0*exp(-h_raph./hs).*(v_ini*exp(0.5*B*hs*(rho_raph - rho_ini)./sind(gamma_ref_250))).^2*S*CD0 - 2650;
    D = a*((-B*rho0.*exp(-h_raph./hs))./sind(gamma_ref_250) - 1/hs) .* exp(-h_raph./hs + (B*hs*(rho0.*exp(-h_raph./hs) - rho_ini))./sind(gamma_ref_250));
    nb_iter = nb_iter + 1;
end
disp(['Valeur de h (250): ', num2str(h_raph), ' en ', num2str(nb_iter),  ' itérations'])
rho_2_250 = rho0 .* exp(-h_raph./hs);
v2_250 = v_ini*exp(0.5*B*hs*(rho_2_250 - rho_ini)/sind(gamma_ref_250));
h2_250 = h_raph;

v_moy_250 = (v2_250 + v1_250)/2;
dist_250 = (h2_250 - h1_250)/sind(-gamma_ref_250);
dt_250 = dist_250 / v_moy_250;
disp(['Valeur de v2 (250): ', num2str(v2_250)])
fprintf('\n')
disp(['Temps limite obtenu: ', num2str(dt_250)])
fprintf('\n')

%% Actionneur translation
Kp =4;
gamma_dot = tf([1], [1/Kp 1]);
% figure()
% step(gamma_dot)


%% Actionneur rotation

Kp = 400;
Kd = 28;

ft = tf([Kp Kd], [1 Kd Kp])
% figure()
% step(ft)

%% Validation
clc 
close all

myfile = fullfile(tempdir,'g_ref_250_rt.mat');
matObj = matfile(myfile,'Writable',true);

  z0 = [v_ini gamma_ini h_ini s_ini teta_ini q_ini 0];
  tspan = [0, 150];

  abs_tol = 1e-10;
  options = odeset('abstol' ,abs_tol);
  [tv, z] = ode45('eqn_dyn', tspan, z0, options);

[~,gamma_ref_250_rt] = cellfun(@(t,x)  eqn_dyn(t,x.'), num2cell(tv), num2cell(z,2),'uni',0);
gamma_ref_250_rt = rad2deg(cell2mat(gamma_ref_250_rt));
% [~,gamma_ref_300_rt] = cellfun(@(t,x)  eqn_dyn(t,x.'), num2cell(tv), num2cell(z,2),'uni',0);
% gamma_ref_300_rt = rad2deg(cell2mat(gamma_ref_300_rt));
%   graph_sans_asser
  graph_asser_250
%   graph_asser_300

%   disp(['Temps > 2650 = ', num2str(z(end,7)), ' s']);

% gamma_ref_250_rt = readmatrix('gamma_ref_250_rt.txt');
% gamma_ref_300_rt = readmatrix('gamma_ref_300_rt.txt');
  







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
hold on
plot(ts, v_RAA)

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


figure()
plot(Y)
hold on
plot(g)