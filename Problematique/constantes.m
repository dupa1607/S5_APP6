m = 50; %50 kg
J = 1.5; %kg/m^2
Rmars = 3397e03; %m
MUmars = 42830e09; %m^3/s^2
S = 0.8; %m^2
d = 0.05; % m
CD0 = 1.2;
CLalpha = 0.8;
CMalpha = -0.07;
CMq = -0.05;
CMdelta = 0.1;

v_ini = 6100; %m/s
gamma_ini = deg2rad(-20.5);%rad
h_ini = 120000; %m
s_ini = deg2rad(0);%rad
teta_ini = deg2rad(-80); %rad
q_ini = deg2rad(0); %rad/s

%conditions finales désirés
vfin = 250; %ou 300 m/s
hfin = 10000; % m

%Contraintes
dt_lim = 40; %s
Pdyn_max = 9500; %N/m^2
Pdyn_lim = 2650;
teta_cmd = 60; %degrés