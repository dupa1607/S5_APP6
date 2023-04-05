close all
[~, index] = min(abs(z(:,3) - hfin));
tfin =  tv(index);

figure()
plot(tv, rad2deg(z(:,2)), "LineWidth",2)
grid minor
title('angle de vol sans asservissment', "FontSize",20)
xlabel("temps(s)", "FontSize",15)
ylabel("Angle de vol gamma (degrés)", "FontSize",15)
xline(tfin, '--')

figure()
plot(z(:,3) / 1000, z(:,1), "LineWidth",2)
grid minor
title('vitesse en fonction de la hauteur sans asservissement', "FontSize",20)
xlabel("hauteur(km)", "FontSize",15)
ylabel("vitesse (m/s)", "FontSize",15)
xline(z(index,3)/1000, '--')
text(z(index,3) / 1000,z(index,1),['\leftarrow v = ',num2str(z(index,1)), ' m/s' ])

figure()
plot(tv, rad2deg(z(:,5)))
hold on
plot(tv, rad2deg(z(:,5)) - rad2deg(z(:,2)))
grid minor
title('angle de vol et d''attaque sans asservissment', "FontSize",20)
xlabel("temps(s)", "FontSize",15)
ylabel("Angle d'attaque alpha (degrés)", "FontSize",15)
xline(tfin, '--')
legend('Angle de tangage teta', 'Angle d''attaque alpha')

figure()
plot(tv, rad2deg(z(:,6)), "LineWidth",2)
grid minor
title('vitesse angulaire de tangage q sans asservissment', "FontSize",20)
xlabel("temps(s)", "FontSize",15)
ylabel("vitesse angulaire de tangage (rad/s)", "FontSize",15)
xline(tfin, '--')

figure()
rho_valid = rho0 * exp(-z(:,3)./hs);
Pdyn_valid = 0.5*rho_valid.*z(:,1).^2;
plot(tv, Pdyn_valid, "LineWidth",1.25)
hold on
Daero_valid = Pdyn_valid*S*CD0;
plot(tv, Daero_valid, "LineWidth",1.25)
[M,I] = max(Pdyn_valid);
text(tv(I),Pdyn_valid(I)+100,['\leftarrow Pdyn max = ',num2str(M), ' N/m^2' ])
grid minor
title('Pression dynamique et trainée aerodynamique sans asservissment', "FontSize",20)
xlabel("temps(s)", "FontSize",15)
ylabel("Pression dynamique et trainée aerodynamique", "FontSize",15)
xline(tfin, '--')
yline(2650, '--')
legend('Pression dynmaique', 'trainée aérodynamique')

figure()
pl = plot(tv, z(:,7), "LineWidth",2);
grid minor
title('Intégrale de Daero', "FontSize",20)
xlabel("temps (s)", "FontSize",15)
ylabel("temps > 2650N (s)", "FontSize",15)
xline(tfin, '--')
dt = datatip(pl,tv(index),z(end,7));




