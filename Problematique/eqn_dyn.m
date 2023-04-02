function f = eqn_dyn(t,x)
vfin = 250;
ctl = 1; % 1 = asservissement, 0 = pas d'asservissement


v = x(1);
gamma = x(2);
h = x(3);
s = x(4);
teta = x(5);
q = x(6);
constantes

Kp_rot = 400;
Kd_rot = 28;
Kp_trans = 4;
rho0 = 0.021403106629152;
hs = 1.101715085089048e+04;
rho = rho0 * exp(-h/hs);
Pdyn = 0.5*rho * v^2;
Daero = Pdyn * S * CD0;
r = Rmars + h;
gr = MUmars/r^2;
alpha = teta - gamma;
B = CD0*S/m;


if ctl == 0;
    delta = 0;
else
    rfin = Rmars + hfin;
    rho_fin = rho0 * exp(-hfin/hs);
    dv_aero = vfin - sqrt(v^2 + 2 * MUmars*(1/rfin - 1/r));
    gamma_ref = asin(0.5 * B * hs * (rho_fin-rho)/(log(1+dv_aero/v)));

    teta_eq = -((-Pdyn * S * CLalpha * gamma)/(v*m) + (v/r - (MUmars)/(v*r^2))*cos(gamma)) / ((Pdyn * S * CLalpha)/(v*m));
    teta_cmd = (Kp_trans * (gamma_ref - gamma))/((Pdyn * S * CLalpha)/(v*m));
    teta_des = teta_cmd + teta_eq;
%     teta_des = min(60, max(-60, teta_des));
    if teta_des > deg2rad(60)
        teta_des = deg2rad(60);
    else 
        if teta_des < deg2rad(-60)
            teta_des = deg2rad(-60);
        end
    end
    
    f_delta = (Pdyn*S*d/J)*(CMalpha*alpha + (d*CMq*q)/(2*v));
    g_delta = (Pdyn*S*d*CMdelta)/J;

    delta_eq = -f_delta/g_delta;
    delta_cmd = (Kp_rot*(teta_des - teta) + Kd_rot*(-q))/g_delta;
    delta = delta_eq + delta_cmd;
    
end

f(1) = -Daero/m - gr*sin(x(2));

Laero = Pdyn*S*CLalpha*alpha;
f(2) = (1/v)*(Laero/m + (v^2/r - gr)*cos(gamma));

f(3) = v*sin(gamma);

f(4) = v/r*cos(gamma);

f(5) = q;

Maero = Pdyn * S * d * (CMalpha * alpha + (d*CMq*q)/(2*v) + CMdelta*delta);
f(6) = (1/J) * Maero;

f(7) = Daero > 2650;


f = f(:);