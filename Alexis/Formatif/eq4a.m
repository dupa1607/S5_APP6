function f = eq4a(t,x)
    Kp = 5; Kd = 8; Ka = 6;
    cmd = 1;
    xdes = 5;

    f(1) = x(2);
    f(2) = x(3);
    f(3) = -1.75*x(3)-1.625*x(2)*abs(x(2))-1.25*x(1)^3;
    if cmd == 1
        f(3) = Kp*(xdes-x(1)) + Kd*(-x(2)) + Ka*(-x(3)); 
    end
    f=f(:);
end