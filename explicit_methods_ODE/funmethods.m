function [euler , heun , ruku] = funmethods(f, dt, f0, size1)
    euler = zeros(1,size1);
    heun = zeros(1,size1);
    ruku = zeros(1,size1);
    euler(1) = f0;
    heun(1) = f0;
    ruku(1) = f0;
     for i = 2:size1
            euler(i) = expEuler(f, dt, euler(i-1));
            heun(i) = expHeun(f, dt, heun(i-1));
            ruku(i) = expRuKu(f, dt, ruku(i-1));
     end
end

function f1 = expEuler(f, dt, f0)
    f1 = f0 + dt * f(f0);
end

function f1 = expHeun(f, dt, f0)
    x = expEuler(f, dt, f0);
    f1 = f0 + dt *  ( 0.5 * ( f(f0) + f(x)) );
end

function f1 = expRuKu(f, dt, f0)
    y  = f0;
    Y1 = f(y);
    Y2 = f(y + dt/2 * Y1);
    Y3 = f(y + dt/2 * Y2);
    Y4 = f(y + dt * Y3);
    f1 = y + dt*(1/6)*(Y1+2*Y2+2*Y3+Y4);
end

