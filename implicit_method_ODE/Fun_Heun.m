function [heun] = Fun_Heun(f, dt, f0, size1)
    heun = zeros(1,size1);
    heun(1) = f0;
     for i = 2:size1
         x = Fun_ExpEuler(f, dt, heun(i-1),2);
         heun(i) = heun(i-1) + ( dt *  ( 0.5 * ( f(heun(i-1)) + f(x(2))) ));
     end
end


