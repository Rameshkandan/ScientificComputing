function [impeuler] = Fun_ImpEuler1(f, dt, f0, size1)

    err = 10e-4;
    impeuler = zeros(1,size1);
    impeuler(1) = f0;
     for i = 2:size1 
         g = @(x) x - (impeuler(i-1) + dt * f(x));
          gprime = @(x) 1 - 7*dt*(1-2*(x/10)); %derivative
         impeuler(i) = Fun_Newton(g,gprime,impeuler(i-1),err);
     end  
end