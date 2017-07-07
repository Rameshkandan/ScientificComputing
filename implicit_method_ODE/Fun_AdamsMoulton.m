function [adams] = Fun_AdamsMoulton(f, dt, f0, size1)    
    adams = zeros(1,size1);
    adams(1) = f0;
    err = 10e-4;
     for i = 2:size1
             x0 = adams(i-1);
             g = @(x) x - (x0 + ( (dt/2)*(f(x0) +f(x) ) )  ); 
              gprime =@(x) 1 - dt/2*(((f(x0+eps)-f(x0))/eps) + ((f(x+eps)-f(x))/eps)); %derivative
             adams(i) = Fun_Newton(g,gprime,adams(i-1),err); 
     end
     
     
end