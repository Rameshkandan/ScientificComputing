function x1 = Fun_Newton(g,gprime,x0,err)   
cnt = 1;
x1 = x0 + 1;
    while abs(x1-x0) >= err
             x0 = x1;
             x1 = x0 - g(x0)/gprime(x0);
             cnt = cnt +1;
             if cnt >100
                 break;
             end  
    end 
         
end
