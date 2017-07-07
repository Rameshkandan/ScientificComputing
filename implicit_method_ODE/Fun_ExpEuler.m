function [euler] = Fun_ExpEuler(f, dt, f0, size1)
    euler = zeros(1,size1);
    euler(1) = f0;
     for i = 2:size1
            euler(i) = euler(i-1) + dt*f(euler(i-1));
     end
end


