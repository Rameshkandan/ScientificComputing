function y = Fun_Linearization2(y0, dt, tend)
t = 0:dt:tend;
size1 = numel(t);
y = zeros(1,size1);
y(1) = y0;
 for i = 2:size1
    X = (7 - (7*y(i-1)/10)) * (dt/2);
    y(i) = ( y(i-1) + (X*y(i-1)) ) / (1-X);
 end
    
    