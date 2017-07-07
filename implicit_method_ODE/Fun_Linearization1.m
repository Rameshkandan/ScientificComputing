function y = Fun_Linearization1(y0, dt, tend)

t = 0:dt:tend;
size1 = numel(t);
y = zeros(1,size1);
y(1) = y0;
 for i = 2:size1
     y(i) = ( y(i-1) + ( 7 * (dt/2) * ( 2 - (y(i-1)/10) ) * y(i-1) ) ) / (1 + (7*y(i-1)/10)*(dt/2)) ; 
 end
    
    