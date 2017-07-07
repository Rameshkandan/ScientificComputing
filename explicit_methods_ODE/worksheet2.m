%Worksheet 2
%Abraham, Lukas, Ramesh
%Group 9
clc 
clear variables
format long

%initiliaze t
t0 = 0;
tmax = 5;
tdelta = [1, 1/2 , 1/4 , 1/8];
deltas = numel(tdelta);
    
%Exact Solution
p = @(t) 10 / ( 1 +  9 * exp(-t) );
%Derivative
p1 = @(p) ( 1 - (p / 10) ) * p;
%initial condition
p0 = p(t0);

for i = 1:deltas
    dt = tdelta(i);
    t = [t0 : dt : tmax];
    size1  = numel(t);
    [euler,heun,ruku] = funmethods(p1, dt, p0, size1);   
    euler1{i} = euler; heun1{i} = heun; ruku1{i} = ruku;
    pfunc{i} = 10 ./ ( 1 +  9 .* exp(-t) );
    tarray{i} = t;
end
%allmeths(methods, timesteps)
allmeths = [euler1; heun1; ruku1; pfunc];

%preallocation:
exacte = zeros(3,deltas);
beste = zeros(3,deltas);
red_error = zeros(3,3);

%Error Calculation
for i = 1:3
    for j = 1:deltas
        %Error compared to Exact Solution:
        exacte(i,j) = compError(allmeths(i,j), tdelta(j), allmeths(4,j), tdelta(j), tmax);
        
        %Error compared to smallest timestep:
        beste(i,j) = compError(allmeths(i,j), tdelta(j), allmeths(i,deltas), tdelta(deltas), tmax);
    end
end

for i = 1:3
    for j = 1:(deltas-1)
        %Error reduced factor 
        red_error(i,j) = exacte(i,j) / exacte(i,j+1);
    end
end

%display:
exacte
beste
red_error


%PLOTS
%Exact solution Graph
figure(1)
plot(tarray{4},pfunc{4});
title('Exact Solution')
legend('Exact Solution')

%Explicit Euler Graph
figure(2)
plot(tarray{4},pfunc{4}, tarray{1},euler1{1},tarray{2},euler1{2},tarray{3},euler1{3},tarray{4},euler1{4});
title('Euler Method')
legend('Exact Solution', 't=1', 't=0.5', 't=0.25','t=0.125')

%Heun Method Graph
figure(3)
plot(tarray{4},pfunc{4}, tarray{1},heun1{1},tarray{2},heun1{2},tarray{3},heun1{3},tarray{4},heun1{4});
title('Heun Method')
legend('Exact Solution', 't=1', 't=0.5', 't=0.25','t=0.125')

%Runge Kutta Method Graph
figure(4)
plot(tarray{4},pfunc{4}, tarray{1},ruku1{1},tarray{2},ruku1{2},tarray{3},ruku1{3},tarray{4},ruku1{4});
title('Runge Kutta Method')
legend('Exact Solution', 't=1', 't=0.5', 't=0.25','t=0.125')


