%Worksheet 3
%Abraham, Lukas, Ramesh
%Group 9
clc 
clear variables
format short

plotGraphs = true;
printResults = true;

%initiliaze t
t0 = 0;
tmax = 5;
exp2 = 5;
for i = 0:exp2
    tdelta(i+1) = 1/(2^i);
end

deltas = numel(tdelta);

%Derivative:
p1 = @(p) ( 7*( 1 - (p / 10) ) * p) ;

%initial condition
p0 = 20;

for i = 1:deltas
    dt = tdelta(i);
    t = [t0 : dt : tmax];
    size1  = numel(t);
    euler = Fun_ExpEuler(p1, dt, p0, size1);  
    impeuler = Fun_ImpEuler1(p1, dt, p0, size1);
    adams = Fun_AdamsMoulton(p1, dt, p0, size1);
    heun = Fun_Heun(p1, dt, p0, size1);
    linAdam1 = Fun_Linearization1(p0, dt, tmax);
    linAdam2 = Fun_Linearization2(p0, dt, tmax);
    euler1{i} = euler; 
    heun1{i} = heun; 
    impeuler1{i} = impeuler;
    adams1{i} = adams;
    lin1{i} = linAdam1;
    lin2{i} = linAdam2;
    pfunc{i} = (200 ./ ( 20 -  (10 * exp(-7*t)) ));
    tarray{i} = t;
end

% Error calculation

%allmeths(methods, timesteps)
allmeths = [euler1; heun1; impeuler1; adams1; lin1; lin2; pfunc];

%preallocation:
err = zeros(6,deltas-1);
red_error = zeros(6,deltas-2);

for i = 1:6
    for j = 2:deltas
        %Error compared to Exact Solution:     
        err(i,j-1) = Fun_Error(allmeths(i,j), tdelta(j), allmeths(7,j), tdelta(j), tmax);
    end
end

for i = 1:6
    for j = 1:(deltas-2)
        %Error reduced factor 
        red_error(i,j) = err(i,j) / err(i,j+1);
    end
end

if plotGraphs == true
    figure(1)
    plot(tarray{6},pfunc{6});
    title('Exact Solution')
    legend('Exact Solution')
    xlim([0,1])
    ylim([0,20])

    figure(2)
    plot(tarray{4},pfunc{4}, tarray{1},euler1{1},tarray{2},euler1{2},tarray{3},euler1{3},tarray{4},euler1{4},tarray{5},euler1{5},tarray{6},euler1{6});
    title(' Explicit Euler Method')
    legend('Exact Solution', 't = 1', 't=0.5', 't=0.25','t=0.125','t=0.0625', 't=0.0313')
    axis([0 5 0 20]);

    figure(3)
    plot(tarray{4},pfunc{4}, tarray{1},heun1{1},tarray{2},heun1{2},tarray{3},heun1{3},tarray{4},heun1{4},tarray{5},heun1{5},tarray{6},heun1{6});
    title('Heun Method')
    legend('Exact Solution', 't = 1', 't=0.5', 't=0.25','t=0.125','t = 0.0625','t=0.0313')
    axis([0 5 0 20]);

      figure(4)
     plot_graph(tarray,impeuler1,pfunc,'Implicit Euler scheme');

      figure(5)
     plot(tarray{4},pfunc{4},tarray{2},adams1{2},'*',tarray{3},adams1{3},tarray{4},adams1{4},tarray{5},adams1{5},tarray{6},adams1{6});
      legend('Exact Solution','t=0.5', 't=0.25','t=0.125','t=0.0625','t=0.0313');
     title('Adams  Method');

      figure(6)
     plot_graph(tarray,lin1,pfunc, 'Adams Moulton with Linearisation 1');

      figure(7)
     plot_graph(tarray,lin2,pfunc, 'Adams Moulton with Linearisation 2');
end
 
if printResults == true
     % Displaying results
    format compact
    % Explicit euler
    fprintf('Calculated errors of all methods: \n \n');
    fprintf('\n Explicit Euler: ');
    disp(err(1,:));
    fprintf('Reduction factor: \t   ')
    disp(red_error(1,:));
    % Heun
    fprintf('\n Method of Heun: ');
    disp(err(2,:));
    fprintf('Reduction factor: \t   ')
    disp(red_error(2,:));
    % Implicit Euler
    fprintf('\n Implicit Euler: ');
    disp(err(3,:));
    fprintf('Reduction factor: \t   ')
    disp(red_error(3,:));
    % Adams-Moulton
    fprintf('\n Adams-Moulton:  ');
    disp(err(4,:));
    fprintf('Reduction factor: \t   ')
    disp(red_error(4,:));
    % Linearization 1
    fprintf('\n Linearization 1:');
    disp(err(5,:));
    fprintf('Reduction factor: \t   ')
    disp(red_error(5,:));
    % Linearization 2
    fprintf('\n Linearization 2:');
    disp(err(6,:));
    fprintf('Reduction factor: \t   ')
    disp(red_error(6,:));
end
