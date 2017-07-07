clear variables
clc
% Create array with initial conditions
count = 0;
count_nx = 0;
count_imp = 0;
plots = 1;
intStability = zeros(4,7);
stability = string(intStability);

tic
for Nx = [3 7 15 31]
    %Initial guess for T_implicit:
    T_imp=[zeros(1,Nx+2); zeros(Nx,1) ones(Nx) zeros(Nx,1) ; zeros(1,Nx+2)];
    count_nx = count_nx + 1;
    [X,Y]=meshgrid(linspace(0,1,Nx+2));
    count_tau = 0;
    for tau = [1/64 1/128 1/256 1/512 1/1024 1/2048 1/4096]
        Ny = Nx;
        count_tau = count_tau + 1;
        count = count + 1;
        % Calculate hx and hy
        hx = 1 / (Nx + 1);
        hy = 1 / (Ny + 1);
    
        % Initialize T Array with initial conditions
        T_exp = ones(Nx*Ny,1);
    
        % Create A matrix
        A = Fun_SparseMat(Nx, Ny);

        % Loop from t = 0 to max time
        for t = 0 : tau : 4/8
           % Exlplicit Euler solution
           T_exp = Fun_ExpEuler(T_exp, A, tau);
           if plots == 1
               if mod(t, (1/8)) == 0 && t ~= 0             
                   %Save Explicit plot to folder
                   %saveFig = figure('visible', 'off');
                   figure(9)
                   Fun_plot(T_exp, Nx, Ny, 'vector', 'surface');
                   view(3)
                   str = sprintf('ExplicitEuler_t=%.5f_Nx=Ny=%.5f_dt=%.5f', t, Nx, tau);
                   str2 = strrep(str,'.','_');
                   title(['Explicit ','Nx=Ny=',num2str(Nx),',tau = ', num2str(tau), ', t= ', num2str(t)])
                   saveas(figure(9),str2,'jpeg')
                   %saveas(saveFig,str2,'jpeg');
                   %close(saveFig);
                   
                   
                   display(t*8);
                   figure(t*8);
                   %hold on
                   subplot(4,7,count);
                   Fun_plot(T_exp, Nx, Ny, 'vector', 'surface');
                   if Nx == 3
                        title(['\tau = ', num2str(tau)])
                   end
                   if tau == 1/64
                       zlabel(['Nx = ', num2str(Nx)]);
                   end
                   %title(['Ex','t= ', num2str(t),',Nx=Ny=',num2str(Nx),',tau =', num2str(tau)])
                   %hold off
               end
           end

           % Implicit Euler solution
            if tau == 1/64
                T_imp = Fun_GaussSeidel(Nx, Ny, tau, T_imp);
                if plots == 1
                    if mod(t, (1/8)) == 0 && t ~= 0 && tau == (1/64)
                        %Save Implicit plot to folder:
                        saveFig = figure('visible', 'off');
                        Fun_plot(T_exp, Nx, Ny, 'vector', 'surface');
                        view(3)
                        str = sprintf('ImplicitEuler_t=%.5f_Nx=Ny=%.5f_dt=%.5f', t, Nx, tau);
                        str2 = strrep(str,'.','_');
                        title(['Implicit ', 'Nx=Ny=',num2str(Nx),',tau = ', num2str(tau), ', t= ', num2str(t)])
                        saveas(saveFig,str2,'jpeg')
                        close(saveFig);
                        
                        %Implicit Subplot, dt=1/64:
                        count_imp = count_imp+1;
                        figure(5);
                        subplot(4,4,count_imp);
                        surf(X,Y,T_imp);             
                        if Nx == 3
                            title(['t = ', num2str(t)])
                        end
                        if t == 1/8
                            zlabel(['Nx = ', num2str(Nx)]);
                        end
                        %title(['Imp ', 'Nx=',num2str(Nx),',t=', num2str(t)])
                    end
                end
            end
               
        end

        % Check for stability via von Neumann stability
        intStability(count_nx, count_tau) = abs(((hx*hy)/2) - tau);
        if tau < ((hx*hy) / 2)
            stability(count_nx, count_tau) = 'stable';
        else
            stability(count_nx, count_tau) = 'unstable';
        end
    end    
end
toc

% Print Solutions
%-----------------------------
display('From equations (1), (2), (3) it can be seen, that if the time tends to infinity the Tempreture decreases to 0 everywhere.');

stability = reshape(stability,[4,7]);
stabilities = array2table(stability,'RowNames', {'3', '7', '15', '31'});
stabilities.Properties.VariableNames = {'dt_1_64', 'dt_1_128', 'dt_1_256', 'dt_1_512', 'dt_1_1024', 'dt_1_2048', 'dt_1_4096'};
display(stabilities);

display('Looking at the solutions shown in the plot the obtained tempreture are indeed decreasing.');
