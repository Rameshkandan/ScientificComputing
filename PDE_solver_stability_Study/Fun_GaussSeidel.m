function [T] = Fun_GaussSeidel (Nx, Ny, dt, T1)
err = 1e-6;
% Precalculation of 1 over N for residual
N1 = 1 / (Nx*Ny);
% Precalculation of h values
hx = 1 / (Nx + 1);
%hy = 1 / (Ny + 1);
% initialisation of T
T = [zeros(1,Nx+2); zeros(Nx,1) zeros(Nx) zeros(Nx,1) ; zeros(1,Nx+2)];
% initialisation of residual to make sure it is bigger than tolerance at start
R = 1111;
count = 0;
%division operator
div = ((4*dt)/(hx*hx)) + 1;
% Loop over T vector until residual is smaller than given tolerance.
while R > err
   for i = 2:Nx+1
    for j = 2:Ny+1
        T(i,j) = ( T1(i,j) + ( (dt/(hx*hx)) * ( T1(i-1,j) + T1(i+1,j) + T1(i,j-1) +T1(i,j+1) ) ) ) / div;
    end 
   end
    %R = sqrt(N1 * sum((T_old - T).^2));
    %R = sqrt(N1 * ((sum(sum(T - T1)))^2) );
    diffMat = T - T1;
    R = sqrt(N1 * sum(sum((diffMat .* diffMat))) );
    T1 = T;
    count = count +1;
    if count > 5000
        T = 0;
        disp('did not converge')
        break;
    end       
end
