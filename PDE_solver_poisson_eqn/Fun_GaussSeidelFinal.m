function [T_GS0] = Fun_GaussSeidelFinal (b, Nx, Ny, N1)
err = 1e-4;
% Precalculation of h values
hx = 1 / (Nx + 1);
hy = 1 / (Ny + 1);
% Precalculation of values of A matrix
a11 = -( 2 / (hx^2) ) - (2 / (hy^2));
a12 = 1 / (hx^2);
a21 = 1 / (hy^2);
% Precalculation of 1 over N for residual
N1 = 1 / (Nx*Ny);
% initialisation of T
T = zeros((Nx+2)*(Ny+2),1);
% initialisation of residual to make sure it is bigger than tolerance at start
R = 1111;
count = 0;
% Loop over T vector until residual is smaller than given tolerance.
while R > err 
    ax = zeros((Nx+2)*(Ny+2),1);
    count = count + 1;
    % Update Temperature of each cell (except for boundary cells).
    for i = (Nx+4) : ((Nx+2)*(Ny+2)-(Nx+2))
        % Check for boundary cells inside T vector
        if mod(i,(Nx+2)) ~= 0 && mod(i,(Nx+2)) ~= 1
            % Compute temperature in current cell
            T(i) = ( T(i-1) + T(i+1) + T(i+(Nx+2)) + T(i-(Nx+2)) - b(i)*(hx^2))*0.25 ;
            
            % Add contributions of current T to ax vector
            ax(i) = ax(i) + a11 * T(i);
            % Contribution to column above (Exclude first columns of all blocks)
            if i ~= (Nx+4) + (floor(i/(Nx+2))-1) * (Nx + 2)
                ax(i-1) = ax(i-1) + a12 * T(i);
            end
            % Contribution to column below (Exclude last column of all blocks)
            if i ~= (2*Nx+3) + (floor(i/(Nx+2))-1) * (Nx + 2)
                ax(i+1) = ax(i+1) + a21 * T(i);
            end
            % Contributen to unit matrix above (Exclude first block)
            if i > (2*Nx + 3)
                ax(i-(Nx+2)) = ax(i-(Nx+2)) + a21 * T(i);
            end
            % Contribution to unit matrix below (Exclude last block)
            if i < ( (Nx+2)*(Ny+2) - (2*Nx+2))
                ax(i+(Nx+2)) = ax(i+(Nx+2)) + a12 * T(i);
            end

        end
    end
    R = sqrt(N1 * sum((b - ax).^2));
end
T_GS0 = reshape(T,[Nx+2 Ny+2]);
