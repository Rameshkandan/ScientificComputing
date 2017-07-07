function [T] = Fun_GaussSeidel0 (b, Nx, Ny, N1)
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
% initialisation of convergence factor to make sure cnvrge is bigger than tolerance
convrge = 1111;
count = 0;
% Loop over T vector until residual is smaller than given tolerance.
while convrge > err && count < 100 
    count = count + 1;
    % Update Temperature of each cell (except for boundary cells).
    for i = (Nx+4) : ((Nx+2)*(Ny+2)-(Nx+2))
        % Check for boundary cells inside T vector
        if mod(i,(Nx+2)) ~= 0 && mod(i,(Nx+2)) ~= 1
            % Compute temperature in current cell
            T(i) = ( T(i-1) + T(i+1) + T(i+(Nx+2)) + T(i-(Nx+2)) - b(i)*(hx^2))*0.25 ;
            
            % Multiply T value with corresponding values of A for residual
            %-------------------------------------------------------------
            % T corresponding to first columns of blocks of A
            if i == (Nx+4) + (floor(i/Nx)-2) * (Nx + 2)
                % Very first column
                if i = Nx+4
                    colSum = T(i) * (a11 + a21 + a12);
                % First column of last block
                elseif i = (Nx+2)*(Ny+2)-(2*Nx+2)
                    colSum = T(i) * (a21 + a11 + a21);
                % First column of every intermediate block
                else
                    colSum = T(i) * (a21 + a11 + a21 + a12);
                end
            % T corresponding to last columns of blocks of A
            elseif i == (2*Nx+3) + (floor(i/Nx)-2) * (Nx + 2)
                % Last column of first block
                if i = (2*Nx + 3)
                    colSum = T(i) * (a12 + a11 + a12);
                % Very last column
                elseif i = (Nx+2)*(Nx+2)
                    colSum = T(i) * (a21 + a12 + a11);
                % Last column of intermediate blocks
                else
                    colSum = T(i) * (a21 + a12 + a11 + a12);
                end
            % T corresponding to columns in between
            else
                % Intermediate columns of first block
                if i < (2*Nx + 3)
                    colSum = T(i) * (a12 + a11 + a21 + a12);
                % Intermediate columns of last block
                elseif i > ( (Nx+2)*(Ny+2) - (2*Nx+2))
                    colSum = T(i) * (a21 + a12 + a11 + a21);
                else
                    colSum = T(i) * (a21 + a12 + a11 + a21 + a12);
                end
            end

            % Add sum of current column to total sum
            totSum = totSum + colSum

        end
    end
    totSum = totSum - sum(b)
    T0 = T(T~=0);
    b0 = b(b~=0);
    if mod(count,Nx) == 0
        convrge = Fun_ResidualNorm (T0, b0, Nx, Ny, 1/(Nx*Ny));
    end
end
