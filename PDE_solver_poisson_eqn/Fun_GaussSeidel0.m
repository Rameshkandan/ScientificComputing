function [T] = Fun_GaussSeidel0 (b, Nx, Ny, N1)
err = 1e-4;
hx = 1 / (Nx + 1);
N1 = 1 / (Nx*Ny);
T = zeros((Nx+2)*(Ny+2),1);
convrge = 1111;
count = 0;
while convrge > err && count < 100 
    count = count + 1;
    for i = (Nx+4) : ((Nx+2)*(Ny+2)-(Nx+2))
        if mod(i,(Nx+2)) ~= 0 && mod(i,(Nx+2)) ~= 1
            % Compute temperature in current cell
            T(i) = ( T(i-1) + T(i+1) + T(i+(Nx+2)) + T(i-(Nx+2)) - b(i)*(hx^2))*0.25 ;
        end
        
    end
    T0 = T(T~=0);
    b0 = b(b~=0);
    convrge = Fun_ResidualNorm (T0, b0, Nx, Ny, 1/(Nx*Ny));
end
