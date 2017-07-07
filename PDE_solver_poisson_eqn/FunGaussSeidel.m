function [Tmesh] = Fun_GaussSeidel (b, Nx, Ny)
err = 1e-4;
T = zeros(Nx*Ny,1);
ij = (-2*((Nx+1)^2)) + (-2*((Ny+1)^2));
ij1 = (Nx+1)^2;
ij2 = (Ny+1)^2;
sum1=zeros(Nx);
convrge = 2;
count = 0;
%while convrge > err
    while convrge > err && count < 100
        count = count +1
        for i = 1:(Nx*Ny)
            if (mod(i,Nx) == 0) && (i > Nx) && (i <((Nx*Ny)-Nx))
                T(i) = (b(i) - T(i-1) - T(i+Nx) - T(i-Nx)) * (-0.25);
            elseif (mod(i,Nx) == 1) && (i > Nx) && (i <((Nx*Ny)-Nx))
                T(i) = (b(i) - T(i+1) - T(i+Nx) - T(i-Nx)) * (-0.25);
            elseif(i == Nx)           
                T(i) = (b(i) - T(i-1) - T(i+Nx) ) * (-0.25);
            elseif (i == 1)
                T(i) = (b(i) - T(i+1) - T(i+Nx) ) * (-0.25);
            elseif (i>1) && (i<Nx)
                T(i) = (b(i) - T(i+1) - T(i-1) - T(i+Nx)) *(-0.25);
            elseif i>((Nx*Ny)-Nx) && i<(Nx*Ny)
                T(i) = (b(i) - T(i+1) - T(i-1) - T(i-Nx)) *(-0.25);
            elseif i==((Nx*Ny)-Nx)
                T(i) = (b(i) - T(i+1) - T(i-Nx)) *(-0.25);
            elseif i==(Nx*Ny)
                T(i) = (b(i) - T(i-1) - T(i-Nx)) *(-0.25);
            else
                T(i) = (b(i) - T(i+1) - T(i-1) - T(i-Nx) - T(i+Nx)) *(-0.25)           
            end
        end
        
    %    convrge = Fun_ResidualNorm1 (T, b, Nx, Ny, 1/(Nx*Ny));

    end
    T1 = reshape(T,[Nx,Ny]);
    Tmesh = [zeros(1,(Nx)+2); zeros((Ny),1), T1, zeros((Ny),1); zeros(1,(Nx)+2)];
    
    %end

    
% while convrge > err
%     for i = 2:(Nx+1)
%                     
% 
%         for j = 2:(Ny+1)
%             T(i,j) = ( b(i) - T(i-1,j) - T(i+1,j)-T(i,j-1)-T(i,j+1) )* (-0.25);
%             if i==j
%                 a = ij;
%             elseif i == Nx + j
%                 a = ij2;
%             elseif j == Ny + i
%                 a = ij1;
%             elseif i == j+1
%                 a = ij1;
%             elseif j == i+1
%                 a = ij2;
%             else
%                 a = 0;
%             end   
%             sum1(i-1) = sum1(i-1) + T(i,j) * a;    
%             
%         end
%         sum2 = sum2 + b(i)- sum1
%         
%     end
%     
%     
%     
% end
