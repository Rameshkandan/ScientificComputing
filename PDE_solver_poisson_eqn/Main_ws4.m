clear variables
clc
%Main File

plotSurface = true % Enable/disable surface plots
plotContour = true % Enable/disable contour plots

i=1;k=1;
e = 1;
err_rdn = zeros(e-1,1);
for Nx = [7 15 31 63 127] 
    Ny=Nx;
    disp(sprintf('\n Nx and Ny are: %.0f %.0f', Nx, Ny));
    disp(sprintf('--------------------------------------- \n'));
    hx = 1 / (Nx + 1);
    hy = 1 / (Ny + 1);
    
    
    %initialize
    f = zeros((Nx*Ny),1);
    f = @(x,y) -2.*(pi^2).*sin(pi*x).*sin(pi*y) ;
    x1 =(1:Nx)./(Nx+1); y1 =(1:Ny)./(Ny+1);
    [X,Y]=meshgrid(x1,y1);                      
    RHS_function=f(X,Y);
    RHS_Mat = [zeros(1,Nx+2); zeros(Nx,1) RHS_function zeros(Nx,1) ; zeros(1,Nx+2) ];
    RHS_Mat = RHS_Mat(:);
    %make use of RHS_Mat : this is the  column vector with the zeros incorporated in the vector 
    RHS_function= RHS_function(:);
    % figure(i);
    % plot_gph(T,Nx,Ny);
    % i =i+1;
    %disp('sparse');
    
    % Calculate results
    % ---------------------------------------------------------------------------------------
    
    % Analytical solution
    for i = 0 : Nx + 1
        for j = 0 : Ny + 1
            T_exact(i+1,j+1) = sin(pi*i*hx)*sin(pi*j*hy);
        end
    end
    
    if Nx<=63
        % Calculation with full A
        disp('Solving with full matrix');
        tic
        A = Fun_FullMat1(Nx,Ny);
        T_full = A\RHS_function;
        toc
        % storage requirement: full matrix A + vector b
        disp(sprintf('Storage requirement: %.0f \n', ((Nx*Ny)^2 + Nx*Ny)))
        
        % sparse A matrix
        disp('Solving with sparse matrix');
        tic
        sparse = Fun_SparseMat(Nx,Ny);
        T_sparse = sparse\RHS_function;
        toc
        nnz = Nx*(Nx + 2*(Nx-1)) + 2*(Nx-1)*Nx; % non-zero entries of A
        A_entires = nnz + nnz + (Nx+1); % values of nnz, row indeces of nnz, cumulative num of entries in column
        b_entries = Nx*Ny; % elements in vector b
        % storage requirement: non-zero entries of matrix A times + vector b
        disp(sprintf('Storage requirement: %.0f \n', (A_entires + b_entries)))
    end
    
    % Gauss-Seidel solver
    disp('Solving with Gauss-Seidel algorithm');
    tic
    T_GS = Fun_GaussSeidelFinal(RHS_Mat, Nx, Ny);
    toc
    % storage requirement: vector b + vector T + vector ax
    disp(sprintf('Storage requirement: %.0f \n', (3*(Nx+2)*(Ny+2))));
    
    % Error estimation for Gauss-Seidel
    err(e) = (1/(Nx*Ny) * sum(sum((T_exact-T_GS).^2,2)))^0.5;
    e = e + 1;
    if (e>=2)
        err_rdn(2:e-1) = err(1:e-2)./err(2:e-1);
    end
    
    if Nx <= 63
        if plotSurface == true
            figure(k);
            plot_gph(T_full,Nx,Ny,'vector', 'surface');
            figure(k+1)
            plot_gph(T_sparse,Nx,Ny, 'vector', 'surface');
            figure(k+2)
            plot_gph(T_GS, Nx, Ny, 'matrix', 'surface');
            k = k+3;
        end
        if plotContour == true
            figure(k);
            plot_gph(T_full,Nx,Ny,'vector', 'contour');
            figure(k+1);
            plot_gph(T_sparse,Nx,Ny, 'vector', 'contour');
            figure(k+2);
            plot_gph(T_GS, Nx, Ny, 'matrix', 'contour');
            k = k+3;
        end
    end
end

Nx = 127;
Ny = Nx;
disp(sprintf('Nx and Ny are: %.0f %.0f \n', Nx, Ny));
% Gauss-Seidel solver
disp('Solving with Gauss-Seidel algorithm');
tic
T_GS = Fun_GaussSeidelFinal(RHS_Mat, Nx, Ny);
toc
% storage requirement: vector b + vector T + vector ax
disp(sprintf('Storage requirement: %.0f \n', (3*(Nx+2)*(Ny+2))));
% Error estimation for Gauss-Seidel
err(e) = (1/(Nx*Ny) * sum(sum((T_exact-T_GS).^2,2)))^0.5;
err_rdn(e) = err(e-1)./err(e);






