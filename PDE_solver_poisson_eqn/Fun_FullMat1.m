function [A] = Fun_FullMat1(Nx,Ny)

hx = (1/(Nx+1));
hy = (1/(Ny+1));
m = (Nx * Ny);

diag0 = diag ( (ones(Nx,1) * ((-2/(hx^2)) + (-2/(hy^2)) ) )  );
diag1 = diag(ones(Nx-1,1)*(1/(hx^2)),1);
diagn1 = diag(ones(Nx-1,1)*(1/(hy^2)),-1);
B = diag0 + diag1 +diagn1;
K = eye(Nx);
K1 = kron(K,B);
diagIx = diag(ones(m-Nx,1)*(1/(hx^2)),Nx);
diagIy = diag(ones(m-Nx,1)*(1/(hy^2)),-1*Nx);
A = K1 + diagIx + diagIy;


