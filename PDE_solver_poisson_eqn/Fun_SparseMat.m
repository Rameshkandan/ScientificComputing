function sparse = Fun_SparseMat(Nx,Ny)
k1 = (-2*((Nx+1)^2)) + (-2*((Ny+1)^2));
k2 = (Nx+1)^2;
k3 = (Ny+1)^2;
idy =speye(Nx);
Super_diag=ones(Nx,1);                       
tri_diagl=spdiags([k1*ones(Nx,1), k2*Super_diag, k3*Super_diag],[0,1,-1],Nx,Nx);   
idy_diag=spdiags([ k3*Super_diag, k2*Super_diag ],[-1,1],Nx,Nx);                   
t1 = kron(idy,tri_diagl);t2 = kron(idy_diag,idy);
sparse=t1+t2; 
