function plot_gph(T, Nx, Ny, form, type)
    
if strcmp(form, 'vector') == true
    Tmesh = reshape(T,[Nx,Ny]);
    Tmesh1 = [zeros(1,(Nx)+2); zeros((Nx),1), Tmesh, zeros((Nx),1); zeros(1,(Nx)+2)];
elseif strcmp(form, 'matrix') == true
    Tmesh1 = T;
end
[X,Y]=meshgrid(linspace(0,1,Nx+2));
if strcmp(type, 'surface') == true
    surf(X,Y,Tmesh1);
elseif strcmp(type, 'contour') == true
    contour(Tmesh1);
end
