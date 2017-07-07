function R = Fun_ResidualNorm1 (T, b, Nx, Ny, N)

% initializing parameters
sum_k = 0;
    
h = 1/(Nx+1);
% looping over the 'virtual' rows of A
for k = 1 : Nx*Ny
sum_m = 0;

    % looping over the 'virtual' columns of A
    for m = 1 : Nx*Ny
        % Select a according to position in A matrix
        if m == k %( (Nx+4)+(k-1)*(Nx+2)+(k-1))
            a = -4;
        elseif m == k+1 && mod(k,Nx)~= 0%(Nx+4)+(k-1)*(Nx+3) - 1 && k ~= Nx+1
            a = 1;
        elseif m == k-1 && mod(k,(Nx))~=1 %(Nx+4)+(k-1)*(Nx+3) + 1 && k ~= Nx+2
            a = 1;
        elseif m == Nx+k
            a = 1;
        elseif k == Nx+m
            a = 1;
        else
            a = 0;
        end   

        % Calculate sum over b (complete row of A
        adding = T(m) * (a/(h^2));
        sum_m = sum_m + adding;
    end
    
    %diff(k,1) = sum_m;
    %diff(k,2) = b(k);
    %diff(k,3) = b(k) - sum_m;
    %disp(diff(1,:));

    % Calculate sum over k (column wise of A)
    sum_k = sum_k + (b(k) - sum_m)^2;
end


% Calculate residual of iteration loop
R = sqrt(N * sum_k);
