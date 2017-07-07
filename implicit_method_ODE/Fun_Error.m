function E = Fun_Error(A1, dt1, A2, dt2, tmax)
a1 = cell2mat(A1);
a2 = cell2mat(A2);
    val_sum = 0;
    for k = 1 : (numel(a1))
        val_sum = val_sum + (a1(k) - a2((dt1/dt2) * (k-1) + 1))^2;
    end
    E = sqrt(dt1/tmax * val_sum);
end