function error_table(lambda, method, phi, f, g1, g2, filename)
    
    M = 5;    
    M_arr = [M];
    err_arr = [];
    
    U = method(phi, f, g1, g2, M, lambda);
    
    
    while M <= 100
        M = M*2;
        M_arr = [M_arr, M];
        
        U2 = method(phi, f, g1, g2, M, lambda);
        temp = U2;
        U2 = U2(1:2:end, 1:4:end);
        
        err_arr = [err_arr, norm(U - U2, Inf)];
        U = temp;
    end
    
    err_arr = [err_arr, nan];
    order_conv = [];
    for i=1:length(err_arr)-1
        if err_arr(i) ~= 0
            order_conv = [order_conv, log2(err_arr(i)/err_arr(i+1))];
        else
            order_conv = [order_conv, nan];
        end
    end
    
    order_conv = [order_conv, nan];
    T = table(M_arr', err_arr', order_conv', 'VariableNames',{'M', 'Error', 'Order of Convergence'});
    writetable(T, ['error_tables/' filename '.txt']);
    
end



