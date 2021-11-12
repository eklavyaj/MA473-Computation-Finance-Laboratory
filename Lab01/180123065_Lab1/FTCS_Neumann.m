function U = FTCS_Neumann(phi, f, g1, g2, M, lambda)

    N = round(M*M/lambda);
    k = 1/N;
    h = 1/M;
    
    U = zeros(M+1, N+1);
    
    Tmax = N*k;
    X = 0:h:M*h;
    T = 0:k:Tmax;
    
    U(1:M+1,1) = f(X);
    U(1, 1:N+1) = g1(T);
    U(M+1, 1:N+1) = g2(T);
        
    for i=2:N+1
        tn = T(i);
        U(1, i) = (1 - 2*lambda)*U(1, i-1) + 2*lambda*U(2, i-1);
        for j=2:M
            xm = X(j);
            U(j, i) = lambda*(U(j-1, i-1) + U(j+1, i-1)) + (1-2*lambda)*(U(j, i-1)) + k*phi(xm, tn);
        end
    end
end