function U = BTCS_Neumann(phi, f, g1, g2, M, lambda)

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
    
    for i= 2:N+1
                
        A = diag((1 + 2* lambda)*ones(1,M+1)) + diag(-lambda*ones(1,M),1) + diag(-lambda*ones(1,M),-1);
        A(1, 1) = 1;
        A(1, 2) = 0;
        A(M+1, M+1) = 1;
        A(M+1, M) = 0;
        
        b = zeros(M+1, 1);
        b(1) = U(1, i-1) + 2*lambda*(U(2, i-1) - U(1, i-1));
        
        b(2:M) = U(2:M, i-1) + (k*phi(X(2:M), T(i)))';
        b(M+1) = U(M+1, i);
        U(:, i) = A\b;
        
    end
end