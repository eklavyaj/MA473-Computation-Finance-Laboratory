function U = CN_Neumann(phi, f, g1, g2, M, lambda)

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
                
        A = diag((1 + lambda)*ones(1,M+1)) + diag((-lambda/2)*ones(1,M),1) + diag((-lambda/2)*ones(1,M),-1);
        A(1, 1) = 1;
        A(1, 2) = 0;
        A(M+1, M+1) = 1/h;
        A(M+1, M) = -1/h;
        
        b = zeros(M+1, 1);
        b(1) = U(1, i);
        b(M+1) = U(M+1, i);
        
        b(2:M) = U(1:M-1, i-1)*lambda/2 + U(2:M, i-1)*(1-lambda) + U(3:M+1, i-1)*lambda/2 + (k*phi(X(2:M), T(i)))';
        
        U(:, i) = A\b;
        
    end
end