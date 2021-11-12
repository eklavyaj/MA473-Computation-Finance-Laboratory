function q3

    close all;
        
    M = 10;
    lambda = 0.5;
    
    N = ceil(M*M/lambda);
    k = 1/N;
    h = 1/M;
    Tmax = N*k;
    X = 0:h:M*h;
    T = 0:k:Tmax;
    
    n = round(0.1/k);
    
    disp("Computing FTCS")
    U1 = FTCS_Neumann(@phi, @f, @g1, @g2, M, lambda);
    LinePlot(X, U1(:, N+1), 'x', 'u(x, T)', 'FTCS Line Plot for u(x,T)', 'q3_FTCS');
    SurfPlot(T, X, U1, 't', 'x', 'u(x, t)', 'FTCS Surface Plot for u(x,t)', 'q3_FTCS');
    error_table(lambda, @FTCS_Neumann, @phi, @f, @g1, @g2, 'q3_FTCS');
    
    disp("Computing BTCS")
    U2 = BTCS_Neumann(@phi, @f, @g1, @g2, M, lambda);
    LinePlot(X, U2(:, N+1), 'x', 'u(x, T)', 'BTCS Line Plot for u(x,T)', 'q3_BTCS');
    SurfPlot(T, X, U2, 't', 'x', 'u(x, t)', 'BTCS Surface Plot for u(x,t)', 'q3_BTCS');
    error_table(lambda, @BTCS_Neumann, @phi, @f, @g1, @g2, 'q3_BTCS');
  
    disp("Computing Crank Nicolson")
    U3 = CN_Neumann(@phi, @f, @g1, @g2, M, lambda);
    LinePlot(X, U3(:, N+1), 'x', 'u(x, T)', 'Crank-Nicolson Line Plot for u(x,T)', 'q3_CN');
    SurfPlot(T, X, U3, 't', 'x', 'u(x, t)', 'Crank-Nicolson Surface Plot for u(x,t)', 'q3_CN');
    error_table(lambda, @CN_Neumann, @phi, @f, @g1, @g2, 'q3_CN');
    
    disp("U at t = 0.1:");
    disp(U2(:,n)');
    
end


function val = phi(x, t)
    val = 0;
end

function val = f(x)
    val = cos(pi*x/2);
end

function val = g1(t)
    val = 0;
end

function val = g2(t)
    val = 0;
end

