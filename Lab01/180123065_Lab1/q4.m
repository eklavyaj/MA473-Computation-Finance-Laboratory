function q4
    
    close all;
    
    M = 2;
    lambda = 0.5;
    
    N = ceil(M*M/lambda);
    k = 1/N;
    h = 0.5/M;
    Tmax = N*k;
    X = 0:h:M*h;
    T = 0:k:Tmax;
    
    disp("Computing FTCS")
    U1 = FTCS_Neumann(@phi, @f, @g1, @g2, M, lambda);
    LinePlot(X, U1(:, N+1), 'x', 'u(x, T)', 'FTCS Line Plot for u(x,T)', 'q4_FTCS');
    SurfPlot(T, X, U1, 't', 'x', 'u(x, t)', 'FTCS Surface Plot for u(x,t)', 'q4_FTCS');
    error_table(lambda, @FTCS_Neumann, @phi, @f, @g1, @g2, 'q4_FTCS');
    
    disp("Computing BTCS")
    U2 = BTCS_Neumann(@phi, @f, @g1, @g2, M, lambda);
    LinePlot(X, U2(:, N+1), 'x', 'u(x, T)', 'BTCS Line Plot for u(x,T)', 'q4_BTCS');
    SurfPlot(T, X, U2, 't', 'x', 'u(x, t)', 'BTCS Surface Plot for u(x,t)', 'q4_BTCS');
    error_table(lambda, @BTCS_Neumann, @phi, @f, @g1, @g2, 'q4_BTCS');
  
    disp("Computing Crank Nicolson")
    U3 = CN_Neumann(@phi, @f, @g1, @g2, M, lambda);
    LinePlot(X, U3(:, N+1), 'x', 'u(x, T)', 'Crank-Nicolson Line Plot for u(x,T)', 'q4_CN');
    SurfPlot(T, X, U3, 't', 'x', 'u(x, t)', 'Crank-Nicolson Surface Plot for u(x,t)', 'q4_CN');
    error_table(lambda, @CN_Neumann, @phi, @f, @g1, @g2, 'q4_CN');
end


function val = phi(x, t)
    val = 0;
end

function val = f(x)
    val = x.*(1-x);
end

function val = g1(t)
    val = 0;
end

function val = g2(t)
    val = t.^2;
end

