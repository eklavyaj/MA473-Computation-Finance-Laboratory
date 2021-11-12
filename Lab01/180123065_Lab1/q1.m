function q1
    close all;
        
    M = 40;
    lambda = 100;
    
    N = round(M*M/lambda);
    k = 1/N;
    h = 1/M;
    Tmax = N*k;
    X = 0:h:M*h;
    T = 0:k:Tmax;
    
    disp("Computing FTCS")
    U1 = FTCS(@phi, @f, @g1, @g2, M, lambda);
    LinePlot(X, U1(:, N+1), 'x', 'u(x, T)', 'FTCS Line Plot for u(x,T)', 'q1_FTCS');
    SurfPlot(T, X, U1, 't', 'x', 'u(x, t)', 'FTCS Surface Plot for u(x,t)', 'q1_FTCS');
    error_table(lambda, @FTCS, @phi, @f, @g1, @g2, 'q1_FTCS');

    
    disp("Computing BTCS")
    U2 = BTCS(@phi, @f, @g1, @g2, M, lambda);
    LinePlot(X, U2(:, N+1), 'x', 'u(x, T)', 'BTCS Line Plot for u(x,T)', 'q1_BTCS');
    SurfPlot(T, X, U2, 't', 'x', 'u(x, t)', 'BTCS Surface Plot for u(x,t)', 'q1_BTCS');
    error_table(lambda, @BTCS, @phi, @f, @g1, @g2, 'q1_BTCS');

  
    disp("Computing Crank Nicolson")
    U3 = CN(@phi, @f, @g1, @g2, M, lambda);
    LinePlot(X, U3(:, N+1), 'x', 'u(x, T)', 'Crank-Nicolson Line Plot for u(x,T)', 'q1_CN');
    SurfPlot(T, X, U3, 't', 'x', 'u(x, t)', 'Crank-Nicolson Surface Plot for u(x,t)', 'q1_CN');
    error_table(lambda, @CN, @phi, @f, @g1, @g2, 'q1_CN');

    
end


function val = phi(x, t)
    val = sin(2*pi*x) .* sin(4*pi*t);
end

function val = f(x)
    val = 0;
end

function val = g1(t)
    val = 0;
end

function val = g2(t)
    val = 0;
end

