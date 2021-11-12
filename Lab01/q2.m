function q2

    close all;
        
    M = 10;
    lambda = 4;
    
    N = ceil(M*M/lambda);
    k = 1/N;
    h = 1/M;
    Tmax = N*k;
    X = 0:h:M*h;
    T = 0:k:Tmax;
    
%     part (a)
    disp("Computing FTCS")
    U1 = FTCS(@phi, @f1, @g1, @g2, M, lambda);
    LinePlot(X, U1(:, N+1), 'x', 'u(x, T)', 'FTCS Line Plot for u(x,T)', 'q2a_FTCS');
    SurfPlot(T, X, U1, 't', 'x', 'u(x, t)', 'FTCS Surface Plot for u(x,t)', 'q2a_FTCS');
    error_table(lambda, @FTCS, @phi, @f1, @g1, @g2, 'q2a_FTCS');
    
    disp("Computing BTCS")
    U2 = BTCS(@phi, @f1, @g1, @g2, M, lambda);
    LinePlot(X, U2(:, N+1), 'x', 'u(x, T)', 'BTCS Line Plot for u(x,T)', 'q2a_BTCS');
    SurfPlot(T, X, U2, 't', 'x', 'u(x, t)', 'BTCS Surface Plot for u(x,t)', 'q2a_BTCS');
    error_table(lambda, @BTCS, @phi, @f1, @g1, @g2, 'q2a_BTCS');
  
    disp("Computing Crank Nicolson")
    U3 = CN(@phi, @f1, @g1, @g2, M, lambda);
    LinePlot(X, U3(:, N+1), 'x', 'u(x, T)', 'Crank-Nicolson Line Plot for u(x,T)', 'q2a_CN');
    SurfPlot(T, X, U3, 't', 'x', 'u(x, t)', 'Crank-Nicolson Surface Plot for u(x,t)', 'q2a_CN');
    error_table(lambda, @CN, @phi, @f1, @g1, @g2, 'q2a_CN');
   
%     part (b)
    disp("Computing FTCS")
    U1 = FTCS(@phi, @f2, @g1, @g2, M, lambda);
    LinePlot(X, U1(:, N+1), 'x', 'u(x, T)', 'FTCS Line Plot for u(x,T)', 'q2b_FTCS');
    SurfPlot(T, X, U1, 't', 'x', 'u(x, t)', 'FTCS Surface Plot for u(x,t)', 'q2b_FTCS');
    error_table(lambda, @FTCS, @phi, @f2, @g1, @g2, 'q2b_FTCS');
    
    disp("Computing BTCS")
    U2 = BTCS(@phi, @f2, @g1, @g2, M, lambda);
    LinePlot(X, U2(:, N+1), 'x', 'u(x, T)', 'BTCS Line Plot for u(x,T)', 'q2b_BTCS');
    SurfPlot(T, X, U2, 't', 'x', 'u(x, t)', 'BTCS Surface Plot for u(x,t)', 'q2b_BTCS');
    error_table(lambda, @BTCS, @phi, @f2, @g1, @g2, 'q2b_BTCS');
  
    disp("Computing Crank Nicolson")
    U3 = CN(@phi, @f2, @g1, @g2, M, lambda);
    LinePlot(X, U3(:, N+1), 'x', 'u(x, T)', 'Crank-Nicolson Line Plot for u(x,T)', 'q2b_CN');
    SurfPlot(T, X, U3, 't', 'x', 'u(x, t)', 'Crank-Nicolson Surface Plot for u(x,t)', 'q2b_CN');
    error_table(lambda, @CN, @phi, @f2, @g1, @g2, 'q2b_CN');
end


function val = phi(x, t)
    val = 0;
end

function val = f1(x)
	val = sin(pi*x);
end

function val = f2(x)
	val = x.*(1-x);
end

function val = g1(t)
    val = 0;
end

function val = g2(t)
    val = 0;
end

