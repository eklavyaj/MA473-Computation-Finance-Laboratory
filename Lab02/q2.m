function q2
	close all; 
    clear;
    
    
	% Option Parameters
	T = 1;
	K = 10;
	r = 0.06;
	sigma = 0.3;
	delta = 0;

	q = 2*r/sigma^2;
	qd = 2*(r-delta)/sigma^2;

	% Computational Parameters
	x_max = 1;
	x_min = -3;

	h = 0.05;
	k = h^2/2;
	M = (x_max - x_min)/h;
	N = ceil((T*sigma^2/2)/k);

	X = x_min:h:x_max;
	tau = 0:k:T*sigma^2/2;

	S = K*exp(X);
	Time = T - 2*tau/sigma^2;

	disp("FTCS");
	U = FTCS(@phi, @f, @g1, @g2, K, q, qd, x_min, x_max, h, k, M, N, X, tau);
    LinePlot(S, U(end, :), U(1, :), 'Option Value at t = 0', 'Option Value at t = T', 'Stock Price', 'u(S, t)', 'FTCS Line Plot', 'Q2_FTCS');
    SurfPlot(S, Time, U, 'Stock Price', 'time', 'u(x, t)', 'FTCS Surface Plot', 'Q2_FTCS');

    disp("BTCS");
	U = BTCS(@phi, @f, @g1, @g2, K, q, qd, x_min, x_max, h, k, M, N, X, tau);
    LinePlot(S, U(end, :), U(1, :), 'Option Value at t = 0', 'Option Value at t = T', 'Stock Price', 'u(S, t)', 'BTCS Line Plot', 'Q2_BTCS');
    SurfPlot(S, Time, U, 'Stock Price', 'time', 'u(x, t)', 'BTCS Surface Plot', 'Q2_BTCS');

    disp("Crank-Nicolson");
	U = CN(@phi, @f, @g1, @g2, K, q, qd, x_min, x_max, h, k, M, N, X, tau);
    LinePlot(S, U(end, :), U(1, :), 'Option Value at t = 0', 'Option Value at t = T', 'Stock Price', 'u(S, t)', 'CN Line Plot', 'Q2_CN');
    SurfPlot(S, Time, U, 'Stock Price', 'time', 'u(x, t)', 'CN Surface Plot', 'Q2_CN');

	
end

function [y] = phi(x, t)
	y = 0;
end

function [y] = f(x, qd)
	temp1 = zeros(size(x));
	temp2 = exp(x*(qd - 1)/2 ) - exp(x*(qd + 1)/2);
	y = max([temp1; temp2]);
end

function [y] = g1(x, t, qd)
	y = exp(x.*(qd - 1)/2 + t.*(qd - 1)^2/4);
end

function [y] = g2(x, t, qd)
	y = 0;
end



