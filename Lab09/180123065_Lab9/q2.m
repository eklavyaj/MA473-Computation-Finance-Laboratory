function q1
	close all; clear;
	% Option Parameters
	T = 1;
	K = 10;
	r = 0.06;
	sig = 0.3;
	delta = 0;

	q = 2*r/sig^2;
	qd = 2*(r-delta)/sig^2;

	% Computational Parameters
	x_max = 1;
	x_min = -5;

	h = 0.05;
	k = h^2/2;
	m = (x_max - x_min)/h;
	n = ceil((T*sig^2/2)/k);

	X = x_min:h:x_max;
	Tau = 0:k:T*sig^2/2;

	S = K*exp(X);
	Time = T - 2*Tau/sig^2;

	Methods = ['Trapezoidal rule with piecewise linear functions'; 'Simpson’s  rule  with piecewise linear functions'];
    i = 0;
	for meth = 1:2
		U = CN(@phi, @f, @g1, @g2, T, K, r, sig, delta, q, qd, x_min, x_max, h, k, m, n, X, Tau, Methods(meth, :));
		LinePlot(S, U(end, :), U(1, :), 'Cost of option at t = 0', 'Cost of option at t = T', 'S', 'u(S, t)', sprintf('Crank-Nicolson using %s method', Methods(meth, :)), ['q2_' num2str(i)])
        SurfPlot(S, Time, U, 'S', 't', 'u(S,t)', sprintf('Crank-Nicolson using %s method', Methods(meth, :)), ['q2_' num2str(i)]);
        i = i+1;
	end
end

function [y] = phi(x, t)
	y = 0;s
end

function [y] = f(x, qd)
	temp1 = zeros(size(x));
	temp2 = exp(x*(qd + 1)/2 ) - exp(x*(qd - 1)/2);
	y = max([temp1; temp2]);
end

function [y] = g1(x, t, qd)
	y = exp(x.*(qd - 1)/2 + t.*(qd - 1)^2/4);
end

function [y] = g2(x, t, qd)
	y = 0;
end
