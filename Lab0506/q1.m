function q1
    close all; clear;

	% Option Parameters
	T = 1;
	K = 10;
	r = 0.25;
	sig = 0.6;
	delta = 0.2;

	q = 2*r/sig^2;
	qd = 2*(r-delta)/sig^2;

    % Computational Parameters
	x_max = 2;
	x_min = -5;

	h = 0.01;
	k = 0.01;
	M = (x_max - x_min)/h;
	N = ceil((T*sig^2/2)/k);
	m_base = M;

	X = x_min:h:x_max;
	Tau = 0:k:T*sig^2/2;

	S = K*exp(X);
    indices = (1 < S) & (S < 30);
	Time = T - 2*Tau/sig^2;

    U = BTCS(@phi, @f,  @g, @g1, @g2, K, q, qd, x_min, x_max, h, k, M, N, X, Tau, @transformation);
    temp = transformation(g(X, Time(end), q, qd), X, Time(end), q, qd, K);
    LinePlot(S(indices), U(end, indices), temp(indices), 'Price of option at t = 0', '(S-K)^+ at t = 0', 'S', 'u(S, t)', 'BTCS with PSOR', 'q1_BTCS');
    SurfPlot(S(indices), Time, U(:, indices), 'S', 't', 'u(S, t)', 'BTCS with PSOR', 'q1_BTCS')
    U1_real = U(end, :);

    U = CN(@phi, @f,  @g, @g1, @g2, K, q, qd, x_min, x_max, h, k, M, N, X, Tau, @transformation);
    temp = transformation(g(X, Time(end), q, qd), X, Time(end), q, qd, K);
    LinePlot(S(indices), U(end, indices), temp(indices), 'Price of option at t = 0', '(S-K)^+ at t = 0', 'S', 'u(S, t)', 'CN with PSOR', 'q1_CN');
    SurfPlot(S(indices), Time, U(:, indices), 'S', 't', 'u(S, t)', 'CN with PSOR', 'q1_CN')
    U2_real = U(end, :);

    Is = [1, 2, 5, 10, 12.5, 25];
	N1 = []; E1 = []; E2 = [];

    for j = 1:6
        h = 0.5/Is(j);
		k = 0.5/Is(j);
		M = (x_max - x_min)/h;
		N = floor((T*sig^2/2)/k);

        N1 = [N1, M];
        X = x_min:h:x_max;
		Tau = 0:k:T*sig^2/2;
        S = K*exp(X);
		Time = T - 2*Tau/sig^2;

        U = BTCS(@phi, @f, @g, @g1, @g2, K, q, qd, x_min, x_max, h, k, M, N, X, Tau, @transformation);
        U1_bad = U(end, :);
        E1 = [E1, max(abs(U1_bad - U1_real(1:(m_base/M):end)))];

        U = CN(@phi, @f,  @g, @g1, @g2, K, q, qd, x_min, x_max, h, k, M, N, X, Tau, @transformation);
        U2_bad = U(end, :);
        E2 = [E2, max(abs(U2_bad - U2_real(1:(m_base/M):end)))];

    end
    SinglePlot(S, abs(U1_bad - U1_real(1:(m_base/M):end)), 'BTCS Error', 'S', 'Error', 'BTCS Error in U for (δx and δτ), and (δx/2 and δτ/2) at t = 0', 'q1_error1');
    SinglePlot(S, abs(U2_bad - U2_real(1:(m_base/M):end)), 'CN Error', 'S', 'Error', 'CN Error in U for (δx and δτ), and (δx/2 and δτ/2) at t = 0', 'q1_error2');
    LinePlot(N1, E1, E2, 'BTCS', 'CN','N','Error', 'Max absolute error vs N', 'q1_abs_err');

end


function y = phi(x, t)
	y = 0;
end

function y = g(x, t, q, qd)
	temp1 = zeros(size(x));
	temp2 = exp(x*(qd + 1)/2 ) - exp(x*(qd - 1)/2);

	y = exp(t*( (qd-1)^2 + 4*q )/4 ) .* max([temp1; temp2]);
end

function y = f(x, qd)
	temp1 = zeros(size(x));
	temp2 = exp(x*(qd + 1)/2 ) - exp(x*(qd - 1)/2);
	y = max([temp1; temp2]);
end

function y = g1(x, t, qd)
	y = 0;
end

function y = g2(x, t, qd)
	y = exp(x.*(qd + 1)/2 + t.*(qd + 1)^2/4);
end

function [y] = transformation(U, X, Tau, q, qd, K)
	y = zeros(size(U));
	if length(Tau) == 1
		for j = 1:length(X)
			y(j) = U(j) * K * exp(-0.5* (qd-1)*X(j) - (0.25*(qd-1)^2 + q)*Tau);
		end
	else
		for i = 1:length(Tau)
			for j = 1:length(X)
				y(i, j) = U(i, j) * K * exp(-0.5* (qd-1)*X(j) - (0.25*(qd-1)^2 + q)*Tau(i));
			end
		end
	end
end