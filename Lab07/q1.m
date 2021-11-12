function q1
    close all;
    clear;
    isTerminal = true;

    T = 0.2;
	K = 100;
	r = 0.05;
	sigma = 0.25;

	% Boundary
	R_min = 0;
	R_max = 1;

	h = 0.01;
	k = 0.001;
	m = (R_max - R_min)/h;
	n = ceil(T/k);

	R = R_min:h:R_max;
	Tau = 0:k:T;
	indices = (0 <= R) & (R < 1);
%{
    H = FTCS(T, r, sigma, R, Tau, h, k, isTerminal, @f, @fa, @fb, @fc);
    LinePlot(R(indices), H(1, indices), H(end, indices), 'H(R, 0)', 'H(R, T)', 'R', 'H(R, t)', 'H(R, t) using FTCS', 'q1_FTCS');
    SurfPlot(R, Tau, H, 'R', 't', 'H(R, t)', 'H(R, t) using FTCS', 'q1_FTCS');
%}
    H = BTCS(T, r, sigma, R, Tau, h, k, isTerminal, @f, @fa, @fb, @fc);
    LinePlot(R(indices), H(1, indices), H(end, indices), 'H(R, 0)', 'H(R, T)', 'R', 'H(R, t)', 'H(R, t) using BTCS', 'q1_BTCS_1');
    SurfPlot(R, Tau, H, 'R', 't', 'H(R, t)', 'H(R, t) using BTCS', 'q1_BTCS_1');

    H = CN(T, r, sigma, R, Tau, h, k, isTerminal, @f, @fa, @fb, @fc);
    LinePlot(R(indices), H(1, indices), H(end, indices), 'H(R, 0)', 'H(R, T)', 'R', 'H(R, t)', 'H(R, t) using Crank-Nicolson', 'q1_CN_1');
    SurfPlot(R, Tau, H, 'R', 't', 'H(R, t)', 'H(R, t) using CN', 'q1_CN_1');

end

function y = f(R, T)
	temp1 = zeros(size(R));
	temp2 = 1 - R./T;
	y = max([temp1; temp2]);
end

function y = fa(sig, R)
	y = sig.^2 .* R.^2 / 2;
end

function y = fb(r, R);
	y = 1 - r*R;
end

function y = fc(r)
	y = 0;
end