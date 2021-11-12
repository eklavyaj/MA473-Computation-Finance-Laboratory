function q1

	close all; 
    clear;
	
    % Terminal Condition flag
	isTerminal = true;

    % Given parameters
	T = 1;
	K = 10;
	r = 0.06;
	sigma = 0.3;
	delta = 0;

	% Boundary
	S_min = 0;
	S_max = 30;

	h = 1;
	k = h^2/70;

	S = S_min:h:S_max;
	Tau = 0:k:T;

    fprintf("Computing FTCS...");
	U = FTCS(K, r, sigma, delta, S, Tau, h, k, isTerminal, @g1, @g2, @f, @fa, @fb, @fc);
    LinePlot(S, U(1,:), U(end,:), 'Cost of option at t = 0', 'Cost of option at t = T', 'S', 'u(S, t)', 'FTCS Line Plot', 'q1_FTCS');
	SurfPlot(S, Tau, U, 'S', 't', 'u(S, t)', 'FTCS Surface Plot', 'q1_FTCS');
	
	Methods = ["Matlab", "Gauss-Seidel", "Jacobi", "SOR"];
    
    for method = Methods
        fprintf(strcat("Using ", method, "\n"));
        fprintf("Computing BTCS...");
		U = BTCS(K, r, sigma, delta, S, Tau, h, k, isTerminal, @g1, @g2, @f, @fa, @fb, @fc, method);
		LinePlot(S, U(1,:), U(end,:), 'Cost of option at t = 0', 'Cost of option at t = T', 'S', 'u(S, t)', ['BTCS Line Plot using ' method], strcat('q1_BTCS_', method));
	    SurfPlot(S, Tau, U, 'S', 't', 'u(S, t)', ['BTCS Line Plot using ' method], strcat('q1_BTCS_', method));
        
        fprintf("Computing Crank-Nicolson...");
        U = CN(K, r, sigma, delta, S, Tau, h, k, isTerminal, @g1, @g2, @f, @fa, @fb, @fc, method);
		LinePlot(S, U(1,:), U(end,:), 'Cost of option at t = 0', 'Cost of option at t = T', 'S', 'u(S, t)', ['Crank-Nicolson Line Plot using ' method], strcat('q1_CN_', method));
	    SurfPlot(S, Tau, U, 'S', 't', 'u(S, t)', ['Crank-Nicolson Line Plot using ' method], strcat('q1_CN_', method));
        
    end
end

function y = f(S, K)
	temp1 = zeros(size(S));
	temp2 = S - K;
	y = max([temp1; temp2]);
end

function y = g1(r, s, Tau, K)
	y = 0;
end

function y = g2(r, s, Tau, K)
	y = s - K*exp(-r.*Tau);
end

function y = fa(sigma, s)
	y = sigma.^2 .* s.^2 / 2;
end

function y = fb(r, delta, s);
	y = (r - delta) .* s;
end

function y = fc(r)
	y = -r;
end
