function q1

	close all; 
    clear;

    % Given parameters
	T = 1;
	r = 0.04;
	sigma = 0.25;
	delta = 0.1;
    q = 1;

	% Boundary
	eta_min = 0;
	eta_max = 1;

    
	h = 0.01;
	k = h/2;
    
	Eta = eta_min:h:eta_max;
	Tau = 0:k:T;

    fprintf("Computing FTCS...");
	U = FTCS(r, sigma, delta, Eta, Tau, h, k, @g1, @g2, @f, @fa, @fb, @fc);
    [S, Time, U] = transform(Eta, Tau, U, q);
    LinePlot(S, U(1,:), U(end,:), 'Cost of option at t = 0', 'Cost of option at t = T', 'S', 'u(S, t)', 'FTCS Line Plot', 'q1_FTCS');
	SurfPlot(S, Time, U, 'S', 't', 'u(S, t)', 'FTCS Surface Plot', 'q1_FTCS');
	
	Methods = ["Matlab", "Gauss-Seidel", "Jacobi", "SOR", "Conjugate-Gradient"];
    
    for method = Methods
        fprintf(strcat("Using ", method, "\n"));
        fprintf("Computing BTCS...");
		U = BTCS(r, sigma, delta, Eta, Tau, h, k, @g1, @g2, @f, @fa, @fb, @fc, method);
        [S, Time, U] = transform(Eta, Tau, U, q);
		LinePlot(S, U(1,:), U(end,:), 'Cost of option at t = 0', 'Cost of option at t = T', 'S', 'u(S, t)', ['BTCS Line Plot using ' method], strcat('q1_BTCS_', method));
	    SurfPlot(S, Time, U, 'S', 't', 'u(S, t)', ['BTCS Line Plot using ' method], strcat('q1_BTCS_', method));
        
        fprintf("Computing Crank-Nicolson...");
        U = CN(r, sigma, delta, Eta, Tau, h, k, @g1, @g2, @f, @fa, @fb, @fc, method);
        [S, Time, U] = transform(Eta, Tau, U, q);
		LinePlot(S, U(1,:), U(end,:), 'Cost of option at t = 0', 'Cost of option at t = T', 'S', 'u(S, t)', ['Crank-Nicolson Line Plot using ' method], strcat('q1_CN_', method));
	    SurfPlot(S, Time, U, 'S', 't', 'u(S, t)', ['Crank-Nicolson Line Plot using ' method], strcat('q1_CN_', method));
        
    end
end

function y = f(Eta)
	temp1 = zeros(size(Eta));
	temp2 = 2*Eta - 1;
	y = max([temp1; temp2]);
end

function y = g1(r, Tau)
	y = f(0) * exp(-r*Tau);
end

function y = g2(delta, Tau)
	y = f(1) * exp(-delta*Tau);
end

function y = fa(sig, eta)
	y = - sig.^2 .* eta.^2 .* (1 - eta.^2) / 2;
end

function y = fb(r, delta, eta)
	y = - (r - delta) .* eta .* (1 - eta);
end

function y = fc(r, delta, eta)
	y = r*(1 - eta) + delta*eta;
end
