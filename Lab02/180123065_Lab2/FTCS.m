function [U] = FTCS(phi, f, g1, g2, K, q, qd, x_min, x_max, h, k, M, N, X, tau)

	lambda = k / h^2;
	U = zeros(N+1, M+1);

	U(1:end, 1) = g1(x_min, tau, qd);
	U(1:end, end) = g2(x_max, tau, qd);
	U(1, 1:end) = f(X, qd);

	for i = 2:N+1
		for j = 2:M
			U(i, j) = lambda*U(i-1,j-1) + (1-2*lambda)*U(i-1,j) + lambda*U(i-1,j+1);
		end
	end

	U = transformation(U, X, tau, q, qd, K);
    
end

function [y] = transformation(U, X, tau, q, qd, K)
	y = zeros(size(U));
	if length(tau) == 1
		for j = 1:length(X)
			y(j) = U(j) * K * exp(-0.5* (qd-1)*X(j) - (0.25*(qd-1)^2 + q)*tau);
		end
	else
		for i = 1:length(tau)
			for j = 1:length(X)
				y(i, j) = U(i, j) * K * exp(-0.5* (qd-1)*X(j) - (0.25*(qd-1)^2 + q)*tau(i));
			end
		end
	end
end
