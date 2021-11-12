
function [U] = CN(phi, f, g1, g2, T, K, r, sigma, delta, q, qd, x_min, x_max, h, k, m, n, X, Tau, method)
	fprintf('\nRunning Crank Nicolson\n');
	fprintf('Using %s method\n\n', method);
	
    lambda = k / h^2;
	U = zeros(n+1, m+1);

	U(1:end, 1) = g1(x_min, Tau, qd);
	U(1:end, end) = g2(x_max, Tau, qd);
	U(1, 1:end) = f(X, qd);

	A = zeros(m-1, m-1);

	for i = 2:n+1
		b = zeros(m+1, 1);

		A(1:m+2:end) = 1 + lambda;
		A(2:m+2:end) = -lambda/2;
		A(m+2:m+2:end) = -lambda/2;

		A(1,1) = 1;
		A(1,2) = 0;
		A(m+1,m+1) = 1;
		A(m+1,m) = 0;

		b(2:m) = U(i-1,1:m-1)*lambda/2 + (1-lambda)*U(i-1,2:m) + U(i-1,3:m+1)*lambda/2;
		b(1) = U(i,1);
		b(end) = U(i,end);

		U(i,:) = (A\b)';
	end

	U = transform(U, X, Tau, q, qd, K);
end


function [y] = transform(U, X, Tau, q, qd, K)
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
