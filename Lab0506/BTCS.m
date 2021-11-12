function [U] = BTCS(phi, f, g, g1, g2, K, q, qd, x_min, x_max, h, k, M, N, X, tau, transformation)

	lambda = k / h^2;
	U = zeros(N+1, M+1);

	U(1:end, 1) = g1(x_min, tau, qd);
	U(1:end, end) = g2(x_max, tau, qd);
	U(1, 1:end) = f(X, qd);

	for i = 2:N+1
		A = zeros(M+1, M+1);
		b = zeros(M+1, 1);

		A(1:M+2:end) = 1 + 2*lambda;
		A(2:M+2:end) = -lambda;
		A(M+2:M+2:end) = -lambda;

		A(1,1) = 1;
		A(1,2) = 0;
		A(M+1,M+1) = 1;
		A(M+1,M) = 0;

		b(2:M) = U(i-1,2:M);
		b(1) = U(i,1);
		b(end) = U(i,end);

		GG = g(X, tau(i), q, qd);
		U(i,:) = psor(A,b - A*GG',1000,1e-5)' + GG;
			
			
	end

	U = transformation(U, X, tau, q, qd, K);
    
end

