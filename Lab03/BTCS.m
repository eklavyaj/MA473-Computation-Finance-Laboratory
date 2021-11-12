function U = BTCS(K, r, sigma, delta, S, Tau, h, k, isTerminal, g1, g2, f, fa, fb, fc, method)

	M = length(S);
	N = length(Tau);
	U = zeros(N, M);

	if isTerminal
		k = -k;
	end

	U(1:end, 1) = g1(r, S(1), Tau, K);
	U(1:end, end) = g2(r, S(end), Tau, K);
	U(1, 1:end) = f(S, K);

	for i = 2:N
		A = zeros(M, M);
		B = zeros(M, 1);

		a = fa(sigma, S);
		b = fb(r, delta, S);
		c = fc(r);

		A(1:M+1:end) = 1 - 2*a*k/h^2 + c*k;
		A(2:M+1:end) = a(2:M)*k/h^2 - b(2:M)*k/(2*h);
		A(M+1:M+1:end) = a(1:M-1)*k/h^2 + b(1:M-1)*k/(2*h);

		A(1,1) = 1;
		A(1,2) = 0;
		A(M,M) = 1;
		A(M,M-1) = 0;

		B(2:M-1) = U(i-1,2:M-1);
		B(1) = U(i,1);
		B(end) = U(i,end);
        
		if method(1) == 'M'
			U(i,:) = (A\B)';
            
		elseif method(1) == 'G'
			U(i,:) = GaussSeidel(A, B, 1000, 1e-5);
            
		elseif method(1) == 'J'
			U(i,:) = Jacobi(A, B, 1000, 1e-5);
            
		else
			U(i,:) = SOR(A, B, 1000, 1e-5);
            
		end			
			
	end

	if isTerminal
		U = U(end:-1:1, :);
	end
end