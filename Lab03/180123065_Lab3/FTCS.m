function U = FTCS(K, r, sigma, delta, S, Tau, h, k, isTerminal, g1, g2, f, fa, fb, fc)

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
		for j = 2:M-1
			a = fa(sigma, S(j));
			b = fb(r, delta, S(j));
			c = fc(r);
			U(i, j) = (-a*k/h^2 + 0.5*b*k/h)*U(i-1,j-1) + (1 + 2*a*k/h^2 - c*k)*U(i-1,j) + (-a*k/h^2 - 0.5*b*k/h)*U(i-1,j+1);
		end
	end

	if isTerminal
		U = U(end:-1:1, :);
	end
end