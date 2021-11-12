function U = FTCS(r, sigma, delta, Eta, Tau, h, k, g1, g2, f, fa, fb, fc)
    
    M = length(Eta);
	N = length(Tau);
	U = zeros(N, M);

	U(1:end, 1) = g1(r, Tau);
	U(1:end, end) = g2(delta, Tau);
	U(1, 1:end) = f(Eta);

    for i = 2:N
		for j = 2:M-1
			aa = fa(sigma, Eta(j));
			bb = fb(r, delta, Eta(j));
			cc = fc(r, delta, Eta(j));
			U(i, j) = (-aa*k/h^2 + 0.5*bb*k/h)*U(i-1,j-1) + (1 + 2*aa*k/h^2 - cc*k)*U(i-1,j) + (-aa*k/h^2 - 0.5*bb*k/h)*U(i-1,j+1);
		end
    end
    
end