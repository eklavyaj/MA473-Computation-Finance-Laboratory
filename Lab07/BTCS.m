function U = BTCS(T, r, sigma, R, Tau, h, k, isTerminal, f, fa, fb, fc)

	M = length(R);
	N = length(Tau);
	U = zeros(N, M);

    if isTerminal
		k = -k;
    end

	U(1, 1:end) = f(R, T);

	for i = 2:N
		A = zeros(M, M);
		B = zeros(M, 1);

		a = fa(sigma, R);
		b = fb(r, R);
		c = fc(r);

		A(1:M+1:end) = 1 - 2*a*k/h^2 + c*k;
		A(2:M+1:end) = a(2:M)*k/h^2 - b(2:M)*k/(2*h);
		A(M+1:M+1:end) = a(1:M-1)*k/h^2 + b(1:M-1)*k/(2*h);

    %{
		A(1,1) = 1 - (3*k)/(2*h);
		A(1,2) = (4*k)/(2*h);
		A(1,3) = (-k)/(2*h);
    %}

		A(1,1) = 1 - k/h;
		A(1,2) = k/h;
        
		A(M,M) = 1;
		A(M,M-1) = 0;

		B(2:M-1) = U(i-1,2:M-1);
		B(1) = U(i-1,1);
		B(end) = 0;

		U(i,:) = (A\B)';
			
	end

	if isTerminal
		U = flipud(U);
	end
end