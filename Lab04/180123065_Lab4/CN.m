function U = CN(r, sigma, delta, Eta, Tau, h, k, g1, g2, f, fa, fb, fc, method)

	M = length(Eta);
	N = length(Tau);
	U = zeros(N, M);

	U(1:end, 1) = g1(r, Tau);
	U(1:end, end) = g2(delta, Tau);
	U(1, 1:end) = f(Eta);

    for i = 2:N
        
		A = zeros(M, M);
		b = zeros(M, 1);
	
		aa = fa(sigma, Eta);
		bb = fb(r, delta, Eta);
		cc = fc(r, delta, Eta);

		A(1:M+1:end) = 2 - 2*aa*k/h^2 + cc*k;
		A(2:M+1:end) = aa(2:M)*k/h^2 - bb(2:M)*k/(2*h);
		A(M+1:M+1:end) = aa(1:M-1)*k/h^2 + bb(1:M-1)*k/(2*h);

		A(1,1) = 1;
		A(1,2) = 0;
		A(M,M) = 1;
		A(M,M-1) = 0;

		b(2:M-1) = (-aa(2:M-1)*k/h^2 + bb(2:M-1)*k/(2*h)) .* U(i-1,1:M-2) + (2 + 2*aa(2:M-1)*k/h^2 - cc(2:M-1)*k) .* U(i-1,2:M-1) + (-aa(2:M-1)*k/h^2 - bb(2:M-1)*k/(2*h)) .* U(i-1,3:M);
		b(1) = U(i,1);
		b(end) = U(i,end);

		if method(1) == 'M'
			U(i,:) = (A\b)';
            
        elseif method(1) == 'J'
			U(i,:) = Jacobi(A, b, 1000, 1e-5);
            
		elseif method(1) == 'G'
			U(i,:) = GaussSeidel(A, b, 1000, 1e-5);
        
        elseif method(1) == 'C'
			U(i,:) = ConjugateGradient(A, b, 1000, 1e-5);
            
        else
			U(i,:) = SOR(A, b, 1000, 1e-5);
            
		end			
		
        
    end
    
end