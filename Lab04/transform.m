function [S, Time, U] = transform(Eta, Tau, U, q)
	S = q*Eta ./ (1 - Eta);
	m = sum(S < 4);
	S = S(1:m);
	U = U(:, 1:m);
	Time = Tau;
	U = flipud(U);
	U = (S + q) .* U;
end