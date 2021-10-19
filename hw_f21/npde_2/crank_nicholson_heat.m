function u = crank_nicholson_heat(f, k, T)
    % Runs a Crank-Nicholson method on the heat equation.
    n = length(f); % dimension of the state space
    Laplace = right_diff(n) * left_diff(n); % discrete Laplacian
    Q = (eye(n) - (k/2) * Laplace) \ (eye(n) + (k/2) * Laplace);
    
    % Q is circulant so we can diagonalize to take powers
    A = dftmtx(n);
    D = A \ Q * A;
    D = diag(diag(D).^(T/k));
    Q = real(A * D / A);
    
    u = Q * f;
end