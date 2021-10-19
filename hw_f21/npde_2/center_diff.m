function A = center_diff(n)
    % Returns the matrix of the operator D_0 acting on a n-dimensional
    % space.
    v = zeros([n 1]);
    h = 2 * pi / (n + 1);
    v(n) = -1 / (2 * h);
    v(2) = 1 / (2 * h);
    A = circulant(v);
end