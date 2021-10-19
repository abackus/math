function A = right_diff(n)
    % Returns the matrix of the operator D_+ acting on a n-dimensional
    % space.
    v = zeros([n 1]);
    h = 2 * pi / (n + 1);
    v(1) = -1 / h;
    v(2) = 1 / h;
    A = circulant(v);
end