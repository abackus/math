N = 2^6;
f = grid2list(pr3_poisson_data(N), N);
u = list2grid(f / graph_laplacian(N), N);

[X,Y] = meshgrid(1/(N+2):1/(N+1):1-1/(N+2),1/(N+2):1/(N+1):1-1/(N+2));
v = (sin(X) .* cos(X)) + (1 - Y).^2;
surf(X,Y, u);
hold on;
surf(X, Y, v);
hold off;

norm = 0;
for i = 2:(N-1)
    for j = 2:(N-1)
        norm = max(norm, abs(u(i, j) - v(i, j)));
    end
end