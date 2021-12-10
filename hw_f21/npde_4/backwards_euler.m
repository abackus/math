N = 20;
k = 1/20;

x = [0:1/(N+2):1];
x = x([2:N+1]);
u = zeros(N); 
for i=1:N
    for j=1:N
        u(i, j) = x(i) * (x(i) - 1) * x(j) * (x(j) - 1);
    end
end

v = grid2list(u, N);
TimeAdvance = inv(eye(N^2) + k * N^2 .* Laplace);
v = v * (TimeAdvance^(1/k));
u = list2grid(v, N);

norm = 0;
for i=1:N
    norm = norm + N^(-2) * v(i);
end
norm = sqrt(norm);