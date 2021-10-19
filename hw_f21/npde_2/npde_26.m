% A script that runs the method of lines and RK4 on a heat equation.

h = 1/100;
T = 1;
k = 1/300;

x = transpose(0:h:2*pi);
f = x - pi;

D = center_diff(length(f));
Laplace = right_diff(length(f)) * left_diff(length(f));

% Q = D - (h^2/6) * D * Laplace + (h^4/30) * D * Laplace^2; % Q_6
Q = D - (h^2/6) * D * Laplace; % Q_4

% save time pre-doing all the matrix multiplication
Q2 = Q^2;
Q3 = Q^3;
Q4 = Q^4;

u = f;
for i=1:round(T/k)
    u = u + k * Q * u + (k^2/2) * Q2 * u + (k^3/6) * Q3 * u + (k^4/24) * Q4 * u;
end

plot(x, f, x, u);